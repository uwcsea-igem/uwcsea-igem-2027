"""Execute a frozen source-phage-only grouped feature/resource pilot.

No host features, phenotype labels, interaction models or target scores enter
this component pilot. Run from the repository root with uv --project.
"""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time

from tipcraft import partition, protein_features

POLICY = Path("software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-fold-feature-policy.json")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def verify(bindings):
    for path, expected in bindings.items():
        if sha(path) != expected:
            raise ValueError(f"Authenticated input/implementation changed: {path}")


def guard_destinations(workdir, output, protected):
    destinations = [Path(workdir).resolve(), Path(output).resolve()]
    if any(Path(p).exists() or Path(p).is_symlink() for p in [workdir, output]):
        raise ValueError("Workdir and summary output must both be new")
    a, b = destinations
    if a == b or a in b.parents or b in a.parents:
        raise ValueError("Workdir and summary output must be separate")
    for destination in destinations:
        for source in protected:
            source = Path(source).resolve()
            if destination == source or source in destination.parents or destination in source.parents:
                raise ValueError("Output overlaps a protected input")
            if source.parent.name == "pilot-v1" and source.parent in destination.parents:
                raise ValueError("Output overlaps historical pilot directory")


def build_plan(policy, group_report, workdir):
    names = group_report["retained_model_table_phages"]
    if len(names) != 93 or names != sorted(set(names)):
        raise ValueError("Expected exactly the frozen 93 source phages")
    groups = group_report["model_phage_groups"]
    assignment = partition.assign_group_folds({groups[n] for n in names}, axis=policy["axis"],
                                               folds=policy["fold_count"], seed=policy["seed"])
    if assignment != policy["group_assignments"] or policy["fold_count"] != 3:
        raise ValueError("Frozen group assignments differ; no automatic regrouping")
    plans = []
    for i, frozen in enumerate(policy["folds"]):
        train = [n for n in names if assignment[groups[n]] != i]
        heldout = [n for n in names if assignment[groups[n]] == i]
        if (frozen != {"fold": i, "training_phages": train, "heldout_phages": heldout,
                       "mapping_control_source": min(train)} or not train or not heldout
                or {groups[n] for n in train} & {groups[n] for n in heldout}):
            raise ValueError("Frozen membership/group boundary differs")
        substitutions = {"mmseqs": policy["mmseqs_executable"], "fold_dir": str(Path(workdir) / f"fold-{i}")}
        commands = [{"argv": [v.format(**substitutions) for v in step["argv"]],
                     "stdin": step["stdin"].format(**substitutions) if step["stdin"] else None,
                     "input_prefixes": [v.format(**substitutions) for v in step["input_prefixes"]],
                     "output_prefixes": [v.format(**substitutions) for v in step["output_prefixes"]]}
                    for step in policy["command_templates"]]
        if len(commands) != 8:
            raise ValueError("Expected eight frozen MMseqs commands per fold")
        plans.append(frozen | {"training_group_count": len({groups[n] for n in train}),
                               "heldout_group_count": len({groups[n] for n in heldout}), "commands": commands})
    if len(plans) != 3 or Counter(n for f in plans for n in f["heldout_phages"]) != Counter(names):
        raise ValueError("Every source phage must be held out exactly once")
    return {"group_assignments": assignment, "folds": plans}


def authenticated(policy_path):
    policy = json.loads(policy_path.read_text())
    bindings = {Path(i["path"]): i["sha256"] for i in policy["inputs"] + policy["implementation_inputs"]}
    bindings[policy_path], bindings[Path(__file__)] = sha(policy_path), sha(Path(__file__))
    bindings[Path(policy["mmseqs_executable"])] = policy["mmseqs_sha256"]
    verify(bindings)
    spec = importlib.util.spec_from_file_location("authenticated_kernel_replay", "scripts/check-protein-feature-kernel.py")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    _, old_policy, files, old_bindings = helper.authenticated_inputs()
    bindings.update(old_bindings)
    roster = json.loads(Path("software/tipcraft/benchmarks/retro-t7-t3/genophi-roster-policy.json").read_text())
    if (roster["group_splits"]["phage_folds"] != policy["fold_count"]
            or roster["group_splits"]["seed"] != policy["seed"] or policy["axis"] != "phage"):
        raise ValueError("Fold settings differ from frozen roster policy")
    groups = json.loads(Path("data/processed/genophi-phage-groups.json").read_text())
    records = json.loads(files["protein-map.json"].read_text())
    if len({r["protein"] for r in records}) != len(records):
        raise ValueError("Duplicate historical protein identifiers")
    records = {r["protein"]: r for r in records if r["role"] == "train"}
    parse = helper.parser()
    helper.checked_proteins(files["training.faa"], records, parse)
    sequences = dict(parse(files["training.faa"].read_bytes()))
    if (sorted({r["genome"] for r in records.values()}) != groups["retained_model_table_phages"]
            or old_policy["training_phages"] != groups["retained_model_table_phages"]):
        raise ValueError("Authenticated proteomes differ from the frozen eligible source set")
    return policy, groups, records, sequences, bindings


def prepare_fold(plan, records, sequences, workdir):
    directory = Path(workdir) / f"fold-{plan['fold']}"
    directory.mkdir()
    train, query, train_map, query_map = {}, {}, {}, {}
    for protein, record in sorted(records.items()):
        genome = record["genome"]
        if genome in plan["training_phages"]:
            train[protein], train_map[protein] = sequences[protein], genome
        elif genome in plan["heldout_phages"]:
            query[protein], query_map[protein] = sequences[protein], genome
        else:
            raise ValueError("Unplanned source protein")
        if genome == plan["mapping_control_source"]:
            copied = f"control{plan['fold']}_{protein}"
            query[copied], query_map[copied] = sequences[protein], "mapping_control"
    if set(train) & set(query) or set(train_map.values()) & set(query_map.values()):
        raise ValueError("Training/query protein or genome identities overlap")
    for name, proteins in [("training.faa", train), ("queries.faa", query)]:
        (directory / name).write_text("".join(f">{p}\n{s}\n" for p, s in proteins.items()))
    dump(directory / "protein-map.json", {"training": train_map, "query": query_map})
    return train, query, train_map, query_map


def artifacts(prefixes):
    result = {}
    for prefix in prefixes:
        path = Path(prefix)
        found = [p for p in path.parent.glob(path.name + "*") if p.is_file()]
        if not found or any(p.is_symlink() for p in found):
            raise ValueError(f"Missing or symlinked command artifact: {prefix}")
        result.update({str(p): sha(p) for p in sorted(found)})
    return result


def run_command(step, index, directory, timeout, bindings, ledger):
    verify(bindings)
    before = artifacts(step["input_prefixes"])
    log = directory / f"step-{index:02d}.log"
    start, usage = time.monotonic(), resource.getrusage(resource.RUSAGE_CHILDREN)
    print(f"START {directory.name} step {index}/8 {step['argv'][1]}", flush=True)
    with log.open("xb") as stream:
        stdin = Path(step["stdin"]).open("rb") if step["stdin"] else None
        try:
            child = subprocess.Popen(step["argv"], stdin=stdin, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                code = -signal.SIGKILL
        finally:
            if stdin:
                stdin.close()
    end = resource.getrusage(resource.RUSAGE_CHILDREN)
    entry = {"command": step["argv"], "stdin": step["stdin"], "input_artifacts": before,
             "output_artifacts": artifacts(step["output_prefixes"]) if code == 0 else {},
             "exit_code": code, "log": str(log), "log_sha256": sha(log),
             "wall_seconds": time.monotonic() - start, "child_user_cpu_seconds": end.ru_utime - usage.ru_utime,
             "child_system_cpu_seconds": end.ru_stime - usage.ru_stime,
             "cumulative_max_child_rss_bytes": end.ru_maxrss * (1 if sys.platform == "darwin" else 1024)}
    ledger.append(entry)
    (directory / "command-ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
    verify({Path(p): h for p, h in before.items()})
    if code:
        raise ValueError(f"Command failed/timed out: {directory.name} step {index}, exit {code}")
    bindings.update({Path(p): h for p, h in entry["output_artifacts"].items()})


def run_fold(plan, data, workdir, policy, bindings):
    train, query, train_map, query_map = data
    directory = Path(workdir) / f"fold-{plan['fold']}"
    lengths, qlengths = {p: len(s) for p, s in train.items()}, {p: len(s) for p, s in query.items()}
    ledger, vocabulary, self_hits, self_diagnostic = [], None, None, None
    local = dict(bindings)
    for index, step in enumerate(plan["commands"], 1):
        run_command(step, index, directory, policy["timeout_seconds_per_command"], local, ledger)
        if index == 5:
            with (directory / "clusters.tsv").open() as stream:
                families = protein_features.read_families(stream, lengths)
            with (directory / "self.tsv").open() as stream:
                self_hits, self_diagnostic = protein_features.best_hits(stream, lengths, lengths, families)
            if set(self_hits) != set(train):
                raise ValueError("Training self-query coverage incomplete")
            vocabulary = protein_features.fit_vocabulary(plan["training_phages"], train_map, families, self_hits)
            dump(directory / "vocabulary.json", vocabulary)
            dump(directory / "pre-query-freeze.json", {"completed_commands": 5, "query_createdb_started": False,
                 "vocabulary_sha256": sha(directory / "vocabulary.json"), "training_proteins": len(train),
                 "training_groups_disjoint_from_heldout": True})
            local.update({directory / n: sha(directory / n) for n in ["vocabulary.json", "pre-query-freeze.json"]})
    with (directory / "query.tsv").open() as stream:
        query_hits, diagnostic = protein_features.best_hits(stream, qlengths, lengths, families)
    names = plan["training_phages"] + plan["heldout_phages"] + ["mapping_control"]
    features, vectors = protein_features.project(vocabulary, names, train_map | query_map, self_hits | query_hits)
    verify(local)
    matrix = io.StringIO(newline="")
    writer = csv.writer(matrix, delimiter="\t", lineterminator="\n")
    writer.writerow(["phage", *features])
    writer.writerows([n, *vectors[n]] for n in names)
    (directory / "feature-matrix.tsv").write_text(matrix.getvalue())
    dump(directory / "best-query.json", query_hits)
    dump(directory / "best-self.json", self_hits)
    copies, assigned = Counter((train_map | query_map).values()), Counter((train_map | query_map)[p] for p in self_hits | query_hits)
    control = bool(any(vectors[plan["mapping_control_source"]]) and vectors["mapping_control"] == vectors[plan["mapping_control_source"]])
    return {"fold": plan["fold"], "status": "passed" if control else "failed_mapping_control",
            "training_phages": plan["training_phages"], "heldout_phages": plan["heldout_phages"],
            "training_group_count": plan["training_group_count"], "heldout_group_count": plan["heldout_group_count"],
            "training_proteins": len(train), "query_proteins": len(query), "features": len(features),
            "clustered_families": vocabulary["clustered_families"],
            "retained_families": len(vocabulary["family_features"]), "all_training_self_queries_represented": True,
            "mapping_control": {"source": plan["mapping_control_source"], "nonzero_exact_match": control},
            "proteomes": [{"phage": n, "role": "train" if n in plan["training_phages"] else "control" if n == "mapping_control" else "heldout",
                           "proteins": copies[n], "assigned_proteins": assigned[n], "assigned_fraction": assigned[n] / copies[n],
                           "retained_feature_count": sum(vectors[n]), "zero_vector": not any(vectors[n])} for n in names],
            "training_hit_diagnostics": self_diagnostic, "query_hit_diagnostics": diagnostic,
            "commands": ledger, "artifacts": {str(directory / n): sha(directory / n) for n in
                ["training.faa", "queries.faa", "protein-map.json", "clusters.tsv", "self.tsv", "query.tsv",
                 "vocabulary.json", "pre-query-freeze.json", "feature-matrix.tsv", "best-query.json", "best-self.json"]}}


def execute(policy_path, workdir, output, plan_only=False):
    policy, groups, records, sequences, bindings = authenticated(policy_path)
    guard_destinations(workdir, output, bindings)
    plan = build_plan(policy, groups, workdir)
    if plan_only:
        return {"status": "planned_only", **plan}
    workdir.mkdir(parents=True)
    data = [prepare_fold(fold, records, sequences, workdir) for fold in plan["folds"]]
    for fold in plan["folds"]:
        for n in ["training.faa", "queries.faa", "protein-map.json"]:
            path = workdir / f"fold-{fold['fold']}" / n
            bindings[path] = sha(path)
    plan["frozen_inputs"] = {str(p): h for p, h in bindings.items()}
    dump(workdir / "frozen-plan.json", plan)
    bindings[workdir / "frozen-plan.json"] = sha(workdir / "frozen-plan.json")
    report = {"schema_version": 1, "status": "running", "purpose": policy["purpose"],
              "frozen_plan_sha256": sha(workdir / "frozen-plan.json"), "authenticated_inputs": plan["frozen_inputs"],
              "execution_context": {"platform": platform.platform(), "architecture": platform.machine(), "python": sys.version.split()[0]},
              "resource_scope": "Serial command timings/CPU; cumulative maximum direct-child RSS, not aggregate workflow memory",
              "software_deviation": policy["software_deviation"], "folds": [], "raw_workdir": str(workdir),
              "phenotype_tables_loaded": False, "host_features_prepared": False, "interaction_model_fitted": False,
              "target_scores_generated": False, "production_run_spec_complete": False,
              "interpretation": policy["still_incomplete"]}
    try:
        for fold, inputs in zip(plan["folds"], data):
            result = run_fold(fold, inputs, workdir, policy, bindings)
            report["folds"].append(result)
            bindings.update({Path(p): h for p, h in result["artifacts"].items()})
            if result["status"] != "passed":
                raise ValueError("Copied nonzero mapping control failed; no retuning")
        verify(bindings)
        report["status"] = "passed"
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        report["status"], report["error"] = "failed", str(exc)
    output.parent.mkdir(parents=True, exist_ok=True)
    dump(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plan-only", action="store_true", help="Authenticate and print the frozen plan; no output writes or MMseqs")
    args = parser.parse_args()
    result = execute(args.policy, args.workdir, args.output, args.plan_only)
    print(json.dumps(result if args.plan_only else {"status": result["status"], "completed_folds": len(result["folds"])}))
    raise SystemExit(1 if result["status"] == "failed" else 0)
