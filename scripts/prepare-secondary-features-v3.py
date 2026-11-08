"""Prepare one authenticated source-panel block/axis with training-only families.

This component accepts no label argument, fits no interaction model, and leaves
final-source/native-query orchestration to a separately frozen future runner.
"""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import platform
import subprocess
import sys

from tipcraft import partition, protein_features
from tipcraft.learning import Sample

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/secondary-feature-policy-v3.json"
CONTROL = "__mapping_control__"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def encoded(value):
    return (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def bound(path):
    return {"path": str(path), "sha256": sha(path)}


def check(binding):
    path = ROOT / binding["path"]
    if sha(path) != binding["sha256"]:
        raise ValueError("Authenticated file differs: " + str(path))
    return path


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def selection(panel, block_id, axis):
    if axis not in {"host", "phage"}:
        raise ValueError("Axis must be host or phage")
    if (panel["status"] != "frozen_secondary_reference_panel" or panel["annotation_join_complete"] is not True
            or panel["ready_for_feature_preparation"] is not True or panel["training_authorized"] is not False):
        raise ValueError("Frozen eligible reference-panel membership/annotation join required")
    samples = [Sample(**row) for row in panel["samples"]]
    rebuilt = partition.build_blocked_folds(samples, panel["host_group_folds"], panel["phage_group_folds"])
    if panel["blocks"] != rebuilt:
        raise ValueError("Panel train/test/embargo blocks differ from label-independent reconstruction")
    matches = [b for b in rebuilt if b["fold_id"] == block_id]
    if len(matches) != 1:
        raise ValueError("Exactly one explicit source-panel block required")
    block = matches[0]
    by_id = {s.row_id: s for s in samples}
    field, group_field = axis + "_id", axis + "_group"
    train = sorted({getattr(by_id[i], field) for i in block["train_ids"]})
    query = sorted({getattr(by_id[i], field) for i in block["test_ids"]})
    tg = {getattr(by_id[i], group_field) for i in block["train_ids"]}
    qg = {getattr(by_id[i], group_field) for i in block["test_ids"]}
    if not train or not query or set(train) & set(query) or tg & qg or CONTROL in train + query:
        raise ValueError("Nonempty disjoint training/query genome groups required")
    genomes = panel[axis + "_genomes"]
    for name in train + query:
        genome = genomes[name]
        if genome["excluded"] is not False or "annotation" not in genome:
            raise ValueError("Excluded or unannotated reference cannot enter feature preparation")
        if {getattr(s, group_field) for s in samples if getattr(s, field) == name} != {genome["group"]}:
            raise ValueError("Genome group differs from Sample membership")
    return {"block_id": block_id, "axis": axis, "training_genomes": train, "query_genomes": query,
            "mapping_control_source": train[0], "training_group_count": len(tg), "query_group_count": len(qg)}


def normalize(panel, selected, parse):
    """Retain per-genome copies/order; neutralize IDs after strict FAA parsing."""
    train, query, records, bindings = {}, {}, [], {}
    genomes = panel[selected["axis"] + "_genomes"]
    for role, names, destination in [("train", selected["training_genomes"], train),
                                     ("query", selected["query_genomes"], query)]:
        for index, genome in enumerate(names):
            annotation = genomes[genome]["annotation"]
            for key in ["FAA", "GFF", "normalized_input", "contig_map", "protein_coordinate_inventory"]:
                path = check(annotation[key])
                bindings[path] = annotation[key]["sha256"]
            proteins = parse((ROOT / annotation["FAA"]["path"]).read_bytes())
            if len(proteins) != annotation["proteins"]:
                raise ValueError("FAA count differs from frozen annotation inventory")
            for number, (source_id, sequence) in enumerate(proteins):
                identifier = f"{selected['axis'][0]}{role[0]}g{index:04d}p{number:07d}"
                destination[identifier] = sequence
                records.append({"protein": identifier, "source_protein_id": source_id, "genome": genome,
                                "role": role, "length": len(sequence), "sequence_sha256": digest(sequence.encode())})
                if role == "train" and genome == selected["mapping_control_source"]:
                    copied = "control_" + identifier
                    query[copied] = sequence
                    records.append({"protein": copied, "source_protein_id": source_id, "genome": CONTROL,
                                    "role": "mapping_control", "copied_training_protein": identifier,
                                    "length": len(sequence), "sequence_sha256": digest(sequence.encode())})
    if not train or not query or set(train) & set(query) or len(records) != len(train) + len(query):
        raise ValueError("Incomplete or overlapping normalized protein inventory")
    return train, query, records, bindings


def command_plan(policy, workdir):
    substitutions = {"mmseqs": str(ROOT / policy["mmseqs_executable"]), "fold_dir": str(workdir)}
    commands = [{"argv": [v.format(**substitutions) for v in step["argv"]],
                 "stdin": step["stdin"].format(**substitutions) if step["stdin"] else None,
                 "input_prefixes": [v.format(**substitutions) for v in step["input_prefixes"]],
                 "output_prefixes": [v.format(**substitutions) for v in step["output_prefixes"]]}
                for step in policy["command_templates"]]
    if len(commands) != 9 or [s["argv"][1] for s in commands] != ["createdb", "cluster", "createtsv", "createdb", "search", "createtsv", "createdb", "search", "createtsv"]:
        raise ValueError("Nine fixed symmetric training/query commands required")
    if (commands[3]["stdin"] != commands[0]["stdin"]
            or commands[4]["argv"][2] == commands[4]["argv"][3]
            or commands[7]["argv"][2] == commands[7]["argv"][3]
            or commands[4]["argv"][3] != commands[7]["argv"][3]):
        raise ValueError("Both projections require independent query prefixes and one training target")
    return commands


def check_control(records, self_hits, query_hits, control_vector, training_vector):
    """Compare mapping outcomes; matched absence is not family identification."""
    comparisons = []
    original_ids, copied_ids = set(), set()
    counts = {"same_assigned_family": 0, "both_unassigned": 0, "assignment_discrepancy": 0}
    for record in records:
        if record["role"] != "mapping_control":
            continue
        original, copied = record["copied_training_protein"], record["protein"]
        if original in original_ids or copied in copied_ids:
            raise ValueError("Duplicate copied-control comparison")
        original_ids.add(original)
        copied_ids.add(copied)
        expected, actual = self_hits.get(original), query_hits.get(copied)
        if expected is None and actual is None:
            state = "both_unassigned"
        elif expected is not None and actual is not None and expected["family"] == actual["family"]:
            state = "same_assigned_family"
        else:
            state = "assignment_discrepancy"
        counts[state] += 1
        comparisons.append({"copied_protein": copied, "training_protein": original,
                            "expected_family": expected["family"] if expected else None,
                            "actual_family": actual["family"] if actual else None,
                            "state": state, "mapping_outcome_equal": state != "assignment_discrepancy"})
    return {"passed": bool(comparisons) and not counts["assignment_discrepancy"] and control_vector == training_vector,
            "control_contract": "symmetric_mapping_outcomes_and_vector_equality_v3",
            "protein_mapping_comparisons": comparisons, **counts,
            "compared_proteins": len(comparisons), "vector_equal": control_vector == training_vector,
            "nonzero_vector": any(control_vector)}


def execute(policy_path, panel_path, expected_panel_sha, block_id, axis, workdir, output, plan_only=False):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if sha(__file__) != policy["helper_sha256"] or policy["workers"] != 1 or policy["timeout_seconds_per_command"] != 1800:
        raise ValueError("Feature helper/resource policy differs")
    bindings = {check(item): item["sha256"] for item in policy["inputs"].values()}
    bindings[policy_path], bindings[Path(__file__)] = digest(policy_bytes), sha(__file__)
    bindings[check({"path": policy["mmseqs_executable"], "sha256": policy["mmseqs_sha256"]})] = policy["mmseqs_sha256"]
    template_policy = json.loads((ROOT / policy["inputs"]["template_policy"]["path"]).read_text())
    if (policy["command_templates"] != template_policy["command_templates"]
            or any(policy[k] != template_policy[k] for k in ["mmseqs_executable", "mmseqs_sha256"])
            or sha(protein_features.__file__) != policy["inputs"]["kernel"]["sha256"]
            or sha(partition.__file__) != policy["inputs"]["partition"]["sha256"]):
        raise ValueError("Pinned templates/kernel/partition differ")
    if not panel_path.exists() or not expected_panel_sha:
        return {"schema_version": 1, "status": "HOLD_missing_frozen_panel", "execution_started": False,
                "panel_exists": panel_path.exists(), "externally_frozen_panel_hash_supplied": bool(expected_panel_sha)}
    bindings[check({"path": str(panel_path), "sha256": expected_panel_sha})] = expected_panel_sha
    panel = json.loads(panel_path.read_text())
    if (panel["helper_sha256"] != policy["inputs"]["panel_helper"]["sha256"]
            or panel["policy_sha256"] != policy["inputs"]["panel_policy"]["sha256"]):
        raise ValueError("Frozen panel producer/policy differs")
    for item in panel["inputs"].values():
        bindings[check(item)] = item["sha256"]
    selected = selection(panel, block_id, axis)
    helper = load("pinned_secondary_feature_ledger", ROOT / policy["inputs"]["ledger_helper"]["path"])
    parser = load("pinned_secondary_feature_parser", ROOT / policy["inputs"]["protein_parser"]["path"]).proteins
    helper.guard_destinations(workdir, output, bindings)
    train, query, records, source_bindings = normalize(panel, selected, parser)
    bindings.update(source_bindings)
    helper.guard_destinations(workdir, output, bindings)
    train_bytes = "".join(f">{p}\n{s}\n" for p, s in train.items()).encode()
    query_bytes = "".join(f">{p}\n{s}\n" for p, s in query.items()).encode()
    map_bytes = encoded(records)
    commands = command_plan(policy, workdir)
    plan = {**selected, "commands": commands, "panel": {"path": str(panel_path), "sha256": expected_panel_sha},
            "source_bindings": [{"path": str(p), "sha256": h} for p, h in sorted(bindings.items())],
            "normalized_sha256": {"training.faa": digest(train_bytes), "queries.faa": digest(query_bytes), "protein-map.json": digest(map_bytes)},
            "training_proteins": len(train), "query_proteins_including_control": len(query),
            "policy_sha256": digest(policy_bytes), "helper_sha256": sha(__file__),
            "labels_used_for_feature_preparation": False, "execution_started": False}
    if plan_only:
        helper.verify(bindings)
        return {"schema_version": 1, "status": "planned_only", **plan}
    workdir.mkdir(parents=True, exist_ok=False)
    for filename, content in [("training.faa", train_bytes), ("queries.faa", query_bytes), ("protein-map.json", map_bytes),
                              ("frozen-plan.json", encoded(plan))]:
        path = workdir / filename
        path.write_bytes(content)
        bindings[path] = sha(path)
    normalized = [bound(workdir / n) for n in ["training.faa", "queries.faa", "protein-map.json"]]
    lengths, qlengths = {p: len(s) for p, s in train.items()}, {p: len(s) for p, s in query.items()}
    tmap = {r["protein"]: r["genome"] for r in records if r["role"] == "train"}
    qmap = {r["protein"]: r["genome"] for r in records if r["role"] != "train"}
    ledger, local = [], dict(bindings)
    result = {"schema_version": 1, "status": "failed_secondary_axis_features", **selected,
              "panel": plan["panel"], "helper_sha256": plan["helper_sha256"], "policy_sha256": plan["policy_sha256"],
              "frozen_plan": bound(workdir / "frozen-plan.json"), "normalized_input_bindings": normalized,
              "source_bindings": plan["source_bindings"], "commands": ledger, "mapping_control_passed": False,
              "execution_started": True, "execution_context": {"architecture": platform.machine(), "platform": platform.platform()},
              "model_fitted_or_scored": False, "new_original_benchmark_coverage": 0,
              "final_source_native_query_mode": "pending_separate_orchestration", "labels_used_for_feature_preparation": False}
    try:
        for index, step in enumerate(commands, 1):
            helper.run_command(step, index, workdir, 1800, local, ledger)
            if index == 6:
                with (workdir / "clusters.tsv").open() as stream:
                    families = protein_features.read_families(stream, lengths)
                with (workdir / "self.tsv").open() as stream:
                    self_hits, self_diagnostics = protein_features.best_hits(stream, lengths, lengths, families)
                vocabulary = protein_features.fit_vocabulary(selected["training_genomes"], tmap, families, self_hits)
                helper.dump(workdir / "vocabulary.json", vocabulary)
                helper.dump(workdir / "best-self.json", self_hits)
                helper.dump(workdir / "pre-query-freeze.json", {"completed_commands": 6, "query_createdb_started": False,
                            "vocabulary": bound(workdir / "vocabulary.json"), "best_self": bound(workdir / "best-self.json"),
                            "training_inputs": [bound(workdir / n) for n in ["training.faa", "clusters.tsv", "self.tsv"]]})
                local.update({workdir / n: sha(workdir / n) for n in ["vocabulary.json", "best-self.json", "pre-query-freeze.json"]})
        with (workdir / "query.tsv").open() as stream:
            hits, query_diagnostics = protein_features.best_hits(stream, qlengths, lengths, families)
        names = selected["training_genomes"] + selected["query_genomes"]
        features, vectors = protein_features.project(vocabulary, names + [CONTROL], tmap | qmap, self_hits | hits)
        control_vector = vectors.pop(CONTROL)
        control = check_control(records, self_hits, hits, control_vector, vectors[selected["mapping_control_source"]])
        passed = control["passed"]
        helper.dump(workdir / "mapping-control-check.json", control)
        helper.dump(workdir / "best-query.json", hits)
        matrix = io.StringIO(newline="")
        writer = csv.writer(matrix, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", *features])
        writer.writerows([name, *vectors[name]] for name in names)
        (workdir / "feature-matrix.tsv").write_text(matrix.getvalue())
        helper.verify(local)
        result.update({"status": "completed_secondary_axis_features" if passed else "failed_mapping_control",
                       "feature_names": features, "vectors": vectors, "vocabulary": bound(workdir / "vocabulary.json"),
                       "pre_query_freeze": bound(workdir / "pre-query-freeze.json"), "mapping_control_passed": passed,
                       "mapping_control": {"source": selected["mapping_control_source"], "vector": control_vector,
                                           "mapping_outcome_and_vector_equality": passed,
                                           "same_assigned_family": control["same_assigned_family"],
                                           "both_unassigned": control["both_unassigned"],
                                           "assignment_discrepancies": control["assignment_discrepancy"], "nonzero_vector": control["nonzero_vector"],
                                           "check": bound(workdir / "mapping-control-check.json"), "training_member": False},
                       "training_projection_mode": "ordinary_external_query_database_v3",
                       "control_contract": control["control_contract"],
                       "training_hit_diagnostics": self_diagnostics, "query_hit_diagnostics": query_diagnostics,
                       "zero_vector_genomes": [n for n in names if not any(vectors[n])],
                       "protein_counts": dict(Counter((tmap | qmap).values())),
                       "ready_for_fitting": passed and bool(features), "biological_interpretation": "Representation/mapping control only; no binding or host-specificity validation"})
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        result["error"] = str(exc)
    # Preserve all completed raw commands on failure; never rerun automatically.
    result["artifacts"] = [bound(p) for p in sorted(workdir.iterdir()) if p.is_file()]
    result["command_logs_retained"] = bool(ledger)
    result["resource_scope"] = "Serial commands; child CPU and cumulative direct-child maximum RSS, not simultaneous aggregate memory"
    output.parent.mkdir(parents=True, exist_ok=True)
    helper.dump(output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--expected-panel-sha256")
    parser.add_argument("--block-id", required=True)
    parser.add_argument("--axis", choices=["host", "phage"], required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    result = execute(args.policy, args.panel, args.expected_panel_sha256, args.block_id, args.axis, args.workdir, args.output, args.plan_only)
    print(json.dumps(result if args.plan_only else {"status": result["status"], "execution_started": result.get("execution_started", False)}))
    return 0 if result["status"] in {"planned_only", "completed_secondary_axis_features"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
