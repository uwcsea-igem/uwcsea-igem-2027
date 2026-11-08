"""Final source-only protein vocabulary with fixed native reference queries.

This mode refits on the full admitted source roster. It is distinct from source
outer-fold validation and accepts neither source-CV performance nor native labels.
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

from tipcraft import partition, protein_features
from tipcraft.learning import Sample

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/native-transfer-feature-policy.json"
CONTROL = "__mapping_control__"
NATIVE = {"host": ["BL21", "BW25113", "MG1655"], "phage": ["T7", "T3"]}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def check(binding):
    path = ROOT / binding["path"]
    if sha(path) != binding["sha256"]:
        raise ValueError("Authenticated input differs: " + str(path))
    return path


def bound(path):
    return {"path": str(path), "sha256": sha(path)}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def selection(panel, design, axis):
    if axis not in NATIVE:
        raise ValueError("Axis must be host or phage")
    if (panel["status"] != "frozen_secondary_reference_panel" or panel["annotation_join_complete"] is not True
            or panel["ready_for_feature_preparation"] is not True or panel["training_authorized"] is not False):
        raise ValueError("Complete authenticated admitted source panel required")
    samples = [Sample(**row) for row in panel["samples"]]
    if partition.build_blocked_folds(samples, panel["host_group_folds"], panel["phage_group_folds"]) != panel["blocks"]:
        raise ValueError("Source membership/partition integrity differs")
    if design["execute_if_source_cv_disappoints"] is not True:
        raise ValueError("Final mode must not condition execution on source-CV performance")
    for kind, names in NATIVE.items():
        if design["native_candidates"][kind + "s"] != names:
            raise ValueError("Fixed native candidate order differs")
    queries = design["native_query_proteomes"]
    expected = {(kind, name) for kind, names in NATIVE.items() for name in names}
    if (len(queries) != 5 or {(r["kind"], r["source_id"]) for r in queries} != expected
            or any(r["query_only"] is not True for r in queries)):
        raise ValueError("All and only five native query-only proteomes required")
    field = axis + "_id"
    training = sorted({getattr(s, field) for s in samples})
    genomes = panel[axis + "_genomes"]
    if training != sorted(n for n, g in genomes.items() if not g["excluded"]):
        raise ValueError("All and only admitted observed source genomes must train final vocabulary")
    for name in training:
        genome = genomes[name]
        if ("annotation" not in genome or genome["excluded"] is not False
                or {getattr(s, axis + "_group") for s in samples if getattr(s, field) == name} != {genome["group"]}):
            raise ValueError("Unannotated or inconsistent source genotype membership")
    query = list(NATIVE[axis])
    if not training or set(training) & set(query) or CONTROL in training + query:
        raise ValueError("Native queries/control must never be source vocabulary members")
    return {"mode": "final_source_native_query", "axis": axis, "training_genomes": training,
            "query_genomes": query, "mapping_control_source": training[0],
            "training_group_count": len({getattr(s, axis + "_group") for s in samples}),
            "source_training_rows": len(samples), "source_cv_performance_used": False}


def normalize_native(panel, selected, design, secondary, parser):
    # Reuse the pure source-normalization stage, including its cloned query
    # control. Native queries have their own actual three-artifact inventory.
    train, query, records, bindings = secondary.normalize(panel, {**selected, "query_genomes": []}, parser)
    inventory = {(r["kind"], r["source_id"]): r for r in design["native_query_proteomes"]}
    for index, genome in enumerate(selected["query_genomes"]):
        row = inventory[selected["axis"], genome]
        for key in ["normalized_DNA", "proteins", "GFF"]:
            path = check(row[key])
            bindings[path] = row[key]["sha256"]
        proteins = parser((ROOT / row["proteins"]["path"]).read_bytes())
        if len(proteins) != row["predicted_proteins"]:
            raise ValueError("Native query protein inventory differs")
        for number, (source_id, sequence) in enumerate(proteins):
            identifier = f"{selected['axis'][0]}ng{index:04d}p{number:07d}"
            if identifier in train or identifier in query:
                raise ValueError("Normalized native protein identifier collision")
            query[identifier] = sequence
            records.append({"protein": identifier, "source_protein_id": source_id, "genome": genome,
                            "role": "native_query", "accession": row["accession"], "annotation_id": row["annotation_id"],
                            "length": len(sequence), "sequence_sha256": secondary.digest(sequence.encode())})
    if len(records) != len(train) + len(query) or set(train) & set(query):
        raise ValueError("Incomplete/disjoint normalized source/native protein inventory required")
    return train, query, records, bindings


def execute(policy_path, panel_path, expected_panel_sha, axis, workdir, output, plan_only=False):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if sha(__file__) != policy["helper_sha256"] or policy["workers"] != 1 or policy["timeout_seconds_per_command"] != 1800:
        raise ValueError("Frozen helper/resource policy differs")
    bindings = {check(item): item["sha256"] for item in policy["inputs"].values()}
    bindings[policy_path], bindings[Path(__file__)] = hashlib.sha256(policy_bytes).hexdigest(), sha(__file__)
    tool = {"path": policy["mmseqs_executable"], "sha256": policy["mmseqs_sha256"]}
    bindings[check(tool)] = tool["sha256"]
    design_path = ROOT / policy["inputs"]["native_design"]["path"]
    design = json.loads(design_path.read_text())
    source_policy = json.loads((ROOT / policy["inputs"]["source_feature_policy"]["path"]).read_text())
    if (any(policy[k] != source_policy[k] for k in ["command_templates", "mmseqs_executable", "mmseqs_sha256"])
            or sha(protein_features.__file__) != policy["inputs"]["kernel"]["sha256"]
            or sha(partition.__file__) != policy["inputs"]["partition"]["sha256"]):
        raise ValueError("Pinned source feature semantics differ")
    if not panel_path.exists() or not expected_panel_sha:
        return {"schema_version": 1, "status": "HOLD_missing_frozen_panel", "mode": "final_source_native_query",
                "execution_started": False, "panel_exists": panel_path.exists(),
                "externally_frozen_panel_hash_supplied": bool(expected_panel_sha)}
    bindings[check({"path": str(panel_path), "sha256": expected_panel_sha})] = expected_panel_sha
    panel = json.loads(panel_path.read_text())
    if (panel["helper_sha256"] != policy["inputs"]["panel_helper"]["sha256"]
            or panel["policy_sha256"] != policy["inputs"]["panel_policy"]["sha256"]):
        raise ValueError("Completed panel producer/policy differs")
    for item in panel["inputs"].values():
        bindings[check(item)] = item["sha256"]
    selected = selection(panel, design, axis)
    secondary = load("pinned_native_source_utilities", ROOT / policy["inputs"]["source_feature_helper"]["path"])
    ledger_helper = load("pinned_native_command_ledger", ROOT / policy["inputs"]["ledger_helper"]["path"])
    parser = load("pinned_native_protein_parser", ROOT / policy["inputs"]["protein_parser"]["path"]).proteins
    ledger_helper.guard_destinations(workdir, output, bindings)
    train, query, records, normalized_sources = normalize_native(panel, selected, design, secondary, parser)
    bindings.update(normalized_sources)
    ledger_helper.guard_destinations(workdir, output, bindings)
    contents = {"training.faa": "".join(f">{p}\n{s}\n" for p, s in train.items()).encode(),
                "queries.faa": "".join(f">{p}\n{s}\n" for p, s in query.items()).encode(),
                "protein-map.json": secondary.encoded(records)}
    commands = secondary.command_plan(policy, workdir)
    plan = {**selected, "panel": {"path": str(panel_path), "sha256": expected_panel_sha},
            "native_design": policy["inputs"]["native_design"], "commands": commands,
            "source_bindings": [{"path": str(p), "sha256": h} for p, h in sorted(bindings.items())],
            "normalized_sha256": {n: secondary.digest(data) for n, data in contents.items()},
            "policy_sha256": secondary.digest(policy_bytes), "helper_sha256": sha(__file__),
            "training_proteins": len(train), "query_proteins_including_control": len(query), "execution_started": False,
            "native_outcome_assets_opened": False, "source_cv_performance_used": False}
    if plan_only:
        ledger_helper.verify(bindings)
        return {"schema_version": 1, "status": "planned_only", **plan}
    workdir.mkdir(parents=True, exist_ok=False)
    for filename, content in {**contents, "frozen-plan.json": secondary.encoded(plan)}.items():
        path = workdir / filename
        path.write_bytes(content)
        bindings[path] = sha(path)
    lengths, qlengths = {p: len(s) for p, s in train.items()}, {p: len(s) for p, s in query.items()}
    tmap = {r["protein"]: r["genome"] for r in records if r["role"] == "train"}
    qmap = {r["protein"]: r["genome"] for r in records if r["role"] != "train"}
    ledger, local = [], dict(bindings)
    report = {"schema_version": 1, "status": "failed_native_transfer_axis_features", **selected,
              "panel": plan["panel"], "native_design": plan["native_design"], "helper_sha256": plan["helper_sha256"],
              "policy_sha256": plan["policy_sha256"], "frozen_plan": bound(workdir / "frozen-plan.json"),
              "normalized_input_bindings": [bound(workdir / n) for n in contents], "source_bindings": plan["source_bindings"],
              "commands": ledger, "mapping_control_passed": False, "execution_started": True,
              "execution_context": {"architecture": platform.machine(), "platform": platform.platform()},
              "native_outcome_assets_opened": False, "source_cv_performance_used": False,
              "model_fitted_or_scored": False, "new_original_benchmark_coverage": 0}
    try:
        for index, command in enumerate(commands, 1):
            ledger_helper.run_command(command, index, workdir, 1800, local, ledger)
            if index == 5:
                with (workdir / "clusters.tsv").open() as stream:
                    families = protein_features.read_families(stream, lengths)
                with (workdir / "self.tsv").open() as stream:
                    self_hits, self_diagnostics = protein_features.best_hits(stream, lengths, lengths, families)
                if set(self_hits) != set(train):
                    raise ValueError("Complete source training self-query coverage required")
                vocabulary = protein_features.fit_vocabulary(selected["training_genomes"], tmap, families, self_hits)
                ledger_helper.dump(workdir / "vocabulary.json", vocabulary)
                ledger_helper.dump(workdir / "best-self.json", self_hits)
                ledger_helper.dump(workdir / "pre-query-freeze.json", {"completed_commands": 5, "query_createdb_started": False,
                                   "vocabulary": bound(workdir / "vocabulary.json"), "best_self": bound(workdir / "best-self.json"),
                                   "training_inputs": [bound(workdir / n) for n in ["training.faa", "clusters.tsv", "self.tsv"]]})
                local.update({workdir / n: sha(workdir / n) for n in ["vocabulary.json", "best-self.json", "pre-query-freeze.json"]})
        with (workdir / "query.tsv").open() as stream:
            hits, query_diagnostics = protein_features.best_hits(stream, qlengths, lengths, families)
        names = selected["training_genomes"] + selected["query_genomes"]
        features, vectors = protein_features.project(vocabulary, names + [CONTROL], tmap | qmap, self_hits | hits)
        control_vector = vectors.pop(CONTROL)
        control = secondary.check_control(records, self_hits, hits, control_vector, vectors[selected["mapping_control_source"]])
        ledger_helper.dump(workdir / "mapping-control-check.json", control)
        ledger_helper.dump(workdir / "best-query.json", hits)
        matrix = io.StringIO(newline="")
        writer = csv.writer(matrix, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", *features])
        writer.writerows([n, *vectors[n]] for n in names)
        (workdir / "feature-matrix.tsv").write_text(matrix.getvalue())
        ledger_helper.verify(local)
        report.update({"status": "completed_native_transfer_axis_features" if control["passed"] else "failed_mapping_control",
                       "feature_names": features, "vectors": vectors, "vocabulary": bound(workdir / "vocabulary.json"),
                       "pre_query_freeze": bound(workdir / "pre-query-freeze.json"), "mapping_control_passed": control["passed"],
                       "mapping_control": {"source": selected["mapping_control_source"], "vector": control_vector,
                                           "family_and_vector_equality": control["passed"], "nonzero_vector": control["nonzero_vector"],
                                           "check": bound(workdir / "mapping-control-check.json"), "training_member": False},
                       "training_hit_diagnostics": self_diagnostics, "query_hit_diagnostics": query_diagnostics,
                       "zero_vector_genomes": [n for n in names if not any(vectors[n])],
                       "protein_counts": dict(Counter((tmap | qmap).values())),
                       "ready_for_fitting": control["passed"] and bool(features),
                       "interpretation": "Full-source vocabulary/native projection only; distinct from held-out source CV, not binding, specificity or assembly validation"})
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        report["error"] = str(exc)
    report["artifacts"] = [bound(p) for p in sorted(workdir.iterdir()) if p.is_file()]
    report["command_logs_retained"] = bool(ledger)
    report["resource_scope"] = "Serial commands; child CPU and cumulative direct-child maximum RSS, not simultaneous aggregate memory"
    output.parent.mkdir(parents=True, exist_ok=True)
    ledger_helper.dump(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--expected-panel-sha256")
    parser.add_argument("--axis", choices=["host", "phage"], required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    result = execute(args.policy, args.panel, args.expected_panel_sha256, args.axis, args.workdir, args.output, args.plan_only)
    print(json.dumps(result if args.plan_only else {"status": result["status"], "execution_started": result.get("execution_started", False)}))
    raise SystemExit(0 if result["status"] in {"planned_only", "completed_native_transfer_axis_features"} else 2)
