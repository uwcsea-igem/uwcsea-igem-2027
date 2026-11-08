"""Refit admitted source observations and score six fixed native pairs.

This stage never parses Ando labels or historical WIsH scores. Native label
evaluation is a separate operation after all three score matrices are saved.
"""
import argparse
from dataclasses import asdict
import importlib.util
from importlib.metadata import version
import json
from pathlib import Path
import sys

from tipcraft import learning

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/native-transfer-fitting-policy-v3.json"
SPEC = importlib.util.spec_from_file_location("pinned_source_fitting", ROOT / "scripts/fit-secondary-block.py")
FIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIT)
HOSTS = ["BL21", "BW25113", "MG1655"]
PHAGES = ["T7", "T3"]


def assemble_native(panel, host, phage, design, host_groups, phage_groups):
    if (panel["status"] != "frozen_secondary_reference_panel" or panel["membership_generated"] is not True
            or panel["annotation_join_complete"] is not True or panel["ready_for_feature_preparation"] is not True):
        raise ValueError("Complete source panel required")
    if design["native_candidates"] != {"hosts": HOSTS, "phages": PHAGES}:
        raise ValueError("Exact original native candidates required")
    source = learning._samples([learning.Sample(**s) for s in panel["samples"]])
    host_map = {}
    for name in HOSTS:
        groups = [gid for gid, group in host_groups["partition"]["source_groups"].items()
                  if "original:" + name in group["members"]]
        if len(groups) != 1:
            raise ValueError("Exactly one measured original host group required")
        host_map[name] = groups[0]
    phage_map = phage_groups["grouping"]["container_groups"]
    target = [learning.Sample("native:" + p + ":" + h, h, p, host_map[h], phage_map["reference:" + p],
                              "Ando et al. (2015); original native plaque observations") for p in PHAGES for h in HOSTS]
    validation = learning.validate_partition(source, target, regime="both_unseen", publication_holdout=True)
    names, vectors = [], {}
    for axis, report, queries in [("host", host, HOSTS), ("phage", phage, PHAGES)]:
        train = sorted({getattr(s, axis + "_id") for s in source})
        if (report["status"] != "completed_native_transfer_axis_features" or report["axis"] != axis
                or report["training_genomes"] != train or report["query_genomes"] != queries
                or report["mapping_control_passed"] is not True):
            raise ValueError("Matching final-source/native feature roster and control required")
        columns = learning._feature_names([axis + ":" + f for f in report["feature_names"]])
        if set(report["vectors"]) != set(train + queries):
            raise ValueError("Exact source/native vector inventory required")
        all_names = train + queries
        matrix = learning._matrix([report["vectors"][n] for n in all_names], len(columns))
        vectors[axis] = dict(zip(all_names, matrix))
        names.extend(columns)
    rows = lambda samples: [vectors["host"][s.host_id] + vectors["phage"][s.phage_id] for s in samples]
    labels = learning._labels([panel["labels"][s.row_id] for s in source], len(source), both_classes=True)
    return {"source": source, "target": target, "X_train": rows(source), "X_query": rows(target),
            "y_train": labels, "feature_names": names, "partition": validation}


def fit_and_score_native(data, settings, workdir, provenance):
    if set(settings) != set(FIT.SCOPES):
        raise ValueError("All three fixed scopes required")
    workdir.mkdir(parents=True, exist_ok=False)
    FIT.write_new(workdir / "plan-before-fitting.json", {
        "mode": "all_admitted_source_native_transfer", "settings": {s: asdict(settings[s]) for s in FIT.SCOPES},
        "source_row_ids": [s.row_id for s in data["source"]],
        "native_pairs": [{"phage": s.phage_id, "host": s.host_id} for s in data["target"]],
        "feature_names": data["feature_names"], "provenance": provenance,
        "native_labels_used": False, "model_selection_using_source_cv_performance": False})
    folders = {scope: workdir / scope for scope in FIT.SCOPES}
    for scope, folder in folders.items():
        learning.fit_bundle(data["X_train"], data["y_train"], data["source"], data["feature_names"],
                            settings[scope], folder, provenance)
    frozen = {s: FIT.sha(p / "manifest.json") for s, p in folders.items()}
    boundary = workdir / "all-scopes-before-native-scoring.json"
    FIT.write_new(boundary, {"scoring_started": False, "manifest_sha256": frozen})
    methods = {}
    for scope, folder in folders.items():
        values = learning.score_bundle(data["X_query"], data["feature_names"], folder,
                                       expected_manifest_sha256=frozen[scope])
        scores = {p: {} for p in PHAGES}
        for row, score in zip(data["target"], values):
            scores[row.phage_id][row.host_id] = score
        if set(scores) != set(PHAGES) or any(set(v) != set(HOSTS) for v in scores.values()):
            raise ValueError("Complete six-native-pair scores required")
        methods[scope] = {"model_id": frozen[scope], "scores": scores}
    if any(FIT.sha(p / "manifest.json") != frozen[s] for s, p in folders.items()):
        raise ValueError("Frozen native inference bundles changed")
    return {"schema_version": 1, "status": "completed_native_transfer_scores", "methods": methods,
            "bundles": {s: FIT.binding(p / "manifest.json") for s, p in folders.items()},
            "pre_scoring_boundary": FIT.binding(boundary), "plan_before_fitting": FIT.binding(workdir / "plan-before-fitting.json"),
            "source_training_rows": len(data["source"]), "native_scored_pairs_per_scope": 6,
            "exchange_construct_pairs_abstained": 4, "native_labels_used_for_fit_or_selection": False,
            "historical_WIsH_scores_used_for_fit_or_selection": False,
            "score_type": "RawFormulaVal; uncalibrated", "partition": data["partition"],
            "label_evaluation_executed": False, "assembly_or_binding_evidence": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    for name in ["panel", "source-evaluation", "host-features", "phage-features"]:
        parser.add_argument("--" + name, type=Path, required=True)
        parser.add_argument("--expected-" + name + "-sha256", required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workdir.exists() or args.output.exists() or args.output.resolve().is_relative_to(args.workdir.resolve()):
        parser.error("Use fresh destinations with public output outside raw workdir")
    policy_bytes = args.policy.read_bytes()
    policy = json.loads(policy_bytes)
    if FIT.sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Native fitting implementation differs")
    FIT.verify_bindings(policy["inputs"])
    if FIT.sha(learning.__file__) != policy["inputs"]["learning"]["sha256"] or FIT.sha(FIT.__file__) != policy["inputs"]["source_fitter"]["sha256"]:
        raise ValueError("Imported native fitting dependency differs")
    reports, inputs = {}, {}
    for name in ["panel", "source_evaluation", "host_features", "phage_features"]:
        path = getattr(args, name)
        FIT.checked({"path": str(path), "sha256": getattr(args, "expected_" + name + "_sha256")})
        inputs[name] = FIT.binding(path)
        reports[name] = json.loads(path.read_text())
        FIT.verify_bindings(reports[name])
    panel = reports["panel"]
    if (panel["helper_sha256"] != policy["inputs"]["panel_helper"]["sha256"]
            or panel["policy_sha256"] != policy["inputs"]["panel_policy"]["sha256"]):
        raise ValueError("Source panel producer differs")
    evaluation = reports["source_evaluation"]
    if (evaluation["status"] != "completed_secondary_matched_evaluation"
            or evaluation["eligible_scoring_coverage"] != 1.0
            or evaluation["source_counts"] != panel["counts"]
            or evaluation["provenance"]["helper_sha256"] != policy["inputs"]["source_evaluator"]["sha256"]
            or evaluation["provenance"]["policy_sha256"] != policy["inputs"]["source_evaluation_policy"]["sha256"]):
        raise ValueError("Complete matched source evaluation required, irrespective of performance")
    source_manifest = json.loads(FIT.checked(evaluation["provenance"]["score_manifest"]).read_text())
    if ((ROOT / source_manifest["panel"]["path"]).resolve() != args.panel.resolve()
            or source_manifest["panel"]["sha256"] != inputs["panel"]["sha256"]):
        raise ValueError("Source evaluation belongs to a different panel")
    native_design = json.loads(FIT.checked(policy["inputs"]["native_design"]).read_text())
    for axis in ["host", "phage"]:
        report = reports[axis + "_features"]
        if ((ROOT / report["panel"]["path"]).resolve() != args.panel.resolve()
                or report["panel"]["sha256"] != inputs["panel"]["sha256"]
                or report["native_design"]["sha256"] != policy["inputs"]["native_design"]["sha256"]
                or report["helper_sha256"] != policy["inputs"]["native_feature_helper"]["sha256"]
                or report["policy_sha256"] != policy["inputs"]["native_feature_policy"]["sha256"]
                or len(report["commands"]) != 9
                or any(type(c["exit_code"]) is not int or c["exit_code"] != 0 for c in report["commands"])):
            raise ValueError("Matching native feature execution required")
    host_groups = json.loads(FIT.checked(panel["inputs"]["host_groups"]).read_text())
    phage_groups = json.loads(FIT.checked(panel["inputs"]["phage_groups"]).read_text())
    data = assemble_native(panel, reports["host_features"], reports["phage_features"], native_design, host_groups, phage_groups)
    design = json.loads(FIT.checked(policy["inputs"]["design"]).read_text())
    settings = FIT.fixed_settings(design, data["feature_names"])
    if sys.version.split()[0] != policy["runtime"]["python"]:
        raise ValueError("Declared locked Python runtime required")
    for name, expected in policy["runtime"]["packages"].items():
        if version(name) != expected:
            raise ValueError("Runtime package differs: " + name)
    provenance = {"inputs": inputs, "helper": FIT.binding(__file__), "policy": FIT.binding(args.policy),
                  "mode": "final_source_native_transfer"}
    try:
        result = fit_and_score_native(data, settings, args.workdir, provenance)
    except Exception as error:
        if args.workdir.exists():
            FIT.write_new(args.workdir / "failed-attempt.json", {"status": "HOLD_native_fit_or_scoring_failure",
                "exception_type": type(error).__name__, "message": str(error), "provenance": provenance,
                "automatic_retry": False, "partial_artifacts": [FIT.binding(p) for p in sorted(args.workdir.rglob("*")) if p.is_file()]})
        raise
    FIT.verify_bindings(inputs)
    FIT.verify_bindings(policy["inputs"])
    if args.policy.read_bytes() != policy_bytes or FIT.sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Executed native fitting source/policy changed")
    result.update({"inputs": inputs, "helper_sha256": policy["helper_sha256"], "policy": FIT.binding(args.policy),
                   "source_evaluation_screen_retained": evaluation["screen"]})
    FIT.write_new(args.output, result)
    print(json.dumps({k: result[k] for k in ["status", "source_training_rows", "native_scored_pairs_per_scope"]}))


if __name__ == "__main__":
    main()
