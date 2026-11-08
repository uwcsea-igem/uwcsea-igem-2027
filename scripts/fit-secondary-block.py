"""Fit three frozen scopes on one authenticated outer training block.

No held-out labels are passed to model fitting or returned with predictions.
Whole panel JSON may contain those labels; this is input-use isolation, not a
claim that a process cannot access its supplied evaluation annotations.
"""
import argparse
from dataclasses import asdict
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys

from tipcraft import learning

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/secondary-fitting-policy.json"
SCOPES = ["paired", "host_only", "phage_only"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binding(path):
    path = Path(path).resolve()
    name = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
    return {"path": str(name), "sha256": sha(path)}


def checked(record):
    path = ROOT / record["path"]
    if sha(path) != record["sha256"]:
        raise ValueError("Input binding differs: " + str(path))
    return path


def verify_bindings(value):
    if isinstance(value, dict):
        if "path" in value and "sha256" in value:
            checked(value)
        for item in value.values():
            verify_bindings(item)
    elif isinstance(value, list):
        for item in value:
            verify_bindings(item)


def verify_implementation(policy):
    expected = policy["inputs"]["learning"]
    if (Path(learning.__file__).resolve() != (ROOT / expected["path"]).resolve()
            or sha(learning.__file__) != expected["sha256"]):
        raise ValueError("Imported fitting implementation differs from pinned source")


def write_new(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def assemble_block(panel, host, phage, block_id):
    """Resolve exact frozen memberships; select labels for training rows only."""
    if (panel["status"] != "frozen_secondary_reference_panel"
            or panel["ready_for_feature_preparation"] is not True
            or panel["annotation_join_complete"] is not True
            or panel["membership_generated"] is not True):
        raise ValueError("Complete authenticated source panel required")
    samples = learning._samples([learning.Sample(**row) for row in panel["samples"]])
    by_id = {row.row_id: row for row in samples}
    found = [block for block in panel["blocks"] if block["fold_id"] == block_id]
    if len(found) != 1:
        raise ValueError("Exactly one declared outer block required")
    block = found[0]
    ids = [block[key] for key in ["train_ids", "test_ids", "embargo_ids"]]
    if (any(not subset or len(set(subset)) != len(subset) for subset in ids)
            or set.union(*(set(s) for s in ids)) != set(by_id)
            or any(set(ids[i]) & set(ids[j]) for i in range(3) for j in range(i))):
        raise ValueError("Block roles must completely and uniquely partition source rows")
    train, test = [[by_id[i] for i in subset] for subset in ids[:2]]
    validation = learning.validate_partition(train, test, regime="both_unseen")
    matrices, names = {}, []
    for axis, report in [("host", host), ("phage", phage)]:
        if (report["status"] != "completed_secondary_axis_features" or report["axis"] != axis
                or report["block_id"] != block_id or report["mapping_control_passed"] is not True):
            raise ValueError("Completed matching axis/block features with control required")
        attribute = axis + "_id"
        train_genomes = sorted({getattr(s, attribute) for s in train})
        query_genomes = sorted({getattr(s, attribute) for s in test})
        if report["training_genomes"] != train_genomes or report["query_genomes"] != query_genomes:
            raise ValueError("Feature training/query genomes differ from block membership")
        features = learning._feature_names([axis + ":" + f for f in report["feature_names"]])
        vectors = report["vectors"]
        if set(vectors) != set(train_genomes + query_genomes):
            raise ValueError("Exact source-only feature vector roster required")
        genomes = sorted(vectors)
        valid = learning._matrix([vectors[g] for g in genomes], len(features))
        matrices[axis] = dict(zip(genomes, valid))
        names.extend(features)
    matrix = lambda rows: [matrices["host"][s.host_id] + matrices["phage"][s.phage_id] for s in rows]
    X_train, X_test = matrix(train), matrix(test)
    y_train = learning._labels([panel["labels"][s.row_id] for s in train], len(train), both_classes=True)
    return {"block_id": block_id, "train": train, "test": test, "X_train": X_train,
            "X_test": X_test, "y_train": y_train, "feature_names": names, "partition": validation}


def fixed_settings(design, names):
    config = design["fit"]
    if (config["scopes"] != SCOPES or config["n_selected_features_rule"] != "min(100,P) after scope filtering; P=0 is HOLD"
            or config["rfe_step_rule"] != "max(1,int((P-K)/10))"
            or any(config[key] is not False for key in ["hyperparameter_search", "early_stopping", "calibration"])):
        raise ValueError("Unsupported or changed prospective fitting design")
    result = {}
    for scope in SCOPES:
        p = len(learning._eligible_indices(names, scope))
        if not p:
            raise ValueError("Empty feature scope; HOLD without fallback")
        k = min(100, p)
        result[scope] = learning.FitSettings(**{name: config[name] for name in
            ["rfe_iterations", "final_iterations", "depth", "learning_rate", "l2_leaf_reg", "seed", "threads", "weight_method"]},
            n_selected_features=k, rfe_step=max(1, int((p-k)/10)), feature_scope=scope)
    return result


def fit_and_score(data, settings, workdir, provenance):
    """Save every scope before scoring any; no hyperparameter selection stage."""
    if set(settings) != set(SCOPES):
        raise ValueError("All three fixed scopes required")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=False)
    plan = {"block_id": data["block_id"], "settings": {s: asdict(settings[s]) for s in SCOPES},
            "train_ids": [s.row_id for s in data["train"]], "test_ids": [s.row_id for s in data["test"]],
            "feature_names": data["feature_names"], "provenance": provenance,
            "test_labels_used_for_fit_or_selection": False, "parameter_search": False}
    write_new(workdir / "plan-before-fitting.json", plan)
    folders = {scope: workdir / scope for scope in SCOPES}
    for scope, folder in folders.items():
        learning.fit_bundle(data["X_train"], data["y_train"], data["train"], data["feature_names"],
                            settings[scope], folder, provenance)
    frozen = {scope: sha(folder / "manifest.json") for scope, folder in folders.items()}
    boundary = workdir / "all-scopes-before-scoring.json"
    write_new(boundary, {"scoring_started": False, "manifest_sha256": frozen})
    methods = {}
    for scope, folder in folders.items():
        values = learning.score_bundle(data["X_test"], data["feature_names"], folder,
                                       expected_manifest_sha256=frozen[scope])
        methods[scope] = {"model_id": frozen[scope], "scores": dict(zip([s.row_id for s in data["test"]], values))}
    if any(sha(folders[s] / "manifest.json") != frozen[s] for s in SCOPES):
        raise ValueError("Frozen bundle manifest changed during scoring")
    return {"status": "completed_secondary_block_scores", "block_id": data["block_id"],
            "methods": methods, "bundles": {s: binding(folders[s] / "manifest.json") for s in SCOPES},
            "pre_scoring_boundary": binding(boundary), "plan_before_fitting": binding(workdir / "plan-before-fitting.json"),
            "partition": data["partition"], "training_rows": len(data["train"]), "test_rows": len(data["test"]),
            "test_labels_used_for_fit_or_selection": False, "test_labels_returned_with_scores": False,
            "score_type": "RawFormulaVal; uncalibrated", "original_Ando_cases_scored": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    for name in ["panel", "host-features", "phage-features"]:
        parser.add_argument("--" + name, type=Path, required=True)
        parser.add_argument("--expected-" + name + "-sha256", required=True)
    parser.add_argument("--block-id", required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workdir.exists() or args.output.exists() or args.output.resolve().is_relative_to(args.workdir.resolve()):
        parser.error("Use fresh destinations, with public output outside raw workdir")
    policy = json.loads(args.policy.read_text())
    if sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Fitting helper differs from policy")
    verify_bindings(policy["inputs"])
    verify_implementation(policy)
    inputs = {}
    reports = {}
    for name in ["panel", "host_features", "phage_features"]:
        path = getattr(args, name)
        bound = {"path": str(path.resolve()), "sha256": getattr(args, "expected_" + name + "_sha256")}
        checked(bound)
        reports[name] = json.loads(path.read_text())
        verify_bindings(reports[name])
        inputs[name] = binding(path)
    panel = reports["panel"]
    if (panel["helper_sha256"] != policy["inputs"]["panel_helper"]["sha256"]
            or panel["policy_sha256"] != policy["inputs"]["panel_policy"]["sha256"]):
        raise ValueError("Panel producer differs from prospective policy")
    for axis in ["host", "phage"]:
        report = reports[axis + "_features"]
        if ((ROOT / report["panel"]["path"]).resolve() != args.panel.resolve()
                or report["panel"]["sha256"] != inputs["panel"]["sha256"]
                or report["helper_sha256"] != policy["inputs"]["feature_helper"]["sha256"]
                or report["policy_sha256"] != policy["inputs"]["feature_policy"]["sha256"]
                or not report["commands"] or any(type(c["exit_code"]) is not int or c["exit_code"] != 0 for c in report["commands"])):
            raise ValueError("Feature execution provenance differs")
    data = assemble_block(panel, reports["host_features"], reports["phage_features"], args.block_id)
    design = json.loads(checked(policy["inputs"]["design"]).read_text())
    settings = fixed_settings(design, data["feature_names"])
    provenance = {"inputs": inputs, "helper": binding(__file__), "policy": binding(args.policy)}
    if sys.version.split()[0] != policy["runtime"]["python"]:
        raise ValueError("Use the declared locked Python runtime")
    for package, expected in policy["runtime"]["packages"].items():
        if version(package) != expected:
            raise ValueError("Runtime package differs: " + package)
    try:
        result = fit_and_score(data, settings, args.workdir, provenance)
    except Exception as error:
        if args.workdir.exists():
            write_new(args.workdir / "failed-attempt.json", {
                "status": "HOLD_fit_or_scoring_failure", "block_id": args.block_id,
                "exception_type": type(error).__name__, "message": str(error), "inputs": inputs,
                "provenance": provenance, "automatic_retry": False,
                "partial_artifacts": [binding(p) for p in sorted(args.workdir.rglob("*")) if p.is_file()]})
        raise
    verify_bindings(policy["inputs"])
    verify_implementation(policy)
    verify_bindings(inputs)
    for report in reports.values():
        verify_bindings(report)
    if sha(__file__) != policy["helper_sha256"] or sha(args.policy) != provenance["policy"]["sha256"]:
        raise ValueError("Executed helper/policy changed")
    result.update({"schema_version": 1, "inputs": inputs, "helper_sha256": sha(__file__), "policy": binding(args.policy)})
    write_new(args.output, result)
    print(json.dumps({key: result[key] for key in ["status", "block_id", "training_rows", "test_rows"]}))


if __name__ == "__main__":
    main()
