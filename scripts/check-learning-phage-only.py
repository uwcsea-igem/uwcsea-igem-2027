"""Check the phage-only ablation with the existing fixed synthetic fixture.

No biological input is accepted. Historical source and results remain unchanged.
"""
import argparse
from dataclasses import replace
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import time

from tipcraft import learning

ROOT = Path(__file__).resolve().parents[1]
LEGACY_REPORT = ROOT / "data/processed/learning-engine-check.json"
LEGACY_SOURCE = ROOT / "data/raw/learning-phage-only-check/source-before-change-v1/learning.py"
LEGACY_REPORT_SHA = "d1b77fab398c93e24f59c33daef2404cdb72f29bfecd5f62b89bc9cb5ef5585a"
LEGACY_SOURCE_SHA = "698024df5a680e24b87ca8d35ebca7c97c0106eefaf9a6c629b78ea895f093f3"
FIXTURE_HELPER = ROOT / "scripts/check-learning-fit.py"
FIXTURE_HELPER_SHA = "400e3c09e7202b692a5c4251150b443eedb212e5de6aa2ef2fdfb055ab184e34"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path):
    return {"path": str(path.resolve().relative_to(ROOT)), "sha256": sha(path)}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def write_new(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def run_check(workdir):
    if workdir.exists():
        raise ValueError("Use a new check directory")
    for path, expected in [(LEGACY_REPORT, LEGACY_REPORT_SHA), (LEGACY_SOURCE, LEGACY_SOURCE_SHA),
                           (FIXTURE_HELPER, FIXTURE_HELPER_SHA)]:
        if sha(path) != expected:
            raise ValueError(f"Historical source/result checksum differs: {path}")
    fixture_helper = load_module("synthetic_fixture_phage_only", FIXTURE_HELPER)
    legacy = load_module("legacy_learning_phage_only", LEGACY_SOURCE)
    spec, X, y, samples = fixture_helper.fixture()
    names = spec["features"]
    settings = replace(learning.FitSettings(**spec["settings"]), feature_scope="phage_only")
    tolerance = spec["absolute_score_comparison_tolerance"]
    perturbed = [[1 - v for v in row[:3]] + row[3:] for row in X]
    old_result = json.loads(LEGACY_REPORT.read_text())
    old_bundles = old_result["bundle_manifests"]
    for item in old_bundles:
        if sha(ROOT / item["path"]) != item["sha256"]:
            raise ValueError("Historical bundle manifest differs")
    inputs = {"component": binding(Path(learning.__file__)), "helper": binding(Path(__file__)),
              "fixture_helper": binding(FIXTURE_HELPER), "fixture_spec": binding(fixture_helper.SPEC),
              "runtime_lock": binding(fixture_helper.LOCK), "legacy_source": binding(LEGACY_SOURCE),
              "legacy_report": binding(LEGACY_REPORT)}
    workdir.mkdir(parents=True, exist_ok=False)
    plan = {"inputs": inputs, "settings": vars(settings), "rows": len(X), "features": names,
            "absolute_score_comparison_tolerance": tolerance,
            "intervention": "Complement every host column, consistently within each synthetic host; phage columns unchanged",
            "expected": ["phage-only selected columns only", "finite nonconstant margins", "save/load equality",
                         "host perturbation invariant at fit and scoring", "legacy paired/host-only bundle score equality"],
            "biological_data_accepted": False, "parameter_search": False}
    write_new(workdir / "plan-before-fitting.json", plan)
    start = time.monotonic()
    bundles = [workdir / "phage-only", workdir / "host-perturbed-fit"]
    provenance = {"purpose": "synthetic-phage-only-check", "biological_data": False,
                  "plan_sha256": sha(workdir / "plan-before-fitting.json")}
    manifests = []
    for matrix, bundle in zip([X, perturbed], bundles):
        manifests.append(learning.fit_bundle(matrix, y, samples, names, settings, bundle, provenance))
    frozen = {folder.name: sha(folder / "manifest.json") for folder in bundles}
    boundary = workdir / "bundles-before-scoring.json"
    write_new(boundary, {"scoring_started": False, "manifest_sha256": frozen,
                         "historical_bundle_manifests": old_bundles})

    def score(matrix, folder):
        return learning.score_bundle(matrix, names, folder, expected_manifest_sha256=frozen[folder.name])

    scores = score(X, bundles[0])
    close = lambda a, b: fixture_helper.close(a, b, tolerance)
    checks = {"selected_columns_phage_only": all(
        m["selected_global_indices"] == [3, 4, 5] and all(n.startswith("phage:") for n in m["selected_feature_names"])
        for m in manifests),
        "finite_nonconstant_margins": all(math.isfinite(v) for v in scores) and max(scores) - min(scores) > tolerance,
        "save_load_scores_equal": close(scores, score(X, bundles[0])),
        "host_perturbation_at_scoring_invariant": close(scores, score(perturbed, bundles[0])),
        "host_perturbation_at_fit_invariant": close(scores, score(X, bundles[1]))}
    for item in old_bundles:
        folder = (ROOT / item["path"]).parent
        kwargs = {"expected_manifest_sha256": item["sha256"]}
        checks["legacy_" + folder.name + "_load_scores_equal"] = close(
            legacy.score_bundle(X, names, folder, **kwargs), learning.score_bundle(X, names, folder, **kwargs))
    checks["input_sources_unchanged"] = all(sha(ROOT / b["path"]) == b["sha256"] for b in inputs.values())
    checks["frozen_bundle_manifests_unchanged"] = all(sha(p / "manifest.json") == frozen[p.name] for p in bundles)
    checks["historical_bundle_manifests_unchanged"] = all(sha(ROOT / b["path"]) == b["sha256"] for b in old_bundles)
    result = {"schema_version": 1, "status": "passed" if all(checks.values()) else "failed",
              "inputs": inputs, "checks": checks, "settings": vars(settings),
              "plan_before_fitting": binding(workdir / "plan-before-fitting.json"),
              "pre_scoring_boundary": binding(boundary),
              "bundle_manifests": [binding(p / "manifest.json") for p in bundles],
              "versions": manifests[0]["versions"], "architecture": platform.machine(),
              "platform": platform.platform(), "wall_seconds": time.monotonic() - start,
              "absolute_score_comparison_tolerance": tolerance,
              "legacy_source_git_revision": "4666b4d60e548847a38d7fd589cebc84cf2c9d83",
              "legacy_artifacts_rewritten": False, "biological_training_executed": False,
              "target_infection_scores_generated": False, "new_method_covered_cases": 0,
              "interpretation": "Synthetic software checks only; phage-only raw margins are not calibrated infection probabilities"}
    write_new(workdir / "check-report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.resolve().is_relative_to(args.workdir.resolve()):
        parser.error("Use a new report outside the check directory")
    result = run_check(args.workdir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_new(args.output, result)
    print(json.dumps({"status": result["status"], "checks": result["checks"]}))
    raise SystemExit(0 if result["status"] == "passed" else 1)
