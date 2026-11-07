"""Exercise the fitting component using a fixed, nonbiological binary fixture.

No genome, published interaction table or target observation can be supplied.
The check establishes software behavior, not host-prediction performance.
"""
import argparse
from dataclasses import replace
import hashlib
import inspect
import json
import math
from pathlib import Path
import platform
import shutil
import sys
import time

from tipcraft import learning


ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "software/tipcraft/experiments/learning/synthetic-check-spec.json"
LOCK = ROOT / "software/tipcraft/experiments/learning/uv.lock"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def close(first, second, tolerance):
    return len(first) == len(second) and all(abs(a - b) <= tolerance for a, b in zip(first, second))


def fixture():
    specification = json.loads(SPEC.read_text())
    features, labels, samples = [], [], []
    for i in range(specification["rows"]):
        b = [(i >> j) & 1 for j in range(5)]
        features.append([b[0], b[1], b[0] ^ b[1], b[2], b[3], b[4]])
        labels.append(b[0] & b[1] & b[2])
        samples.append(learning.Sample(f"synthetic_row_{i}", f"synthetic_host_{i % 4}",
                                       f"synthetic_phage_{i // 4}", f"host_group_{i % 4}",
                                       f"phage_group_{i // 4}", "synthetic_publication"))
    return specification, features, labels, samples


def run_check(workdir):
    if workdir.exists():
        raise ValueError("Use a new check directory")
    workdir.mkdir(parents=True)
    specification, features, labels, samples = fixture()
    settings = learning.FitSettings(**specification["settings"])
    names = specification["features"]
    tolerance = specification["absolute_score_comparison_tolerance"]
    source_hash = digest(Path(learning.__file__).read_bytes())
    start = time.monotonic()
    provenance = {"purpose": "synthetic-engine-check", "real_biological_data": False,
                  "check_spec_sha256": digest(SPEC.read_bytes()), "runtime_lock_sha256": digest(LOCK.read_bytes())}
    paired = workdir / "paired"
    renamed = workdir / "renamed"
    host = workdir / "host-only"
    learning.fit_bundle(features, labels, samples, names, settings, paired, provenance=provenance)
    renamed_samples = [replace(s, row_id=f"anonymous_row_{i}", host_id=f"anonymous_h_{i % 4}",
                               phage_id=f"anonymous_p_{i // 4}", host_group=f"anonymous_hg_{i % 4}",
                               phage_group=f"anonymous_pg_{i // 4}") for i, s in enumerate(samples)]
    learning.fit_bundle(features, labels, renamed_samples, names, settings, renamed, provenance=provenance)
    learning.fit_bundle(features, labels, samples, names, replace(settings, feature_scope="host_only"),
                        host, provenance=provenance)
    bundles = (paired, renamed, host)
    frozen_hashes = {folder.name: digest((folder / "manifest.json").read_bytes()) for folder in bundles}
    boundary = workdir / "bundles-before-scoring.json"
    with boundary.open("x") as stream:
        json.dump({"scoring_started": False, "manifest_sha256": frozen_hashes}, stream, indent=2)
        stream.write("\n")
    original_scores = learning.score_bundle(features, names, paired, expected_manifest_sha256=frozen_hashes[paired.name])
    renamed_scores = learning.score_bundle(features, names, renamed, expected_manifest_sha256=frozen_hashes[renamed.name])
    host_scores = learning.score_bundle(features, names, host, expected_manifest_sha256=frozen_hashes[host.name])
    perturbed = [row[:3] + [1 - value for value in row[3:]] for row in features]
    perturbed_scores = learning.score_bundle(perturbed, names, host, expected_manifest_sha256=frozen_hashes[host.name])
    checks = {
        "scores_finite_and_nonconstant": all(math.isfinite(value) for value in original_scores)
            and max(original_scores) - min(original_scores) > tolerance,
        "saved_bundle_reload_scores_match": close(original_scores, learning.score_bundle(features, names, paired, expected_manifest_sha256=frozen_hashes[paired.name]), tolerance),
        "bijective_metadata_renaming_invariant": close(original_scores, renamed_scores, tolerance),
        "host_only_ignores_phage_columns": close(host_scores, perturbed_scores, tolerance),
        "fit_api_has_no_evaluation_labels_or_early_stopping_input": not any(
            name in inspect.signature(learning.fit_bundle).parameters for name in ("eval_set", "y_test", "y_validation", "early_stopping_rounds")),
    }
    unrelated = learning.Sample("heldout", "new_host", "new_phage", "new_hg", "new_pg", "new_publication")
    learning.validate_partition(samples, [unrelated], regime="both_unseen", publication_holdout=True)
    for key, heldout in [("shared_host_group_rejected", replace(unrelated, host_group=samples[0].host_group)),
                         ("shared_phage_group_rejected", replace(unrelated, phage_group=samples[0].phage_group))]:
        try:
            learning.validate_partition(samples, [heldout], regime="both_unseen")
        except ValueError:
            checks[key] = True
        else:
            checks[key] = False
    try:
        learning.score_bundle(features, names[::-1], paired, expected_manifest_sha256=frozen_hashes[paired.name])
    except ValueError:
        checks["reordered_feature_schema_rejected"] = True
    else:
        checks["reordered_feature_schema_rejected"] = False
    for filename in ("manifest.json", "model.cbm", "fit-spec.json", "selected-features.json"):
        copied = workdir / ("tampered-" + filename.replace(".", "-"))
        shutil.copytree(paired, copied)
        with (copied / filename).open("ab") as stream:
            stream.write(b"\n")
        try:
            learning.score_bundle(features, names, copied, expected_manifest_sha256=frozen_hashes[paired.name])
        except ValueError:
            checks[filename + "_tampering_rejected"] = True
        else:
            checks[filename + "_tampering_rejected"] = False
    if source_hash != digest(Path(learning.__file__).read_bytes()):
        raise ValueError("Fitting component changed during the check")
    if any(digest((folder / "manifest.json").read_bytes()) != frozen_hashes[folder.name] for folder in bundles):
        raise ValueError("A fitted bundle changed after its pre-scoring hash was frozen")
    result = {"schema_version": 1, "purpose": specification["purpose"],
              "status": "passed" if all(checks.values()) else "failed", "checks": checks,
              "spec_sha256": digest(SPEC.read_bytes()), "runtime_lock_sha256": digest(LOCK.read_bytes()),
              "component_sha256": source_hash, "runner_sha256": digest(Path(__file__).read_bytes()),
              "fixture_sha256": digest(json.dumps({"X": features, "y": labels, "features": names}, separators=(",", ":")).encode()),
              "samples": len(samples), "input_features": len(names), "settings": specification["settings"],
              "absolute_score_comparison_tolerance": tolerance,
              "platform": platform.platform(), "architecture": platform.machine(), "python": sys.version.split()[0],
              "wall_seconds": time.monotonic() - start,
              "pre_scoring_boundary": {"path": str(boundary), "sha256": digest(boundary.read_bytes())},
              "bundle_manifests": [{"path": str(folder / "manifest.json"), "sha256": frozen_hashes[folder.name]}
                                   for folder in bundles],
              "biological_training_executed": False, "biological_validation_executed": False,
              "target_infection_scores_generated": False, "new_method_covered_cases": 0,
              "interpretation": "Actual fitting-engine integrity checks on synthetic binary data; no biological accuracy or full-run readiness established"}
    (workdir / "check-report.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.resolve().is_relative_to(args.workdir.resolve()):
        parser.error("Use a new report path outside the check directory")
    result = run_check(args.workdir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "checks": result["checks"]}))
    raise SystemExit(0 if result["status"] == "passed" else 1)
