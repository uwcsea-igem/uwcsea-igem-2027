"""Evaluate all nine frozen source blocks with matched scope coverage."""
import argparse
import hashlib
import json
from pathlib import Path

from tipcraft import blocked_ranking, learning, partition

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/secondary-evaluation-policy.json"
SCOPES = {"paired", "host_only", "phage_only"}
BLOCKS = {f"host{i}-phage{j}" for i in range(3) for j in range(3)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked(binding):
    path = ROOT / binding["path"]
    if sha(path) != binding["sha256"]:
        raise ValueError("Bound artifact differs: " + str(path))
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


def verify_bundle(binding, scope, report, policy):
    path = checked(binding)
    bundle = json.loads(path.read_text())
    expected = {"fit-spec.json", "fit-spec.sha256", "model.cbm", "provenance.json", "selected-features.json"}
    if bundle["status"] != "complete" or bundle["feature_scope"] != scope or set(bundle["artifacts"]) != expected:
        raise ValueError("Complete scope-specific bundle required")
    for filename, record in bundle["artifacts"].items():
        checked({"path": str(path.parent / filename), "sha256": record["sha256"]})
    spec = json.loads((path.parent / "fit-spec.json").read_text())
    provenance = json.loads((path.parent / "provenance.json").read_text())
    if (spec["implementation_sha256"] != policy["inputs"]["learning"]["sha256"]
            or spec["training_rows"] != report["training_rows"]
            or spec["settings"]["feature_scope"] != scope
            or provenance["inputs"] != report["inputs"]
            or provenance["helper"]["sha256"] != policy["inputs"]["fitting_helper"]["sha256"]
            or provenance["policy"]["sha256"] != policy["inputs"]["fitting_policy"]["sha256"]):
        raise ValueError("Bundle fitting provenance differs from score receipt")


def summarize(panel, reports):
    """No incomplete-run metric or favorable-common-subset fallback."""
    if panel["status"] != "frozen_secondary_reference_panel" or panel["membership_generated"] is not True:
        raise ValueError("Frozen complete panel required")
    samples = learning._samples([learning.Sample(**s) for s in panel["samples"]])
    blocks = partition.build_blocked_folds(samples, panel["host_group_folds"], panel["phage_group_folds"])
    if blocks != panel["blocks"] or {b["fold_id"] for b in blocks} != BLOCKS:
        raise ValueError("Exact nine frozen both-unseen blocks required")
    by_block = {}
    for report in reports:
        name = report["block_id"]
        if (name in by_block or name not in BLOCKS or report["status"] != "completed_secondary_block_scores"
                or set(report["methods"]) != SCOPES or set(report["bundles"]) != SCOPES
                or report["test_labels_used_for_fit_or_selection"] is not False
                or report["test_labels_returned_with_scores"] is not False):
            raise ValueError("Unique completed block scores with all three scopes required")
        by_block[name] = report
    if set(by_block) != BLOCKS:
        raise ValueError("All nine block reports required; no partial-run evaluation")
    methods = {scope: {} for scope in sorted(SCOPES)}
    rows = {s.row_id: s for s in samples}
    for block in blocks:
        name = block["fold_id"]
        report = by_block[name]
        if report["training_rows"] != len(block["train_ids"]) or report["test_rows"] != len(block["test_ids"]):
            raise ValueError("Score receipt row counts differ from frozen block")
        for scope in SCOPES:
            value = report["methods"][scope]
            if value["model_id"] != report["bundles"][scope]["sha256"]:
                raise ValueError("Scores must identify the frozen scope bundle")
            if set(value["scores"]) != set(block["test_ids"]):
                raise ValueError("Each scope must cover exactly the complete test block")
            methods[scope][name] = value
        # Same-input controls are evaluated within one fitted model, never across folds.
        for scope, attribute in [("host_only", "host_id"), ("phage_only", "phage_id")]:
            seen = {}
            for row_id, score in report["methods"][scope]["scores"].items():
                genome = getattr(rows[row_id], attribute)
                if genome in seen and seen[genome] != score:
                    raise ValueError("Scope attribution control violated: " + scope)
                seen[genome] = score
    labels = {s.row_id: panel["labels"][s.row_id] for s in samples}
    ranking = blocked_ranking.evaluate_blocked_rankings(samples, labels,
        {b["fold_id"]: b["test_ids"] for b in blocks}, methods)
    counts = panel["counts"]
    if (counts["eligible_observed_rows"] != len(samples)
            or counts["eligible_observed_rows"] + counts["excluded_observed_rows"] != counts["source_observed_rows"]):
        raise ValueError("Source and evaluated observation denominators disagree")
    macro = ranking["group_macro_concordance"]
    informative = ranking["informative_units"] > 0
    signal = informative and macro["paired"] > .5 and macro["paired"] > macro["host_only"]
    return {"schema_version": 1, "status": "completed_secondary_matched_evaluation",
            "screen": ("positive_secondary_interaction_attribution_screen" if signal else
                       "no_positive_secondary_interaction_attribution_screen" if informative else "HOLD_no_informative_units"),
            "ranking": ranking, "predictions": methods, "source_counts": counts,
            "eligible_scoring_coverage": 1.0,
            "scored_fraction_of_full_source_observations": len(samples) / counts["source_observed_rows"],
            "all_source_unknown_pairs_retained": counts["source_unknown_pairs"],
            "group_differences": {g: {"paired_minus_host_only": v-ranking["per_group"]["host_only"][g],
                                       "paired_minus_tied": v-.5}
                                  for g, v in ranking["per_group"]["paired"].items()},
            "native_Ando_cases_scored": 0, "original_benchmark_replaced": False,
            "interpretation": "Descriptive retrospective spot-clearance screen. No independent-replicate, significance, reliable generalization, receptor causality, native plaque or assembly inference."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--score-manifest", type=Path, required=True)
    parser.add_argument("--expected-score-manifest-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a fresh evaluation output")
    policy_bytes = args.policy.read_bytes()
    policy = json.loads(policy_bytes)
    if sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Evaluator differs from prospective policy")
    verify_bindings(policy["inputs"])
    for module, name in [(blocked_ranking, "ranking"), (partition, "partition"), (learning, "learning")]:
        if sha(module.__file__) != policy["inputs"][name]["sha256"]:
            raise ValueError("Imported evaluator implementation differs: " + name)
    manifest_binding = {"path": str(args.score_manifest), "sha256": args.expected_score_manifest_sha256}
    manifest = json.loads(checked(manifest_binding).read_text())
    if len(manifest["scores"]) != 9:
        raise ValueError("Nine authenticated score receipts required")
    verify_bindings(manifest)
    panel_path = checked(manifest["panel"])
    panel = json.loads(panel_path.read_text())
    if (panel["helper_sha256"] != policy["inputs"]["panel_helper"]["sha256"]
            or panel["policy_sha256"] != policy["inputs"]["panel_policy"]["sha256"]):
        raise ValueError("Panel producer differs from policy")
    reports = []
    paths = set()
    for item in manifest["scores"]:
        path = checked(item)
        if path.resolve() in paths:
            raise ValueError("Duplicate score receipt path")
        paths.add(path.resolve())
        report = json.loads(path.read_text())
        if (report["helper_sha256"] != policy["inputs"]["fitting_helper"]["sha256"]
                or report["policy"]["sha256"] != policy["inputs"]["fitting_policy"]["sha256"]
                or report["inputs"]["panel"]["sha256"] != manifest["panel"]["sha256"]
                or (ROOT / report["inputs"]["panel"]["path"]).resolve() != panel_path.resolve()):
            raise ValueError("Score producer/policy/panel differs")
        verify_bindings(report)
        boundary = json.loads(checked(report["pre_scoring_boundary"]).read_text())
        expected = {scope: rec["sha256"] for scope, rec in report["bundles"].items()}
        if boundary != {"scoring_started": False, "manifest_sha256": expected}:
            raise ValueError("All three bundle hashes must be frozen before scoring")
        for scope, item in report["bundles"].items():
            verify_bundle(item, scope, report, policy)
        reports.append(report)
    result = summarize(panel, reports)
    verify_bindings(manifest)
    verify_bindings(policy["inputs"])
    if args.policy.read_bytes() != policy_bytes or sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Executed evaluator/policy changed")
    result["provenance"] = {"score_manifest": manifest_binding, "policy_sha256": hashlib.sha256(policy_bytes).hexdigest(),
                            "helper_sha256": policy["helper_sha256"]}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "screen": result["screen"], "source_counts": result["source_counts"]}))


if __name__ == "__main__":
    main()
