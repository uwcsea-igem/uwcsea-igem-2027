"""Evaluate authenticated saved native scores after source-only fitting.

Original native labels and WIsH values are opened only after score, bundle,
pre-score boundary and same-panel source-evaluation authentication succeeds.
"""
import argparse
from copy import deepcopy
import hashlib
import json
import math
from numbers import Real
from pathlib import Path

from tipcraft import benchmark, retrospective

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/native-transfer-evaluation-policy.json"
SCOPES = ["paired", "host_only", "phage_only"]
PHAGES, HOSTS = ["T7", "T3"], ["BL21", "BW25113", "MG1655"]
BLOCKS = {f"host{i}-phage{j}" for i in range(3) for j in range(3)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked(binding):
    path = ROOT / binding["path"]
    if sha(path) != binding["sha256"]:
        raise ValueError("Authenticated artifact differs: " + str(path))
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


def read(binding):
    return json.loads(checked(binding).read_text())


def same_binding(first, second):
    return first["sha256"] == second["sha256"] and (ROOT / first["path"]).resolve() == (ROOT / second["path"]).resolve()


def score_matrix(scores):
    if not isinstance(scores, dict) or set(scores) != set(PHAGES):
        raise ValueError("Exactly T7/T3 native score rows required")
    for hosts in scores.values():
        if not isinstance(hosts, dict) or set(hosts) != set(HOSTS):
            raise ValueError("Exactly three original native hosts per phage required")
        if any(isinstance(v, bool) or not isinstance(v, Real) or not math.isfinite(v) for v in hosts.values()):
            raise ValueError("All six raw scores must be finite numeric values")


def verify_bundle(binding, scope, receipt, plan, policy):
    path = checked(binding)
    manifest = json.loads(path.read_text())
    files = {"fit-spec.json", "fit-spec.sha256", "model.cbm", "provenance.json", "selected-features.json"}
    if (manifest["status"] != "complete" or manifest["feature_scope"] != scope or set(manifest["artifacts"]) != files):
        raise ValueError("Complete scope-specific native bundle required")
    for name, item in manifest["artifacts"].items():
        checked({"path": str(path.parent / name), "sha256": item["sha256"]})
    spec_path = path.parent / "fit-spec.json"
    spec = json.loads(spec_path.read_text())
    provenance_path = path.parent / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    selected = json.loads((path.parent / "selected-features.json").read_text())
    if ((path.parent / "fit-spec.sha256").read_text() != sha(spec_path) + "\n"
            or spec["implementation_sha256"] != policy["inputs"]["learning"]["sha256"]
            or spec["training_rows"] != receipt["source_training_rows"]
            or spec["settings"] != plan["settings"][scope] or spec["settings"]["feature_scope"] != scope
            or spec["feature_names"] != plan["feature_names"]
            or spec["provenance_sha256"] != sha(provenance_path)
            or provenance != plan["provenance"] or provenance["inputs"] != receipt["inputs"]
            or provenance["helper"]["sha256"] != policy["inputs"]["native_fitter"]["sha256"]
            or provenance["policy"] != receipt["policy"]
            or selected["global_indices"] != manifest["selected_global_indices"]
            or selected["names"] != manifest["selected_feature_names"]):
        raise ValueError("Native bundle specification/provenance differs from saved-score receipt")


def authenticate(receipt_binding, policy):
    receipt = read(receipt_binding)
    if (receipt["status"] != "completed_native_transfer_scores" or set(receipt["methods"]) != set(SCOPES)
            or set(receipt["bundles"]) != set(SCOPES) or receipt["native_scored_pairs_per_scope"] != 6
            or receipt["exchange_construct_pairs_abstained"] != 4
            or receipt["native_labels_used_for_fit_or_selection"] is not False
            or receipt["historical_WIsH_scores_used_for_fit_or_selection"] is not False
            or receipt["label_evaluation_executed"] is not False
            or receipt["helper_sha256"] != policy["inputs"]["native_fitter"]["sha256"]
            or receipt["policy"]["sha256"] != policy["inputs"]["native_fitting_policy"]["sha256"]):
        raise ValueError("Completed six-pair source-only native score receipt required")
    verify_bindings(receipt)
    expected = {s: receipt["bundles"][s]["sha256"] for s in SCOPES}
    if read(receipt["pre_scoring_boundary"]) != {"scoring_started": False, "manifest_sha256": expected}:
        raise ValueError("All three manifests must be frozen before native scoring")
    plan = read(receipt["plan_before_fitting"])
    pairs = [{"phage": p, "host": h} for p in PHAGES for h in HOSTS]
    if (plan["mode"] != "all_admitted_source_native_transfer" or plan["native_pairs"] != pairs
            or plan["native_labels_used"] is not False or plan["model_selection_using_source_cv_performance"] is not False
            or set(plan["settings"]) != set(SCOPES)):
        raise ValueError("Fixed source-only native plan and candidate order required")
    for scope in SCOPES:
        method = receipt["methods"][scope]
        if method["model_id"] != expected[scope]:
            raise ValueError("Saved native scores do not identify their frozen model")
        score_matrix(method["scores"])
        verify_bundle(receipt["bundles"][scope], scope, receipt, plan, policy)
    panel = read(receipt["inputs"]["panel"])
    if (panel["status"] != "frozen_secondary_reference_panel" or panel["membership_generated"] is not True
            or len(panel["samples"]) != receipt["source_training_rows"]
            or plan["source_row_ids"] != [s["row_id"] for s in panel["samples"]]):
        raise ValueError("All admitted observed source rows must enter final refit")
    source = read(receipt["inputs"]["source_evaluation"])
    if (source["status"] != "completed_secondary_matched_evaluation" or source["eligible_scoring_coverage"] != 1.0
            or source["source_counts"] != panel["counts"]
            or source["provenance"]["helper_sha256"] != policy["inputs"]["source_evaluator"]["sha256"]
            or source["provenance"]["policy_sha256"] != policy["inputs"]["source_evaluation_policy"]["sha256"]
            or receipt["source_evaluation_screen_retained"] != source["screen"]):
        raise ValueError("Complete same-panel source evaluation must be retained, including negative results")
    manifest = read(source["provenance"]["score_manifest"])
    if not same_binding(manifest["panel"], receipt["inputs"]["panel"]) or len(manifest["scores"]) != 9:
        raise ValueError("Source evaluation must bind all nine scores and the same source panel")
    verify_bindings(manifest)
    blocks, paths = set(), set()
    for item in manifest["scores"]:
        path = checked(item).resolve()
        score = read(item)
        name = score["block_id"]
        if (path in paths or name in blocks or name not in BLOCKS or score["status"] != "completed_secondary_block_scores"
                or not same_binding(score["inputs"]["panel"], manifest["panel"])):
            raise ValueError("Unique completed source blocks for the same panel required")
        paths.add(path)
        blocks.add(name)
        for scope in SCOPES:
            if source["predictions"][scope][name] != score["methods"][scope]:
                raise ValueError("Persisted source evaluation predictions differ from score receipts")
    if blocks != BLOCKS:
        raise ValueError("Complete nine-block source execution required")
    if receipt["partition"]["regime"] != "both_unseen" or receipt["partition"]["publication_holdout"] is not True:
        raise ValueError("Both-unseen source/native groups and publication holdout required")
    if any(receipt["partition"]["overlaps"][k] for k in ["host_id", "phage_id", "host_group", "phage_group", "publication"]):
        raise ValueError("Native source partition overlap prohibited")
    phage = read(receipt["inputs"]["phage_features"])
    if (phage["status"] != "completed_native_transfer_axis_features" or phage["axis"] != "phage"
            or phage["query_genomes"] != PHAGES or not same_binding(phage["panel"], receipt["inputs"]["panel"])
            or phage["native_design"]["sha256"] != policy["inputs"]["native_design"]["sha256"]
            or phage["helper_sha256"] != policy["inputs"]["native_feature_helper"]["sha256"]
            or phage["policy_sha256"] != policy["inputs"]["native_feature_policy"]["sha256"]
            or phage["mapping_control_passed"] is not True or len(phage["commands"]) != 8
            or any(type(c["exit_code"]) is not int or c["exit_code"] != 0 for c in phage["commands"])):
        raise ValueError("Authenticated exact native phage feature report required")
    verify_bindings(phage)
    for p in PHAGES:
        vector = phage["vectors"][p]
        if (len(vector) != len(phage["feature_names"]) or any(type(v) is not int or v not in (0, 1) for v in vector)):
            raise ValueError("Native phage vectors must be complete binary vectors")
    for h in HOSTS:
        if receipt["methods"]["host_only"]["scores"]["T7"][h] != receipt["methods"]["host_only"]["scores"]["T3"][h]:
            raise ValueError("Host-only attribution control varies by phage")
    if any(len(set(receipt["methods"]["phage_only"]["scores"][p].values())) != 1 for p in PHAGES):
        raise ValueError("Phage-only attribution control varies by host")
    return receipt, source, phage


def evaluate(cases, methods, wish, source, phage):
    """Post-score descriptive evaluation; no model loading, fitting or tuning."""
    expected = {(p, h) for p in PHAGES for h in HOSTS} | {
        (p, h) for p in ["T7_T3_gp17", "T3_T7_gp17"] for h in ["BW25113", "MG1655"]}
    if (len(cases) != 10 or len({c["case_id"] for c in cases}) != 10
            or {(c["phage"], c["host"]) for c in cases} != expected
            or any(c["label"] not in {"0", "1"} for c in cases)):
        raise ValueError("Exact ten original observations required; no favorable subset")
    if (sum(c["label"] == "1" for c in cases if c["phage"] == "T7") != 3
            or sum(c["label"] == "1" for c in cases if c["phage"] == "T3") != 1):
        raise ValueError("Frozen all-positive T7 / single-positive T3 endpoint required")
    if set(methods) != set(SCOPES):
        raise ValueError("All three saved scope score matrices required")
    matrices = {scope: deepcopy(methods[scope]["scores"]) for scope in SCOPES}
    matrices["WIsH"] = deepcopy(wish["scores"])
    matrices["tied"] = {p: {h: 0.0 for h in HOSTS} for p in PHAGES}
    reports, primary, permutations, failures = {}, {}, {}, {}
    for method, scores in matrices.items():
        score_matrix(scores)
        inference = {"method": method, "scores": scores}
        report = retrospective.evaluate_scores(cases, inference)
        report.pop("gate")  # Historical helper hardcodes a WIsH-specific gate.
        report["secondary_top1"] = report.pop("primary")
        for row in report["cases"]:
            row["raw_score"] = row.pop("mean_log_likelihood")
        values = [p["pair_concordance"] for p in report["per_phage"].values() if p["positive_negative_pairs"]]
        primary[method] = sum(values) / len(values) if values else None
        reports[method] = report
        placements = []
        for host in HOSTS:
            altered = [dict(c, label=str(int(c["host"] == host))) if c["phage"] == "T3" else dict(c) for c in cases]
            value = retrospective.evaluate_scores(altered, inference)["per_phage"]["T3"]["pair_concordance"]
            placements.append({"positive_host": host, "pair_concordance": value})
        observed = report["per_phage"]["T3"]["pair_concordance"]
        permutations[method] = {"single_positive_placements": placements,
                                "exact_upper_tail_fraction": sum(p["pair_concordance"] >= observed for p in placements)/3,
                                "best_possible_upper_tail_fraction": 1/3}
        contrasts = []
        for p in PHAGES:
            positive = [c for c in cases if c["phage"] == p and c["label"] == "1"]
            negative = [c for c in cases if c["phage"] == p and c["label"] == "0"]
            for a in positive:
                for b in negative:
                    x, y = scores[p][a["host"]], scores[p][b["host"]]
                    contrasts.append({"phage": p, "positive_host": a["host"], "negative_host": b["host"],
                                      "positive_score": x, "negative_score": y,
                                      "outcome": "concordant" if x > y else "tied" if x == y else "misordered"})
        failures[method] = [c for c in contrasts if c["outcome"] != "concordant"]
    paired = primary["paired"]
    screen = paired > .5 and paired > primary["WIsH"]
    source_positive = source["screen"] == "positive_secondary_interaction_attribution_screen"
    equal = phage["vectors"]["T7"] == phage["vectors"]["T3"]
    return {"schema_version": 1, "status": "completed_native_transfer_evaluation",
            "primary": {"endpoint": "Mean within-phage positive-negative concordance across informative native phages; ties 0.5",
                        "informative_phages": ["T3"], "scores": primary,
                        "strict_native_screen_passed": screen,
                        "screen": "positive_exploratory_native_ranking_screen" if screen else "no_positive_native_ranking_screen"},
            "coverage": {"scored_cases_per_method": 6, "original_cases": 10, "fraction": .6,
                         "abstained_case_ids": reports["paired"]["abstained_case_ids"], "complete_ten_case_gate_passed": False},
            "raw_scores": matrices, "methods": reports, "failed_or_tied_positive_negative_comparisons": failures,
            "T7_all_three_documented_positive": True, "T7_ranking_informative": False,
            "uniform_random_top_host_expected_documented_positive_fraction": 2/3,
            "T3_post_score_label_permutations": permutations,
            "attribution": {"paired_minus_host_only": paired-primary["host_only"],
                            "paired_minus_phage_only": paired-primary["phage_only"],
                            "host_only_matches_paired_concordance": primary["host_only"] == paired,
                            "T7_T3_full_phage_feature_vectors_equal": equal,
                            "paired_native_score_rows_equal": matrices["paired"]["T7"] == matrices["paired"]["T3"],
                            "phage_vector_equality_is_new_admission_gate": False,
                            "interpretation": "Equal phage vectors prevent assigning phage-dependent input effects; whole-proteome discrimination does not establish recognition-head representation."},
            "source_evaluation": {"screen": source["screen"], "group_macro_concordance": source["ranking"]["group_macro_concordance"],
                                  "informative_units": source["ranking"]["informative_units"], "eligible_scoring_coverage": source["eligible_scoring_coverage"],
                                  "source_counts": source["source_counts"], "negative_source_result_overridden_by_native": False},
            "decision": "HOLD_source_no_informative_units" if source["screen"] == "HOLD_no_informative_units" else
                        "HOLD_exploratory_signal_only" if source_positive and screen else
                        "NO_GO_source_interaction_attribution_signal" if not source_positive else "NO_GO_native_ranking_screen",
            "full_original_benchmark_gate": "HOLD_four_exchange_constructs_abstained",
            "saved_score_numeric_replay_executed": False, "new_fitting_or_inference_executed": False,
            "actual_anonymous_native_replay_verified": False,
            "independence": {"source_native_both_unseen_and_publication_holdout_required": True,
                             "target_labels_historical_WIsH_or_donor_lookup_used_for_fit_or_selection": False,
                             "basis": "Authenticated pinned fitter/bundle provenance; no claim of blinded method selection",
                             "remaining_inference_control": "Actual separate anonymous native replay is not verified here"},
            "interpretation": "Exploratory retrospective diagnostic after labels were known during method selection. Public references do not authenticate experimental stocks. T3 alone is informative; best possible three-placement permutation fraction is 1/3. Native improvement cannot override poor source validation or establish generalization, recognition-head binding, assembly or productive infection. Negative labels mean no plaques detected under reported conditions."}


def execute(policy_path, score_path, expected_sha):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Native evaluator differs from frozen policy")
    verify_bindings(policy["inputs"])
    if (sha(retrospective.__file__) != policy["inputs"]["retrospective"]["sha256"]
            or sha(benchmark.__file__) != policy["inputs"]["benchmark"]["sha256"]):
        raise ValueError("Imported native evaluation implementation differs")
    binding = {"path": str(score_path), "sha256": expected_sha}
    receipt, source, phage = authenticate(binding, policy)  # Boundary before original-label access.
    design = read(policy["inputs"]["native_design"])
    for item in policy["outcome_inputs"].values():
        if item not in design["bound_inputs"]:
            raise ValueError("Outcome authority differs from frozen native design")
    cases = benchmark.load_cases(checked(policy["outcome_inputs"]["cases"]))
    wish = read(policy["outcome_inputs"]["wish"])
    result = evaluate(cases, receipt["methods"], wish, source, phage)
    verify_bindings(receipt)
    verify_bindings(policy["inputs"])
    verify_bindings(policy["outcome_inputs"])
    checked(binding)
    if policy_path.read_bytes() != policy_bytes or sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Executed evaluation source/policy changed")
    result["provenance"] = {"score_receipt": binding, "policy_sha256": hashlib.sha256(policy_bytes).hexdigest(),
                            "helper_sha256": policy["helper_sha256"], "inputs": policy["inputs"],
                            "outcome_inputs_opened_only_after_score_authentication": policy["outcome_inputs"]}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--expected-scores-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a fresh native evaluation output")
    result = execute(args.policy, args.scores, args.expected_scores_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "primary": result["primary"], "decision": result["decision"]}))
