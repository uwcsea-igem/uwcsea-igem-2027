"""Synthetic saved-score/provenance fixtures only; no model execution."""
from copy import deepcopy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("native_transfer_evaluation", ROOT / "scripts/evaluate-native-transfer.py")
EVAL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EVAL)


def cases():
    rows = [{"case_id": p + "__" + h, "phage": p, "host": h,
             "label": str(int(p == "T7" or h == "BL21")), "source": "synthetic", "locator": "fixture", "interpretation": "synthetic"}
            for p in EVAL.PHAGES for h in EVAL.HOSTS]
    rows += [{"case_id": p + "__" + h, "phage": p, "host": h, "label": "0",
              "source": "synthetic", "locator": "fixture", "interpretation": "synthetic"}
             for p in ["T7_T3_gp17", "T3_T7_gp17"] for h in ["BW25113", "MG1655"]]
    return rows


def matrices():
    return {"paired": {"T7": dict(zip(EVAL.HOSTS, [3., 2., 1.])), "T3": dict(zip(EVAL.HOSTS, [3., 2., 1.]))},
            "host_only": {p: dict(zip(EVAL.HOSTS, [2., 3., 1.])) for p in EVAL.PHAGES},
            "phage_only": {p: {h: float(i) for h in EVAL.HOSTS} for i, p in enumerate(EVAL.PHAGES)}}


def source():
    return {"status": "completed_secondary_matched_evaluation", "screen": "no_positive_secondary_interaction_attribution_screen",
            "eligible_scoring_coverage": 1.0, "source_counts": {"eligible_observed_rows": 1},
            "ranking": {"group_macro_concordance": {"paired": .5, "host_only": .5}, "informative_units": 1}}


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")
    return {"path": str(path), "sha256": EVAL.sha(path)}


def receipt_fixture(directory):
    """Fake CBM bytes are authenticated artifacts, never loaded as models."""
    policy = {"helper_sha256": EVAL.sha(EVAL.__file__), "inputs": {}, "outcome_inputs": {}}
    for key in ["native_fitter", "native_fitting_policy", "learning", "source_evaluator", "source_evaluation_policy",
                "native_feature_helper", "native_feature_policy"]:
        policy["inputs"][key] = write(directory / (key + ".json"), {"synthetic": key})
    policy["inputs"]["retrospective"] = {"path": EVAL.retrospective.__file__, "sha256": EVAL.sha(EVAL.retrospective.__file__)}
    policy["inputs"]["benchmark"] = {"path": EVAL.benchmark.__file__, "sha256": EVAL.sha(EVAL.benchmark.__file__)}
    case_path = directory / "synthetic-cases.csv"
    with case_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases()[0]))
        writer.writeheader()
        writer.writerows(cases())
    policy["outcome_inputs"] = {"cases": {"path": str(case_path), "sha256": EVAL.sha(case_path)},
        "wish": write(directory / "synthetic-wish.json", {"scores": {p: dict(zip(EVAL.HOSTS, [1., 3., 2.])) for p in EVAL.PHAGES}})}
    policy["inputs"]["native_design"] = write(directory / "native-design.json", {"bound_inputs": list(policy["outcome_inputs"].values())})
    panel = {"status": "frozen_secondary_reference_panel", "membership_generated": True,
             "samples": [{"row_id": "source:one"}], "counts": source()["source_counts"]}
    panel_binding = write(directory / "panel.json", panel)
    score_bindings, predictions = [], {s: {} for s in EVAL.SCOPES}
    for block in sorted(EVAL.BLOCKS):
        methods = {s: {"model_id": "synthetic-model", "scores": {"synthetic-row": .25}} for s in EVAL.SCOPES}
        for scope in EVAL.SCOPES:
            predictions[scope][block] = methods[scope]
        score_bindings.append(write(directory / (block + ".json"), {"block_id": block,
            "status": "completed_secondary_block_scores", "inputs": {"panel": panel_binding}, "methods": methods}))
    manifest = write(directory / "source-score-manifest.json", {"panel": panel_binding, "scores": score_bindings})
    evaluation = source() | {"predictions": predictions, "provenance": {"score_manifest": manifest,
        "helper_sha256": policy["inputs"]["source_evaluator"]["sha256"], "policy_sha256": policy["inputs"]["source_evaluation_policy"]["sha256"]}}
    source_binding = write(directory / "source-evaluation.json", evaluation)
    phage = {"status": "completed_native_transfer_axis_features", "axis": "phage", "query_genomes": EVAL.PHAGES,
             "panel": panel_binding, "native_design": policy["inputs"]["native_design"],
             "feature_names": ["pc_0"], "vectors": {"T7": [1], "T3": [1]},
             "helper_sha256": policy["inputs"]["native_feature_helper"]["sha256"],
             "policy_sha256": policy["inputs"]["native_feature_policy"]["sha256"],
             "mapping_control_passed": True, "commands": [{"exit_code": 0} for _ in range(8)]}
    phage_binding = write(directory / "phage-features.json", phage)
    inputs = {"panel": panel_binding, "source_evaluation": source_binding, "phage_features": phage_binding}
    provenance = {"inputs": inputs, "helper": policy["inputs"]["native_fitter"], "policy": policy["inputs"]["native_fitting_policy"]}
    settings = {scope: {"feature_scope": scope} for scope in EVAL.SCOPES}
    plan = {"mode": "all_admitted_source_native_transfer", "settings": settings, "feature_names": ["host:pc_0", "phage:pc_0"],
            "native_pairs": [{"phage": p, "host": h} for p in EVAL.PHAGES for h in EVAL.HOSTS],
            "source_row_ids": ["source:one"], "native_labels_used": False,
            "model_selection_using_source_cv_performance": False, "provenance": provenance}
    bundles = {}
    for scope in EVAL.SCOPES:
        folder = directory / scope
        folder.mkdir()
        provenance_binding = write(folder / "provenance.json", provenance)
        spec_binding = write(folder / "fit-spec.json", {"implementation_sha256": policy["inputs"]["learning"]["sha256"],
            "training_rows": 1, "settings": settings[scope], "feature_names": plan["feature_names"], "provenance_sha256": provenance_binding["sha256"]})
        (folder / "fit-spec.sha256").write_text(spec_binding["sha256"] + "\n")
        selected = {"global_indices": [1 if scope == "phage_only" else 0], "names": ["phage:pc_0" if scope == "phage_only" else "host:pc_0"]}
        write(folder / "selected-features.json", selected)
        (folder / "model.cbm").write_bytes(b"SYNTHETIC ARTIFACT; NEVER LOAD")
        artifacts = {p.name: {"sha256": EVAL.sha(p)} for p in folder.iterdir()}
        bundles[scope] = write(folder / "manifest.json", {"status": "complete", "feature_scope": scope, "artifacts": artifacts,
            "selected_global_indices": selected["global_indices"], "selected_feature_names": selected["names"]})
    boundary = write(directory / "pre-scoring.json", {"scoring_started": False, "manifest_sha256": {s: b["sha256"] for s, b in bundles.items()}})
    receipt = {"status": "completed_native_transfer_scores", "methods": {s: {"model_id": bundles[s]["sha256"], "scores": m} for s, m in matrices().items()},
        "bundles": bundles, "native_scored_pairs_per_scope": 6, "exchange_construct_pairs_abstained": 4,
        "native_labels_used_for_fit_or_selection": False, "historical_WIsH_scores_used_for_fit_or_selection": False,
        "label_evaluation_executed": False, "helper_sha256": policy["inputs"]["native_fitter"]["sha256"],
        "policy": policy["inputs"]["native_fitting_policy"], "pre_scoring_boundary": boundary,
        "plan_before_fitting": write(directory / "plan.json", plan), "inputs": inputs, "source_training_rows": 1,
        "source_evaluation_screen_retained": evaluation["screen"],
        "partition": {"regime": "both_unseen", "publication_holdout": True,
                      "overlaps": {k: [] for k in ["host_id", "phage_id", "host_group", "phage_group", "publication"]}}}
    receipt_binding = write(directory / "native-scores.json", receipt)
    return policy, receipt_binding, receipt


class NativeTransferEvaluationTests(unittest.TestCase):
    def test_source_without_informative_units_remains_hold_despite_native_win(self):
        methods = {s: {"scores": m} for s, m in matrices().items()}
        wish = {"scores": {p: dict(zip(EVAL.HOSTS, [1., 3., 2.])) for p in EVAL.PHAGES}}
        no_contrast = source()
        no_contrast["screen"] = "HOLD_no_informative_units"
        no_contrast["ranking"] = {"group_macro_concordance": {"paired": None, "host_only": None}, "informative_units": 0}
        result = EVAL.evaluate(cases(), methods, wish, no_contrast, {"vectors": {"T7": [1], "T3": [1]}})
        self.assertTrue(result["primary"]["strict_native_screen_passed"])
        self.assertEqual(result["decision"], "HOLD_source_no_informative_units")
        self.assertEqual(result["source_evaluation"]["screen"], "HOLD_no_informative_units")

    def test_concordance_primary_coverage_failures_and_negative_source_retained(self):
        scores = matrices()
        methods = {s: {"scores": m} for s, m in scores.items()}
        wish = {"scores": {p: dict(zip(EVAL.HOSTS, [1., 3., 2.])) for p in EVAL.PHAGES}}
        original = deepcopy(methods)
        result = EVAL.evaluate(cases(), methods, wish, source(), {"vectors": {"T7": [1], "T3": [1]}})
        self.assertEqual(result["primary"]["scores"], {"paired": 1., "host_only": .5, "phage_only": .5, "WIsH": 0., "tied": .5})
        self.assertTrue(result["primary"]["strict_native_screen_passed"])
        self.assertEqual(result["decision"], "NO_GO_source_interaction_attribution_signal")
        self.assertEqual(result["coverage"]["fraction"], .6)
        self.assertEqual(len(result["coverage"]["abstained_case_ids"]), 4)
        self.assertFalse(result["T7_ranking_informative"])
        self.assertEqual(result["T3_post_score_label_permutations"]["paired"]["exact_upper_tail_fraction"], 1/3)
        self.assertEqual(len(result["failed_or_tied_positive_negative_comparisons"]["WIsH"]), 2)
        self.assertTrue(result["attribution"]["T7_T3_full_phage_feature_vectors_equal"])
        self.assertFalse(result["attribution"]["phage_vector_equality_is_new_admission_gate"])
        self.assertEqual(methods, original)

    def test_incomplete_nonfinite_or_different_case_subset_rejected(self):
        wish = {"scores": matrices()["paired"]}
        methods = {s: {"scores": m} for s, m in matrices().items()}
        for mutate in [lambda m: m["paired"]["scores"]["T3"].pop("MG1655"),
                       lambda m: m["paired"]["scores"]["T3"].update(BL21=float("nan")),
                       lambda m: m["paired"]["scores"].update(extra={"BL21": 1})]:
            changed = deepcopy(methods)
            mutate(changed)
            with self.assertRaises(ValueError):
                EVAL.evaluate(cases(), changed, wish, source(), {"vectors": {"T7": [], "T3": []}})
        with self.assertRaises(ValueError):
            EVAL.evaluate(cases()[:-1], methods, wish, source(), {"vectors": {"T7": [], "T3": []}})

    def test_bundle_boundary_same_panel_and_source_receipt_authentication(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            policy, binding, receipt = receipt_fixture(directory)
            actual, evaluation, phage = EVAL.authenticate(binding, policy)
            self.assertEqual(evaluation["screen"], "no_positive_secondary_interaction_attribution_screen")
            self.assertEqual(actual["methods"], receipt["methods"])
            model = directory / "paired/model.cbm"
            model.write_bytes(b"tampered")
            with self.assertRaisesRegex(ValueError, "artifact"):
                EVAL.authenticate(binding, policy)

    def test_missing_pre_score_scope_and_source_panel_substitution_fail(self):
        for defect in ["boundary", "panel"]:
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                policy, binding, receipt = receipt_fixture(directory)
                if defect == "boundary":
                    boundary = {"scoring_started": False, "manifest_sha256": {"paired": receipt["bundles"]["paired"]["sha256"]}}
                    receipt["pre_scoring_boundary"] = write(directory / "pre-scoring.json", boundary)
                else:
                    manifest_path = directory / "source-score-manifest.json"
                    manifest = json.loads(manifest_path.read_text())
                    manifest["panel"] = write(directory / "different-panel.json", {"different": True})
                    new_manifest_binding = write(manifest_path, manifest)
                    # All changed bindings remain internally consistent so the same-panel guard is exercised.
                    evaluation_path = directory / "source-evaluation.json"
                    evaluation = json.loads(evaluation_path.read_text())
                    evaluation["provenance"]["score_manifest"] = new_manifest_binding
                    receipt["inputs"]["source_evaluation"] = write(evaluation_path, evaluation)
                    # Other provenance intentionally remains frozen and will also reject substitution.
                binding = write(directory / "native-scores.json", receipt)
                with self.assertRaises(ValueError):
                    EVAL.authenticate(binding, policy)

    def test_bad_saved_score_hash_prevents_any_original_label_access(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            policy, binding, receipt = receipt_fixture(directory)
            policy["outcome_inputs"] = {"cases": {"path": "DO-NOT-OPEN.csv", "sha256": "x"}}
            policy_path = directory / "policy.json"
            write(policy_path, policy)
            with patch.object(EVAL.benchmark, "load_cases") as reader, self.assertRaises(ValueError):
                EVAL.execute(policy_path, Path(binding["path"]), "0" * 64)
            reader.assert_not_called()

    def test_full_synthetic_evaluation_and_strict_wish_comparison(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            policy, binding, receipt = receipt_fixture(directory)
            policy_path = directory / "policy.json"
            write(policy_path, policy)
            result = EVAL.execute(policy_path, Path(binding["path"]), binding["sha256"])
            self.assertEqual(result["status"], "completed_native_transfer_evaluation")
            self.assertEqual(result["coverage"]["scored_cases_per_method"], 6)
            self.assertFalse(result["new_fitting_or_inference_executed"])
        methods = {s: {"scores": m} for s, m in matrices().items()}
        result = EVAL.evaluate(cases(), methods, {"scores": matrices()["paired"]}, source(), {"vectors": {"T7": [0], "T3": [1]}})
        self.assertFalse(result["primary"]["strict_native_screen_passed"])  # Equality with WIsH is insufficient.


if __name__ == "__main__":
    unittest.main()
