"""Complete matched evaluation on synthetic observations and score receipts."""
from copy import deepcopy
from dataclasses import asdict
import importlib.util
from pathlib import Path
import hashlib
import json
import tempfile
import unittest

from tipcraft.learning import Sample
from tipcraft.partition import build_blocked_folds

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("secondary_evaluation", ROOT / "scripts/evaluate-secondary-scores.py")
EVAL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EVAL)


def fixture(phages_per_fold=1):
    samples, hf, pf, labels = [], {}, {}, {}
    for p in range(3 * phages_per_fold):
        pf[f"pg{p}"] = p % 3
        for h in range(3):
            for v in range(2):
                hf[f"hg{h}{v}"] = h
                row = Sample(f"r{p}{h}{v}", f"h{h}{v}", f"p{p}", f"hg{h}{v}", f"pg{p}", "synthetic")
                samples.append(row)
                labels[row.row_id] = 1-v
    blocks = build_blocked_folds(samples, hf, pf)
    panel = {"status": "frozen_secondary_reference_panel", "membership_generated": True,
             "samples": [asdict(s) for s in samples], "blocks": blocks,
             "host_group_folds": hf, "phage_group_folds": pf, "labels": labels,
             "counts": {"eligible_observed_rows": len(samples), "source_observed_rows": len(samples)+2,
                        "excluded_observed_rows": 2, "source_unknown_pairs": 3}}
    reports = []
    for i, block in enumerate(blocks):
        methods, bundles = {}, {}
        for j, scope in enumerate(sorted(EVAL.SCOPES)):
            model = f"{i*3+j:064x}"
            methods[scope] = {"model_id": model, "scores": {r: float(labels[r]) if scope == "paired" else 0.
                                                           for r in block["test_ids"]}}
            bundles[scope] = {"path": "synthetic-unused", "sha256": model}
        reports.append({"status": "completed_secondary_block_scores", "block_id": block["fold_id"],
                        "methods": methods, "bundles": bundles, "training_rows": len(block["train_ids"]),
                        "test_rows": len(block["test_ids"]), "test_labels_used_for_fit_or_selection": False,
                        "test_labels_returned_with_scores": False})
    return panel, reports


class SecondaryEvaluationTests(unittest.TestCase):
    def test_model_artifact_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            report = {"inputs": {}, "training_rows": 8}
            policy = {"inputs": {name: {"sha256": name} for name in ["learning", "fitting_helper", "fitting_policy"]}}
            contents = {"fit-spec.json": json.dumps({"implementation_sha256": "learning", "training_rows": 8,
                                                    "settings": {"feature_scope": "paired"}}),
                        "provenance.json": json.dumps({"inputs": {}, "helper": {"sha256": "fitting_helper"},
                                                       "policy": {"sha256": "fitting_policy"}}),
                        "fit-spec.sha256": "synthetic checksum fixture", "model.cbm": "synthetic artifact fixture",
                        "selected-features.json": "{}"}
            for name, value in contents.items():
                (directory / name).write_text(value)
            artifacts = {name: {"sha256": EVAL.sha(directory / name)} for name in contents}
            path = directory / "manifest.json"
            path.write_text(json.dumps({"status": "complete", "feature_scope": "paired", "artifacts": artifacts}))
            binding = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            EVAL.verify_bundle(binding, "paired", report, policy)
            (directory / "model.cbm").write_text("modified synthetic artifact")
            with self.assertRaisesRegex(ValueError, "artifact differs"):
                EVAL.verify_bundle(binding, "paired", report, policy)

    def test_full_matched_report_has_explicit_source_denominator(self):
        report = EVAL.summarize(*fixture())
        self.assertEqual(report["screen"], "positive_secondary_interaction_attribution_screen")
        self.assertEqual(report["ranking"]["group_macro_concordance"],
                         {"tied": .5, "paired": 1., "host_only": .5, "phage_only": .5})
        self.assertEqual(report["ranking"]["informative_units"], 9)
        self.assertEqual(report["scored_fraction_of_full_source_observations"], .9)
        self.assertEqual(report["eligible_scoring_coverage"], 1.)
        self.assertEqual(report["native_Ando_cases_scored"], 0)

    def test_equal_host_only_performance_gives_no_interaction_screen(self):
        panel, reports = fixture()
        for report in reports:
            report["methods"]["host_only"]["scores"] = deepcopy(report["methods"]["paired"]["scores"])
        result = EVAL.summarize(panel, reports)
        self.assertEqual(result["screen"], "no_positive_secondary_interaction_attribution_screen")

    def test_no_informative_units_is_hold_not_zero_or_success(self):
        panel, reports = fixture()
        panel["labels"] = {row: 1 for row in panel["labels"]}
        result = EVAL.summarize(panel, reports)
        self.assertEqual(result["screen"], "HOLD_no_informative_units")
        self.assertIsNone(result["ranking"]["group_macro_concordance"]["paired"])

    def test_host_only_scores_cannot_depend_on_phage_within_one_model(self):
        panel, reports = fixture(2)
        hosts = {s["row_id"]: s["host_id"] for s in panel["samples"]}
        scores = reports[0]["methods"]["host_only"]["scores"]
        first = next(iter(scores))
        self.assertEqual(sum(hosts[r] == hosts[first] for r in scores), 2)
        scores[first] = 1.
        with self.assertRaisesRegex(ValueError, "host_only"):
            EVAL.summarize(panel, reports)

    def test_incomplete_duplicate_or_mismatched_receipts_fail(self):
        for mutate in [lambda p, r: r.pop(), lambda p, r: r.append(deepcopy(r[0])),
                       lambda p, r: r[0]["methods"].pop("host_only"),
                       lambda p, r: r[0]["methods"]["paired"]["scores"].pop(next(iter(r[0]["methods"]["paired"]["scores"]))),
                       lambda p, r: r[0]["methods"]["paired"].update(model_id="wrong"),
                       lambda p, r: r[0].update(status="HOLD"),
                       lambda p, r: r[0].update(test_rows=999)]:
            panel, reports = fixture()
            mutate(panel, reports)
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                EVAL.summarize(panel, reports)

    def test_nonfinite_or_host_dependent_phage_only_scores_fail(self):
        for score in [float("nan"), float("inf"), 1.]:
            panel, reports = fixture()
            first = reports[0]["methods"]["phage_only"]["scores"]
            first[next(iter(first))] = score
            with self.subTest(score=score), self.assertRaises(ValueError):
                EVAL.summarize(panel, reports)

    def test_frozen_block_edits_and_false_source_counts_fail(self):
        panel, reports = fixture()
        panel["blocks"][0]["test_ids"].reverse()
        with self.assertRaisesRegex(ValueError, "frozen"):
            EVAL.summarize(panel, reports)
        panel, reports = fixture()
        panel["counts"]["source_observed_rows"] = 18
        with self.assertRaisesRegex(ValueError, "denominators"):
            EVAL.summarize(panel, reports)


if __name__ == "__main__":
    unittest.main()
