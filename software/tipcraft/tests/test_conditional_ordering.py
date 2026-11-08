from copy import deepcopy
import math
import unittest

from tipcraft.conditional_ordering import diagnose_conditional_ordering
from tipcraft.learning import Sample


def fixture():
    rows = [Sample(h + p, h, p, "g" + h, "g" + p, "synthetic")
            for p in ["p1", "p2"] for h in ["A", "B"]]
    predictions = {"block": {"model_id": "model1", "scores": {"Ap1": 1., "Bp1": 0., "Ap2": 11., "Bp2": 10.}}}
    features = {"block": ["host:f0", "phage:f1"]}
    return rows, predictions, features


class ConditionalOrderingTests(unittest.TestCase):
    def test_additive_counterexample_has_no_conditional_ordering(self):
        result = diagnose_conditional_ordering(*fixture())
        block = result["blocks"]["block"]
        self.assertEqual((block["selected_host_columns"], block["selected_phage_columns"]), (1, 1))
        self.assertEqual(block["eligible_host_pairs"], 1)
        self.assertEqual(block["host_pair_comparison_coverage"], 1.)
        pair = block["eligible_pairs"][0]
        self.assertEqual(pair["order_states"], [1])
        self.assertEqual(pair["difference_range"], 0.)
        self.assertFalse(pair["state_change"])
        self.assertFalse(pair["strict_reversal"])

    def test_strict_reversal_and_tie_change_are_distinct(self):
        for score, states, strict in [(-1., [-1, 1], True), (0., [0, 1], False)]:
            rows, predictions, features = fixture()
            predictions["block"]["scores"].update(Ap2=score, Bp2=0.)
            result = diagnose_conditional_ordering(rows, predictions, features)["blocks"]["block"]
            self.assertEqual(result["eligible_pairs"][0]["order_states"], states)
            self.assertEqual(result["host_pairs_with_state_change"], 1)
            self.assertEqual(result["host_pairs_with_strict_reversal"], int(strict))

    def test_incomplete_rectangles_are_reported_without_imputation(self):
        rows, predictions, features = fixture()
        rows = [r for r in rows if r.row_id != "Bp2"]
        predictions["block"]["scores"].pop("Bp2")
        block = diagnose_conditional_ordering(rows, predictions, features)["blocks"]["block"]
        self.assertEqual(block["eligible_host_pairs"], 0)
        self.assertEqual(block["host_pairs_with_one_jointly_observed_phage"], 1)
        self.assertEqual(block["unassessed_pairs"][0]["jointly_observed_phages"], ["p1"])
        self.assertEqual(block["host_pair_comparison_coverage"], 0.)

    def test_host_pairs_in_different_models_are_never_pooled(self):
        rows, _, _ = fixture()
        predictions = {"b1": {"model_id": "m1", "scores": {"Ap1": 1., "Bp1": 0.}},
                       "b2": {"model_id": "m2", "scores": {"Ap2": -1., "Bp2": 0.}}}
        result = diagnose_conditional_ordering(rows, predictions, {b: ["host:f"] for b in predictions})
        self.assertFalse(result["cross_model_comparisons"])
        self.assertTrue(all(b["eligible_host_pairs"] == 0 for b in result["blocks"].values()))

    def test_zero_common_phages_and_single_host_coverage(self):
        rows, _, _ = fixture()
        rows = [r for r in rows if r.row_id in {"Ap1", "Bp2"}]
        predictions = {"b": {"model_id": "m", "scores": {r.row_id: 0. for r in rows}}}
        block = diagnose_conditional_ordering(rows, predictions, {"b": ["phage:f"]})["blocks"]["b"]
        self.assertEqual(block["host_pairs_with_no_jointly_observed_phage"], 1)
        one = diagnose_conditional_ordering(rows[:1], {"b": {"model_id": "m", "scores": {"Ap1": 0.}}},
                                           {"b": ["phage:f"]})["blocks"]["b"]
        self.assertEqual(one["possible_host_pairs"], 0)
        self.assertIsNone(one["host_pair_comparison_coverage"])

    def test_score_variation_without_order_change_is_exposed(self):
        rows, predictions, features = fixture()
        predictions["block"]["scores"]["Ap2"] = 12.
        pair = diagnose_conditional_ordering(rows, predictions, features)["blocks"]["block"]["eligible_pairs"][0]
        self.assertEqual((pair["minimum_difference"], pair["maximum_difference"], pair["difference_range"]), (1., 2., 1.))
        self.assertFalse(pair["state_change"])

    def test_missing_duplicate_unknown_rows_and_sample_pairs_fail(self):
        for mutation in ["missing", "duplicate", "unknown", "samples"]:
            rows, predictions, features = deepcopy(fixture())
            if mutation == "missing":
                predictions["block"]["scores"].pop("Ap1")
            elif mutation == "duplicate":
                predictions["other"] = {"model_id": "m2", "scores": {"Ap1": 0.}}
                features["other"] = ["host:f"]
            elif mutation == "unknown":
                predictions["block"]["scores"]["unknown"] = 0.
            else:
                rows.append(rows[0])
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                diagnose_conditional_ordering(rows, predictions, features)

    def test_nonfinite_boolean_and_multiple_model_declarations_fail(self):
        for value in [math.nan, math.inf, -math.inf, True, "0"]:
            rows, predictions, features = fixture()
            predictions["block"]["scores"]["Ap1"] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                diagnose_conditional_ordering(rows, predictions, features)
        rows, predictions, features = fixture()
        predictions["block"]["model_id"] = ["m1", "m2"]
        with self.assertRaises(ValueError):
            diagnose_conditional_ordering(rows, predictions, features)

    def test_feature_roster_and_duplicate_columns_fail(self):
        for features in [{}, {"block": ["host:f", "host:f"]}, {"block": ["label"]}]:
            with self.subTest(features=features), self.assertRaises(ValueError):
                rows, predictions, _ = fixture()
                diagnose_conditional_ordering(rows, predictions, features)

    def test_derived_overflow_fails_and_input_order_does_not_change_result(self):
        rows, predictions, features = fixture()
        expected = diagnose_conditional_ordering(rows, predictions, features)
        reversed_scores = dict(reversed(list(predictions["block"]["scores"].items())))
        self.assertEqual(expected, diagnose_conditional_ordering(list(reversed(rows)),
            {"block": {"model_id": "model1", "scores": reversed_scores}}, features))
        predictions["block"]["scores"].update(Ap1=1e308, Bp1=-1e308)
        with self.assertRaisesRegex(ValueError, "not finite"):
            diagnose_conditional_ordering(rows, predictions, features)


if __name__ == "__main__":
    unittest.main()
