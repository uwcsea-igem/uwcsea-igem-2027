import copy
import unittest

from tipcraft.blocked_ranking import evaluate_blocked_rankings
from tipcraft.learning import Sample


class BlockedRankingTests(unittest.TestCase):
    def fixture(self):
        rows = [Sample(str(i), f"h{i}", "p", f"hg{i}", "pg", "source") for i in range(4)]
        labels = {"0": 1, "1": 0, "2": 1, "3": 0}
        blocks = {"a": ["0", "1"], "b": ["2", "3"]}
        methods = {"paired": {"a": {"model_id": "fit-a", "scores": {"0": 2., "1": 1.}},
                              "b": {"model_id": "fit-b", "scores": {"2": -98., "3": -99.}}}}
        return rows, labels, blocks, methods

    def test_cross_model_offsets_never_create_host_comparisons(self):
        args = self.fixture()
        report = evaluate_blocked_rankings(*args)
        self.assertEqual(report["group_macro_concordance"], {"tied": .5, "paired": 1.})
        self.assertEqual(sum(u["positive_negative_pairs"] for u in report["units"]), 2)
        self.assertFalse(report["cross_model_host_comparisons"])
        args[3]["paired"]["b"]["scores"] = {"2": 10002., "3": 10001.}
        self.assertEqual(evaluate_blocked_rankings(*args)["group_macro_concordance"],
                         report["group_macro_concordance"])

    def test_ties_and_uninformative_units(self):
        rows, labels, blocks, methods = self.fixture()
        methods["paired"]["a"]["scores"] = {"0": 0., "1": 0.}
        labels["2"] = 0
        report = evaluate_blocked_rankings(rows, labels, blocks, methods)
        self.assertEqual(report["informative_units"], 1)
        self.assertEqual(report["informative_phages"], 1)
        self.assertEqual(report["informative_phage_groups"], 1)
        self.assertEqual(report["group_macro_concordance"]["paired"], .5)
        labels["0"] = 0
        self.assertIsNone(evaluate_blocked_rankings(rows, labels, blocks, methods)
                          ["group_macro_concordance"]["paired"])

    def test_group_macro_does_not_weight_large_groups_more(self):
        rows, labels, blocks, methods = self.fixture()
        # Two perfect queries share one group; a failing query has its own group.
        for phage, group, success in [("q", "pg", True), ("r", "other", False)]:
            ids = [phage + "0", phage + "1"]
            rows.extend(Sample(i, i, phage, i, group, "source") for i in ids)
            labels.update(dict(zip(ids, [1, 0])))
            blocks[phage] = ids
            methods["paired"][phage] = {"model_id": phage, "scores": dict(zip(ids, [1, 0] if success else [0, 1]))}
        report = evaluate_blocked_rankings(rows, labels, blocks, methods)
        self.assertEqual(report["per_group"]["paired"], {"other": 0., "pg": 1.})
        self.assertEqual(report["group_macro_concordance"]["paired"], .5)

    def test_missing_extra_and_nonfinite_scores_fail(self):
        for value in [float("nan"), float("inf"), True, "1"]:
            rows, labels, blocks, methods = self.fixture()
            methods["paired"]["a"]["scores"]["0"] = value
            with self.assertRaises(ValueError):
                evaluate_blocked_rankings(rows, labels, blocks, methods)
        rows, labels, blocks, methods = self.fixture()
        methods["host_only"] = copy.deepcopy(methods["paired"])
        del methods["host_only"]["a"]["scores"]["0"]
        with self.assertRaisesRegex(ValueError, "same rows"):
            evaluate_blocked_rankings(rows, labels, blocks, methods)

    def test_block_partition_must_be_complete_and_disjoint(self):
        for replacement in [["0", "2"], ["0", "0"], ["unknown"], ["0"]]:
            rows, labels, blocks, methods = self.fixture()
            blocks["a"] = replacement
            with self.assertRaises(ValueError):
                evaluate_blocked_rankings(rows, labels, blocks, methods)


if __name__ == "__main__":
    unittest.main()
