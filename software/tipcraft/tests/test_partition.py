from collections import Counter
from dataclasses import replace
import unittest

from tipcraft.learning import Sample
from tipcraft.partition import assign_group_folds, build_blocked_folds


def samples():
    return [Sample(f"r{h}{p}", f"h{h}", f"p{p}", f"hg{h}", f"pg{p}", "same-publication")
            for h in range(3) for p in range(3)]


class PartitionTests(unittest.TestCase):
    def test_group_assignment_is_order_independent_and_group_balanced(self):
        names = [f"group{i}" for i in range(11)]
        first = assign_group_folds(names, axis="host", folds=3, seed=2027)
        self.assertEqual(first, assign_group_folds(names[::-1] + names, axis="host", folds=3, seed=2027))
        counts = Counter(first.values())
        self.assertLessEqual(max(counts.values()) - min(counts.values()), 1)

    def test_invalid_or_insufficient_group_assignment_never_relaxes(self):
        cases = [(["a", "b"], {"axis": "host", "folds": 3, "seed": 0}),
                 (["a", "b"], {"axis": "row", "folds": 2, "seed": 0}),
                 (["a", "b"], {"axis": "host", "folds": True, "seed": 0}),
                 (["a", "b"], {"axis": "host", "folds": 2, "seed": -1}),
                 (["a", " "], {"axis": "host", "folds": 2, "seed": 0})]
        for names, kwargs in cases:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                assign_group_folds(names, **kwargs)

    def test_crossed_blocks_cover_each_pair_once_and_embargo_mixed_quadrants(self):
        rows = samples()
        hf, pf = {f"hg{i}": i for i in range(3)}, {f"pg{i}": i for i in range(3)}
        blocks = build_blocked_folds(rows, hf, pf)
        self.assertEqual(len(blocks), 9)
        self.assertEqual(Counter(i for b in blocks for i in b["test_ids"]), Counter(s.row_id for s in rows))
        block = blocks[0]
        self.assertEqual(block["test_ids"], ["r00"])
        self.assertEqual(block["train_ids"], ["r11", "r12", "r21", "r22"])
        self.assertEqual(block["embargo_ids"], ["r01", "r02", "r10", "r20"])
        self.assertEqual(block["validation"]["publication_overlap_disclosed"], ["same-publication"])
        self.assertEqual(blocks, build_blocked_folds(rows[::-1], hf, pf))

    def test_one_relatedness_group_is_never_split_across_train_and_test(self):
        rows = samples() + [Sample(f"extra{p}", "h0-relative", f"p{p}", "hg0", f"pg{p}", "same-publication") for p in range(3)]
        blocks = build_blocked_folds(rows, {f"hg{i}": i for i in range(3)}, {f"pg{i}": i for i in range(3)})
        by_id = {s.row_id: s for s in rows}
        for block in blocks:
            for axis in ("host_group", "phage_group"):
                train = {getattr(by_id[i], axis) for i in block["train_ids"]}
                test = {getattr(by_id[i], axis) for i in block["test_ids"]}
                self.assertTrue(train.isdisjoint(test))

    def test_empty_blocks_bad_mappings_or_inconsistent_groups_are_rejected(self):
        rows = samples()
        hf, pf = {f"hg{i}": i for i in range(3)}, {f"pg{i}": i for i in range(3)}
        for hosts, phages in [({"hg0": 0}, pf), (hf | {"unused": 0}, pf),
                              ({"hg0": 0, "hg1": 2, "hg2": 2}, pf),
                              (hf, {"pg0": True, "pg1": 1, "pg2": 2})]:
            with self.subTest(hosts=hosts, phages=phages), self.assertRaises(ValueError):
                build_blocked_folds(rows, hosts, phages)
        with self.assertRaises(ValueError):
            build_blocked_folds([rows[0], rows[4], rows[8]], hf, pf)
        with self.assertRaises(ValueError):
            build_blocked_folds(rows + [replace(rows[0], row_id="new", phage_id="new", host_group="wrong")], hf, pf)
