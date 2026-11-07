import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("host_overlap", ROOT / "scripts/audit-genophi-host-overlap.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
POLICY = json.loads((ROOT / "software/tipcraft/benchmarks/retro-t7-t3/genophi-overlap-policy.json").read_text())["host_exclusion"]


class HostOverlapTests(unittest.TestCase):
    def test_near_target_requires_reciprocal_aligned_fraction(self):
        self.assertEqual(AUDIT.classify_pair(99.8, 98, 98, POLICY), "near_target_exclude")
        self.assertEqual(AUDIT.classify_pair(99.8, 40, 99, POLICY), "outside_declared_near_target_rule")

    def test_boundary_is_excluded_pending_review(self):
        self.assertEqual(AUDIT.classify_pair(99.48, 99, 99, POLICY), "boundary_exclude_pending_review")
        self.assertEqual(AUDIT.classify_pair(99.8, 89.8, 99, POLICY), "boundary_exclude_pending_review")

    def test_missing_comparison_and_named_target_are_not_retained(self):
        result = AUDIT.summarize([], {"x": "host"}, {"q": "BL21"}, POLICY)
        self.assertEqual(result["excluded_hosts"], ["host"])
        row = {"Ref_file": "x", "Query_file": "q", "ANI": "90", "Align_fraction_ref": "95", "Align_fraction_query": "95"}
        result = AUDIT.summarize([row], {"x": "BL21"}, {"q": "BL21"}, POLICY)
        self.assertEqual(result["excluded_hosts"], ["BL21"])

    def test_bad_values_unknown_ids_and_duplicates_are_rejected(self):
        for value in [float("nan"), float("inf"), -1, 101]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                AUDIT.classify_pair(value, 90, 90, POLICY)
        row = {"Ref_file": "x", "Query_file": "q", "ANI": "98", "Align_fraction_ref": "95", "Align_fraction_query": "95"}
        with self.assertRaises(ValueError):
            AUDIT.summarize([row], {}, {"q": "BL21"}, POLICY)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            AUDIT.summarize([row, row], {"x": "host"}, {"q": "BL21"}, POLICY)
