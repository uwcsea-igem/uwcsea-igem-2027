"""Keep source test records while removing related fitting records, without labels."""
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "role_overlap", Path(__file__).resolve().parents[3] / "scripts/audit-receptor-role-overlap.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def node(name, role):
    return {"container_key": name, "cohort_member": True, "published_model_role": role}


class RoleOverlapTests(unittest.TestCase):
    def test_all_test_members_remain_and_training_relatives_are_quarantined(self):
        nodes = [node("train1", "training"), node("train2", "training"), node("test", "testing"), node("mutant", "mutant")]
        result = AUDIT.partition(nodes, {"train1": "g1", "train2": "g2", "test": "g1", "mutant": "g1"}, set())
        self.assertEqual(result["buckets"]["candidate_testing"], ["test"])
        self.assertEqual(result["buckets"]["candidate_training"], ["train2"])
        self.assertEqual(result["buckets"]["training_group_quarantine"], ["train1"])
        self.assertTrue(result["all_cohort_records_accounted_once"])

    def test_target_exclusion_precedes_published_role_partition(self):
        result = AUDIT.partition([node("target", "training"), node("test", "testing")],
                                 {"target": "g", "test": "g"}, {"target", "test"})
        self.assertEqual(result["counts"]["target_related_cohort_exclusions"], 2)
        self.assertEqual(result["counts"]["candidate_testing"], 0)

    def test_duplicate_missing_group_and_unknown_role_fail(self):
        for nodes, groups in (([node("a", "training")]*2, {"a": "g"}),
                              ([node("a", "training")], {}), ([node("a", "unknown")], {"a": "g"})):
            with self.subTest(nodes=nodes), self.assertRaises(ValueError):
                AUDIT.partition(nodes, groups, set())


if __name__ == "__main__":
    unittest.main()
