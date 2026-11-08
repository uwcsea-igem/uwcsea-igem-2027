"""Keep absent assignments distinct from successful family identification."""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("symmetric_panel_probe", ROOT / "scripts/probe-mmseqs-symmetric-panel.py")
PROBE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PROBE)


class SymmetricPanelProbeTests(unittest.TestCase):
    def setUp(self):
        self.controls = [{"copied_training_protein": f"t{i}", "protein": f"c{i}"} for i in range(4)]

    def test_family_equality_absence_and_discrepancies_are_separate(self):
        result = PROBE.compare_assignments(self.controls,
            {"t0": {"family": "f0"}, "t2": {"family": "f1"}, "t3": {"family": "f2"}},
            {"c0": {"family": "f0"}, "c3": {"family": "f3"}})
        self.assertEqual((result["same_assigned_family"], result["both_unassigned"], result["assignment_discrepancies"]), (1, 1, 2))
        self.assertEqual(result["comparisons"][1]["state"], "both_unassigned")
        self.assertIsNone(result["comparisons"][1]["training_family"])

    def test_all_absent_does_not_report_successful_family_assignment(self):
        result = PROBE.compare_assignments(self.controls, {}, {})
        self.assertEqual(result["same_assigned_family"], 0)
        self.assertEqual(result["both_unassigned"], 4)
        self.assertEqual(result["assignment_discrepancies"], 0)

    def test_invalid_rosters_or_out_of_roster_hits_fail(self):
        for controls, left, right in [([], {}, {}), (self.controls * 2, {}, {}),
                                     (self.controls, {"extra": {"family": "f"}}, {})]:
            with self.subTest(controls=controls), self.assertRaises(ValueError):
                PROBE.compare_assignments(controls, left, right)


if __name__ == "__main__":
    unittest.main()
