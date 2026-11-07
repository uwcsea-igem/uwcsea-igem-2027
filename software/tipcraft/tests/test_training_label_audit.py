import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("label_audit", Path(__file__).resolve().parents[3] / "scripts/audit-genophi-labels.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class TrainingLabelAuditTests(unittest.TestCase):
    def test_blank_source_is_not_an_observed_negative(self):
        result = AUDIT.reconcile([{"strain": "A1", "phage": "P", "interaction": "0"}],
                                 [{"bacteria": "A-1", "P": ""}])
        self.assertEqual(result["eligible_label_rows_before_sequence_exclusions"], 0)
        self.assertEqual(len(result["source_blank_pairs"]), 1)

    def test_nonzero_assay_is_binarized_but_conflicts_are_visible(self):
        result = AUDIT.reconcile([{"strain": "A1", "phage": "P", "interaction": "1"},
                                  {"strain": "A1", "phage": "Q", "interaction": "1"}],
                                 [{"bacteria": "A-1", "P": "2.0", "Q": "0.0"}])
        self.assertEqual(result["concordant_nonblank_pairs"], 1)
        self.assertEqual(result["conflicting_nonblank_pairs"][0]["phage"], "Q")
        self.assertEqual(result["label_curation_status"], "source_conflicts_require_review")

    def test_ambiguous_names_and_duplicate_pairs_are_rejected(self):
        row = {"strain": "A1", "phage": "P", "interaction": "0"}
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            AUDIT.reconcile([row], [{"bacteria": "A-1", "P": "0.0"}, {"bacteria": "A_1", "P": "0.0"}])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            AUDIT.reconcile([row, row], [{"bacteria": "A-1", "P": "0.0"}])

    def test_missing_phage_and_unrecognized_assay_value_are_not_zero(self):
        row = {"strain": "A1", "phage": "P", "interaction": "0"}
        for source in [{"bacteria": "A-1"}, {"bacteria": "A-1", "P": "nan"}, {"bacteria": "A-1", "P": "-1"}]:
            with self.subTest(source=source), self.assertRaises(ValueError):
                AUDIT.reconcile([row], [source])
