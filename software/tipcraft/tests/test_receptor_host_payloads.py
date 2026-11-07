"""Host-source inventory must preserve record boundaries and unresolved formats."""
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "host_payload_audit", Path(__file__).resolve().parents[3] / "scripts/audit-receptor-host-payloads.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class HostPayloadTests(unittest.TestCase):
    def test_distinct_records_and_source_case_are_preserved(self):
        records = AUDIT.fasta(b">chr source\nacgt\n>plasmid other\nNNAA\n", set("ACGTN"))
        self.assertEqual(records, [("chr source", "acgt"), ("plasmid other", "NNAA")])

    def test_duplicate_empty_and_invalid_sequences_fail(self):
        for raw in (b">chr a\nACG\n>chr b\nTTA\n", b">chr\n", b">chr\nACX\n", b">\nACG\n"):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                AUDIT.fasta(raw, set("ACGTN"))

    def test_tabular_protein_text_is_not_guessed_to_be_fasta(self):
        with self.assertRaises(ValueError):
            AUDIT.fasta(b"gene\tprotein\nx\tMELK\n", set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"))


if __name__ == "__main__":
    unittest.main()
