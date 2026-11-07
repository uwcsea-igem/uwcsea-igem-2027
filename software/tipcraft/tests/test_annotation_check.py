import importlib.util
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location("annotation_check", Path(__file__).resolve().parents[3] / "scripts/check-genophi-annotation.py")
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class AnnotationCheckTests(unittest.TestCase):
    def test_proteome_comparison_ignores_headers_order_and_terminal_stop_but_preserves_copies(self):
        first = CHECK.protein_multiset(b">a\nMACT*\n>b\nMACT\n>c\nGHI\n")
        second = CHECK.protein_multiset(b">z\nGHI*\n>y\nMA\nCT*\n>x\nMACT*\n")
        self.assertTrue(CHECK.compare_multisets(first, second)["exact_multiset_match"])
        self.assertEqual(CHECK.multiset_digest(first), CHECK.multiset_digest(second))
        fewer = CHECK.protein_multiset(b">a\nMACT\n>b\nGHI\n")
        result = CHECK.compare_multisets(first, fewer)
        self.assertFalse(result["exact_multiset_match"])
        self.assertEqual(result["predicted_only_copies"], 1)

    def test_invalid_or_empty_protein_records_are_rejected(self):
        for value in (b">a\nMA*CT\n", b">a\nMACT**\n", b">a\n*\n", b">a\n", b">\nMACT\n", b"MACT\n", b">a\nMA\n>a\nCT\n"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                CHECK.protein_multiset(value)

    def test_host_header_normalization_retains_contig_boundaries_and_rejects_empty_contigs(self):
        value = CHECK.records_with_duplicate_headers(b">same\naCG\n>same\nTTA\n")
        self.assertEqual(value, [("c000000", "ACG"), ("c000001", "TTA")])
        with self.assertRaises(ValueError):
            CHECK.records_with_duplicate_headers(b">a\n>b\nACGT\n")

    def test_changed_cached_inputs_fail_checksums(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            path.write_bytes(b"changed")
            with self.assertRaises(ValueError):
                CHECK.checked(path, CHECK.digest(b"original"))
