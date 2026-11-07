import hashlib
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("phage_lengths", Path(__file__).resolve().parents[3] / "scripts/audit-genophi-phage-lengths.py")
LENGTHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LENGTHS)


class PhageSourceLengthTests(unittest.TestCase):
    def test_terminal_overlap_and_n_runs_are_measured_without_changes(self):
        data = b">source\nnnacnn\n"
        before = hashlib.sha256(data).hexdigest()
        result = LENGTHS.inspect_sequence(data)
        self.assertEqual(result["exact_prefix_suffix_match_lengths"], [1, 2, 6])
        self.assertEqual((result["leading_n_bases"], result["trailing_n_bases"]), (2, 2))
        self.assertEqual(result["total_n_bases"], 4)
        self.assertEqual(result["sequence_sha256"], hashlib.sha256(b"NNACNN").hexdigest())
        self.assertEqual(hashlib.sha256(data).hexdigest(), before)

    def test_expanded_comparison_stops_at_the_frozen_limit(self):
        result = LENGTHS.inspect_sequence(b">source\n" + b"A" * 130 + b"\n")
        self.assertEqual(result["exact_prefix_suffix_match_lengths"], list(range(1, 121)))

    def test_contigs_empty_records_and_invalid_bases_fail(self):
        for data in [b">a\nACGT\n>b\nACGT\n", b">empty\n>a\nACGT\n", b">a\nACGTZ\n", b"ACGT\n"]:
            with self.subTest(data=data), self.assertRaises(ValueError):
                LENGTHS.inspect_sequence(data)

    def test_git_blob_authentication_rejects_changed_source(self):
        data = b">a\nACGT\n"
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        entry = {"type": "blob", "sha": blob, "size": len(data)}
        self.assertEqual(LENGTHS.authenticate(data, entry), blob)
        with self.assertRaises(ValueError):
            LENGTHS.authenticate(data.replace(b"ACGT", b"ACGA"), entry)
