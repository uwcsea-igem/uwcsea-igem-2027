"""Synthetic checks for the standalone GPL adaptation; no aligner is run."""
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[3] / "scripts/third_party/viridic_directional.py"
SPEC = importlib.util.spec_from_file_location("viridic_directional", SCRIPT)
VIRIDIC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VIRIDIC)


def hsp(query="AAAAAA", subject=None, qstart=1, sstart=1, qlen=20, slen=20,
        qid="A", sid="B", reverse_subject=False, reverse_query=False):
    subject = query if subject is None else subject
    qend = qstart + len(query.replace("-", "")) - 1
    send = sstart + len(subject.replace("-", "")) - 1
    if reverse_subject:
        sstart, send = send, sstart
    if reverse_query:
        qstart, qend = qend, qstart
    return VIRIDIC.HSP(qid, sid, 1e-8, 40.0, qlen, slen, qstart, qend, sstart, send,
                       query, subject, sum(a == b and a != "-" for a, b in zip(query, subject)),
                       query.count("-") + subject.count("-"))


def tsv(row):
    return "\t".join(str(value) for value in asdict(row).values())


class IntergenomicTests(unittest.TestCase):
    def test_no_hits_and_reciprocal_genome_normalization(self):
        rows = [hsp(qlen=10, slen=20), hsp("AAAA", qid="B", sid="A", qlen=20, slen=10)]
        pairs = VIRIDIC.calculate_similarities(rows, {"A": 10, "B": 20, "C": 5})
        self.assertEqual(len(pairs), 3)
        self.assertAlmostEqual(pairs[0]["similarity_percent_unrounded"], 100 * 10 / 30)
        self.assertEqual(pairs[0]["similarity_percent"], 33.333)
        self.assertEqual(pairs[1]["similarity_percent"], 0)
        self.assertEqual(pairs[2]["similarity_percent"], 0)
        self.assertEqual(VIRIDIC.directional_nident([]).nident, 0)

    def test_no_overlap_sums_directional_matches(self):
        result = VIRIDIC.directional_nident([hsp("AAA"), hsp("CCC", qstart=5)])
        self.assertEqual((result.nident, result.original_aligned_query_bases), (6, 6))

    def test_contained_and_duplicate_partial_hsps_are_removed(self):
        result = VIRIDIC.directional_nident([hsp("AAAAAAAA"), hsp("AAAA", qstart=3)])
        self.assertEqual((result.nident, result.original_aligned_query_bases), (8, 8))
        duplicate = hsp(qstart=2)
        result = VIRIDIC.directional_nident([duplicate, duplicate, duplicate])
        self.assertEqual((result.nident, result.original_aligned_query_bases), (6, 6))

    def test_equal_range_ties_keep_input_order_not_highest_score(self):
        exact = hsp(qstart=2)
        mismatch = hsp(subject="AAACCC", qstart=2)
        self.assertEqual(VIRIDIC.directional_nident([mismatch, exact]).nident, 3)
        self.assertEqual(VIRIDIC.directional_nident([exact, mismatch]).nident, 6)

    def test_range_covered_by_union_is_removed_before_overlap_trimming(self):
        rows = [hsp(qstart=1), hsp(qstart=5), hsp(qstart=3)]
        result = VIRIDIC.directional_nident(rows)
        self.assertEqual(result.nident, 10)
        self.assertEqual(result.original_aligned_query_bases, 12)

    def test_partial_overlap_trims_identical_columns_and_retains_original_spans(self):
        rows = [hsp("ACGTAC", qstart=1), hsp("TACGTA", qstart=4)]
        result = VIRIDIC.directional_nident(rows)
        self.assertEqual((result.nident, result.original_aligned_query_bases), (9, 12))
        rows[0] = hsp("ACGTAC", "ACGAAT", qstart=1)
        # Only one of the three overlap columns is identical in the first HSP.
        self.assertEqual(VIRIDIC.directional_nident(rows).nident, 9)

    def test_overlap_trimming_uses_gapped_aligned_strings(self):
        rows = [hsp("ACG-T-A", "ACGGTTA"), hsp("TACC", qstart=4)]
        result = VIRIDIC.directional_nident(rows)
        self.assertEqual((result.nident, result.original_aligned_query_bases), (7, 9))

    def test_upstream_unanchored_suffix_regex_is_preserved(self):
        # The first regex match is earlier than the actual aligned-query tail.
        rows = [hsp("AAAAA-", "AAA-AA"), hsp("AAAA", qstart=4)]
        self.assertEqual(VIRIDIC.directional_nident(rows).nident, 6)

    def test_reverse_subject_and_single_reverse_query_are_supported(self):
        row = hsp("ACGT", reverse_subject=True)
        parsed = VIRIDIC.parse_blast_tsv(tsv(row), {"A": 20, "B": 20})
        self.assertGreater(parsed[0].sstart, parsed[0].send)
        self.assertEqual(VIRIDIC.directional_nident(parsed).nident, 4)
        row = hsp("ACGT", reverse_query=True)
        self.assertEqual(VIRIDIC.directional_nident([row]).nident, 4)
        with self.assertRaisesRegex(VIRIDIC.UnsupportedUpstreamCase, "reverse-query"):
            VIRIDIC.directional_nident([row, hsp(qstart=6)])

    def test_full_query_branch_and_ambiguous_sentinel_are_explicit(self):
        full = hsp("A" * 20)
        result = VIRIDIC.directional_nident([hsp(qstart=5), full])
        self.assertEqual((result.nident, result.original_aligned_query_bases), (20, 20))
        with self.assertRaisesRegex(VIRIDIC.UnsupportedUpstreamCase, "100000"):
            VIRIDIC.directional_nident([full, full])

    def test_parser_preserves_duplicates_and_rejects_invalid_rows(self):
        row = tsv(hsp())
        self.assertEqual(len(VIRIDIC.parse_blast_tsv(row + "\n" + row)), 2)
        alterations = {2: "nan", 3: "inf", 4: "0", 6: "0", 7: "21", 10: "AAA",
                       11: "AAZAAA", 12: "7", 13: "1"}
        for index, value in alterations.items():
            columns = row.split("\t")
            columns[index] = value
            with self.subTest(index=index), self.assertRaises(ValueError):
                VIRIDIC.parse_blast_tsv("\t".join(columns))
        with self.assertRaisesRegex(ValueError, "14"):
            VIRIDIC.parse_blast_tsv("A\tB")
        with self.assertRaisesRegex(ValueError, "declared"):
            VIRIDIC.parse_blast_tsv(row, {"A": 21, "B": 20})
        with self.assertRaisesRegex(ValueError, "Inconsistent"):
            VIRIDIC.parse_blast_tsv(row + "\n" + tsv(hsp(qlen=21)))

    def test_standalone_cli_marks_adaptation_and_preserves_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            blast, lengths, output = root / "blast.tsv", root / "lengths.json", root / "report.json"
            blast.write_text("")
            lengths.write_text(json.dumps({"A": 20, "B": 20}))
            args = [sys.executable, str(SCRIPT), "--blast-tsv", str(blast),
                    "--genome-lengths", str(lengths), "--output", str(output)]
            subprocess.run(args, check=True, capture_output=True)
            report = json.loads(output.read_text())
            self.assertFalse(report["upstream_equivalence_verified"])
            self.assertEqual(report["license"], "GPL-3.0-only")
            self.assertEqual(report["pairs"][0]["similarity_percent"], 0)
            bad = subprocess.run(args[:-1] + [str(blast)], capture_output=True)
            self.assertNotEqual(bad.returncode, 0)
            self.assertEqual(blast.read_text(), "")


if __name__ == "__main__":
    unittest.main()
