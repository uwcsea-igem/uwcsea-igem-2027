import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from tipcraft.benchmark import load_cases
from tipcraft.retrospective import evaluate_scores, leakage_audit
from tipcraft.wish import parse_matrix, run_scores, sequence

ROOT = Path(__file__).resolve().parents[3]
CASES = ROOT / "software/tipcraft/benchmarks/retro-t7-t3/cases.csv"
NAMED = ROOT / "data/processed/wish/named/scores.json"
ANONYMOUS = ROOT / "data/processed/wish/anonymous/scores.json"
MANIFEST = ROOT / "software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json"

class WishTests(unittest.TestCase):
    def setUp(self):
        self.cases = load_cases(CASES)
        self.scores = json.loads(NAMED.read_text())

    def test_real_result_has_failure_and_abstentions(self):
        report = evaluate_scores(self.cases, self.scores)
        self.assertEqual(report["coverage"], 0.6)
        self.assertEqual(len(report["abstained_case_ids"]), 4)
        self.assertEqual(report["primary"]["top1_documented_positive_rate"], 0.5)
        self.assertEqual(report["per_phage"]["T3"]["pair_concordance"], 0)
        self.assertIsNone(report["per_phage"]["T7"]["pair_concordance"])
        self.assertEqual(report["forced_choice_sensitivity"]["balanced_accuracy"], 0.375)
        self.assertEqual(report["always_positive_on_same_subset"]["balanced_accuracy"], 0.5)

    def test_renamed_upstream_run_matches_and_labels_do_not_mutate_scores(self):
        frozen = copy.deepcopy(self.scores)
        audit = leakage_audit(self.cases, self.scores, json.loads(ANONYMOUS.read_text()))
        self.assertTrue(audit["name_and_fasta_header_invariance"])
        self.assertEqual(audit["t3_exact_permutation_upper_tail_fraction"], 1)
        self.assertFalse(audit["prospective_blind_validation"])
        self.assertEqual(self.scores, frozen)

    def test_name_change_failure_is_not_hidden(self):
        altered = json.loads(ANONYMOUS.read_text())
        altered["scores"]["T3"]["BL21"] += 0.01
        audit = leakage_audit(self.cases, self.scores, altered)
        self.assertFalse(audit["name_and_fasta_header_invariance"])
        self.assertEqual(audit["status"], "failed")
        self.assertTrue(audit["conclusion"].startswith("FAILED"))

    def test_reusing_named_run_or_mismatched_provenance_is_rejected(self):
        with self.assertRaises(ValueError):
            leakage_audit(self.cases, self.scores, self.scores)
        altered = json.loads(ANONYMOUS.read_text())
        altered["alpha"] = 17
        with self.assertRaises(ValueError):
            leakage_audit(self.cases, self.scores, altered)

    def test_ties_are_not_broken_by_host_name(self):
        scores = copy.deepcopy(self.scores)
        scores["scores"]["T3"] = {h: -1.0 for h in ("BL21", "BW25113", "MG1655")}
        report = evaluate_scores(self.cases, scores)
        self.assertAlmostEqual(report["per_phage"]["T3"]["top1_documented_positive"], 1/3)
        self.assertEqual(report["per_phage"]["T3"]["pair_concordance"], 0.5)

    def test_invalid_matrix_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "matrix"
            for text in ["", "T3\nBL21 nan\n", "T3\nBL21 -1\nBL21 -2\n", "T3 T3\nBL21 -1 -2\n"]:
                path.write_text(text)
                with self.subTest(text=text), self.assertRaises(ValueError):
                    parse_matrix(path)

    def test_fasta_ambiguity_and_multiple_records_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "seq.fna"
            for text in [">x\nACNT", ">x\nACGT\n>y\nACGT", "ACGT"]:
                path.write_text(text)
                with self.subTest(text=text), self.assertRaises(ValueError):
                    sequence(path)

    def test_wrong_genome_stops_before_external_execution(self):
        manifest = json.loads(MANIFEST.read_text())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path / "hosts").mkdir()
            first = next(iter(manifest["hosts"]))
            (path / "hosts" / (first + ".fna")).write_text(">x\nACGT\n")
            with patch("tipcraft.wish.subprocess.run") as external, self.assertRaises(ValueError):
                run_scores(MANIFEST, path, path / "missing", path / "out")
            external.assert_not_called()

    def test_output_cannot_overwrite_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            manifest = work / "scores.json"
            original = MANIFEST.read_bytes()
            manifest.write_bytes(original)
            with patch("tipcraft.wish.subprocess.run") as external, self.assertRaises(ValueError):
                run_scores(manifest, work / "genomes", work / "missing", work)
            external.assert_not_called()
            self.assertEqual(manifest.read_bytes(), original)

    def test_empty_and_nonfinite_scores_are_rejected(self):
        for scores in [{}, {"T3": {"BL21": float("nan")}}]:
            with self.subTest(scores=scores), self.assertRaises(ValueError):
                evaluate_scores(self.cases, {"scores": scores})

if __name__ == "__main__":
    unittest.main()
