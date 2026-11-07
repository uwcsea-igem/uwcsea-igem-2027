import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from tipcraft.assembly import AssemblyEvidence, summarize
from tipcraft.benchmark import evaluate, load_cases, load_predictions
from tipcraft.rank import classification_metrics

CASES = Path(__file__).resolve().parents[1] / "benchmarks/retro-t7-t3/cases.csv"

class BenchmarkTests(unittest.TestCase):
    def test_baseline_and_no_claim_of_prediction(self):
        report = evaluate(load_cases(CASES))
        self.assertEqual(report["cases"], 10)
        self.assertEqual(report["baseline"]["accuracy"], 0.6)
        self.assertEqual(report["baseline"]["balanced_accuracy"], 0.5)
        self.assertEqual(report["model"]["status"], "not_evaluated")
        self.assertEqual(report["assembly"]["unknown"], 10)

    def test_partial_predictions_are_explicit(self):
        report = evaluate(load_cases(CASES), {"T7__BL21": 0.9})
        self.assertEqual(report["model"]["prediction_coverage"], 0.1)
        self.assertIsNone(report["model"]["balanced_accuracy"])
        self.assertTrue(report["gate"].startswith("HOLD"))
        with self.assertRaises(ValueError):
            evaluate(load_cases(CASES), {"unknown": 0.2})

    def test_metric_edge_cases(self):
        for labels, scores in [([], []), ([1], []), ([2], [0.2]), ([1], [float("nan")]), ([0], [1.1])]:
            with self.subTest(labels=labels, scores=scores), self.assertRaises(ValueError):
                classification_metrics(labels, scores)
        self.assertEqual(classification_metrics([0, 1], [0.2, 0.5])["balanced_accuracy"], 1)

    def test_unknown_is_not_negative_assembly(self):
        report = summarize([AssemblyEvidence(), AssemblyEvidence(False, "reference", "direct assay")])
        self.assertEqual(report["unknown"], 1)
        self.assertEqual(report["unsupported"], 1)
        self.assertEqual(report["direct_evidence_coverage"], 0.5)
        with self.assertRaises(ValueError):
            AssemblyEvidence(True)

    def test_prediction_file_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "predictions.csv"
            for body in ["", "case_id,score,method\nx,nan,m\n", "case_id,score,method\nx,0.5,m\nx,0.7,m\n", "case_id,score,method\nx,0.5,\n", "case_id,score,method\nx,0.5\n"]:
                p.write_text(body)
                with self.subTest(body=body), self.assertRaises(ValueError):
                    load_predictions(p)

    def test_duplicate_case_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "cases.csv"
            rows = load_cases(CASES)
            with p.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows([rows[0], rows[0]])
            with self.assertRaises(ValueError):
                load_cases(p)

    def test_cli_writes_reproducible_report_and_preserves_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "report.json"
            cmd = [sys.executable, "-m", "tipcraft", "benchmark", "--cases", str(CASES)]
            done = subprocess.run(cmd + ["--output", str(out)], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(done.stdout), json.loads(out.read_text()))
            before = CASES.read_bytes()
            bad = subprocess.run(cmd + ["--output", str(CASES)], capture_output=True, text=True)
            self.assertEqual(bad.returncode, 2)
            self.assertEqual(before, CASES.read_bytes())

if __name__ == "__main__":
    unittest.main()
