import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("phage_features", ROOT / "scripts/check-genophi-phage-features.py")
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)
gp_spec = importlib.util.spec_from_file_location("gp17_mapping", ROOT / "scripts/check-genophi-gp17-feature-mapping.py")
gp17 = importlib.util.module_from_spec(gp_spec)
gp_spec.loader.exec_module(gp17)


class FeaturePilotTests(unittest.TestCase):
    def test_normalization_keeps_copies_x_and_rejects_internal_or_double_stops(self):
        self.assertEqual(pilot.proteins(b">a\nMXX*\n>b\nMXX*\n"), [("a", "MXX"), ("b", "MXX")])
        for data in [b">a\nMA**\n", b">a\nMA*MX\n"]:
            with self.assertRaises(ValueError):
                pilot.proteins(data)
        with self.assertRaises(ValueError):
            pilot.proteins(b">a\nMA\n>a\nMX\n")

    def test_best_hits_float32_ties_keep_raw_first_and_report_cross_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.tsv"
            path.write_text("q\tb\t40\t0.800000001\t1e-60\t0\t2\t3\t0\t2\t3\n"
                            "q\ta\t99\t0.800000002\t2e-60\t0\t2\t3\t0\t2\t3\n")
            hits, summary = pilot.best_hits(path, {"q": "MAA"}, {"a": "MAA", "b": "MAA"},
                                            {"a": "fa", "b": "fb"})
            self.assertEqual(hits["q"]["target"], "b")
            self.assertEqual(summary["cross_family_best_hit_ties"][0]["families"], ["fa", "fb"])

    def test_hits_reject_nontraining_targets_nonfinite_and_bad_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.tsv"
            for line in ["q\tx\t1\t0.8\t1e-4\t0\t2\t3\t0\t2\t3\n",
                         "q\ta\t1\tnan\t1e-4\t0\t2\t3\t0\t2\t3\n",
                         "q\ta\t1\t0.8\t1e-4\t0\t3\t3\t0\t2\t3\n"]:
                path.write_text(line)
                with self.assertRaises(ValueError):
                    pilot.best_hits(path, {"q": "MAA"}, {"a": "MAA"}, {"a": "fa"})

    def test_family_coverage_requires_every_training_protein_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "clusters.tsv"
            for text in ["a\ta\n", "a\ta\na\ta\n", "a\ta\nx\tb\n"]:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    pilot.read_families(path, {"a": "MA", "b": "MA"})
            path.write_text("a\ta\na\tb\n")
            self.assertEqual(pilot.read_families(path, {"a": "MA", "b": "MA"}), {"a": "a", "b": "a"})

    def test_occurrence_train_only_consolidation_or_and_explicit_zero_rows(self):
        hits = {"a": {"family": "f1"}, "b": {"family": "f1"},
                "c": {"family": "f2"}, "d": {"family": "f2"}, "e": {"family": "singleton"}}
        owners = {"a": "g1", "b": "g2", "c": "g1", "d": "g2", "e": "g1"}
        fmap, occurrence, summary = pilot.consolidate(["g1", "g2", "g3"], owners, hits)
        self.assertEqual(fmap, {"f1": "pc_0", "f2": "pc_0"})
        self.assertEqual(occurrence["f1"], (1, 1, 0))
        features, rows = pilot.project(["target", "zero"], {"q1": "target", "q2": "target"},
                                      {"q1": {"family": "f1"}, "q2": {"family": "f2"}}, fmap)
        self.assertEqual(rows, {"target": [1], "zero": [0]})
        self.assertEqual(summary["features_after_exact_occurrence_consolidation"], 1)

    def test_corrupt_input_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "protein.faa"
            path.write_bytes(b">p\nMA\n")
            with self.assertRaises(ValueError):
                pilot.checked(path, "0" * 64)

    def test_replay_keeps_external_context_separate_from_current_machine(self):
        original = {"platform": "original-system", "architecture": "x86_64", "rss_scope": "original children"}
        with patch.object(pilot.platform, "machine", return_value="arm64"), patch.object(pilot.platform, "platform", return_value="current-system"):
            result = pilot.contexts(original)
        self.assertEqual(result["external_execution_context"]["architecture"], "x86_64")
        self.assertEqual(result["external_execution_context"]["platform"], "original-system")
        self.assertEqual(result["postprocessing_context"]["architecture"], "arm64")
        self.assertEqual(original, {"platform": "original-system", "architecture": "x86_64", "rss_scope": "original children"})

    def test_posthoc_gp17_uses_exact_hash_and_requires_raw_best_hit_consistency(self):
        reference = {"protein_id": "ref", "protein_sequence_sha256": "abc"}
        protein_map = [{"genome": "T7", "role": "target", "sequence_sha256": "abc", "protein": "q", "length": 10},
                       {"genome": "T7", "role": "target", "sequence_sha256": "different", "protein": "q2", "length": 10}]
        result = gp17.summarize_reference("T7", reference, protein_map, {}, {}, gp17.Counter())
        self.assertEqual(result["exact_reference_copies_in_target_proteome"], 1)
        self.assertEqual(result["status"], "no_hit_returned")
        self.assertIsNone(result["copies"][0]["retained_feature"])
        with self.assertRaises(ValueError):
            gp17.summarize_reference("T7", reference, protein_map, {}, {}, gp17.Counter({"q": 1}))


if __name__ == "__main__":
    unittest.main()
