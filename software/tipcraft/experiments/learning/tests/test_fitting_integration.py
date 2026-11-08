import importlib.util
from dataclasses import replace
import hashlib
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[5]
SPEC = importlib.util.spec_from_file_location("fit_check", ROOT / "scripts/check-learning-fit.py")
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class FittingIntegrationTests(unittest.TestCase):
    def test_phage_only_fit_load_and_host_invariance_on_synthetic_data(self):
        learning = CHECK.learning
        spec, X, y, samples = CHECK.fixture()
        perturbed = [[1 - v for v in row[:3]] + row[3:] for row in X]
        cfg = replace(learning.FitSettings(**spec["settings"]), feature_scope="phage_only")
        with tempfile.TemporaryDirectory() as temporary:
            bundles = [Path(temporary) / "original", Path(temporary) / "host-perturbed"]
            for matrix, bundle in zip([X, perturbed], bundles):
                manifest = learning.fit_bundle(matrix, y, samples, spec["features"], cfg, bundle)
                self.assertEqual(manifest["selected_global_indices"], [3, 4, 5])
                self.assertEqual(manifest["feature_scope"], "phage_only")
            # Freeze both manifests before any scoring.
            hashes = [hashlib.sha256((p / "manifest.json").read_bytes()).hexdigest() for p in bundles]
            scores = learning.score_bundle(X, spec["features"], bundles[0], expected_manifest_sha256=hashes[0])
            for matrix, bundle, expected in [(X, bundles[0], hashes[0]),
                                             (perturbed, bundles[0], hashes[0]), (X, bundles[1], hashes[1])]:
                actual = learning.score_bundle(matrix, spec["features"], bundle, expected_manifest_sha256=expected)
                self.assertTrue(CHECK.close(scores, actual, spec["absolute_score_comparison_tolerance"]))
            self.assertGreater(max(scores) - min(scores), spec["absolute_score_comparison_tolerance"])
            with self.assertRaisesRegex(ValueError, "ordered feature schema"):
                learning.score_bundle([row[3:] for row in X], spec["features"][3:], bundles[0],
                                      expected_manifest_sha256=hashes[0])

    def test_real_fitting_and_isolation_controls_on_synthetic_data(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = CHECK.run_check(Path(temporary) / "synthetic-check")
        self.assertEqual(result["status"], "passed")
        self.assertTrue(all(result["checks"].values()), result["checks"])
        self.assertFalse(result["biological_training_executed"])
        self.assertFalse(result["target_infection_scores_generated"])
        self.assertEqual(result["new_method_covered_cases"], 0)
