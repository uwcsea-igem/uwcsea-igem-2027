import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[5]
SPEC = importlib.util.spec_from_file_location("fit_check", ROOT / "scripts/check-learning-fit.py")
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class FittingIntegrationTests(unittest.TestCase):
    def test_real_fitting_and_isolation_controls_on_synthetic_data(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = CHECK.run_check(Path(temporary) / "synthetic-check")
        self.assertEqual(result["status"], "passed")
        self.assertTrue(all(result["checks"].values()), result["checks"])
        self.assertFalse(result["biological_training_executed"])
        self.assertFalse(result["target_infection_scores_generated"])
        self.assertEqual(result["new_method_covered_cases"], 0)
