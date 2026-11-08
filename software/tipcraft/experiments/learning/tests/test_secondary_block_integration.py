"""Actual three-scope block fitting on a synthetic fixture, never biology."""
from dataclasses import replace
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[5]


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FIT = load("secondary_block_integration", "scripts/fit-secondary-block.py")
CHECK = load("secondary_block_fixture", "scripts/check-learning-fit.py")


class SecondaryBlockIntegrationTests(unittest.TestCase):
    def test_three_frozen_scopes_score_same_unseen_rows(self):
        spec, X, y, train = CHECK.fixture()
        test = [FIT.learning.Sample(f"test{i}{j}", f"newh{i}", f"newp{j}", f"newhg{i}", f"newpg{j}", "synthetic")
                for i in range(2) for j in range(2)]
        host = [[0, 1, 1], [1, 0, 1]]
        phage = [[0, 0, 1], [1, 1, 0]]
        data = {"block_id": "synthetic", "train": train, "test": test, "X_train": X, "y_train": y,
                "X_test": [h + p for h in host for p in phage], "feature_names": spec["features"],
                "partition": FIT.learning.validate_partition(train, test)}
        base = FIT.learning.FitSettings(**spec["settings"])
        settings = {s: replace(base, feature_scope=s) for s in FIT.SCOPES}
        with tempfile.TemporaryDirectory() as directory:
            report = FIT.fit_and_score(data, settings, Path(directory) / "block", {"synthetic": True})
            self.assertEqual(report["status"], "completed_secondary_block_scores")
            self.assertFalse(report["test_labels_used_for_fit_or_selection"])
            boundary = json.loads(Path(report["pre_scoring_boundary"]["path"]).read_text())
            self.assertFalse(boundary["scoring_started"])
            for name, predictions in report["methods"].items():
                self.assertEqual(set(predictions["scores"]), {s.row_id for s in test})
                self.assertTrue(all(math.isfinite(v) for v in predictions["scores"].values()))
                self.assertEqual(predictions["model_id"], boundary["manifest_sha256"][name])
                self.assertEqual(FIT.sha(report["bundles"][name]["path"]), predictions["model_id"])
            h = report["methods"]["host_only"]["scores"]
            p = report["methods"]["phage_only"]["scores"]
            self.assertEqual(h["test00"], h["test01"])
            self.assertEqual(h["test10"], h["test11"])
            self.assertEqual(p["test00"], p["test10"])
            self.assertEqual(p["test01"], p["test11"])


if __name__ == "__main__":
    unittest.main()
