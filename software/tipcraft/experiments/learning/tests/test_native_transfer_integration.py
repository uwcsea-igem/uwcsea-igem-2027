"""Actual final-mode fitting on synthetic source data, never biology."""
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


NATIVE = load("native_fit_integration", "scripts/fit-native-transfer.py")
FIXTURE = load("native_assembly_fixture", "software/tipcraft/tests/test_native_transfer_fitting.py")
CHECK = load("native_learning_fixture", "scripts/check-learning-fit.py")


class NativeTransferIntegrationTests(unittest.TestCase):
    def test_final_mode_saves_three_bundles_before_six_pair_scores(self):
        data = NATIVE.assemble_native(*FIXTURE.fixture())
        spec, _, _, _ = CHECK.fixture()
        base = NATIVE.learning.FitSettings(**spec["settings"])
        settings = {s: replace(base, n_selected_features=2, feature_scope=s) for s in NATIVE.FIT.SCOPES}
        with tempfile.TemporaryDirectory() as directory:
            result = NATIVE.fit_and_score_native(data, settings, Path(directory) / "final", {"synthetic": True})
            self.assertEqual(result["status"], "completed_native_transfer_scores")
            self.assertEqual(result["native_scored_pairs_per_scope"], 6)
            self.assertEqual(result["exchange_construct_pairs_abstained"], 4)
            self.assertFalse(result["label_evaluation_executed"])
            boundary = json.loads(Path(result["pre_scoring_boundary"]["path"]).read_text())
            self.assertFalse(boundary["scoring_started"])
            for scope, method in result["methods"].items():
                self.assertEqual(method["model_id"], boundary["manifest_sha256"][scope])
                self.assertEqual(sum(len(v) for v in method["scores"].values()), 6)
                self.assertTrue(all(math.isfinite(v) for hosts in method["scores"].values() for v in hosts.values()))
            host_only = result["methods"]["host_only"]["scores"]
            self.assertEqual(host_only["T7"], host_only["T3"])
            for hosts in result["methods"]["phage_only"]["scores"].values():
                self.assertEqual(len(set(hosts.values())), 1)


if __name__ == "__main__":
    unittest.main()
