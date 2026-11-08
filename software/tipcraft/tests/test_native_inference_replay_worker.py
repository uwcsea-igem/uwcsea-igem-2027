"""Technical packet checks; no biological data or accuracy assertions."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("native_replay", ROOT / "scripts/replay-native-inference-worker.py")
WORKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WORKER)


class NativeReplayWorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.packet = {"schema_version": 1, "row_ids": [f"q{i:03d}" for i in range(6)],
                       "feature_names": ["host:f", "phage:f"],
                       "matrix": [[0, 1], [1, 0], [1, 1], [0, 0], [1, 0], [0, 1]], "bundles": {}}
        for scope in WORKER.SCOPES:
            directory = Path(self.temp.name) / scope
            directory.mkdir()
            path = directory / "manifest.json"
            path.write_text(json.dumps({"feature_scope": scope}))
            self.packet["bundles"][scope] = {"path": str(path), "sha256": WORKER.sha(path)}

    def scorer(self, matrix, names, directory, *, expected_manifest_sha256):
        self.assertEqual(names, self.packet["feature_names"])
        self.assertEqual(WORKER.sha(directory / "manifest.json"), expected_manifest_sha256)
        return [float(row[0] - row[1]) for row in matrix]

    def test_scores_remain_attached_to_packet_order(self):
        result = WORKER.score_packet(self.packet, self.scorer)
        order = [4, 0, 5, 2, 1, 3]
        changed = deepcopy(self.packet)
        changed["matrix"] = [self.packet["matrix"][i] for i in order]
        replay = WORKER.score_packet(changed, self.scorer)
        for scope in WORKER.SCOPES:
            original = list(result["methods"][scope]["scores"].values())
            self.assertEqual(list(replay["methods"][scope]["scores"].values()), [original[i] for i in order])
        self.assertFalse(result["model_refitted"])

    def test_labels_wrong_coverage_and_invalid_features_fail_before_scoring(self):
        mutations = [lambda p: p.update(labels=[0] * 6),
                     lambda p: p.update(schema_version=True),
                     lambda p: p["row_ids"].reverse(),
                     lambda p: p["matrix"].pop(),
                     lambda p: p["matrix"][0].__setitem__(0, 2),
                     lambda p: p["bundles"].pop("paired")]
        for mutate in mutations:
            packet = deepcopy(self.packet)
            mutate(packet)
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                WORKER.score_packet(packet, lambda *a, **k: self.fail("Invalid packet reached scorer"))

    def test_tampering_scope_mismatch_and_mutation_during_scoring_fail(self):
        packet = deepcopy(self.packet)
        packet["bundles"]["paired"]["sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            WORKER.score_packet(packet, self.scorer)
        packet = deepcopy(self.packet)
        packet["bundles"]["paired"] = packet["bundles"]["host_only"]
        with self.assertRaises(ValueError):
            WORKER.score_packet(packet, self.scorer)
        def mutate_manifest(matrix, names, directory, **kwargs):
            (directory / "manifest.json").write_text("changed")
            return [0.] * 6
        with self.assertRaises(ValueError):
            WORKER.score_packet(self.packet, mutate_manifest)

    def test_missing_or_nonfinite_scores_fail(self):
        for scores in [[0.] * 5, [float("nan")] * 6, [True] * 6]:
            with self.subTest(scores=scores), self.assertRaises(ValueError):
                WORKER.score_packet(self.packet, lambda *a, **k: scores)


if __name__ == "__main__":
    unittest.main()
