"""Distinct proteomes can collide after training-pattern union projection."""
import importlib.util
from pathlib import Path
import unittest

from tipcraft.kmer_features import fit_vocabulary, project_proteomes

SPEC = importlib.util.spec_from_file_location(
    "host_feature_capacity", Path(__file__).resolve().parents[3] / "scripts/check-receptor-host-features.py")
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class HostFeatureCapacityTests(unittest.TestCase):
    def test_real_input_failure_mode_can_be_reproduced_without_phenotypes(self):
        training = {"a": {"a1": "AAAAA", "a2": "AAAAA", "c1": "CCCCC", "c2": "CCCCC"},
                    "b": {"a1": "AAAAA", "a2": "AAAAA", "d1": "DDDDD", "d2": "DDDDD"}}
        vocabulary = fit_vocabulary(training)
        self.assertEqual(len(vocabulary.features), 3)
        original = vocabulary.features
        targets = {"x": {"p": "AAAAACCCCCDDDDD"},
                   "y": {"p": "AAAAA", "q": "CCCCC", "r": "DDDDD"}}
        projected = project_proteomes(vocabulary, targets)
        self.assertNotEqual(targets["x"], targets["y"])
        self.assertEqual(projected["matrix"], [[1, 1, 1], [1, 1, 1]])
        self.assertEqual(CHECK.collisions(projected), [{"sample_ids": ["x", "y"], "vector": [1, 1, 1]}])
        self.assertEqual(vocabulary.features, original)

    def test_zero_vector_is_not_confused_with_a_nonzero_collision(self):
        vocabulary = fit_vocabulary({"a": {"p": "AAAAA", "q": "AAAAA"}})
        projected = project_proteomes(vocabulary, {"present": {"p": "AAAAA"}, "absent": {"p": "CCCCC"}})
        self.assertEqual(projected["zero_vector_sample_ids"], ["absent"])
        self.assertEqual(CHECK.collisions(projected), [])


if __name__ == "__main__":
    unittest.main()
