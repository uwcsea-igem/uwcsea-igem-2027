from dataclasses import FrozenInstanceError, replace
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tipcraft.learning import (FitSettings, Sample, fit_bundle, score_bundle,
                               training_weights, validate_partition)


def sample(row="r1", host="h1", phage="p1", hg="hg1", pg="pg1", publication="paper1"):
    return Sample(row, host, phage, hg, pg, publication)


def settings(**changes):
    return replace(FitSettings(rfe_iterations=5, final_iterations=7, depth=2,
                              learning_rate=0.1, l2_leaf_reg=3.0, n_selected_features=1,
                              rfe_step=1, seed=42, threads=1), **changes)


class PartitionTests(unittest.TestCase):
    def test_metadata_and_settings_are_immutable(self):
        with self.assertRaises(FrozenInstanceError):
            sample().host_id = "changed"
        with self.assertRaises(FrozenInstanceError):
            settings().seed = 0

    def test_both_unseen_blocks_ids_and_transitive_relatedness_groups(self):
        train = [sample()]
        valid = [sample("r2", "h2", "p2", "hg2", "pg2")]
        result = validate_partition(train, valid)
        self.assertFalse(result["random_fallback_used"])
        self.assertEqual(result["publication_overlap_disclosed"], ["paper1"])
        for test in [sample("r2", "h2", "p2", "hg1", "pg2"),
                     sample("r2", "h2", "p2", "hg2", "pg1"),
                     sample("r2", "h1", "p2", "hg1", "pg2"),
                     sample("r2", "h2", "p1", "hg2", "pg1")]:
            with self.subTest(test=test), self.assertRaises(ValueError):
                validate_partition(train, [test])

    def test_single_axis_regimes_disclose_allowed_overlap(self):
        train = [sample()]
        cold_host = [sample("r2", "h2", "p1", "hg2", "pg1")]
        cold_phage = [sample("r2", "h1", "p2", "hg1", "pg2")]
        self.assertEqual(validate_partition(train, cold_host, "cold_host")["overlaps"]["phage_id"], ["p1"])
        self.assertEqual(validate_partition(train, cold_phage, "cold_phage")["overlaps"]["host_id"], ["h1"])
        for regime, test in [("cold_host", cold_phage), ("cold_phage", cold_host)]:
            with self.assertRaises(ValueError):
                validate_partition(train, test, regime)

    def test_publication_holdout_is_explicit_and_enforced(self):
        test = [sample("r2", "h2", "p2", "hg2", "pg2")]
        with self.assertRaises(ValueError):
            validate_partition([sample()], test, publication_holdout=True)
        test = [replace(test[0], publication="paper2")]
        self.assertEqual(validate_partition([sample()], test, publication_holdout=True)["overlaps"]["publication"], [])

    def test_no_fallback_for_missing_invalid_or_conflicting_metadata(self):
        bad_partitions = [([], [sample()]), ([sample()], []),
                          ([sample()], [sample("r1", "h2", "p2", "hg2", "pg2")]),
                          ([sample()], [sample("r2", "h1", "p1", "hg1", "pg1")]),
                          ([sample()], [sample("r2", "h1", "p2", "inconsistent", "pg2")]),
                          ([sample()], [sample("r2", "h2", "p2", "", "pg2")])]
        for train, test in bad_partitions:
            with self.subTest(train=train, test=test), self.assertRaises(ValueError):
                validate_partition(train, test)
        with self.assertRaises(ValueError):
            validate_partition([sample()], [sample("r2", "h2", "p2", "hg2", "pg2")], "random")


class WeightTests(unittest.TestCase):
    def test_exact_training_only_formula_and_negative_unit_weights(self):
        rows = [sample(f"r{i}", f"h{i}") for i in range(21)]
        weights = training_weights([1] + [0] * 20, rows, "per_phage_log10")
        self.assertEqual(weights[0], 2.0)  # 1+log10(20/(1+1)).
        self.assertEqual(weights[1:], [1.0] * 20)
        self.assertEqual(training_weights([1] + [0] * 20, rows, "none"), [1.0] * 21)

    def test_all_positive_and_all_negative_phages_have_unit_weights(self):
        rows = [sample("r1", "h1", "p1"), sample("r2", "h2", "p1"),
                sample("r3", "h1", "p2", pg="pg2"), sample("r4", "h2", "p2", pg="pg2")]
        self.assertEqual(training_weights([1, 1, 0, 0], rows, "per_phage_log10"), [1.0] * 4)
        with self.assertRaises(ValueError):
            training_weights([1], rows, "per_phage_log10")


class FitValidationTests(unittest.TestCase):
    def test_genotype_vectors_cannot_encode_interaction_specific_information(self):
        rows = [sample("r1", "h1", "p1", "hg1", "pg1"),
                sample("r2", "h1", "p2", "hg1", "pg2"),
                sample("r3", "h2", "p1", "hg2", "pg1"),
                sample("r4", "h2", "p2", "hg2", "pg2")]
        # Valid vectors would be [[0,0],[0,1],[1,0],[1,1]].
        host_inconsistent = [[0, 0], [1, 1], [1, 0], [1, 1]]
        phage_inconsistent = [[0, 0], [0, 1], [1, 1], [1, 1]]
        with tempfile.TemporaryDirectory() as tmp:
            for scope in ["paired", "host_only"]:
                for matrix, prefix in [(host_inconsistent, "host:"), (phage_inconsistent, "phage:")]:
                    with self.subTest(scope=scope, prefix=prefix), self.assertRaisesRegex(ValueError, f"Inconsistent {prefix}"):
                        fit_bundle(matrix, [0, 1, 1, 0], rows, ["host:f1", "phage:f2"],
                                   settings(feature_scope=scope), Path(tmp) / "bundle")
                self.assertFalse((Path(tmp) / "bundle").exists())

    def test_invalid_inputs_fail_before_optional_ml_imports_or_output_creation(self):
        rows = [sample(), sample("r2", "h2", "p2", "hg2", "pg2")]
        names = ["host:f1", "phage:f2"]
        cases = [([[1, 0]], [0, 1], rows, names, settings()),
                 ([[1, 0], [0, 1]], [1, 1], rows, names, settings()),
                 ([[1, 0], [0, float("nan")]], [0, 1], rows, names, settings()),
                 ([[1, 0], [0, 2]], [0, 1], rows, names, settings()),
                 ([[1, 0], [0]], [0, 1], rows, names, settings()),
                 ([[1, 0], [0, 1]], [0, 1], rows, ["host:f1", "host:f1"], settings()),
                 ([[1, 0], [0, 1]], [0, 1], rows, ["row_id", "host:f1"], settings()),
                 ([[1, 0], [0, 1]], [0, 1], rows, ["host:host_id", "phage:f2"], settings()),
                 ([[1, 0], [0, 1]], [0, 1], [rows[0], rows[0]], names, settings()),
                 ([[1, 0], [0, 1]], [0, 1], rows, names, settings(n_selected_features=3)),
                 ([[1, 0], [0, 1]], [0, 1], rows, names, settings(feature_scope="host_only", n_selected_features=2)),
                 ([[1, 0], [0, 1]], [0, 1], rows, names, settings(rfe_step=0)),
                 ([[1, 0], [0, 1]], [0, 1], rows, names, settings(learning_rate=float("inf"))),
                 ([[1, 0], [0, 1]], [0, 1], rows, names, settings(weight_method="lookup"))]
        import builtins
        original_import = builtins.__import__

        def guard(name, *args, **kwargs):
            if name.split(".")[0] in {"numpy", "sklearn", "catboost"}:
                raise AssertionError("Optional ML imported before validation")
            return original_import(name, *args, **kwargs)

        with tempfile.TemporaryDirectory() as temp, patch("builtins.__import__", side_effect=guard):
            out = Path(temp) / "new-bundle"
            for X, y, samples, feature_names, cfg in cases:
                with self.subTest(X=X, settings=cfg), self.assertRaises(ValueError):
                    fit_bundle(X, y, samples, feature_names, cfg, out)
                self.assertFalse(out.exists())

    def test_existing_output_dir_is_never_overwritten(self):
        rows = [sample(), sample("r2", "h2", "p2", "hg2", "pg2")]
        with tempfile.TemporaryDirectory() as temp:
            marker = Path(temp) / "marker"
            marker.write_text("retain")
            with self.assertRaises(ValueError):
                fit_bundle([[0], [1]], [0, 1], rows, ["host:f1"], settings(), Path(temp))
            self.assertEqual(marker.read_text(), "retain")

    def test_manifest_tampering_fails_against_external_hash_before_artifact_or_ml_loading(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            original = json.dumps({"schema_version": 1, "status": "complete", "artifacts": {}}).encode()
            expected = hashlib.sha256(original).hexdigest()
            path.write_bytes(original.replace(b'"complete"', b'"tampered"'))
            # No model/spec files or ML dependencies are needed to detect this.
            with self.assertRaisesRegex(ValueError, "manifest checksum"):
                score_bundle([[0]], ["host:f1"], tmp, expected_manifest_sha256=expected)
            with self.assertRaises(TypeError):
                score_bundle([[0]], ["host:f1"], tmp)


if __name__ == "__main__":
    unittest.main()
