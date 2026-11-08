"""Block assembly checks without ML dependencies or biological inputs."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("secondary_fitting", ROOT / "scripts/fit-secondary-block.py")
FIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIT)


def fixture():
    rows = [{"row_id": str(i), "host_id": f"h{i//2}", "phage_id": f"p{i%2}",
             "host_group": f"hg{i//2}", "phage_group": f"pg{i%2}", "publication": "synthetic"}
            for i in range(4)]
    rows += [{"row_id": "test", "host_id": "ht", "phage_id": "pt", "host_group": "hgt",
              "phage_group": "pgt", "publication": "synthetic"},
             {"row_id": "embargo", "host_id": "ht", "phage_id": "p0", "host_group": "hgt",
              "phage_group": "pg0", "publication": "synthetic"}]
    panel = {"status": "frozen_secondary_reference_panel", "ready_for_feature_preparation": True,
             "annotation_join_complete": True, "membership_generated": True,
             "samples": rows, "blocks": [{"fold_id": "block", "train_ids": [str(i) for i in range(4)],
                                           "test_ids": ["test"], "embargo_ids": ["embargo"]}],
             "labels": {"0": 0, "1": 1, "2": 1, "3": 0, "test": "not-read", "embargo": "not-read"}}
    reports = []
    for axis in ["host", "phage"]:
        prefix = "h" if axis == "host" else "p"
        reports.append({"status": "completed_secondary_axis_features", "axis": axis, "block_id": "block",
                        "mapping_control_passed": True, "training_genomes": [prefix + "0", prefix + "1"],
                        "query_genomes": [prefix + "t"], "feature_names": ["f0", "f1"],
                        "vectors": {prefix + "0": [0, 1], prefix + "1": [1, 0], prefix + "t": [1, 1]}})
    return panel, *reports, "block"


class SecondaryFittingTests(unittest.TestCase):
    def test_executed_learning_module_must_be_the_bound_source(self):
        policy = {"inputs": {"learning": FIT.binding(FIT.learning.__file__)}}
        FIT.verify_implementation(policy)
        with patch.object(FIT.learning, "__file__", str(ROOT / "scripts/fit-secondary-block.py")):
            with self.assertRaisesRegex(ValueError, "Imported fitting implementation"):
                FIT.verify_implementation(policy)
        policy["inputs"]["learning"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "Imported fitting implementation"):
            FIT.verify_implementation(policy)

    def test_only_training_labels_enter_fit_inputs(self):
        args = fixture()
        data = FIT.assemble_block(*args)
        self.assertEqual(data["y_train"], [0, 1, 1, 0])
        self.assertEqual(data["feature_names"], ["host:f0", "host:f1", "phage:f0", "phage:f1"])
        args[0]["labels"]["test"] = {"arbitrary": "unused"}
        args[0]["labels"]["embargo"] = float("nan")
        self.assertEqual(data, FIT.assemble_block(*args))
        self.assertNotIn("y_test", data)

    def test_incomplete_or_mismatched_features_fail(self):
        for mutate in [lambda a: a[0].update(status="HOLD"),
                       lambda a: a[1].update(mapping_control_passed=False),
                       lambda a: a[2].update(axis="host"),
                       lambda a: a[1]["training_genomes"].append("ht"),
                       lambda a: a[1]["vectors"].update(control=[1, 1]),
                       lambda a: a[2]["vectors"].update(pt=[1, 2]),
                       lambda a: a[0]["blocks"][0]["train_ids"].append("test")]:
            args = fixture()
            mutate(args)
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                FIT.assemble_block(*args)

    def test_both_unseen_group_boundary_is_required(self):
        args = fixture()
        args[0]["samples"][4]["host_group"] = "hg0"
        args[0]["samples"][5]["host_group"] = "hg0"
        with self.assertRaisesRegex(ValueError, "protected axes"):
            FIT.assemble_block(*args)

    def test_exact_fixed_scope_feature_count_rules(self):
        design = json.loads((ROOT / "software/tipcraft/benchmarks/receptor-resource/secondary-model-design-policy.json").read_text())
        names = ["host:f" + str(i) for i in range(105)] + ["phage:f" + str(i) for i in range(7)]
        result = FIT.fixed_settings(design, names)
        self.assertEqual(result["paired"].n_selected_features, 100)
        self.assertEqual(result["phage_only"].n_selected_features, 7)
        self.assertEqual(result["host_only"].rfe_step, 1)
        self.assertTrue(all(s.final_iterations == 500 and s.seed == 2027 for s in result.values()))
        with self.assertRaisesRegex(ValueError, "Empty feature"):
            FIT.fixed_settings(design, ["host:f0"])
        altered = deepcopy(design)
        altered["fit"]["early_stopping"] = True
        with self.assertRaises(ValueError):
            FIT.fixed_settings(altered, names)


if __name__ == "__main__":
    unittest.main()
