import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tipcraft import partition


SPEC = importlib.util.spec_from_file_location("phage_fold_features", Path(__file__).resolve().parents[3] / "scripts/check-genophi-phage-fold-features.py")
FOLD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FOLD)


def fixture():
    names = [f"p{i:03d}" for i in range(93)]
    groups = {name: f"group{i % 6}" for i, name in enumerate(names)}
    assignment = partition.assign_group_folds(set(groups.values()), axis="phage", folds=3, seed=2027)
    folds = []
    for i in range(3):
        train = [n for n in names if assignment[groups[n]] != i]
        test = [n for n in names if assignment[groups[n]] == i]
        folds.append({"fold": i, "training_phages": train, "heldout_phages": test, "mapping_control_source": train[0]})
    step = {"argv": ["{mmseqs}", "fake", "{fold_dir}/input"], "stdin": None,
            "input_prefixes": ["{fold_dir}/input"], "output_prefixes": ["{fold_dir}/output"]}
    return {"retained_model_table_phages": names, "model_phage_groups": groups}, {
        "axis": "phage", "fold_count": 3, "seed": 2027, "group_assignments": assignment,
        "folds": folds, "command_templates": [step] * 8, "mmseqs_executable": "unused"}


def alignment(query, target):
    return f"{query}\t{target}\t10\t1\t0\t0\t2\t3\t0\t2\t3\n"


class PhageFoldFeatureTests(unittest.TestCase):
    def test_frozen_plan_holds_groups_out_once_without_targets_or_optimization(self):
        groups, policy = fixture()
        result = FOLD.build_plan(policy, groups, Path("new-work"))
        heldout = []
        for fold in result["folds"]:
            train, test = fold["training_phages"], fold["heldout_phages"]
            self.assertFalse(set(train) & set(test))
            self.assertFalse({groups["model_phage_groups"][n] for n in train}
                             & {groups["model_phage_groups"][n] for n in test})
            self.assertEqual(fold["mapping_control_source"], min(train))
            self.assertEqual(len(fold["commands"]), 8)
            heldout.extend(test)
        self.assertEqual(sorted(heldout), groups["retained_model_table_phages"])
        self.assertEqual(result, FOLD.build_plan(copy.deepcopy(policy), copy.deepcopy(groups), Path("new-work")))
        damaged = copy.deepcopy(policy)
        damaged["group_assignments"][next(iter(damaged["group_assignments"]))] = 99
        with self.assertRaises(ValueError):
            FOLD.build_plan(damaged, groups, Path("new-work"))

    def test_vocab_frozen_before_query_createdb_and_heldout_rare_family_stays_absent(self):
        self._run_mock_fold(control=True)

    def test_empty_control_is_failed_even_while_heldout_zero_vectors_are_valid_observations(self):
        self._run_mock_fold(control=False)

    def _run_mock_fold(self, control):
        train = {"a1": "MAA", "a2": "MAA", "rare": "MXX"}
        query = {"q": "MXX", "control_a1": "MAA", "control_rare": "MXX"}
        tmap, qmap = {"a1": "A", "a2": "B", "rare": "A"}, {"q": "H", "control_a1": "mapping_control", "control_rare": "mapping_control"}
        plan = {"fold": 0, "training_phages": ["A", "B"], "heldout_phages": ["H"],
                "mapping_control_source": "A", "training_group_count": 2, "heldout_group_count": 1,
                "commands": [{} for _ in range(8)]}
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp) / "fold-0"
            directory.mkdir()
            for name in ["training.faa", "queries.faa", "protein-map.json"]:
                (directory / name).write_text("fixture")
            seen = []

            def fake_command(step, index, directory, timeout, bindings, ledger):
                seen.append(index)
                if index == 5:
                    (directory / "clusters.tsv").write_text("a1\ta1\na1\ta2\nrare\trare\n")
                    (directory / "self.tsv").write_text(alignment("a1", "a1") + alignment("a2", "a2") + alignment("rare", "rare"))
                if index == 6:
                    frozen = json.loads((directory / "pre-query-freeze.json").read_text())
                    self.assertEqual(frozen["completed_commands"], 5)
                    self.assertFalse(frozen["query_createdb_started"])
                    self.assertEqual(frozen["vocabulary_sha256"], FOLD.sha(directory / "vocabulary.json"))
                    vocabulary = json.loads((directory / "vocabulary.json").read_text())
                    self.assertEqual(vocabulary["family_features"], {"a1": "pc_0"})
                if index == 8:
                    text = alignment("q", "rare")
                    if control:
                        text += alignment("control_a1", "a1") + alignment("control_rare", "rare")
                    (directory / "query.tsv").write_text(text)

            with patch.object(FOLD, "run_command", side_effect=fake_command), patch.object(FOLD.subprocess, "Popen", side_effect=AssertionError("No external execution in tests")):
                result = FOLD.run_fold(plan, (train, query, tmap, qmap), Path(temp), {"timeout_seconds_per_command": 600}, {})
            self.assertEqual(seen, list(range(1, 9)))
            self.assertEqual(result["status"], "passed" if control else "failed_mapping_control")
            self.assertEqual(result["features"], 1)
            heldout = next(p for p in result["proteomes"] if p["phage"] == "H")
            self.assertTrue(heldout["zero_vector"])
            self.assertEqual(heldout["assigned_fraction"], 1.0)

    def test_no_aliasing_existing_inputs_outputs_or_changed_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "input.json"
            source.write_text("source")
            expected = FOLD.sha(source)
            with self.assertRaises(ValueError):
                FOLD.guard_destinations(root / "new", source, [source])
            with self.assertRaises(ValueError):
                FOLD.guard_destinations(root / "new", root / "new" / "output.json", [source])
            source.write_text("changed")
            with self.assertRaises(ValueError):
                FOLD.verify({source: expected})


if __name__ == "__main__":
    unittest.main()
