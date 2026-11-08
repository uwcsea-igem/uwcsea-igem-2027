"""Synthetic membership/planning guards; no external alignments or model scores."""
from copy import deepcopy
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tipcraft.learning import Sample
from tipcraft.partition import assign_group_folds, build_blocked_folds

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("secondary_features", ROOT / "scripts/prepare-secondary-features.py")
FEATURES = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FEATURES)
PARSER_SPEC = importlib.util.spec_from_file_location("strict_secondary_protein_parser", ROOT / "scripts/check-genophi-phage-features.py")
PARSER = importlib.util.module_from_spec(PARSER_SPEC)
PARSER_SPEC.loader.exec_module(PARSER)


def panel():
    samples = [Sample(f"r{h}{p}", f"h{h}", f"p{p}", f"hg{h}", f"pg{p}", "synthetic") for h in range(3) for p in range(3)]
    hf = assign_group_folds([s.host_group for s in samples], axis="host", folds=3, seed=2027)
    pf = assign_group_folds([s.phage_group for s in samples], axis="phage", folds=3, seed=2027)
    result = {"status": "frozen_secondary_reference_panel", "annotation_join_complete": True,
              "ready_for_feature_preparation": True, "training_authorized": False,
              "samples": [asdict(s) for s in samples], "host_group_folds": hf, "phage_group_folds": pf,
              "blocks": build_blocked_folds(samples, hf, pf), "labels": {s.row_id: 0 for s in samples}}
    for axis, prefix in [("host", "h"), ("phage", "p")]:
        result[axis + "_genomes"] = {f"{prefix}{i}": {"group": f"{prefix}g{i}", "excluded": False, "annotation": {}} for i in range(3)}
    return result


class SecondaryFeatureTests(unittest.TestCase):
    def test_zero_mapping_control_requires_complete_family_and_vector_equality(self):
        records = [{"role": "mapping_control", "protein": "copy", "copied_training_protein": "original"}]
        expected = {"original": {"family": "singleton-not-retained"}}
        actual = {"copy": {"family": "singleton-not-retained"}}
        result = FEATURES.check_control(records, expected, actual, [0], [0])
        self.assertTrue(result["passed"])
        self.assertFalse(result["nonzero_vector"])
        self.assertFalse(FEATURES.check_control(records, expected, {}, [0], [0])["passed"])
        self.assertFalse(FEATURES.check_control(records, expected, {"copy": {"family": "other"}}, [0], [0])["passed"])
        self.assertFalse(FEATURES.check_control(records, expected, actual, [0], [1])["passed"])

    def test_selection_is_exact_block_axis_membership_and_ignores_labels(self):
        source = panel()
        for axis in ["host", "phage"]:
            selected = FEATURES.selection(source, "host0-phage0", axis)
            self.assertEqual(len(selected["training_genomes"]), 2)
            self.assertEqual(len(selected["query_genomes"]), 1)
            self.assertFalse(set(selected["training_genomes"]) & set(selected["query_genomes"]))
            changed = deepcopy(source)
            changed["labels"] = {r: 1 for r in source["labels"]}
            self.assertEqual(selected, FEATURES.selection(changed, "host0-phage0", axis))

    def test_modified_roles_exclusions_groups_and_unknown_modes_fail(self):
        for mutate in [lambda p: p["blocks"][0]["train_ids"].append(p["blocks"][0]["embargo_ids"][0]),
                       lambda p: p.update(status="HOLD"),
                       lambda p: p.update(training_authorized=True),
                       lambda p: p["host_genomes"]["h0"].update(excluded=True),
                       lambda p: p["host_genomes"]["h0"].update(group="wrong")]:
            source = panel()
            mutate(source)
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                FEATURES.selection(source, "host0-phage0", "host")
        for block, axis in [("final-source", "host"), ("host0-phage0", "paired")]:
            with self.assertRaises(ValueError):
                FEATURES.selection(panel(), block, axis)

    def test_eight_command_plan_keeps_query_db_out_of_training(self):
        policy = json.loads((ROOT / "software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-fold-feature-policy.json").read_text())
        commands = FEATURES.command_plan(policy, Path("new-synthetic-workdir"))
        self.assertEqual(len(commands), 8)
        self.assertEqual(commands[0]["stdin"], "new-synthetic-workdir/training.faa")
        self.assertEqual(commands[5]["stdin"], "new-synthetic-workdir/queries.faa")
        self.assertTrue(all("query_db" not in " ".join(c["argv"]) for c in commands[:5]))
        for index in [1, 3, 6]:
            argv = commands[index]["argv"]
            self.assertEqual(argv[argv.index("-c")+1], "0.8")
            self.assertEqual(argv[argv.index("--min-seq-id")+1], "0.4")
            self.assertEqual(argv[argv.index("-s")+1], "7.5")

    def test_normalization_authenticates_sources_and_preserves_copies_and_neutral_ids(self):
        source = panel()
        selected = FEATURES.selection(source, "host0-phage0", "host")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for name in selected["training_genomes"] + selected["query_genomes"]:
                path = directory / (name + ".faa")
                path.write_text(">reused_id description\nACDX*\n>copy_id\nACDX\n")
                binding = {"path": str(path), "sha256": FEATURES.sha(path)}
                source["host_genomes"][name]["annotation"] = {k: binding for k in ["FAA", "GFF", "normalized_input", "contig_map", "protein_coordinate_inventory"]} | {"proteins": 2}
            train, query, records, bindings = FEATURES.normalize(source, selected, PARSER.proteins)
            self.assertEqual(len(train), 4)
            self.assertEqual(len(query), 4)  # Two held-out proteins and two cloned control proteins.
            self.assertEqual(len(set(r["protein"] for r in records)), 8)
            self.assertTrue(all(s == "ACDX" for s in list(train.values()) + list(query.values())))
            self.assertFalse(set(train) & set(query))
            self.assertEqual(sum(r["role"] == "mapping_control" for r in records), 2)
            first = selected["training_genomes"][0]
            Path(source["host_genomes"][first]["annotation"]["FAA"]["path"]).write_text(">changed\nAAA\n")
            with self.assertRaisesRegex(ValueError, "Authenticated"):
                FEATURES.normalize(source, selected, PARSER.proteins)

    def test_missing_frozen_panel_starts_no_tools_or_outputs(self):
        # Mock only policy authentication: a missing future panel is a prerequisite HOLD.
        template = ROOT / "software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-fold-feature-policy.json"
        old = json.loads(template.read_text())
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            policy = {"helper_sha256": FEATURES.sha(FEATURES.__file__), "workers": 1,
                      "timeout_seconds_per_command": 1800, "mmseqs_executable": old["mmseqs_executable"],
                      "mmseqs_sha256": old["mmseqs_sha256"], "command_templates": old["command_templates"],
                      "inputs": {"template_policy": {"path": str(template), "sha256": FEATURES.sha(template)},
                                 "kernel": {"path": FEATURES.protein_features.__file__, "sha256": FEATURES.sha(FEATURES.protein_features.__file__)},
                                 "partition": {"path": FEATURES.partition.__file__, "sha256": FEATURES.sha(FEATURES.partition.__file__)}}}
            path = directory / "policy.json"
            path.write_text(json.dumps(policy))
            with patch.object(FEATURES, "check", side_effect=lambda b: ROOT / b["path"]), patch.object(FEATURES, "load") as loader:
                result = FEATURES.execute(path, directory / "missing-panel.json", None, "host0-phage0", "host", directory / "raw", directory / "report.json")
            self.assertEqual(result["status"], "HOLD_missing_frozen_panel")
            self.assertFalse(result["execution_started"])
            loader.assert_not_called()
            self.assertFalse((directory / "raw").exists())
            self.assertFalse((directory / "report.json").exists())

    def test_mocked_orchestration_freezes_vocabulary_before_query_and_keeps_zero_rows(self):
        # Stand-ins produce tiny synthetic HSPs, never execute MMseqs or generate model scores.
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = panel()
            selected = FEATURES.selection(source, "host0-phage0", "host")
            for name in selected["training_genomes"] + selected["query_genomes"]:
                path = directory / (name + ".faa")
                path.write_text(">one\nACDX*\n>copy\nACDX*\n")
                binding = {"path": str(path), "sha256": FEATURES.sha(path)}
                source["host_genomes"][name]["annotation"] = {k: binding for k in ["FAA", "GFF", "normalized_input", "contig_map", "protein_coordinate_inventory"]} | {"proteins": 2}
            bind = lambda p: {"path": str(p), "sha256": FEATURES.sha(p)}
            dummy = directory / "dummy-bound-producer"
            dummy.write_text("synthetic fixture; never executed")
            source.update(helper_sha256=FEATURES.sha(dummy), policy_sha256=FEATURES.sha(dummy), inputs={})
            panel_path = directory / "panel.json"
            panel_path.write_bytes(FEATURES.encoded(source))
            old = json.loads((ROOT / "software/tipcraft/benchmarks/retro-t7-t3/genophi-phage-fold-feature-policy.json").read_text())
            old.update(mmseqs_executable=str(dummy), mmseqs_sha256=FEATURES.sha(dummy))
            template = directory / "templates.json"
            template.write_bytes(FEATURES.encoded(old))
            policy = {"helper_sha256": FEATURES.sha(FEATURES.__file__), "workers": 1, "timeout_seconds_per_command": 1800,
                      "mmseqs_executable": str(dummy), "mmseqs_sha256": FEATURES.sha(dummy), "command_templates": old["command_templates"],
                      "inputs": {"template_policy": bind(template), "kernel": bind(Path(FEATURES.protein_features.__file__)),
                                 "partition": bind(Path(FEATURES.partition.__file__)), "panel_helper": bind(dummy), "panel_policy": bind(dummy),
                                 "ledger_helper": bind(ROOT / "scripts/check-genophi-phage-fold-features.py"),
                                 "protein_parser": bind(ROOT / "scripts/check-genophi-phage-features.py")}}
            policy_path = directory / "policy.json"
            policy_path.write_bytes(FEATURES.encoded(policy))
            original_loader = FEATURES.load
            events = []

            def fake_command(step, index, workdir, timeout, bindings, ledger):
                self.assertEqual(timeout, 1800)
                self.assertTrue((workdir / "frozen-plan.json").exists())
                train = dict(PARSER.proteins((workdir / "training.faa").read_bytes()))
                representative = next(iter(train))
                if index == 6:
                    freeze = json.loads((workdir / "pre-query-freeze.json").read_text())
                    self.assertEqual(freeze["completed_commands"], 5)
                    self.assertFalse(freeze["query_createdb_started"])
                    self.assertEqual(freeze["vocabulary"]["sha256"], FEATURES.sha(workdir / "vocabulary.json"))
                for prefix in step["output_prefixes"]:
                    Path(prefix).write_text("synthetic-command-artifact\n")
                if index == 3:
                    (workdir / "clusters.tsv").write_text("".join(f"{representative}\t{p}\n" for p in train))
                if index in [5, 8]:
                    proteins = train if index == 5 else dict(PARSER.proteins((workdir / "queries.faa").read_bytes()))
                    rows = [f"{p}\t{representative}\t10\t1\t0\t0\t3\t4\t0\t3\t4\n"
                            for p in proteins if index == 5 or p.startswith("control_")]
                    (workdir / ("self.tsv" if index == 5 else "query.tsv")).write_text("".join(rows))
                events.append(index)
                ledger.append({"command": step["argv"], "exit_code": 0, "synthetic_test_standin": True})

            def loader(name, path):
                module = original_loader(name, path)
                if name == "pinned_secondary_feature_ledger":
                    module.run_command = fake_command
                return module

            with patch.object(FEATURES, "load", side_effect=loader):
                planned = FEATURES.execute(policy_path, panel_path, FEATURES.sha(panel_path), "host0-phage0", "host",
                                           directory / "raw", directory / "report.json", plan_only=True)
                self.assertEqual(planned["status"], "planned_only")
                self.assertFalse(planned["execution_started"])
                self.assertFalse((directory / "raw").exists())
                self.assertFalse((directory / "report.json").exists())
                self.assertEqual(events, [])
                result = FEATURES.execute(policy_path, panel_path, FEATURES.sha(panel_path), "host0-phage0", "host",
                                          directory / "raw", directory / "report.json")
            self.assertEqual(events, list(range(1, 9)))
            self.assertEqual(result["status"], "completed_secondary_axis_features")
            self.assertTrue(result["mapping_control_passed"])
            self.assertEqual(set(result["vectors"]), set(selected["training_genomes"] + selected["query_genomes"]))
            self.assertEqual(result["zero_vector_genomes"], selected["query_genomes"])
            self.assertNotIn(FEATURES.CONTROL, result["vectors"])


if __name__ == "__main__":
    unittest.main()
