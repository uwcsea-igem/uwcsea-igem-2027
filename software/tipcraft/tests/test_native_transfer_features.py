"""Synthetic metadata and external-command stand-ins only; no biological scores."""
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


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FEATURES = load("native_transfer_features", "scripts/prepare-native-transfer-features.py")
SECONDARY = load("native_transfer_source_utilities_test", "scripts/prepare-secondary-features.py")
PARSER = load("native_transfer_strict_parser_test", "scripts/check-genophi-phage-features.py")


def metadata():
    samples = [Sample(f"r{h}{p}", f"h{h}", f"p{p}", f"hg{h}", f"pg{p}", "synthetic") for h in range(3) for p in range(3)]
    hf = assign_group_folds([s.host_group for s in samples], axis="host", folds=3, seed=2027)
    pf = assign_group_folds([s.phage_group for s in samples], axis="phage", folds=3, seed=2027)
    panel = {"status": "frozen_secondary_reference_panel", "annotation_join_complete": True, "ready_for_feature_preparation": True,
             "training_authorized": False, "samples": [asdict(s) for s in samples], "host_group_folds": hf, "phage_group_folds": pf,
             "blocks": build_blocked_folds(samples, hf, pf), "labels": {s.row_id: 0 for s in samples}}
    for axis in ["host", "phage"]:
        prefix = axis[0]
        panel[axis + "_genomes"] = {f"{prefix}{i}": {"group": f"{prefix}g{i}", "excluded": False, "annotation": {}} for i in range(3)}
    design = {"execute_if_source_cv_disappoints": True, "native_candidates": {"hosts": FEATURES.NATIVE["host"], "phages": FEATURES.NATIVE["phage"]},
              "native_query_proteomes": [{"kind": kind, "source_id": name, "query_only": True} for kind, names in FEATURES.NATIVE.items() for name in names]}
    return panel, design


def fixture(directory):
    panel, design = metadata()
    bind = lambda path: {"path": str(path), "sha256": FEATURES.sha(path)}
    for kind in ["host", "phage"]:
        for name, genome in panel[kind + "_genomes"].items():
            path = directory / (name + ".faa")
            path.write_text(">raw_one\nACDX*\n>raw_copy\nACDX\n")
            genome["annotation"] = {k: bind(path) for k in ["FAA", "GFF", "normalized_input", "contig_map", "protein_coordinate_inventory"]} | {"proteins": 2}
    for row in design["native_query_proteomes"]:
        path = directory / (row["source_id"] + ".faa")
        path.write_text(">native_one\nACDX*\n>native_copy\nACDX\n")
        row.update({k: bind(path) for k in ["normalized_DNA", "proteins", "GFF"]})
        row.update(predicted_proteins=2, accession="synthetic-" + row["source_id"], annotation_id="synthetic-" + row["source_id"])
    return panel, design


class NativeTransferFeatureTests(unittest.TestCase):
    def test_all_source_genomes_and_fixed_query_order_without_label_or_cv_selection(self):
        panel, design = metadata()
        for axis in ["host", "phage"]:
            selected = FEATURES.selection(panel, design, axis)
            self.assertEqual(selected["training_genomes"], sorted(panel[axis + "_genomes"]))
            self.assertEqual(selected["query_genomes"], FEATURES.NATIVE[axis])
            self.assertEqual(selected["source_training_rows"], 9)
            changed = deepcopy(panel)
            changed["labels"] = {key: 1 for key in panel["labels"]}
            changed["source_cv_status"] = "negative_result_retained"
            self.assertEqual(selected, FEATURES.selection(changed, design, axis))

    def test_malformed_roster_query_roles_and_cv_gating_fail(self):
        for mutation in [lambda p, d: d.update(execute_if_source_cv_disappoints=False),
                         lambda p, d: d["native_query_proteomes"][0].update(query_only=False),
                         lambda p, d: d["native_query_proteomes"].pop(),
                         lambda p, d: d["native_candidates"].update(phages=["T3", "T7"]),
                         lambda p, d: p["host_genomes"].update(extra={"excluded": False}),
                         lambda p, d: p["host_genomes"]["h0"].update(excluded=True),
                         lambda p, d: p["blocks"][0]["train_ids"].pop()]:
            panel, design = metadata()
            mutation(panel, design)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                FEATURES.selection(panel, design, "host")

    def test_native_normalization_uses_actual_three_bindings_and_query_only_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            panel, design = fixture(Path(temporary))
            selected = FEATURES.selection(panel, design, "phage")
            train, query, records, bindings = FEATURES.normalize_native(panel, selected, design, SECONDARY, PARSER.proteins)
            self.assertEqual(len(train), 6)
            self.assertEqual(len(query), 6)  # Two copied controls plus four native query proteins.
            self.assertFalse(set(train) & set(query))
            self.assertEqual({r["genome"] for r in records if r["role"] == "native_query"}, {"T7", "T3"})
            self.assertFalse(any(r["genome"] in {"T7", "T3"} for r in records if r["role"] == "train"))
            self.assertEqual(len(records), len(train) + len(query))
            self.assertTrue(all(seq == "ACDX" for seq in list(train.values()) + list(query.values())))
            row = next(r for r in design["native_query_proteomes"] if r["source_id"] == "T7")
            Path(row["GFF"]["path"]).write_text("corrupted")
            with self.assertRaisesRegex(ValueError, "Authenticated"):
                FEATURES.normalize_native(panel, selected, design, SECONDARY, PARSER.proteins)

    def test_mocked_eight_steps_freeze_source_vocabulary_before_native_queries(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            panel, design = fixture(directory)
            bind = lambda p: {"path": str(p), "sha256": FEATURES.sha(p)}
            dummy = directory / "dummy-bound-producer"
            dummy.write_text("synthetic stand-in; never executed")
            panel.update(helper_sha256=FEATURES.sha(dummy), policy_sha256=FEATURES.sha(dummy), inputs={})
            panel_path = directory / "panel.json"
            panel_path.write_bytes(SECONDARY.encoded(panel))
            design_path = directory / "native-design.json"
            # Outcome assets are deliberately nonexistent and must not be opened.
            design["bound_inputs"] = [{"path": "nonexistent-outcomes.csv", "sha256": "not-an-input"}]
            design_path.write_bytes(SECONDARY.encoded(design))
            source_policy = json.loads((ROOT / "software/tipcraft/benchmarks/receptor-resource/secondary-feature-policy.json").read_text())
            source_policy.update(mmseqs_executable=str(dummy), mmseqs_sha256=FEATURES.sha(dummy))
            source_policy_path = directory / "source-feature-policy.json"
            source_policy_path.write_bytes(SECONDARY.encoded(source_policy))
            policy = {"helper_sha256": FEATURES.sha(FEATURES.__file__), "workers": 1, "timeout_seconds_per_command": 1800,
                      "mmseqs_executable": str(dummy), "mmseqs_sha256": FEATURES.sha(dummy), "command_templates": source_policy["command_templates"],
                      "inputs": {"native_design": bind(design_path), "source_feature_policy": bind(source_policy_path),
                                 "source_feature_helper": bind(ROOT / "scripts/prepare-secondary-features.py"),
                                 "kernel": bind(Path(FEATURES.protein_features.__file__)), "partition": bind(Path(FEATURES.partition.__file__)),
                                 "panel_helper": bind(dummy), "panel_policy": bind(dummy),
                                 "ledger_helper": bind(ROOT / "scripts/check-genophi-phage-fold-features.py"),
                                 "protein_parser": bind(ROOT / "scripts/check-genophi-phage-features.py")}}
            policy_path = directory / "policy.json"
            policy_path.write_bytes(SECONDARY.encoded(policy))
            loader, events = FEATURES.load, []

            def command(step, index, workdir, timeout, bindings, ledger):
                self.assertEqual(timeout, 1800)
                plan = json.loads((workdir / "frozen-plan.json").read_text())
                self.assertEqual(plan["mode"], "final_source_native_query")
                train = dict(PARSER.proteins((workdir / "training.faa").read_bytes()))
                representative = next(iter(train))
                if index == 6:
                    freeze = json.loads((workdir / "pre-query-freeze.json").read_text())
                    self.assertEqual(freeze["completed_commands"], 5)
                    self.assertFalse(freeze["query_createdb_started"])
                    self.assertEqual(freeze["vocabulary"]["sha256"], FEATURES.sha(workdir / "vocabulary.json"))
                for prefix in step["output_prefixes"]:
                    Path(prefix).write_text("synthetic raw artifact")
                if index == 3:
                    (workdir / "clusters.tsv").write_text("".join(f"{representative}\t{p}\n" for p in train))
                if index in [5, 8]:
                    proteins = train if index == 5 else dict(PARSER.proteins((workdir / "queries.faa").read_bytes()))
                    (workdir / ("self.tsv" if index == 5 else "query.tsv")).write_text("".join(
                        f"{p}\t{representative}\t10\t1\t0\t0\t3\t4\t0\t3\t4\n" for p in proteins if index == 5 or p.startswith("control_")))
                events.append(index)
                ledger.append({"command": step["argv"], "exit_code": 0, "synthetic_test_standin": True})

            def patched_loader(name, path):
                module = loader(name, path)
                if name == "pinned_native_command_ledger":
                    module.run_command = command
                return module

            with patch.object(FEATURES, "load", side_effect=patched_loader):
                missing = FEATURES.execute(policy_path, directory / "missing-panel.json", None, "phage", directory / "raw", directory / "report.json", True)
                self.assertEqual(missing["status"], "HOLD_missing_frozen_panel")
                planned = FEATURES.execute(policy_path, panel_path, FEATURES.sha(panel_path), "phage", directory / "raw", directory / "report.json", True)
                self.assertEqual(planned["status"], "planned_only")
                self.assertEqual(events, [])
                self.assertFalse((directory / "raw").exists())
                result = FEATURES.execute(policy_path, panel_path, FEATURES.sha(panel_path), "phage", directory / "raw", directory / "report.json")
            self.assertEqual(events, list(range(1, 9)))
            self.assertEqual(result["status"], "completed_native_transfer_axis_features")
            self.assertEqual(result["query_genomes"], ["T7", "T3"])
            self.assertEqual(set(result["vectors"]), {"p0", "p1", "p2", "T7", "T3"})
            self.assertEqual(result["zero_vector_genomes"], ["T7", "T3"])
            self.assertTrue(result["mapping_control_passed"])
            self.assertNotIn(FEATURES.CONTROL, result["vectors"])
            self.assertFalse(result["native_outcome_assets_opened"])


if __name__ == "__main__":
    unittest.main()
