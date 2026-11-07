"""Small source-curation fixtures; no training, target labels or ANI execution."""
from argparse import Namespace
from collections import Counter
from contextlib import ExitStack
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "training_roster", Path(__file__).resolve().parents[3] / "scripts/build-genophi-training-roster.py")
ROSTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ROSTER)
FASTA = b">public input\nACGTACGT\n"


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return path


def host(source, file_id):
    return {"kind": "host", "source_id": source, "deposit_file_id": file_id,
            "path": f"hosts/{file_id}.fasta", "expected_bytes": len(FASTA),
            "expected_md5": hashlib.md5(FASTA).hexdigest()}


def simple_context():
    rows = [{"strain": "H1", "phage": "P1", "interaction": "0"},
            {"strain": "H2", "phage": "P1", "interaction": "0"},
            {"strain": "H3", "phage": "P2", "interaction": "1"}]
    core = {"groups": {"related": {"source_members": ["H-1", "bridge", "H-2"]},
                       "other": {"source_members": ["H-3"]}},
            "training_host_groups": {"H1": "related", "H2": "related", "H3": "other"}}
    names = {"H1": {"original_source_id": "H-1"}, "H2": {"original_source_id": "H-2"},
             "H3": {"original_source_id": "H-3"}}
    labels = {"source_blank_pairs": [{"strain": "H2", "phage": "P1"}], "conflicting_nonblank_pairs": []}
    return rows, core, names, labels


class TrainingRosterTests(unittest.TestCase):
    def test_exact_original_source_correspondence_never_accepts_alias_candidates(self):
        rows, _, _, _ = simple_context()
        audit = {"crosschecked_host_correspondence": {"H1": "H-1", "H2": "H-2", "H3": "H-3"}}
        source_map = {"hosts": {r["strain"]: {"unaccepted_alias_candidates": ["H-1"]} for r in rows},
                      "phages": {"P1": {}, "P2": {}}}
        inventory = [host("H-1", 1), host("H-2", 2), host("H-3", 3)]
        mapping = ROSTER.crosswalk(rows, audit, source_map, inventory)
        self.assertEqual(mapping["H1"]["deposit_file_id"], 1)
        self.assertFalse(mapping["H1"]["stock_identity_verified"])
        inventory[0]["source_id"] = "H1"  # Punctuation-only similarity is insufficient here.
        with self.assertRaisesRegex(ValueError, "exactly and uniquely"):
            ROSTER.crosswalk(rows, audit, source_map, inventory)

    def test_full_component_exclusion_propagates_and_keeps_all_overlapping_reasons(self):
        rows, core, names, labels = simple_context()
        exclusions = ROSTER.propagate_host_exclusions(core, {"H-1": ["near_target"]})
        self.assertEqual(exclusions["excluded_source_ids"], ["H-1", "H-2", "bridge"])
        curated = ROSTER.curate_rows(rows, labels, names, core, exclusions, {"P1": "pg1", "P2": "pg2"}, {"pg1"})
        second = next(r for r in curated if r["host_id"] == "H2")
        self.assertEqual(second["exclusion_reasons"],
                         ["original_source_empty", "excluded_host_core_group", "excluded_phage_component"])
        summary = ROSTER.label_summary(curated)
        self.assertEqual((summary["excluded_candidate_rows"], summary["retained_candidate_rows"]), (2, 1))
        self.assertEqual(summary["reason_counts_nonexclusive"],
                         {"excluded_host_core_group": 2, "excluded_phage_component": 2, "original_source_empty": 1})
        self.assertEqual([r["phage_id"] for r in curated if r["retained"]], ["P2"])

    def test_conflicts_unmapped_inputs_and_duplicate_rows_cannot_be_retained(self):
        rows, core, names, labels = simple_context()
        labels["conflicting_nonblank_pairs"] = [{"strain": "H3", "phage": "P2"}]
        exclusions = ROSTER.propagate_host_exclusions(core, {})
        curated = ROSTER.curate_rows(rows, labels, names, core, exclusions, {}, set())
        self.assertFalse(any(r["retained"] for r in curated))
        self.assertIn("nonblank_label_conflict", next(r for r in curated if r["host_id"] == "H3")["exclusion_reasons"])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            ROSTER.curate_rows(rows + rows[:1], labels, names, core, exclusions, {}, set())
        core["training_host_groups"]["H1"] = "other"
        with self.assertRaisesRegex(ValueError, "membership"):
            ROSTER.curate_rows(rows, labels, names, core, exclusions, {}, set())

    def test_seed_outside_full_core_groups_is_refused(self):
        _, core, _, _ = simple_context()
        with self.assertRaisesRegex(ValueError, "absent"):
            ROSTER.propagate_host_exclusions(core, {"unmapped": ["near_target"]})

    def test_host_verifier_distinguishes_missing_inputs_from_corrupt_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            genomes = Path(directory)
            first, second = host("H-1", 1), host("H-2", 2)
            path = genomes / first["path"]
            path.parent.mkdir()
            path.write_bytes(FASTA)
            prepare = ROSTER.module("prepare-genophi-genomes.py")
            hashes, missing = ROSTER.verify_hosts([first, second], genomes, prepare)
            self.assertEqual(set(hashes), {"H-1"})
            self.assertEqual(missing, ["hosts/2.fasta"])
            path.write_bytes(FASTA.replace(b"ACGT", b"ACGN"))
            with self.assertRaisesRegex(ValueError, "MD5"):
                ROSTER.verify_hosts([first, second], genomes, prepare)

    def test_changed_frozen_evidence_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "pinned.json"
            source.write_text("original")
            policy = {"inputs": [{"path": str(source), "sha256": ROSTER.sha(source)}],
                      "group_splits": {"regime": "both_unseen", "host_folds": 3, "phage_folds": 3, "seed": 2027},
                      "production_run_spec_complete": False}
            ROSTER.pinned_inputs(policy)
            source.write_text("changed")
            with self.assertRaisesRegex(ValueError, "Frozen input hash"):
                ROSTER.pinned_inputs(policy)

    def test_single_class_diagnostics_do_not_change_deterministic_partitions(self):
        rows = [{"row_id": f"h{h}p{p}", "host_id": f"h{h}", "phage_id": f"p{p}",
                 "host_group": f"hg{h}", "phage_group": f"pg{p}", "publication": "same_source",
                 "label": 0, "retained": True} for h in range(3) for p in range(3)]
        policy = {"group_splits": {"host_folds": 3, "phage_folds": 3, "seed": 2027}}
        first = ROSTER.partitions(rows, policy)
        self.assertEqual(first, ROSTER.partitions(list(reversed(rows)), policy))
        self.assertEqual(len(first["folds"]), 9)
        self.assertTrue(all(b["single_class_training"] and b["single_class_test"]
                            and b["single_class_fit_must_abstain"] for b in first["folds"]))
        self.assertEqual(Counter(i for b in first["folds"] for i in b["test_ids"]), Counter(r["row_id"] for r in rows))
        with self.assertRaisesRegex(ValueError, "Fewer relatedness groups"):
            ROSTER.partitions([r for r in rows if r["host_id"] != "h2"], policy)

    def test_partial_host_cache_writes_hold_only_and_never_calls_partitions(self):
        with tempfile.TemporaryDirectory() as directory, ExitStack() as stack:
            root = Path(directory)
            model = root / "model.csv"
            model.write_text("strain,phage,interaction\nH1,P1,0\nH2,P1,1\n")
            provenance = {"source_commit": "source", "source_tree_sha1": "tree", "reviewed_host_file_inventory_sha256": "inventory"}
            label = {"crosschecked_host_correspondence": {"H1": "H-1", "H2": "H-2"},
                     "source_blank_pairs": [], "conflicting_nonblank_pairs": [],
                     "concordant_nonblank_pairs": 2, "genophi_table_sha256": ROSTER.sha(model)}
            core = {"groups": {"g1": {"source_members": ["H-1"]}, "g2": {"source_members": ["H-2"]}},
                    "training_host_groups": {"H1": "g1", "H2": "g2"}}
            source_map = {**provenance, "genophi_table_sha256": ROSTER.sha(model), "hosts": {"H1": {}, "H2": {}}, "phages": {"P1": {}}}
            values = {"genophi-label-audit.json": label, "genophi-host-core-groups.json": core,
                      "genophi-training-source-map.json": source_map, "genophi-phage-inputs.json": {},
                      "genophi-overlap-policy.json": {"host_exclusion": {"named_target_identifiers": []}}}
            paths = [write_json(root / "evidence" / name, value) for name, value in values.items()]
            policy = write_json(root / "policy.json", {"inputs": [{"path": str(p), "sha256": ROSTER.sha(p)} for p in paths],
                "group_splits": {"regime": "both_unseen", "host_folds": 3, "phage_folds": 3, "seed": 2027},
                "production_run_spec_complete": False})
            inventory = [host("H-1", 1), host("H-2", 2)]
            genomes = root / "genomes"
            (genomes / "hosts").mkdir(parents=True)
            (genomes / "hosts/1.fasta").write_bytes(FASTA)
            prepare = ROSTER.module("prepare-genophi-genomes.py")
            stack.enter_context(patch.object(prepare, "load_inventory", return_value=(inventory, provenance)))
            def modules(name):
                return {"prepare-genophi-genomes.py": prepare,
                        "audit-genophi-labels.py": SimpleNamespace(run=lambda *a: label),
                        "group-genophi-hosts.py": SimpleNamespace(run=lambda *a: core)}[name]
            stack.enter_context(patch.object(ROSTER, "module", side_effect=modules))
            stack.enter_context(patch.object(ROSTER, "MODEL_TABLE", model))
            stack.enter_context(patch.object(ROSTER, "verify_phages", return_value=({}, [])))
            stack.enter_context(patch.object(ROSTER, "verify_phage_groups", return_value={"script_sha256": "bound",
                "groups": {"pg1": {"excluded": False}}, "model_phage_groups": {"P1": "pg1"}}))
            write_json(root / "phages.json", {})
            split = stack.enter_context(patch.object(ROSTER, "partitions", side_effect=AssertionError("No partial admission")))
            output = root / "output"
            argv = ["build-genophi-training-roster.py", "--policy", str(policy), "--metadata-cache", str(root / "metadata"),
                    "--genomes", str(genomes), "--target-genomes", str(root / "targets"),
                    "--phage-groups", str(root / "phages.json"), "--output-dir", str(output)]
            stack.enter_context(patch.object(sys, "argv", argv))
            self.assertEqual(ROSTER.main(), 1)
            split.assert_not_called()
            report = json.loads((output / "report.json").read_text())
            self.assertEqual(report["status"], "HOLD")
            self.assertEqual(report["host_inputs"]["verified"], 1)
            self.assertEqual(report["host_inputs"]["missing_paths"], ["hosts/2.fasta"])
            self.assertFalse(report["final_roster_written"])
            self.assertFalse(report["production_run_spec_complete"])
            self.assertTrue(report["candidate_row_summary"]["provisional"])
            self.assertEqual([p.name for p in output.iterdir()], ["report.json"])

    def test_host_audit_recomputes_missing_and_near_flags_and_rejects_report_tampering(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            audit_dir, genomes, targets = root / "audit", root / "genomes", root / "targets"
            audit_dir.mkdir()
            inventory = [host("H-1", 1), host("H-2", 2)]
            for item in inventory:
                path = genomes / item["path"]
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(FASTA)
            target = targets / "hosts/Target.fna"
            target.parent.mkdir(parents=True)
            target.write_bytes(FASTA)
            manifest = write_json(root / "targets.json", {"hosts": {"Target": "accession"},
                "genome_sources": {"Target": {"file_sha256": ROSTER.sha(target), "length": 8}}})
            exclusion = {"ani_percent": 99.5, "minimum_af_in_each_genome_percent": 90.0,
                         "ani_boundary_buffer_percent": 0.05, "af_boundary_buffer_percent": 0.5,
                         "named_target_identifiers": []}
            params = {"compression": 125, "marker_compression": 1000, "screen_percent": 80,
                      "minimum_reported_af_percent": 15, "threads": 4}
            policy = write_json(root / "policy.json", {"host_parameters": params, "host_exclusion": exclusion})
            refs = {str(genomes / item["path"]): item["source_id"] for item in inventory}
            queries = {str(target): "Target"}
            (audit_dir / "references.txt").write_text("\n".join(refs) + "\n")
            (audit_dir / "queries.txt").write_text(str(target) + "\n")
            matrix = audit_dir / "skani.tsv"
            raw_rows = [{"Ref_file": str(genomes / "hosts/1.fasta"), "Query_file": str(target),
                         "ANI": "99.7", "Align_fraction_ref": "95", "Align_fraction_query": "95"}]
            with matrix.open("w") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(raw_rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(raw_rows)
            executable = root / "skani"
            executable.write_bytes(b"synthetic tool, never executed")
            audit = ROSTER.module("audit-genophi-host-overlap.py")
            hashes = {item["source_id"]: ROSTER.sha(genomes / item["path"]) for item in inventory}
            command = [str(executable), "dist", "--rl", str(audit_dir / "references.txt"),
                       "--ql", str(audit_dir / "queries.txt"), "-c", "125", "-m", "1000", "-s", "80",
                       "--min-af", "15", "-t", "4", "-o", str(matrix)]
            report = {"method": "skani 0.4.0", "executable_sha256": ROSTER.sha(executable),
                      "wrapper_sha256": ROSTER.sha(Path(audit.__file__)), "policy_sha256": ROSTER.sha(policy),
                      "target_manifest_sha256": ROSTER.sha(manifest), "source_provenance": {},
                      "input_file_sha256": dict(hashes), "target_file_sha256": {"Target": ROSTER.sha(target)},
                      "matrix_sha256": ROSTER.sha(matrix), "command": command,
                      "no_phenotype_labels_loaded": True, "infection_or_assembly_scores_produced": False,
                      **audit.summarize(raw_rows, refs, queries, exclusion)}
            write_json(audit_dir / "report.json", report)
            with patch.object(ROSTER, "SKANI_SHA256", ROSTER.sha(executable)):
                seeds, bindings = ROSTER.verify_host_audit(audit_dir, inventory, {}, hashes, genomes, targets, manifest, policy)
                self.assertIn("Target:near_target_exclude", seeds["H-1"])
                self.assertIn("Target:unreported_exclude_pending_review", seeds["H-2"])
                self.assertEqual((bindings["reported_pairs"], bindings["expected_pairs"]), (1, 2))
                report["hosts"]["H-2"]["exclude_from_training"] = False
                write_json(audit_dir / "report.json", report)
                with self.assertRaisesRegex(ValueError, "recomputed raw matrix"):
                    ROSTER.verify_host_audit(audit_dir, inventory, {}, hashes, genomes, targets, manifest, policy)
                report["input_file_sha256"].pop("H-2")
                write_json(audit_dir / "report.json", report)
                with self.assertRaisesRegex(ValueError, "input_file_sha256"):
                    ROSTER.verify_host_audit(audit_dir, inventory, {}, hashes, genomes, targets, manifest, policy)


if __name__ == "__main__":
    unittest.main()
