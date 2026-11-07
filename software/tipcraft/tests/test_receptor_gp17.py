"""Test DNA-only near-hit triage, completion controls and failure retention."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"


def module(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / name)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


AUDIT = module("audit-receptor-gp17-overlap.py")
LOGIC = module("audit-genophi-phage-overlap.py")
THRESHOLDS = {"gp17_near_match": {"identity_percent": 90.0, "minimum_query_coverage_percent": 80.0},
              "broader_gp17_family_diagnostic": {"identity_percent": 40.0, "minimum_query_coverage_percent": 80.0}}


def entry(root, primary="P", container="c000", protein_error=True):
    raw = ("synthetic deposited container " + primary).encode()
    original = root / (container + ".gbk")
    original.write_bytes(raw)
    pairs = [(container + "r000", "ACGTRN"), (container + "r001", "TTAC")]
    content = "".join(f">{name}\n{seq}\n" for name, seq in pairs).encode()
    nucleotide = root / (container + ".fna")
    nucleotide.write_bytes(content)
    return {"source": {"primary_id": primary, "container_id": container, "cohort_member": True,
                       "published_model_role": "training", "declared_size_bytes": len(raw),
                       "git_blob_sha1": hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()},
            "container_path": str(original), "status": "annotation_error" if protein_error else "inventoried",
            "protein_status": "annotation_error" if protein_error else "available",
            "protein_errors": [{"error": "synthetic translation issue"}] if protein_error else [],
            "nucleotide_status": "available", "outputs": {"nucleotide_fasta": {"path": str(nucleotide),
                "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content), "sequences": len(pairs)}},
            "records": [{"record_id": f"original-{i}", "neutral_nucleotide_id": name,
                         "length": len(seq), "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest()}
                        for i, (name, seq) in enumerate(pairs)]}


def hit(query, subject, identity=100, start=1, end=100, frame=1):
    return [query, subject, str(identity), str(start), str(end), "100", "1", str(3 * (end-start+1)),
            str(frame), "200", "1e-30"]


def genomes():
    return [{"id": n, "source_id": source, "role": role, "retained_bases": 3000}
            for n, source, role in [("t000", "T7", "target_reference"), ("t001", "T3", "target_reference"),
                                    ("c000r000", "P", "external_collection"), ("c000r001", "P", "external_collection"),
                                    ("c001r000", "B", "external_collection"), ("c002r000", "split", "external_collection")]]


class ReceptorGp17Tests(unittest.TestCase):
    def test_multirecord_dna_is_preserved_despite_protein_annotation_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            source = entry(Path(directory))
            bindings = {}
            rows, records, sequences = AUDIT.source_inputs([source], bindings)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["nucleotide_status"], "authenticated_available")
            self.assertEqual(rows[0]["protein_status"], "annotation_error")
            self.assertEqual([g["source_id"] for g in records], ["P", "P"])
            self.assertEqual(sequences, [("c000r000", "ACGTRN"), ("c000r001", "TTAC")])
            self.assertEqual(len(bindings), 2)

    def test_corrupt_or_missing_dna_retains_failure_row_without_partial_records(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            valid, bad = entry(root), entry(root, "bad", "c001")
            Path(bad["outputs"]["nucleotide_fasta"]["path"]).write_text("corrupt")
            rows, records, _ = AUDIT.source_inputs([valid, bad], {})
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[1]["nucleotide_status"], "unassessed_input_error")
            self.assertEqual(rows[1]["nucleotide_ids"], [])
            self.assertEqual({g["source_id"] for g in records}, {"P"})
            self.assertIn("differs", rows[1]["nucleotide_error"])

    def test_changed_original_blob_and_changed_record_hash_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = entry(root)
            bad = copy.deepcopy(original)
            bad["source"]["git_blob_sha1"] = "0" * 40
            rows, records, _ = AUDIT.source_inputs([bad], {})
            self.assertEqual(records, [])
            self.assertIn("Original container", rows[0]["nucleotide_error"])
            bad = copy.deepcopy(original)
            bad["records"][1]["sequence_sha256"] = "0" * 64
            rows, records, _ = AUDIT.source_inputs([bad], {})
            self.assertEqual(records, [])
            self.assertIn("bases differ", rows[0]["nucleotide_error"])

    def test_thresholds_are_inclusive_single_alignment_and_exclude_whole_container(self):
        rows = [hit("t000_gp17", "t000"), hit("t001_gp17", "t001"),
                hit("t000_gp17", "c000r001", 90, 1, 80), hit("t001_gp17", "c000r000", 91, 1, 80),
                hit("t000_gp17", "c001r000", 40, 1, 80),
                hit("t000_gp17", "c002r000", 100, 1, 50), hit("t000_gp17", "c002r000", 100, 51, 100)]
        matches = LOGIC.gp17_matches(rows, genomes(), THRESHOLDS, {"t000_gp17": 100, "t001_gp17": 100})
        self.assertEqual(matches["near_match_excluded_source_ids"], ["P"])
        self.assertEqual({r["source_id"] for r in matches["broader_family_matches"]}, {"P", "B"})
        candidates = [{"primary_id": name, "nucleotide_status": "authenticated_available", "admitted": False}
                      for name in ("P", "B", "split")]
        result = AUDIT.mark_candidates(candidates, matches)
        self.assertEqual(result[0]["direct_gp17_disposition"], "direct_near_match_exclude")
        self.assertEqual(result[2]["direct_gp17_disposition"], "no_reported_hit_meeting_frozen_thresholds")
        self.assertTrue(all(r["whole_genome_or_group_clearance"] is False and r["admitted"] is False for r in result))

    def test_missing_original_reference_selfcontrol_and_invalid_frame_are_rejected(self):
        rows = [hit("t000_gp17", "t000")]
        with self.assertRaisesRegex(ValueError, "self controls incomplete"):
            LOGIC.gp17_matches(rows, genomes(), THRESHOLDS, {"t000_gp17": 100, "t001_gp17": 100})
        with self.assertRaisesRegex(ValueError, "frame"):
            LOGIC.gp17_matches([hit("t000_gp17", "t000", frame=0)], genomes(), THRESHOLDS,
                               {"t000_gp17": 100, "t001_gp17": 100})

    def test_exact_perquery_completion_requires_counts_and_terminal_record(self):
        rows = [hit("t000_gp17", "t000"), hit("t001_gp17", "t001")]
        text = "\n".join("# Query: " + r[0] + "\n# 1 hits found\n" + "\t".join(r) for r in rows)
        text += "\n# BLAST processed 2 queries\n"
        parsed, completion = LOGIC.completed_blast(text, {"t000_gp17", "t001_gp17"}, 11)
        self.assertEqual(parsed, rows)
        self.assertEqual(completion["completed_queries"], 2)
        for broken in (text.replace("# BLAST processed 2 queries", ""), text.replace("# 1 hits found", "# 2 hits found", 1)):
            with self.assertRaises(ValueError):
                LOGIC.completed_blast(broken, {"t000_gp17", "t001_gp17"}, 11)

    def test_database_count_requires_every_source_and_control_record(self):
        self.assertEqual(AUDIT.database_count("Adding sequences from FASTA; added 270 sequences in 1 seconds.", 270), 270)
        for text in ("", "Adding sequences from FASTA; added 269 sequences", "Adding sequences from FASTA; added 270 sequences\nAdding sequences from FASTA; added 270 sequences"):
            with self.assertRaises(ValueError):
                AUDIT.database_count(text, 270)

    def test_fasta_rejects_contig_joining_empty_or_duplicate_records(self):
        self.assertEqual(AUDIT.fasta(b">a\nAC\n>b\nGT\n", AUDIT.DNA), [("a", "AC"), ("b", "GT")])
        for text in (b">a\n>b\nAC\n", b">a\nAC\n>a\nGT\n", b"AC\n", b">a\nACZ\n"):
            with self.assertRaises(ValueError):
                AUDIT.fasta(text, AUDIT.DNA)

    def test_incomplete_input_roster_is_reported_and_never_executes_commands(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first, second = entry(root), entry(root, "missing", "c001")
            second["nucleotide_status"] = "unavailable"
            policy_path = root / "policy.json"
            policy_path.write_text("{}")
            bindings = {policy_path: AUDIT.file_sha(policy_path), Path(AUDIT.__file__): AUDIT.file_sha(Path(AUDIT.__file__))}
            policy = {"purpose": "synthetic incomplete roster", "expected_source_nucleotide_records": 4}
            targets = [{"id": "t000", "source_id": "T7"}, {"id": "t001", "source_id": "T3"}]
            with patch.object(AUDIT, "authenticated", return_value=(policy, [first, second], {"manifest.json": "{}", "reference-gp17.faa": b""}, bindings, LOGIC, {})), patch.object(
                    AUDIT, "reference_inputs", return_value=(targets, [], {"t000_gp17": 100, "t001_gp17": 100})), patch.object(AUDIT, "command") as run:
                report = AUDIT.execute(policy_path, root / "work", root / "report.json")
            run.assert_not_called()
            self.assertEqual(report["status"], "HOLD_incomplete_nucleotide_inputs")
            self.assertEqual(len(report["candidate_rows"]), 2)
            self.assertEqual(report["direct_near_hit_excluded_primary_ids"], [])
            self.assertFalse(report["query_completion_verified"])
            self.assertFalse(report["all_rows_admitted"])

    def test_outputs_must_be_new_and_separate_from_every_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            source.write_text("input")
            for work, output in ((root, root / "report.json"), (root / "work", root / "work/report.json"),
                                 (root / "work", source), (source / "child", root / "report.json")):
                with self.subTest(work=work, output=output), self.assertRaises(ValueError):
                    AUDIT.guard(work, output, {source: "hash"})


if __name__ == "__main__":
    unittest.main()
