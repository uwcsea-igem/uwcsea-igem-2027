"""Adversarial completion tests for replay only; never rerun an alignment."""
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"


def module(filename):
    spec = importlib.util.spec_from_file_location(filename, SCRIPTS / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


V1 = module("align-receptor-cohort.py")
REPLAY = module("replay-receptor-alignments.py")


def hsp(query, subject=None):
    return "\t".join([query, subject or query, "0", "20", "4", "4", "1", "4", "1", "4", "ACGT", "ACGT", "4", "0"]) + "\n"


def batch(names, terminal=None):
    result = []
    for name in names:
        result += ["# BLASTN 2.17.0+\n", "# Query: " + name + "\n", "# 1 hits found\n", hsp(name)]
    result += [f"# BLAST processed {len(names) if terminal is None else terminal} queries\n"]
    return result


class ReceptorReplayTests(unittest.TestCase):
    def test_original_67_plus_32_batch_structure_preserves_all_hsp_rows(self):
        names = [f"q{i:03d}" for i in range(99)]
        lines = batch(names[:67]) + batch(names[67:])
        sink = io.StringIO()
        result = REPLAY.validate_batches(iter(lines), dict.fromkeys(names, 4), sink, V1.validate_hsp)
        self.assertEqual(result["completed_batch_query_counts"], [67, 32])
        self.assertEqual(result["completed_queries"], 99)
        self.assertEqual(sink.getvalue(), "".join(hsp(n) for n in names))
        self.assertEqual(len(result["all_record_self_controls_passed"]), 99)

    def test_full270_global_coverage_can_span_multiple_variable_batches(self):
        names = [f"q{i:03d}" for i in range(270)]
        lines = batch(names[:67]) + batch(names[67:134]) + batch(names[134:201]) + batch(names[201:])
        result = REPLAY.validate_batches(iter(lines), dict.fromkeys(names, 4), io.StringIO(), V1.validate_hsp)
        self.assertEqual(result["completed_batch_query_counts"], [67, 67, 67, 69])
        self.assertEqual(result["reported_hsp_rows"], 270)

    def test_single_batch_remains_valid_without_any_hsp_semantic_change(self):
        lines = batch(["a", "b"])
        old, new = io.StringIO(), io.StringIO()
        first = V1.validate_blast(iter(lines), {"a": 4, "b": 4}, old)
        second = REPLAY.validate_batches(iter(lines), {"a": 4, "b": 4}, new, V1.validate_hsp)
        self.assertEqual(old.getvalue(), new.getvalue())
        self.assertEqual(first["reported_rows_by_query"], second["reported_rows_by_query"])

    def test_globally_repeated_query_cannot_be_hidden_in_a_new_batch(self):
        with self.assertRaisesRegex(ValueError, "globally repeated"):
            REPLAY.validate_batches(iter(batch(["a"]) + batch(["a"])), {"a": 4, "b": 4}, io.StringIO(), V1.validate_hsp)

    def test_per_batch_count_must_match_even_if_global_declared_sum_is_correct(self):
        lines = batch(["a"], terminal=2) + batch(["b"], terminal=0)
        with self.assertRaisesRegex(ValueError, "Batch completion"):
            REPLAY.validate_batches(iter(lines), {"a": 4, "b": 4}, io.StringIO(), V1.validate_hsp)

    def test_repeated_terminal_or_empty_batch_is_not_accepted(self):
        for lines in (batch(["a"]) + ["# BLAST processed 0 queries\n"], ["# BLAST processed 0 queries\n"]):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                REPLAY.validate_batches(iter(lines), {"a": 4}, io.StringIO(), V1.validate_hsp)

    def test_terminal_must_follow_last_batch_and_be_last_nonblank_record(self):
        for lines in (batch(["a"]) + batch(["b"])[:-1], batch(["a"]) + ["# trailing uncompleted header\n"]):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                REPLAY.validate_batches(iter(lines), {"a": 4, "b": 4}, io.StringIO(), V1.validate_hsp)
        result = REPLAY.validate_batches(iter(batch(["a"]) + ["\n", "\n"]), {"a": 4}, io.StringIO(), V1.validate_hsp)
        self.assertEqual(result["completed_queries"], 1)

    def test_misreported_hsp_count_or_missing_declared_hit_record_is_rejected(self):
        first = batch(["a"])
        first[2] = "# 2 hits found\n"
        second = batch(["a"])
        del second[2]
        for lines in (first, second):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                REPLAY.validate_batches(iter(lines), {"a": 4}, io.StringIO(), V1.validate_hsp)

    def test_data_after_terminal_and_unknown_queries_are_not_repaired(self):
        for lines in (batch(["a"]) + [hsp("a")], batch(["unknown"])):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                REPLAY.validate_batches(iter(lines), {"a": 4}, io.StringIO(), V1.validate_hsp)

    def test_missing_self_control_is_rejected_despite_complete_query_batches(self):
        lines = batch(["a"]) + batch(["b"])
        lines[-2] = hsp("b", "a")
        with self.assertRaisesRegex(ValueError, "self controls"):
            REPLAY.validate_batches(iter(lines), {"a": 4, "b": 4}, io.StringIO(), V1.validate_hsp)

    def test_original_hsp_validator_still_rejects_invalid_coordinates(self):
        lines = batch(["a"])
        fields = lines[3].strip().split("\t")
        fields[6] = "0"
        lines[3] = "\t".join(fields) + "\n"
        with self.assertRaisesRegex(ValueError, "coordinates"):
            REPLAY.validate_batches(iter(lines), {"a": 4}, io.StringIO(), V1.validate_hsp)

    def test_read_only_aggregation_replay_rejects_numeric_reordering(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, numeric, mapping_file = root / "raw.tsv", root / "numeric.tsv", root / "map.json"
            raw.write_text("".join(batch(["a"]) + batch(["b"])))
            with raw.open() as source, numeric.open("w") as target:
                completion = REPLAY.validate_batches(source, {"a": 4, "b": 4}, target, V1.validate_hsp)
            mapping = {"records": [{"id": "a", "retained_bases": 4}, {"id": "b", "retained_bases": 4}]}
            mapping_file.write_text(json.dumps(mapping))
            binding = lambda p: {"path": str(p), "sha256": REPLAY.sha(p)}
            report = {"status": "validated_all_record_nucleotide_hsps", "search_rerun": False,
                      "original_execution": {"helper_snapshot": binding(Path(V1.__file__))},
                      "record_manifest": binding(mapping_file), "raw_hsps": binding(raw), "numeric_hsps": binding(numeric),
                      "nucleotide_queries": 2, "query_completion": completion}
            self.assertEqual(REPLAY.validate_for_aggregation(report, mapping), completion)
            numeric.write_text("".join(reversed(numeric.read_text().splitlines(keepends=True))))
            report["numeric_hsps"] = binding(numeric)
            with self.assertRaisesRegex(ValueError, "raw row/order"):
                REPLAY.validate_for_aggregation(report, mapping)

    def test_full_protocol_replay_preserves_failed_receipt_and_never_runs_search(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            names = [f"q{i:03d}" for i in range(270)]
            raw = root / "blastn.tsv"
            raw.write_text("".join(batch(names[:67]) + batch(names[67:])))
            partial = root / "old-partial.tsv"
            partial.write_text("original failed-parser artifact\n")
            mapping = root / "record-manifest.json"
            mapping.write_text(json.dumps({"records": [{"id": name, "retained_bases": 4} for name in names]}))
            subjects = root / "subjects.fna"
            subjects.write_text("".join(f">{n}\nACGT\n" for n in names))
            original_policy = root / "original-policy.json"
            original_policy.write_text('{"synthetic_protocol_fixture":true}')
            binding = lambda p: {"path": str(p), "sha256": REPLAY.sha(p)}
            commands = [{"argv": ["not-executed-makeblastdb"], "exit_code": 0},
                        {"argv": ["not-executed-blastn"], "exit_code": 0}]
            plan = root / "frozen-plan.json"
            plan.write_text(json.dumps({"commands": [c["argv"] for c in commands],
                "policy_sha256": REPLAY.sha(original_policy), "wrapper_sha256": REPLAY.sha(Path(V1.__file__)),
                "record_manifest_sha256": REPLAY.sha(mapping), "subjects_sha256": REPLAY.sha(subjects)}))
            original = root / "original-report.json"
            original.write_text(json.dumps({"status": "failed_or_incomplete_hsp_audit",
                "error": "Query appears after terminal completion", "database_records_verified": 270,
                "policy_sha256": REPLAY.sha(original_policy), "wrapper_sha256": REPLAY.sha(Path(V1.__file__)),
                "commands": commands, "frozen_plan": binding(plan), "record_manifest": binding(mapping),
                "subjects_fasta": binding(subjects), "authenticated_inputs": [binding(mapping), binding(subjects)],
                "raw_artifacts": [binding(raw), binding(partial)]}))
            original_bytes = original.read_bytes()
            policy = root / "replay-policy.json"
            policy.write_text(json.dumps({"replay_helper_sha256": REPLAY.sha(Path(REPLAY.__file__)),
                "search_rerun": False, "purpose": "synthetic protocol replay", "original_report": str(original),
                "original_helper": str(Path(V1.__file__)), "original_policy": str(original_policy),
                "inputs": [binding(p) for p in (original, original_policy, Path(V1.__file__), raw, mapping, plan, subjects)]}))
            with patch("subprocess.Popen", side_effect=AssertionError("No executable may run during protocol replay")):
                result = REPLAY.execute(policy, root / "replay", root / "replay-report.json")
            self.assertEqual(result["status"], "validated_all_record_nucleotide_hsps")
            self.assertEqual(result["query_completion"]["completed_batch_query_counts"], [67, 203])
            self.assertFalse(result["search_rerun"])
            self.assertEqual(original.read_bytes(), original_bytes)
            self.assertEqual(partial.read_text(), "original failed-parser artifact\n")
            self.assertEqual(REPLAY.sha(Path(result["original_execution"]["helper_snapshot"]["path"])), REPLAY.sha(Path(V1.__file__)))


if __name__ == "__main__":
    unittest.main()
