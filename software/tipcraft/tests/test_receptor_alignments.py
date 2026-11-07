"""Exercise streaming HSP integrity checks without running alignments."""
import copy
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import MagicMock, patch


SPEC = importlib.util.spec_from_file_location(
    "receptor_alignments", Path(__file__).resolve().parents[3] / "scripts/align-receptor-cohort.py")
ALIGN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ALIGN)


def hsp(query, subject, reverse=False):
    return [query, subject, "0", "20", "4", "4", "1", "4", "4" if reverse else "1",
            "1" if reverse else "4", "ACGT", "ACGT", "4", "0"]


def transcript(blocks, terminal=2):
    for name, rows in blocks:
        yield "# Query: " + name + "\n"
        yield f"# {len(rows)} hits found\n"
        for row in rows:
            yield "\t".join(row) + "\n"
    if terminal is not None:
        yield f"# BLAST processed {terminal} queries\n"


class ReceptorAlignmentTests(unittest.TestCase):
    def test_all_queries_selfcontrols_and_raw_order_are_preserved(self):
        rows = [hsp("a", "a"), hsp("b", "a", True), hsp("b", "b")]
        sink = io.StringIO()
        result = ALIGN.validate_blast(transcript([("a", rows[:1]), ("b", rows[1:])]), {"a": 4, "b": 4}, sink)
        self.assertEqual(result["completed_queries"], 2)
        self.assertEqual(result["reported_hsp_rows"], 3)
        self.assertEqual(result["all_record_self_controls_passed"], ["a", "b"])
        self.assertEqual(sink.getvalue(), "".join("\t".join(row) + "\n" for row in rows))
        self.assertFalse(result["coordinate_offsets_applied"])

    def test_missing_query_terminal_or_selfcontrol_is_refused(self):
        both = [("a", [hsp("a", "a")]), ("b", [hsp("b", "b")])]
        cases = [transcript(both, None), transcript(both, 1), transcript(both[:1]),
                 transcript([("a", [hsp("a", "a")]), ("b", [hsp("b", "a")])]),
                 transcript([("a", [hsp("a", "a")]), ("b", [])])]
        for lines in cases:
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                ALIGN.validate_blast(lines, {"a": 4, "b": 4}, io.StringIO())

    def test_repeated_unknown_or_empty_query_blocks_are_refused(self):
        for lines in (transcript([("a", []), ("a", [])]), transcript([("unknown", [])]), ["# Query: \n"]):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                ALIGN.validate_blast(lines, {"a": 4, "b": 4}, io.StringIO())

    def test_declared_hsp_count_must_equal_actual_rows(self):
        text = list(transcript([("a", [hsp("a", "a")]), ("b", [hsp("b", "b")])]))
        text[1] = "# 2 hits found\n"
        with self.assertRaisesRegex(ValueError, "Incomplete query counts"):
            ALIGN.validate_blast(iter(text), {"a": 4, "b": 4}, io.StringIO())

    def test_gapped_alignment_spans_and_reverse_subjects_are_validated(self):
        ALIGN.validate_hsp(hsp("a", "b", True), {"a": 4, "b": 4})
        row = ["a", "b", "1e-10", "10", "4", "4", "1", "3", "1", "4", "A-CG", "ATCG", "3", "1"]
        ALIGN.validate_hsp(row, {"a": 4, "b": 4})
        wrong = copy.copy(row)
        wrong[-1] = "0"
        with self.assertRaisesRegex(ValueError, "alignment strings"):
            ALIGN.validate_hsp(wrong, {"a": 4, "b": 4})

    def test_ambiguous_nucleotides_are_preserved_without_redefining_identity(self):
        row = hsp("a", "a")
        row[10:12] = ["ACRN", "ACRN"]
        row[12] = "2"  # The aggregation step, not this parser, checks N/ambiguity identity semantics.
        ALIGN.validate_hsp(row, {"a": 4})
        sink = io.StringIO()
        ALIGN.validate_blast(transcript([("a", [row])], 1), {"a": 4}, sink)
        self.assertIn("ACRN\tACRN\t2", sink.getvalue())

    def test_invalid_lengths_fields_coordinates_or_alphabet_are_refused(self):
        original = hsp("a", "a")
        cases = [original[:-1]]
        for index, value in [(0, "unknown"), (2, "NaN"), (3, "-1"), (4, "5"), (6, "0"),
                             (10, "ACXG"), (11, "ACG"), (12, "5"), (13, "1")]:
            row = copy.copy(original)
            row[index] = value
            cases.append(row)
        for row in cases:
            with self.subTest(row=row), self.assertRaises(ValueError):
                ALIGN.validate_hsp(row, {"a": 4})

    def test_streaming_validation_never_collects_hsp_strings(self):
        class Sink:
            count = 0
            def write(self, text):
                self.count += 1
        def lines():
            yield "# Query: a\n"
            yield "# 10000 hits found\n"
            for _ in range(10000):
                yield "\t".join(hsp("a", "a")) + "\n"
            yield "# BLAST processed 1 queries\n"
        sink = Sink()
        result = ALIGN.validate_blast(lines(), {"a": 4}, sink)
        self.assertEqual((sink.count, result["reported_hsp_rows"]), (10000, 10000))

    def test_role_qualified_container_keys_separate_source_and_reference_names(self):
        source = [{"id": "c000r000", "source_id": "T7", "role": "external_collection", "retained_bases": 5, "source_record_id": "deposited1"},
                  {"id": "c000r001", "source_id": "T7", "role": "external_collection", "retained_bases": 3, "source_record_id": "deposited2"},
                  {"id": "t000", "source_id": "T7", "role": "target_reference", "retained_bases": 4, "accession": "NC_001604.1"}]
        records = ALIGN.record_map(source)
        self.assertEqual([r["container_key"] for r in records], ["source:T7", "source:T7", "reference:T7"])
        self.assertEqual([r["retained_bases"] for r in records], [5, 3, 4])
        self.assertEqual(records[2]["source_record_id"], "NC_001604.1")

    def test_changed_parameters_are_refused_before_any_tool_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            policy = Path(directory) / "policy.json"
            policy.write_text(json.dumps({"wrapper_sha256": ALIGN.sha(Path(ALIGN.__file__)),
                "blastn_parameters": ["-word_size", "8"], "blastn_format": "7 " + " ".join(ALIGN.FIELDS)}))
            with patch.object(ALIGN.subprocess, "check_output") as tool, self.assertRaisesRegex(ValueError, "parameters differ"):
                ALIGN.authenticate(policy)
            tool.assert_not_called()

    def test_poll_timeout_keeps_same_process_without_killing_or_restarting(self):
        with tempfile.TemporaryDirectory() as directory:
            child = MagicMock()
            child.pid = 123
            child.wait.side_effect = [subprocess.TimeoutExpired("blastn", 30), 0]
            with patch.object(ALIGN.subprocess, "Popen", return_value=child) as launch, patch("builtins.print"):
                result = ALIGN.run_command(["synthetic-blastn"], Path(directory), 2, SimpleNamespace(verify=lambda x: None), {})
            launch.assert_called_once()
            self.assertEqual(child.wait.call_count, 2)
            child.kill.assert_not_called()
            child.terminate.assert_not_called()
            self.assertEqual(result["exit_code"], 0)

    def test_aggregation_replay_compares_every_numeric_row_and_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, numeric, mapping = root / "raw.tsv", root / "numeric.tsv", root / "map.json"
            raw.write_text("".join(transcript([("a", [hsp("a", "a")]), ("b", [hsp("b", "b")])])))
            with raw.open() as source, numeric.open("w") as sink:
                completion = ALIGN.validate_blast(source, {"a": 4, "b": 4}, sink)
            record_map = {"records": [{"id": "a", "retained_bases": 4}, {"id": "b", "retained_bases": 4}]}
            mapping.write_text(json.dumps(record_map))
            binding = lambda p: {"path": str(p), "sha256": ALIGN.sha(p)}
            report = {"status": "validated_all_record_nucleotide_hsps", "record_manifest": binding(mapping),
                      "raw_hsps": binding(raw), "numeric_hsps": binding(numeric), "nucleotide_queries": 2,
                      "query_completion": completion}
            self.assertEqual(ALIGN.validate_for_aggregation(report, record_map), completion)
            rows = numeric.read_text().splitlines()
            numeric.write_text("\n".join(reversed(rows)) + "\n")
            report["numeric_hsps"] = binding(numeric)  # Even a self-consistent new hash cannot hide reordered rows.
            with self.assertRaisesRegex(ValueError, "raw row/order"):
                ALIGN.validate_for_aggregation(report, record_map)


if __name__ == "__main__":
    unittest.main()
