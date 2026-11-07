"""Test literal parsing and exact name correspondence without biological scoring."""
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "receptor_resource", Path(__file__).resolve().parents[3] / "scripts/audit-receptor-resource.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def table_bytes(header, rows, delimiter):
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter=delimiter, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue().encode()


def metadata_row(name, **changes):
    fields = dict.fromkeys(AUDIT.METADATA_HEADER, "")
    fields.update({"Phage": name, **changes})
    return [fields[key] for key in AUDIT.METADATA_HEADER]


def metadata_records(rows):
    return AUDIT.read_table(table_bytes(AUDIT.METADATA_HEADER, rows, "\t"), "\t", AUDIT.METADATA_HEADER)


def measurement_records(rows):
    return AUDIT.read_table(table_bytes(AUDIT.MEASUREMENT_HEADER, rows, ","), ",", AUDIT.MEASUREMENT_HEADER)


def fixture_policy(root, measurements, metadata):
    root.mkdir(parents=True, exist_ok=True)
    files = []
    for name, data, delimiter in [(AUDIT.MEASUREMENTS, measurements, ","), (AUDIT.METADATA, metadata, "\t")]:
        (root / name).write_bytes(data)
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        files.append({"name": name, "bytes": len(data), "git_blob_sha1": blob, "delimiter": delimiter})
    policy = root / "policy.json"
    policy.write_text(json.dumps({"purpose": "synthetic source audit", "repository": "example/source",
                                  "commit": "frozen", "files": files}))
    return policy


class ReceptorResourceTests(unittest.TestCase):
    def test_nonnegative_decimal_precision_and_raw_strings_survive_without_labels(self):
        for raw in ("0", "-0", "  +1.000000000000000000001  ", "1e-1200", ".125"):
            with self.subTest(raw=raw):
                parsed = AUDIT.measurement(raw)
                self.assertEqual(parsed["raw"], raw)
                self.assertEqual(parsed["kind"], "numeric")
                self.assertNotIn("label", parsed)
        self.assertEqual(AUDIT.measurement("1.000000000000000000001")["value_decimal"], "1.000000000000000000001")
        self.assertEqual(AUDIT.measurement("1e-1200")["value_decimal"], "1E-1200")

    def test_explicit_inequalities_are_bounds_and_never_substituted_with_zero(self):
        for raw, operator in (("<0.001", "<"), (" >= 1e-20", ">="), ("≤0", "≤"), ("≥2", "≥")):
            with self.subTest(raw=raw):
                parsed = AUDIT.measurement(raw)
                self.assertEqual(parsed["kind"], "bounded")
                self.assertEqual(parsed["operator"], operator)
                self.assertNotIn("value_decimal", parsed)
                self.assertEqual(parsed["raw"], raw)
        self.assertEqual(AUDIT.measurement("<-1")["reason"], "negative_bound")
        self.assertEqual(AUDIT.measurement("<Infinity")["reason"], "nonfinite_bound")

    def test_negative_nonfinite_and_uninterpreted_literals_are_not_linear_eop(self):
        for raw in ("-0.918884254", "NaN", "Infinity", "-Infinity", "NA", "n.d.",
                    "0.1 ± 0.02", "1_000", "1,000", "1e", "=0", "<", "< 0.1 extra"):
            with self.subTest(raw=raw):
                parsed = AUDIT.measurement(raw)
                self.assertEqual(parsed["kind"], "unparsed")
                self.assertEqual(parsed["raw"], raw)
                self.assertNotIn("value_decimal", parsed)
        self.assertEqual(AUDIT.measurement("-0.918884254")["reason"], "negative_numeric")

    def test_blank_and_missing_fields_remain_distinct_and_unspecified(self):
        self.assertEqual(AUDIT.measurement(None), {"raw": None, "kind": "blank", "reason": "missing_field"})
        for raw in ("", " \t "):
            parsed = AUDIT.measurement(raw)
            self.assertEqual(parsed["kind"], "blank")
            self.assertEqual(parsed["raw"], raw)
            self.assertNotIn("value_decimal", parsed)

    def test_schema_retains_quoted_cells_physical_lines_and_ragged_rows(self):
        data = b'phage,Genotype,EOP,EOP (complemented)\nP,"line one\nline two",0,\nQ,g,\nR,g,1,2,extra\n'
        rows = AUDIT.read_table(data, ",", AUDIT.MEASUREMENT_HEADER)
        self.assertEqual([r["source_row"] for r in rows], [2, 3, 4])
        self.assertEqual((rows[0]["source_line_start"], rows[0]["source_line_end"]), (2, 3))
        self.assertEqual(rows[0]["fields"]["Genotype"], "line one\nline two")
        self.assertEqual(rows[1]["missing_fields"], ["EOP (complemented)"])
        self.assertIsNone(rows[1]["fields"]["EOP (complemented)"])
        self.assertEqual(rows[2]["extra_cells"], ["extra"])
        with self.assertRaisesRegex(ValueError, "header"):
            AUDIT.read_table(b"phage,wrong\nP,0\n", ",", AUDIT.MEASUREMENT_HEADER)

    def test_primary_id_precedes_alias_and_only_complete_exact_alias_cells_match(self):
        rows = metadata_records([metadata_row("P", **{"Alternative name": "Alias", "Full name": "Full phage P"}),
                                 metadata_row("Other", **{"Alternative name": "P", "Full name": "A; B"})])
        primary, aliases = AUDIT.indexes(rows)
        match = AUDIT.match_phage("P", primary, aliases)
        self.assertEqual(match["rule"], "exact_primary_id")
        self.assertEqual(match["candidates"][0]["metadata_phage"], "P")
        for name in ("Alias", "Full phage P", "A; B"):
            self.assertEqual(AUDIT.match_phage(name, primary, aliases)["status"], "linked")
        for name in ("alias", " Alias", "A", "B", "Full phage"):
            self.assertEqual(AUDIT.match_phage(name, primary, aliases)["status"], "unlinked")
        self.assertFalse(match["stock_identity_verified"])

    def test_ambiguous_primary_or_cross_column_alias_is_never_selected(self):
        rows = metadata_records([metadata_row("P", **{"Alternative name": "shared"}),
                                 metadata_row("Q", **{"Full name": "shared"}), metadata_row("P")])
        primary, aliases = AUDIT.indexes(rows)
        self.assertEqual(AUDIT.match_phage("P", primary, aliases)["status"], "ambiguous")
        self.assertEqual(AUDIT.match_phage("shared", primary, aliases)["status"], "ambiguous")
        rows = metadata_records([metadata_row("P", **{"Alternative name": "same", "Full name": "same"})])
        match = AUDIT.match_phage("same", *AUDIT.indexes(rows))
        self.assertEqual(match["status"], "linked")
        self.assertEqual(match["candidates"][0]["matched_fields"], ["Alternative name", "Full name"])

    def test_duplicates_and_original_complement_contexts_are_preserved_without_averaging(self):
        rows = measurement_records([["P", "Δgene", "-1", ""], ["P", "Δgene", "0", "<0.01"],
                                    ["P", "other", "NA", "1"]])
        summary, records, details = AUDIT.audit_tables(rows, metadata_records([metadata_row("P")]))
        self.assertEqual((summary["measurement_rows"], summary["duplicate_pair_groups"], summary["duplicate_excess_rows"]), (3, 1, 1))
        self.assertEqual(details["duplicate_measurement_pairs"][0]["source_rows"], [2, 3])
        self.assertEqual(records[0]["duplicate_source_rows"], [2, 3])
        self.assertEqual(summary["contexts"]["EOP"]["unparsed"], 2)
        self.assertEqual(summary["contexts"]["EOP (complemented)"]["blank"], 1)
        self.assertEqual(summary["contexts"]["EOP (complemented)"]["bounded"], 1)
        self.assertEqual(records[0]["fields"]["EOP"], "-1")
        self.assertTrue(all(r["admitted"] is False for r in records))

    def test_completeness_is_literal_only_and_name_flags_do_not_promote_stock_identity(self):
        metadata = metadata_records([metadata_row("T7", **{"Genome accesion": "NC_001604.1", "Reference (doi)": "doi-literal"}),
                                     metadata_row("Other", **{"Alternative name": "T3"})])
        rows = measurement_records([["T7", "g", "0", ""], ["T3", "g", "1", ""], ["missing", "g", "1", ""]])
        summary, records, _ = AUDIT.audit_tables(rows, metadata)
        self.assertEqual(summary["linkage_status_counts"], {"linked": 2, "unlinked": 1})
        self.assertEqual(summary["unlinked_or_ambiguous_source_rows"], [4])
        self.assertEqual(summary["metadata_field_completeness"]["Genome accesion"], {"nonblank_literal": 1, "blank_or_absent": 1})
        self.assertEqual(summary["linked_rows_with_nonblank_source_reference_literal"], 1)
        self.assertEqual({r["target_name"] for r in summary["target_name_overlap_warnings"]}, {"T7", "T3"})
        self.assertTrue(all(r["phage_mapping"]["stock_identity_verified"] is False for r in records))

    def test_authenticated_execution_binds_sources_and_writes_only_nonadmitted_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = table_bytes(AUDIT.MEASUREMENT_HEADER, [["P", "g", "-1", "0"]], ",")
            meta = table_bytes(AUDIT.METADATA_HEADER, [metadata_row("P")], "\t")
            policy = fixture_policy(root / "cache", data, meta)
            workdir, output = root / "cache/audit-v1", root / "processed/report.json"
            report = AUDIT.execute(root / "cache", policy, workdir, output)
            self.assertEqual(report["measurement_rows"], 1)
            self.assertEqual(report["policy_sha256"], hashlib.sha256(policy.read_bytes()).hexdigest())
            self.assertFalse(report["measurement_scale_verified"])
            self.assertFalse(report["all_rows_admitted"])
            self.assertFalse(report["training_or_scoring_executed"])
            self.assertFalse(report["original_benchmark_modified"])
            for artifact in report["artifacts"]:
                self.assertEqual(hashlib.sha256(Path(artifact["path"]).read_bytes()).hexdigest(), artifact["sha256"])
            self.assertEqual(json.loads((workdir / "literal-measurements.json").read_text())[0]["measurements"]["EOP"]["raw"], "-1")
            with self.assertRaisesRegex(ValueError, "new"):
                AUDIT.execute(root / "cache", policy, workdir, output)

    def test_changed_source_bytes_are_refused_before_any_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = table_bytes(AUDIT.MEASUREMENT_HEADER, [["P", "g", "1", ""]], ",")
            meta = table_bytes(AUDIT.METADATA_HEADER, [metadata_row("P")], "\t")
            policy = fixture_policy(root / "cache", data, meta)
            (root / "cache" / AUDIT.MEASUREMENTS).write_bytes(data.replace(b",1,", b",2,"))
            with self.assertRaisesRegex(ValueError, "Source bytes differ"):
                AUDIT.execute(root / "cache", policy, root / "work", root / "report.json")
            self.assertFalse((root / "work").exists())
            self.assertFalse((root / "report.json").exists())

    def test_output_locations_must_be_fresh_separate_and_not_overlap_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.csv"
            source.write_text("input")
            for workdir, output in ((root, root / "report.json"), (root / "work", root / "work/report.json"),
                                    (root / "work", source), (source / "child", root / "report.json")):
                with self.subTest(workdir=workdir, output=output), self.assertRaises(ValueError):
                    AUDIT.guard(workdir, output, {source: "hash"})


if __name__ == "__main__":
    unittest.main()
