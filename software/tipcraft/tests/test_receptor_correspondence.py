"""Test exact serialization joins and direct OOXML parsing without predictions."""
from decimal import Decimal
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile


SPEC = importlib.util.spec_from_file_location(
    "receptor_correspondence", Path(__file__).resolve().parents[3] / "scripts/compare-receptor-measurements.py")
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)
WORKBOOK_SPEC = {"sheet": "Dataset S4", "header_row": 2, "columns": {"A": "Phage", "B": "Strain", "C": "EOP"}}


def workbook_bytes(target="worksheets/sheet4.xml", external=False, body=None):
    stream = io.BytesIO()
    mode = ' TargetMode="External"' if external else ""
    if body is None:
        body = ('<row r="3"><c r="A3" t="s"><v>3</v></c>'
                '<c r="B3" t="inlineStr"><is><t>Δx-pCA24N</t></is></c>'
                '<c r="C3"><v>0.1</v></c></row>'
                '<row r="4"><c r="A4" t="s"><v>3</v></c>'
                '<c r="B4" t="inlineStr"><is><t>Δx-pCA24N(x)</t></is></c>'
                '<c r="C4"><v>1.0</v></c></row>'
                '<row r="5"><c r="C5" s="2"/></row>')
    with ZipFile(stream, "w") as z:
        z.writestr("xl/workbook.xml", '<workbook xmlns="' + COMPARE.NS["m"] + '" xmlns:r="' + COMPARE.NS["r"] + '"><sheets><sheet name="Dataset S4" sheetId="4" r:id="rId4"/></sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId4" Target="' + target + '"' + mode + '/></Relationships>')
        z.writestr("xl/sharedStrings.xml", '<sst xmlns="' + COMPARE.NS["m"] + '"><si><t>Phage</t></si><si><t>Strain</t></si><si><t>EOP</t></si><si><r><t>Phage</t></r><r><t>1</t></r></si></sst>')
        z.writestr("xl/worksheets/sheet4.xml", '<worksheet xmlns="' + COMPARE.NS["m"] + '"><sheetData><row r="2"><c r="A2" t="s"><v>0</v></c><c r="B2" t="s"><v>1</v></c><c r="C2" t="s"><v>2</v></c></row>' + body + '</sheetData></worksheet>')
    return stream.getvalue()


def source_row(number, phage, strain, value):
    return {"row": number, "cells": {"A": phage, "B": strain, "C": value}}


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False))
    return path


def fixture_inputs(root):
    workbook = root / "datasets-v1.xlsx"
    workbook.write_bytes(workbook_bytes())
    literal = write_json(root / "literal-workbook-rows.json", [source_row(3, "Phage1", "Δx-pCA24N", "0.1"),
                                                             source_row(4, "Phage1", "Δx-pCA24N(x)", "1.0")])
    browser = root / "KEIO_EOP_reformatted.csv"
    browser.write_text("phage,Genotype,EOP,EOP (complemented)\nPhage1,Δx,0.1000000,1\n")
    original_policy = write_json(root / "source-audit-policy.json", {"synthetic": True})
    prior_artifact = write_json(root / "original-literals.json", {"synthetic": "preserved"})
    auditor = Path(COMPARE.__file__).with_name("audit-receptor-resource.py")
    binding = lambda p: {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "bytes": p.stat().st_size}
    source_audit = write_json(root / "receptor-resource-audit.json", {
        "policy_sha256": binding(original_policy)["sha256"], "admitted": False,
        "script_sha256": binding(auditor)["sha256"], "source_files": [binding(browser)],
        "artifacts": [binding(prior_artifact)]})
    policy = write_json(root / "comparison-policy.json", {"purpose": "synthetic correspondence", "workbook": WORKBOOK_SPEC,
        "inputs": [{"path": str(p), "sha256": binding(p)["sha256"]}
                   for p in (workbook, literal, browser, source_audit, original_policy)],
        "decimal_comparison": {"absolute_tolerance": "1e-8"}})
    return policy, browser, prior_artifact


class ReceptorCorrespondenceTests(unittest.TestCase):
    def test_signed_decimal_exactness_and_fixed_tolerance_boundary(self):
        tolerance = Decimal("1e-8")
        exact = COMPARE.compare_value("-0.1", "-0.100000000", tolerance)
        self.assertEqual(exact["status"], "decimal_exact")
        self.assertFalse(exact["literal_string_equal"])
        self.assertEqual(COMPARE.compare_value("0", "1e-8", tolerance)["status"], "within_fixed_tolerance")
        self.assertEqual(COMPARE.compare_value("0", "0.0000000100000000000000000000001", tolerance)["status"], "value_mismatch")
        tiny = COMPARE.compare_value("1e-1200", "0", tolerance)
        self.assertEqual(tiny["status"], "within_fixed_tolerance")
        self.assertNotEqual(Decimal(tiny["absolute_difference_decimal"]), 0)

    def test_missing_bounds_nonfinite_and_nonstandard_values_are_not_coerced(self):
        for raw in (None, "", "<0.1", "NaN", "Infinity", "1_000", "unknown"):
            with self.subTest(raw=raw):
                result = COMPARE.compare_value(raw, "0", Decimal("1e-8"))
                self.assertEqual(result["status"], "noncomparable")
                self.assertEqual(result["source_raw"], raw)
                self.assertNotIn("absolute_difference_decimal", result)

    def test_direct_ooxml_reads_shared_inline_and_numeric_literals(self):
        rows, details = COMPARE.parse_workbook(workbook_bytes(), WORKBOOK_SPEC)
        self.assertEqual([{k: r[k] for k in ("row", "cells")} for r in rows],
                         [source_row(3, "Phage1", "Δx-pCA24N", "0.1"), source_row(4, "Phage1", "Δx-pCA24N(x)", "1.0")])
        self.assertEqual(details["empty_formatted_rows_after_header"], [5])
        self.assertEqual(details["formula_source_rows"], [])

    def test_ooxml_preserves_formula_cached_value_without_evaluation(self):
        body = '<row r="3"><c r="A3" t="s"><v>3</v></c><c r="B3" t="inlineStr"><is><t>g-pCA24N</t></is></c><c r="C3"><f>1+2</f><v>3</v></c></row>'
        rows, details = COMPARE.parse_workbook(workbook_bytes(body=body), WORKBOOK_SPEC)
        self.assertEqual(rows[0]["cells"]["C"], "3")
        self.assertEqual(rows[0]["formulas"], {"C": "1+2"})
        self.assertEqual(details["formula_source_rows"], [3])

    def test_ooxml_rejects_external_and_escaping_relationships(self):
        for target, external in (("https://example.invalid/sheet.xml", True), ("../../outside.xml", False)):
            with self.subTest(target=target), self.assertRaises(ValueError):
                COMPARE.parse_workbook(workbook_bytes(target=target, external=external), WORKBOOK_SPEC)

    def test_formatter_keys_are_exact_and_unknown_suffixes_are_retained(self):
        rows = [source_row(3, "P", "Δg::Y-marker-pCA24N", "-1"),
                source_row(4, "P", "Δg-pCA24N(g)", "0"), source_row(5, "P", "Δg-unknown", "1")]
        records = COMPARE.source_records(rows)
        self.assertEqual([r["context"] for r in records], ["EOP", "EOP (complemented)", None])
        self.assertEqual(records[0]["derived_genotype"], "Δg")
        self.assertEqual(records[2]["key_status"], "unrecognized_suffix")
        self.assertEqual(records[2]["value_raw"], "1")

    def test_both_browser_contexts_and_source_rows_are_paired_without_transforms(self):
        source = COMPARE.source_records([source_row(3, "P", "g-pCA24N", "-0.91888425400542195"),
                                         source_row(4, "P", "g-pCA24N(g)", "0")])
        browser = COMPARE.browser_records(b"phage,Genotype,EOP,EOP (complemented)\nP,g,-0.918884254,0.0\n")
        summary, details = COMPARE.correspond(source, browser, Decimal("1e-8"))
        self.assertTrue(summary["all_source_and_browser_contexts_correspond_within_fixed_tolerance"])
        self.assertEqual(summary["comparison_counts"]["decimal_exact"], 1)
        self.assertEqual(summary["comparison_counts"]["within_fixed_tolerance"], 1)
        self.assertEqual((len(details["source_records"]), len(details["browser_context_records"])), (2, 2))

    def test_duplicates_on_either_side_stay_ambiguous_without_collapse(self):
        source = COMPARE.source_records([source_row(3, "P", "g-pCA24N", "1"), source_row(4, "P", "g-pCA24N", "1")])
        browser = COMPARE.browser_records(b"phage,Genotype,EOP,EOP (complemented)\nP,g,1,0\nP,g,1,0\n")
        summary, details = COMPARE.correspond(source, browser, Decimal("1e-8"))
        self.assertEqual(summary["paired_unique_keys"], 0)
        self.assertEqual(len(summary["ambiguous_keys"]), 2)
        self.assertEqual(len(details["source_records"]), 2)
        self.assertEqual(len(details["browser_context_records"]), 4)
        self.assertTrue(all(r["mapping_status"] == "ambiguous_key" for r in details["source_records"]))

    def test_unrecognized_and_unpaired_keys_are_reported_without_aliases(self):
        source = COMPARE.source_records([source_row(3, "P", "g-pCA24N", "1"), source_row(4, "P", "g-other", "2")])
        browser = COMPARE.browser_records(b"phage,Genotype,EOP,EOP (complemented)\np,g,1,2\n")
        summary, _ = COMPARE.correspond(source, browser, Decimal("1e-8"))
        self.assertEqual(summary["source_unpaired_rows"], [3])
        self.assertEqual(summary["unrecognized_suffix_source_rows"], [4])
        self.assertEqual(len(summary["browser_unpaired_contexts"]), 2)
        self.assertFalse(summary["all_source_and_browser_contexts_correspond_within_fixed_tolerance"])

    def test_execution_authenticates_prior_audit_and_retains_all_rows_unadmitted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            policy, _, original = fixture_inputs(root)
            before = original.read_bytes()
            result = COMPARE.execute(policy, root / "work", root / "report.json")
            self.assertEqual((result["source_rows"], result["browser_context_records"]), (2, 2))
            self.assertEqual(result["comparison_counts"]["decimal_exact"], 2)
            self.assertTrue(result["direct_ooxml_extraction_matches_prior_literal_artifact"])
            self.assertFalse(result["all_rows_admitted"])
            self.assertFalse(result["measurement_scale_verified"])
            self.assertEqual(original.read_bytes(), before)
            artifact = result["artifact"]
            self.assertEqual(hashlib.sha256(Path(artifact["path"]).read_bytes()).hexdigest(), artifact["sha256"])

    def test_changed_authenticated_input_is_refused_without_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            policy, browser, _ = fixture_inputs(root)
            browser.write_text("changed")
            with self.assertRaisesRegex(ValueError, "changed frozen input"):
                COMPARE.execute(policy, root / "work", root / "report.json")
            self.assertFalse((root / "work").exists())

    def test_input_change_during_comparison_is_refused_before_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            policy, browser, _ = fixture_inputs(root)
            original = COMPARE.correspond
            def changed(*args):
                result = original(*args)
                browser.write_text("changed during parsing")
                return result
            with patch.object(COMPARE, "correspond", side_effect=changed), self.assertRaisesRegex(ValueError, "changed during comparison"):
                COMPARE.execute(policy, root / "work", root / "report.json")
            self.assertFalse((root / "work").exists())

    def test_existing_outputs_and_input_overlap_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            source.write_text("input")
            for work, output in ((root, root / "out.json"), (root / "work", root / "work/out.json"),
                                 (root / "work", source), (source / "child", root / "out.json")):
                with self.subTest(work=work, output=output), self.assertRaises(ValueError):
                    COMPARE.guard(work, output, {source: "hash"})


if __name__ == "__main__":
    unittest.main()
