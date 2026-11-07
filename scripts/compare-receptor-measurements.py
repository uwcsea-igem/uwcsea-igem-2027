"""Check frozen Dataset S4/browser serialization correspondence, not prediction."""
import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal, InvalidOperation, localcontext
import hashlib
import io
import json
from pathlib import Path
import posixpath
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

POLICY = Path("software/tipcraft/benchmarks/receptor-resource/measurement-correspondence-policy.json")
CONTEXTS = ("EOP", "EOP (complemented)")
NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z", re.ASCII)
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def parse_workbook(data, spec):
    """Read saved OOXML literals directly; never evaluate a spreadsheet formula."""
    with ZipFile(io.BytesIO(data)) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = [s for s in workbook.findall("m:sheets/m:sheet", NS) if s.get("name") == spec["sheet"]]
        if len(sheets) != 1:
            raise ValueError("Expected one exact source worksheet name")
        relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        matches = [r for r in relations if r.get("Id") == sheets[0].get("{" + NS["r"] + "}id")]
        if len(matches) != 1 or matches[0].get("TargetMode") == "External":
            raise ValueError("Missing, ambiguous or external worksheet relationship")
        target = matches[0].get("Target", "")
        part = posixpath.normpath(target.lstrip("/") if target.startswith("/") else posixpath.join("xl", target))
        if not part.startswith("xl/worksheets/"):
            raise ValueError("Worksheet relationship escapes the worksheet directory")
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ["".join(t.text or "" for t in item.findall(".//m:t", NS))
                      for item in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("m:si", NS)]
        sheet = ET.fromstring(archive.read(part))
    rows, header, empty_rows, seen = [], None, [], set()
    for row in sheet.findall("m:sheetData/m:row", NS):
        number = int(row.attrib["r"])
        if number in seen:
            raise ValueError("Repeated source worksheet row index")
        seen.add(number)
        cells, types, formulas = {}, {}, {}
        for cell in row.findall("m:c", NS):
            address = cell.get("r", "")
            match = re.fullmatch(r"([A-Z]+)(\d+)", address)
            if match is None or int(match[2]) != number:
                raise ValueError("Worksheet cell/row coordinates disagree")
            column = match[1]
            if column not in spec["columns"]:
                continue
            if column in types:
                raise ValueError("Duplicate worksheet cell address")
            kind, value = cell.get("t", "n"), cell.find("m:v", NS)
            if kind == "inlineStr":
                literal = "".join(t.text or "" for t in cell.findall("m:is//m:t", NS))
            elif kind == "s":
                if value is None or value.text is None or not value.text.isdecimal() or int(value.text) >= len(shared):
                    raise ValueError("Invalid shared-string reference")
                literal = shared[int(value.text)]
            else:
                literal = value.text if value is not None else None
            types[column] = kind
            if literal is not None:
                cells[column] = literal
            formula = cell.find("m:f", NS)
            if formula is not None:
                formulas[column] = formula.text or ""
        if number == spec["header_row"]:
            header = cells
        elif number > spec["header_row"]:
            if any(str(value).strip() for value in cells.values()):
                rows.append({"row": number, "cells": cells, "cell_types": types, "formulas": formulas})
            else:
                empty_rows.append(number)
    if header != spec["columns"]:
        raise ValueError("Exact frozen worksheet header differs")
    return rows, {"worksheet_part": part, "empty_formatted_rows_after_header": empty_rows,
                  "formula_source_rows": [r["row"] for r in rows if r["formulas"]]}


def source_records(rows):
    records = []
    for row in rows:
        phage, strain, value = (row["cells"].get(c) for c in ("A", "B", "C"))
        genotype = strain.split("-pCA")[0].split("::Y")[0] if strain is not None else None
        context = "EOP" if strain is not None and strain.endswith("24N") else (
            "EOP (complemented)" if strain is not None and strain.endswith(")") else None)
        status = "recognized" if context is not None else "unrecognized_suffix"
        if not phage or not genotype:
            status = "missing_key"
        records.append({**row, "phage": phage, "strain_raw": strain, "derived_genotype": genotype,
                        "context": context, "value_raw": value, "key_status": status})
    return records


def browser_records(data):
    reader = csv.reader(io.StringIO(data.decode("utf-8-sig"), newline=""), strict=True)
    if next(reader, None) != ["phage", "Genotype", *CONTEXTS]:
        raise ValueError("Browser schema differs from the frozen source audit")
    records = []
    for row, cells in enumerate(reader, 2):
        if len(cells) != 4:
            raise ValueError("Browser row has missing or extra fields")
        for index, context in enumerate(CONTEXTS, 2):
            records.append({"source_row": row, "phage": cells[0], "genotype": cells[1], "context": context,
                            "value_raw": cells[index], "key_status": "recognized" if cells[0] and cells[1] else "missing_key"})
    return records


def compare_value(source, browser, tolerance):
    result = {"source_raw": source, "browser_raw": browser, "literal_string_equal": source == browser}
    try:
        a, b = Decimal(source), Decimal(browser)
    except (InvalidOperation, TypeError, ValueError):
        return {**result, "status": "noncomparable", "reason": "missing_or_uninterpreted_literal"}
    if not a.is_finite() or not b.is_finite():
        return {**result, "status": "noncomparable", "reason": "nonfinite_literal"}
    if NUMBER.fullmatch(source.strip()) is None or NUMBER.fullmatch(browser.strip()) is None:
        return {**result, "status": "noncomparable", "reason": "nonstandard_decimal_literal"}
    # Precision derives only from the input literal lengths/exponents, not outcomes.
    with localcontext() as context:
        lowest = min(a.as_tuple().exponent, b.as_tuple().exponent, tolerance.as_tuple().exponent)
        context.prec = max(28, max(a.adjusted(), b.adjusted(), tolerance.adjusted()) - lowest + 3)
        difference = abs(a - b)
        within = difference <= tolerance
    return {**result, "status": "decimal_exact" if a == b else "within_fixed_tolerance" if within else "value_mismatch",
            "decimal_exact_equal": a == b, "within_absolute_tolerance_including_exact": within,
            "absolute_difference_decimal": str(difference)}


def correspond(source, browser, tolerance):
    si, bi = defaultdict(list), defaultdict(list)
    for index, row in enumerate(source):
        if row["key_status"] == "recognized":
            si[(row["phage"], row["derived_genotype"], row["context"])].append(index)
    for index, row in enumerate(browser):
        if row["key_status"] == "recognized":
            bi[(row["phage"], row["genotype"], row["context"])].append(index)
    links, ambiguous = [], []
    for key in sorted(set(si) | set(bi)):
        left, right = si.get(key, []), bi.get(key, [])
        if len(left) > 1 or len(right) > 1:
            status = "ambiguous_key"
            ambiguous.append({"key": list(key), "source_rows": [source[i]["row"] for i in left],
                              "browser_rows": [browser[i]["source_row"] for i in right]})
        elif not left or not right:
            status = "unpaired"
        else:
            status = "paired"
            value = compare_value(source[left[0]]["value_raw"], browser[right[0]]["value_raw"], tolerance)
            links.append({"source_row": source[left[0]]["row"], "browser_row": browser[right[0]]["source_row"],
                          "phage": key[0], "genotype": key[1], "context": key[2], **value})
        for index in left:
            source[index]["mapping_status"] = status
            source[index]["browser_candidate_rows"] = [browser[i]["source_row"] for i in right]
        for index in right:
            browser[index]["mapping_status"] = status
            browser[index]["source_candidate_rows"] = [source[i]["row"] for i in left]
    for row in source + browser:
        row.setdefault("mapping_status", row["key_status"])
    counts = Counter(link["status"] for link in links)
    complete = (len(links) == len(source) == len(browser) and all(
        link["status"] in ("decimal_exact", "within_fixed_tolerance") for link in links))
    summary = {
        "source_rows": len(source), "browser_context_records": len(browser), "paired_unique_keys": len(links),
        "comparison_counts": {k: counts[k] for k in ("decimal_exact", "within_fixed_tolerance", "value_mismatch", "noncomparable")},
        "within_absolute_tolerance_including_exact": counts["decimal_exact"] + counts["within_fixed_tolerance"],
        "literal_string_equal_pairs": sum(link["literal_string_equal"] for link in links),
        "source_context_counts": dict(Counter(r["context"] or "unrecognized" for r in source)),
        "browser_context_counts": dict(Counter(r["context"] for r in browser)),
        "source_unpaired_rows": [r["row"] for r in source if r["mapping_status"] == "unpaired"],
        "browser_unpaired_contexts": [{"source_row": r["source_row"], "context": r["context"]} for r in browser if r["mapping_status"] == "unpaired"],
        "unrecognized_suffix_source_rows": [r["row"] for r in source if r["key_status"] == "unrecognized_suffix"],
        "missing_key_source_rows": [r["row"] for r in source if r["key_status"] == "missing_key"],
        "missing_key_browser_contexts": [{"source_row": r["source_row"], "context": r["context"]} for r in browser if r["key_status"] == "missing_key"],
        "ambiguous_keys": ambiguous,
        "value_mismatch_source_rows": [r["source_row"] for r in links if r["status"] == "value_mismatch"],
        "all_source_and_browser_contexts_correspond_within_fixed_tolerance": complete,
    }
    return summary, {"source_records": source, "browser_context_records": browser, "unique_key_comparisons": links}


def authenticate(policy_path):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    raw, paths = {}, {}
    bindings = {policy_path: sha(policy_bytes), Path(__file__): sha(Path(__file__).read_bytes())}
    for item in policy["inputs"]:
        path = Path(item["path"])
        data = path.read_bytes()
        if path.name in raw or sha(data) != item["sha256"]:
            raise ValueError("Duplicate or changed frozen input")
        raw[path.name], paths[path.name], bindings[path] = data, path, item["sha256"]
    audit = json.loads(raw["receptor-resource-audit.json"])
    if audit["policy_sha256"] != sha(raw["source-audit-policy.json"]) or audit["admitted"] is not False:
        raise ValueError("Initial source-audit policy/admission binding differs")
    auditor = Path(__file__).with_name("audit-receptor-resource.py")
    bindings[auditor] = audit["script_sha256"]
    for item in audit["source_files"] + audit["artifacts"]:
        path = Path(item["path"])
        data = path.read_bytes()
        if sha(data) != item["sha256"] or len(data) != item["bytes"]:
            raise ValueError("Initial source-audit artifact/input differs")
        bindings[path] = item["sha256"]
    browser = [i for i in audit["source_files"] if Path(i["path"]).name == "KEIO_EOP_reformatted.csv"]
    if len(browser) != 1 or browser[0]["sha256"] != sha(raw["KEIO_EOP_reformatted.csv"]):
        raise ValueError("Browser table is not the original audited source")
    if policy["decimal_comparison"]["absolute_tolerance"] != "1e-8":
        raise ValueError("Only the frozen absolute serialization tolerance is supported")
    for path, expected in bindings.items():
        if sha(path.read_bytes()) != expected:
            raise ValueError("Initial audit implementation/input hash differs")
    return policy, raw, bindings


def guard(workdir, output, inputs):
    if any(p.exists() or p.is_symlink() for p in (workdir, output)):
        raise ValueError("Workdir and output must be new")
    a, b = workdir.resolve(), output.resolve()
    if a.is_relative_to(b) or b.is_relative_to(a):
        raise ValueError("Workdir and report must be separate")
    for destination in (a, b):
        if any(destination.is_relative_to(p.resolve()) or p.resolve().is_relative_to(destination) for p in inputs):
            raise ValueError("Output overlaps an authenticated input")


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def execute(policy_path, workdir, output):
    policy, raw, bindings = authenticate(policy_path)
    guard(workdir, output, bindings)
    rows, workbook = parse_workbook(raw["datasets-v1.xlsx"], policy["workbook"])
    if [{"row": row["row"], "cells": row["cells"]} for row in rows] != json.loads(raw["literal-workbook-rows.json"]):
        raise ValueError("Direct OOXML extraction differs from the authenticated prior literal rows")
    summary, details = correspond(source_records(rows), browser_records(raw["KEIO_EOP_reformatted.csv"]), Decimal("1e-8"))
    for path, expected in bindings.items():
        if sha(path.read_bytes()) != expected:
            raise ValueError("Authenticated input/implementation changed during comparison")
    workdir.mkdir(parents=True, exist_ok=False)
    artifact = workdir / "row-correspondence.json"
    dump(artifact, details)
    report = {"schema_version": 1, "status": "literal_correspondence_checked", "purpose": policy["purpose"],
              "policy_sha256": bindings[policy_path], "script_sha256": bindings[Path(__file__)],
              "authenticated_inputs": [{"path": str(p), "sha256": h} for p, h in bindings.items()],
              "direct_ooxml_extraction_matches_prior_literal_artifact": True, **workbook, **summary,
              "absolute_tolerance_decimal": "1E-8", "comparison_rule": "Signed finite Decimal literals; exact equality and fixed absolute serialization tolerance; no transform",
              "artifact": {"path": str(artifact), "sha256": sha(artifact.read_bytes()), "bytes": artifact.stat().st_size},
              "initial_source_audit_modified": False, "measurement_scale_verified": False,
              "censoring_or_detection_limit_verified": False, "stock_identity_verified": False,
              "all_rows_admitted": False, "admitted_rows": 0, "training_or_scoring_executed": False,
              "predictions_or_performance_produced": False, "original_benchmark_modified": False,
              "interpretation": "Serialization/key correspondence only. All source rows and browser contexts are retained. No infection labels, inverse logarithms, averaging, assay-denominator inference or endpoint admission are produced."}
    output.parent.mkdir(parents=True, exist_ok=True)
    dump(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = execute(args.policy, args.workdir, args.output)
    print(json.dumps({"status": result["status"], "source_rows": result["source_rows"],
                      "browser_context_records": result["browser_context_records"],
                      "comparison_counts": result["comparison_counts"],
                      "all_correspond_within_tolerance": result["all_source_and_browser_contexts_correspond_within_fixed_tolerance"]}))
