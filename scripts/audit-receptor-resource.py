"""Authenticate and inventory literal receptor-resource measurements.

No interpretation of the measurement scale, assay denominator, bacterial
background, infection threshold or model performance is made by this audit.
"""
import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
from pathlib import Path
import re


POLICY = Path("software/tipcraft/benchmarks/receptor-resource/source-audit-policy.json")
MEASUREMENTS = "KEIO_EOP_reformatted.csv"
METADATA = "Table_S1_updated.tsv"
CONTEXTS = ["EOP", "EOP (complemented)"]
MEASUREMENT_HEADER = ["phage", "Genotype", *CONTEXTS]
METADATA_HEADER = [
    "Phage", "Alternative name", "Full name", "Production host strain", "Lifestyle",
    "Genome accesion", "Genome size (bp)", "Morphotype", "Class", "Family", "Subfamily", "Genus",
    "BW25113 susceptibility", "BL21 susceptibility", "BW25113 receptor", "BW25113 LPS sugar",
    "BL21 receptor", "Receptor-binding protein", "LPS sugar-binding protein", "AF3 model",
    "Predictive models set", "Reference (doi)", "Environmental source (for in-house phages)", "Notes",
]
NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z", re.ASCII)
BOUND = re.compile(r"(<=|>=|<|>|≤|≥)\s*(.+)\Z")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def measurement(raw):
    """Classify literals using Decimal; never substitute, threshold or transform."""
    if raw is None or not raw.strip():
        return {"raw": raw, "kind": "blank", "reason": "missing_field" if raw is None else "empty_or_whitespace"}
    text = raw.strip()
    bound = BOUND.fullmatch(text)
    token = bound[2] if bound else text
    try:
        value = Decimal(token)
    except (InvalidOperation, ValueError):
        return {"raw": raw, "kind": "unparsed", "reason": "uninterpreted_literal"}
    if not value.is_finite():
        return {"raw": raw, "kind": "unparsed", "reason": "nonfinite_bound" if bound else "nonfinite_numeric"}
    if NUMBER.fullmatch(token) is None:
        return {"raw": raw, "kind": "unparsed", "reason": "nonstandard_decimal_literal"}
    if value < 0:
        return {"raw": raw, "kind": "unparsed", "reason": "negative_bound" if bound else "negative_numeric"}
    if bound:
        return {"raw": raw, "kind": "bounded", "operator": bound[1], "bound_decimal": str(value)}
    return {"raw": raw, "kind": "numeric", "value_decimal": str(value)}


def read_table(data, delimiter, expected_header):
    reader = csv.reader(io.StringIO(data.decode("utf-8-sig"), newline=""), delimiter=delimiter, strict=True)
    if next(reader, None) != expected_header:
        raise ValueError("Source header differs from the reviewed schema")
    records = []
    previous_line = reader.line_num
    for ordinal, cells in enumerate(reader, 2):
        records.append({"source_row": ordinal, "source_line_start": previous_line + 1,
                        "source_line_end": reader.line_num,
                        "fields": {name: cells[i] if i < len(cells) else None for i, name in enumerate(expected_header)},
                        "extra_cells": cells[len(expected_header):],
                        "missing_fields": expected_header[len(cells):]})
        previous_line = reader.line_num
    return records


def indexes(metadata):
    primary, other = defaultdict(list), defaultdict(dict)
    for row in metadata:
        fields = row["fields"]
        if fields["Phage"]:
            primary[fields["Phage"]].append(row)
        for field in ("Alternative name", "Full name"):
            value = fields[field]
            if value is not None and value.strip():
                other[value].setdefault(row["source_row"], {"record": row, "fields": []})["fields"].append(field)
    return primary, other


def match_phage(name, primary, other):
    direct = primary.get(name, [])
    if direct:
        candidates = [{"record": row, "fields": ["Phage"]} for row in direct]
        rule = "exact_primary_id"
    else:
        candidates = list(other.get(name, {}).values())
        rule = "unique_exact_complete_name_cell"
    result = {"status": "linked" if len(candidates) == 1 else "ambiguous" if candidates else "unlinked",
              "rule": rule, "stock_identity_verified": False,
              "candidates": [{"metadata_source_row": c["record"]["source_row"],
                              "metadata_phage": c["record"]["fields"]["Phage"],
                              "matched_fields": c["fields"]} for c in candidates]}
    if len(candidates) == 1:
        result["genome_accession_raw"] = candidates[0]["record"]["fields"]["Genome accesion"]
        result["source_reference_raw"] = candidates[0]["record"]["fields"]["Reference (doi)"]
    return result


def nonblank(value):
    return value is not None and bool(value.strip())


def audit_tables(rows, metadata):
    primary, other = indexes(metadata)
    repeated = defaultdict(list)
    for row in rows:
        repeated[(row["fields"]["phage"], row["fields"]["Genotype"])].append(row["source_row"])
    duplicates = [{"phage": key[0], "Genotype": key[1], "source_rows": numbers}
                  for key, numbers in repeated.items() if len(numbers) > 1]
    records = []
    for row in rows:
        fields = row["fields"]
        numbers = repeated[(fields["phage"], fields["Genotype"])]
        records.append({**row, "measurements": {name: measurement(fields[name]) for name in CONTEXTS},
                        "phage_mapping": match_phage(fields["phage"], primary, other),
                        "duplicate_phage_genotype": len(numbers) > 1,
                        "duplicate_source_rows": numbers if len(numbers) > 1 else [], "admitted": False})
    contexts = {}
    for context in CONTEXTS:
        kinds = Counter(r["measurements"][context]["kind"] for r in records)
        contexts[context] = {"rows": len(records), **{kind: kinds[kind] for kind in ("numeric", "bounded", "blank", "unparsed")},
                             "unparsed_reasons": dict(sorted(Counter(r["measurements"][context]["reason"] for r in records
                                                      if r["measurements"][context]["kind"] == "unparsed").items()))}
    warnings = []
    for target in ("T7", "T3"):
        named_rows = [r["source_row"] for r in records if (r["fields"]["phage"] or "").strip().casefold() == target.casefold()]
        metadata_matches = [{"metadata_source_row": r["source_row"], "field": field, "raw": r["fields"][field]}
                            for r in metadata for field in ("Phage", "Alternative name", "Full name")
                            if (r["fields"][field] or "").strip().casefold() == target.casefold()]
        if named_rows or metadata_matches:
            warnings.append({"target_name": target, "measurement_source_rows": named_rows,
                             "metadata_name_matches": metadata_matches,
                             "interpretation": "Name-overlap warning only; no stock identity, homology clearance or prediction"})
    completeness = {}
    for field in ("Phage", "Genome accesion", "Reference (doi)", "Production host strain"):
        present = sum(nonblank(row["fields"][field]) for row in metadata)
        completeness[field] = {"nonblank_literal": present, "blank_or_absent": len(metadata) - present}
    summary = {
        "measurement_rows": len(records), "distinct_literal_phage_ids": len({r["fields"]["phage"] for r in records}),
        "distinct_literal_genotypes": len({r["fields"]["Genotype"] for r in records}),
        "distinct_literal_phage_genotype_pairs": len(repeated), "contexts": contexts,
        "duplicate_pair_groups": len(duplicates), "rows_in_duplicate_groups": sum(len(d["source_rows"]) for d in duplicates),
        "duplicate_excess_rows": sum(len(d["source_rows"]) - 1 for d in duplicates),
        "linkage_status_counts": dict(sorted(Counter(r["phage_mapping"]["status"] for r in records).items())),
        "unlinked_or_ambiguous_source_rows": [r["source_row"] for r in records if r["phage_mapping"]["status"] != "linked"],
        "linked_rows_with_nonblank_genome_accession_literal": sum(nonblank(r["phage_mapping"].get("genome_accession_raw")) for r in records),
        "linked_rows_with_nonblank_source_reference_literal": sum(nonblank(r["phage_mapping"].get("source_reference_raw")) for r in records),
        "metadata_rows": len(metadata), "distinct_nonblank_metadata_primary_ids": len(primary),
        "metadata_duplicate_primary_ids": {name: [r["source_row"] for r in values] for name, values in primary.items() if len(values) > 1},
        "metadata_field_completeness": completeness,
        "measurement_schema_issue_rows": [r["source_row"] for r in records if r["extra_cells"] or r["missing_fields"]],
        "metadata_schema_issue_rows": [r["source_row"] for r in metadata if r["extra_cells"] or r["missing_fields"]],
        "target_name_normalization": "strip surrounding whitespace and casefold; complete cell must equal T7 or T3",
        "target_name_overlap_warnings": warnings,
    }
    return summary, records, {"records": metadata, "duplicate_measurement_pairs": duplicates}


def authenticated(cache, policy_path):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if set(item["name"] for item in policy["files"]) != {MEASUREMENTS, METADATA} or len(policy["files"]) != 2:
        raise ValueError("Expected the two frozen source files")
    raw, bindings = {}, {policy_path: sha(policy_bytes), Path(__file__): sha(Path(__file__).read_bytes())}
    sources = []
    for item in policy["files"]:
        path = cache / item["name"]
        data = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if len(data) != item["bytes"] or blob != item["git_blob_sha1"]:
            raise ValueError(f"Source bytes differ from frozen policy: {item['name']}")
        raw[item["name"]] = data
        bindings[path] = sha(data)
        sources.append({"path": str(path), "bytes": len(data), "sha256": sha(data), "git_blob_sha1": blob,
                        "source_url": f"https://github.com/{policy['repository']}/blob/{policy['commit']}/Ecoli_phages/data/{item['name']}"})
    return policy, raw, bindings, sources


def guard(workdir, output, protected):
    if any(path.exists() or path.is_symlink() for path in (workdir, output)):
        raise ValueError("Workdir and report output must both be new")
    a, b = workdir.resolve(), output.resolve()
    if a.is_relative_to(b) or b.is_relative_to(a):
        raise ValueError("Workdir and report must be separate")
    for destination in (a, b):
        for source in protected:
            source = source.resolve()
            if destination.is_relative_to(source) or source.is_relative_to(destination):
                raise ValueError("Output overlaps an authenticated input")


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def execute(cache, policy_path, workdir, output):
    policy, raw, bindings, sources = authenticated(cache, policy_path)
    guard(workdir, output, bindings)
    specs = {item["name"]: item for item in policy["files"]}
    rows = read_table(raw[MEASUREMENTS], specs[MEASUREMENTS]["delimiter"], MEASUREMENT_HEADER)
    metadata = read_table(raw[METADATA], specs[METADATA]["delimiter"], METADATA_HEADER)
    summary, measurements, mapping = audit_tables(rows, metadata)
    for path, expected in bindings.items():
        if sha(path.read_bytes()) != expected:
            raise ValueError("Authenticated input or audit implementation changed during parsing")
    workdir.mkdir(parents=True, exist_ok=False)
    literal_path, metadata_path = workdir / "literal-measurements.json", workdir / "metadata-and-duplicates.json"
    dump(literal_path, measurements)
    dump(metadata_path, mapping)
    report = {"schema_version": 1, "purpose": policy["purpose"], "status": "source_literals_audited_not_admitted",
              "policy_path": str(policy_path), "policy_sha256": bindings[policy_path],
              "script_sha256": bindings[Path(__file__)], "source_repository": policy["repository"],
              "source_commit": policy["commit"], "source_files": sources, **summary,
              "artifacts": [{"path": str(path), "sha256": sha(path.read_bytes()), "bytes": path.stat().st_size}
                            for path in (literal_path, metadata_path)],
              "stock_identity_verified": False, "admitted": False, "admitted_rows": 0,
              "all_rows_admitted": False, "measurement_scale_verified": False,
              "accessions_or_source_references_externally_verified": False,
              "measurement_scale_normalizer_conditions_or_genotype_background_inferred": False,
              "binding_entry_or_assembly_inferred": False, "training_or_scoring_executed": False,
              "predictions_or_performance_produced": False, "original_benchmark_modified": False,
              "interpretation": "Literal source audit only. Negative numeric strings remain uninterpreted under the frozen nonnegative rule; no log-scale interpretation is assumed. Reference/accession completeness means nonblank text, not verification. Complementation blanks are unspecified, not failed rescue. Repeated rows are retained without averaging or claims of independent replication."}
    output.parent.mkdir(parents=True, exist_ok=True)
    dump(output, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = execute(args.cache, args.policy, args.workdir, args.output)
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": report["status"], "measurement_rows": report["measurement_rows"],
                      "metadata_rows": report["metadata_rows"], "contexts": report["contexts"], "admitted": False}))


if __name__ == "__main__":
    main()
