"""Cross-check published training labels against the original experiment table.

Names are matched by a fixed punctuation-normalization rule, never phenotype
similarity. Blank source cells are excluded pending source clarification.
This is data curation, not a scoring run or sequence-identity test.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re

SOURCE_COMMIT = "99656d8c020fcd1dbbc30a02ee10274a07655562"
SOURCE_BLOB = "25a3f894c5a3832e81b42c8de60c12ffb9429569"
MODEL_TABLE_SHA256 = "c4795d786a1e13ac998e263dbf4aaf95fc0f94c52b842d02b6f7fc6e99ccfef2"


def reconcile(model_rows: list[dict], source_rows: list[dict]) -> dict:
    normalize = lambda s: re.sub(r"[^a-z0-9]", "", s.casefold())
    source_index = {}
    for row in source_rows:
        key = normalize(row["bacteria"])
        if not key or key in source_index:
            raise ValueError("Ambiguous original-source host identifiers")
        source_index[key] = row
    if not source_index or not model_rows:
        raise ValueError("Both source tables must contain records")
    names, missing, conflicts, seen = {}, [], [], set()
    concordant = 0
    for row in model_rows:
        host, phage, label = row["strain"], row["phage"], row["interaction"]
        if (host, phage) in seen or label not in {"0", "1"}:
            raise ValueError("Invalid or duplicate model-table interaction")
        seen.add((host, phage))
        if normalize(host) not in source_index:
            raise ValueError(f"No original-source host ID for {host}")
        source = source_index[normalize(host)]
        if phage not in source:
            raise ValueError(f"No original-source phage ID for {phage}")
        names[host] = source["bacteria"]
        value = source[phage]
        record = {"strain": host, "source_strain": source["bacteria"], "phage": phage,
                  "imported_label": int(label)}
        if value == "":
            missing.append(record)
            continue
        if value not in {"0.0", "1.0", "2.0", "3.0", "4.0"}:
            raise ValueError(f"Unexpected original assay value for {host}/{phage}")
        expected = int(float(value) > 0)
        if expected != int(label):
            conflicts.append({**record, "source_value": value, "expected_binary": expected})
        else:
            concordant += 1
    # Different model IDs must not be silently merged onto the same source host.
    if len(set(names.values())) != len(names):
        raise ValueError("Multiple model host IDs map to one source host")
    return {
        "model_rows": len(model_rows), "source_hosts": len(source_rows),
        "crosschecked_host_correspondence": dict(sorted(names.items())),
        "concordant_nonblank_pairs": concordant,
        "source_blank_pairs": sorted(missing, key=lambda r: (r["strain"], r["phage"])),
        "conflicting_nonblank_pairs": sorted(conflicts, key=lambda r: (r["strain"], r["phage"])),
        "eligible_label_rows_before_sequence_exclusions": concordant,
        "label_curation_status": "consistent_nonblank_subset" if not conflicts else "source_conflicts_require_review",
    }


def run(source: Path, model_table: Path) -> dict:
    original = source.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(original)).encode() + b"\0" + original).hexdigest()
    if blob != SOURCE_BLOB:
        raise ValueError("Original interaction table differs from pinned source")
    model = model_table.read_bytes()
    if hashlib.sha256(model).hexdigest() != MODEL_TABLE_SHA256:
        raise ValueError("GenoPHI table differs from reviewed source")
    result = reconcile(list(csv.DictReader(io.StringIO(model.decode()))),
                       list(csv.DictReader(io.StringIO(original.decode()), delimiter=";")))
    return {
        "schema_version": 1, "reviewed_date": "2026-11-07",
        "source_commit": SOURCE_COMMIT,
        "original_table_url": f"https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/{SOURCE_COMMIT}/data/interactions/interaction_matrix.csv",
        "original_table_git_blob_sha1": blob,
        "original_table_sha256": hashlib.sha256(original).hexdigest(),
        "genophi_table_sha256": MODEL_TABLE_SHA256,
        "name_matching_rule": "unique case-insensitive alphanumeric identifier; no phenotype-based matching",
        "binary_rule": "original nonblank MLC value > 0; consistent with published method",
        "source_blank_policy": "exclude from future fitting/calibration/evaluation pending source clarification; never infer negative from absence",
        **result,
        "sequence_identity_verified": False,
        "training_or_prediction_executed": False,
        "interpretation": "Cross-source curation of external training data. Agreement corroborates naming correspondence, not stock identity or model performance. Blank-cell semantics remain unresolved.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--model-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() in {args.source.resolve(), args.model_table.resolve()}:
        parser.error("Output cannot overwrite source tables")
    report = run(args.source, args.model_table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
