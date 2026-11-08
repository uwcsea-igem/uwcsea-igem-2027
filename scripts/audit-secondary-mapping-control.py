"""Post-hoc read-only diagnosis of one preserved copied-proteome failure.

No external searches, phenotype inputs, model scores or control changes.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = "htg0000p0000695"
COPY = "control_" + ORIGINAL


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def checked(binding):
    path = ROOT / binding["path"]
    if sha(path) != binding["sha256"]:
        raise ValueError("Artifact checksum differs: " + str(path))
    return path


def f32(value):
    return struct.unpack("f", struct.pack("f", float(value)))[0]


def focal_rows(path, query):
    rows = []
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith(query + "\t"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 11:
                    raise ValueError("Expected original 11-field row")
                rows.append({"raw_row": number, "fields": fields})
    if not rows:
        raise ValueError("Missing focal raw query rows")
    return rows


def extract(path, wanted):
    sequences, current, chunks = {}, None, []
    with path.open() as stream:
        for line in stream:
            if line.startswith(">"):
                if current in wanted:
                    if current in sequences:
                        raise ValueError("Duplicate normalized protein ID")
                    sequences[current] = "".join(chunks)
                current, chunks = line[1:].split()[0], []
            elif current in wanted:
                chunks.append(line.rstrip("\r\n"))
        if current in wanted:
            sequences[current] = "".join(chunks)
    return sequences


def audit(report_path, expected_sha):
    report_binding = {"path": str(report_path), "sha256": expected_sha}
    report = json.loads(checked(report_binding).read_text())
    if report["status"] != "failed_mapping_control" or report["mapping_control_passed"] is not False:
        raise ValueError("Require original preserved failed control report")
    input_bindings = [report_binding, report["frozen_plan"]]
    checked(report["frozen_plan"])
    names = ["training.faa", "queries.faa", "protein-map.json", "self.tsv", "query.tsv",
             "best-self.json", "best-query.json", "vocabulary.json", "mapping-control-check.json"]
    files = {}
    for name in names:
        matches = [b for b in report["artifacts"] if Path(b["path"]).name == name]
        if len(matches) != 1:
            raise ValueError("Expected exactly one bound input artifact: " + name)
        files[name] = checked(matches[0])
        input_bindings.append(matches[0])
    source = json.loads(files["best-self.json"].read_text())[ORIGINAL]
    clone = json.loads(files["best-query.json"].read_text())[COPY]
    vocabulary = json.loads(files["vocabulary.json"].read_text())
    mapping_rows = json.loads(files["protein-map.json"].read_text())
    mapping = {r["protein"]: r for r in mapping_rows}
    if len(mapping) != len(mapping_rows):
        raise ValueError("Duplicate protein map IDs")
    ids = {ORIGINAL, COPY, source["target"], clone["target"], source["family"], clone["family"]}
    sequences = extract(files["training.faa"], ids) | extract(files["queries.faa"], ids)
    if set(sequences) != ids:
        raise ValueError("Missing focal protein sequence")
    sequence_metadata = {}
    for identifier, sequence in sequences.items():
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        row = mapping[identifier]
        if len(sequence) != row["length"] or digest != row["sequence_sha256"]:
            raise ValueError("Actual normalized sequence differs from frozen map")
        sequence_metadata[identifier] = {"length": len(sequence), "sequence_sha256": digest,
                                         "genome": row["genome"], "role": row["role"],
                                         "source_protein_id": row["source_protein_id"],
                                         "equal_to_focal_sequence": sequence == sequences[ORIGINAL]}
    rows = {"self": focal_rows(files["self.tsv"], ORIGINAL), "clone": focal_rows(files["query.tsv"], COPY)}
    for name, hit in [("self", source), ("clone", clone)]:
        selected = min(rows[name], key=lambda r: (-f32(r["fields"][3]), f32(r["fields"][4])))
        fields = selected["fields"]
        if (selected["raw_row"] != hit["raw_row"] or fields[1] != hit["target"]
                or vocabulary["training_protein_families"][fields[1]] != hit["family"]
                or f32(fields[3]) != hit["identity_float32"] or f32(fields[4]) != hit["evalue_float32"]):
            raise ValueError("Raw-first float32 hit selection does not replay saved best hit")
    by_target = {name: {r["fields"][1]: r for r in values} for name, values in rows.items()}
    if any(len(by_target[n]) != len(rows[n]) for n in rows):
        raise ValueError("Duplicate focal target rows; require explicit multiset audit")
    a, b = by_target["self"], by_target["clone"]
    common = set(a) & set(b)
    differences = [{"target": t, "self_row": a[t], "clone_row": b[t],
                    "changed_field_indices": [i for i in range(1, 11) if a[t]["fields"][i] != b[t]["fields"][i]]}
                   for t in sorted(common) if a[t]["fields"][1:] != b[t]["fields"][1:]]
    tie_metadata = {}
    for name, values in rows.items():
        best_key = min((-f32(r["fields"][3]), f32(r["fields"][4])) for r in values)
        tied = [r for r in values if (-f32(r["fields"][3]), f32(r["fields"][4])) == best_key]
        tie_metadata[name] = {"best_identity_float32": -best_key[0], "best_evalue_float32": best_key[1],
                              "best_key_rows": len(tied), "families_at_best_key": sorted({vocabulary["training_protein_families"][r["fields"][1]] for r in tied}),
                              "raw_first_row": tied[0]["raw_row"], "raw_first_target": tied[0]["fields"][1]}
    genome = mapping[ORIGINAL]["genome"]
    original_vector, clone_vector = report["vectors"][genome], report["mapping_control"]["vector"]
    differences_in_vector = [{"column": i, "feature": f, "source": original_vector[i], "clone": clone_vector[i]}
                             for i, f in enumerate(report["feature_names"]) if original_vector[i] != clone_vector[i]]
    control = json.loads(files["mapping-control-check.json"].read_text())
    failures = [r for r in control["protein_family_comparisons"] if not r["passed"]]
    output = {"schema_version": 1, "status": "completed_post_hoc_artifact_diagnostic_failure_preserved",
              "helper_sha256": sha(__file__), "input_bindings": input_bindings,
              "focal_original": ORIGINAL, "focal_clone": COPY, "sequence_metadata": sequence_metadata,
              "original_clone_sequence_equal": sequences[ORIGINAL] == sequences[COPY],
              "selected_hits": {"self": source, "clone": clone}, "raw_focal_rows": rows,
              "candidate_comparison": {"self_rows": len(rows["self"]), "clone_rows": len(rows["clone"]),
                  "same_target_order": [r["fields"][1] for r in rows["self"]] == [r["fields"][1] for r in rows["clone"]],
                  "self_only_targets": sorted(set(a)-set(b)), "clone_only_targets": sorted(set(b)-set(a)),
                  "identical_nonquery_hit_tuples": len(common)-len(differences), "changed_tuples": differences},
              "float32_best_key_ties": tie_metadata,
              "family_features": {f: {"feature": vocabulary["family_features"].get(f),
                                       "training_occurrence_pattern": vocabulary["retained_family_patterns"].get(f)}
                                  for f in [source["family"], clone["family"]]},
              "vector_differences": differences_in_vector,
              "copied_control_proteins": len(control["protein_family_comparisons"]), "failed_copied_proteins": failures,
              "original_failed_control_status_unchanged": True, "new_search_parameter_or_control_changes": False,
              "phenotype_labels_models_or_target_scores_used": False,
              "cause_supported_by_artifacts": "Identical input sequence and identical candidate target inventory/order, but exact-self row has reported identity1.00 in training and0.923 for the clone. Training therefore selects that1.00 row; clone uses raw-first among equal float32 keys and selects a different pre-existing family, changing a retained feature bit.",
              "upstream_code_level_explanation": "Not established by this artifact-only audit; inspect pinned MMseqs self-hit/identity-estimation code and effective defaults before proposing any separately frozen diagnostic.",
              "audit_scope": "Selected normalized protein, HSP and vocabulary artifacts; not a complete independent preprocessing replay, no biological score or relaxed control."}
    for binding in input_bindings:
        checked(binding)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a fresh diagnostic report")
    result = audit(args.report, args.expected_report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "original_clone_sequence_equal": result["original_clone_sequence_equal"],
                      "candidate_comparison": {k: v for k, v in result["candidate_comparison"].items() if k != "changed_tuples"},
                      "float32_best_key_ties": result["float32_best_key_ties"], "vector_differences": result["vector_differences"]}))
