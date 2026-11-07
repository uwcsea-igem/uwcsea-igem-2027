"""Post hoc audit of exact reference gp17 proteins in completed feature mapping.

This endpoint was added after inspecting the whole-proteome feature pilot.
It runs no searches, changes no vocabulary, and makes no binding inference.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


MANIFEST = "data/processed/genophi-phage-inputs.json"
MANIFEST_SHA256 = "d0f66b3bba527aa8cc32b327ed7f15107fbcd0943faefa6b1fda2496f528185c"
FEATURE_REPORT = "data/processed/genophi-phage-feature-pilot-verified.json"
FEATURE_REPORT_SHA256 = "b8a4ea98682d64fe1bc6cb84bf5976fc578170d41745be7a6859cacf9e0bd7e8"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def checked(path, expected):
    data = Path(path).read_bytes()
    if digest(data) != expected:
        raise ValueError(f"Checksum differs: {path}")
    return data


def summarize_reference(name, reference, protein_map, best_query, family_features, raw_counts):
    matches = [p for p in protein_map if p["genome"] == name and p["role"] == "target"
               and p["sequence_sha256"] == reference["protein_sequence_sha256"]]
    copies = []
    for protein in matches:
        query = protein["protein"]
        hit = best_query.get(query)
        if bool(raw_counts[query]) != bool(hit):
            raise ValueError("Raw query hits and best-hit mapping disagree")
        family = hit["family"] if hit else None
        copies.append({"query_id": query, "length": protein["length"],
                       "raw_hit_records": raw_counts[query], "best_hit": hit,
                       "assigned_family": family, "retained_feature": family_features.get(family)})
    return {"target": name, "reference_protein_accession": reference["protein_id"],
            "reference_sequence_sha256": reference["protein_sequence_sha256"],
            "exact_reference_copies_in_target_proteome": len(matches), "copies": copies,
            "any_hit_returned_under_frozen_settings": any(c["raw_hit_records"] for c in copies),
            "any_retained_feature": any(c["retained_feature"] is not None for c in copies),
            "status": "exact_reference_not_found" if not copies else
                      "no_hit_returned" if not any(c["raw_hit_records"] for c in copies) else "hit_returned"}


def run():
    manifest = json.loads(checked(MANIFEST, MANIFEST_SHA256))
    report = json.loads(checked(FEATURE_REPORT, FEATURE_REPORT_SHA256))
    artifacts = {Path(a["path"]).name: a for a in report["artifacts"]}
    loaded, bindings = {}, []
    for name in ["protein-map.json", "best-query.json", "feature-map.json", "vocabulary-freeze.json", "query.tsv"]:
        spec = artifacts[name]
        data = checked(spec["path"], spec["sha256"])
        loaded[name] = data if name.endswith(".tsv") else json.loads(data)
        bindings.append(spec)
    frozen_map = [a for a in loaded["vocabulary-freeze.json"]["artifacts"] if Path(a["path"]).name == "feature-map.json"]
    if len(frozen_map) != 1 or frozen_map[0] != artifacts["feature-map.json"]:
        raise ValueError("Feature map does not match pre-target vocabulary freeze")
    counts = Counter()
    for line in loaded["query.tsv"].decode().splitlines():
        fields = line.split("\t")
        if len(fields) != 11:
            raise ValueError("Expected 11-field raw query records")
        counts[fields[0]] += 1
    family_features = {item["family"]: item["feature"] for item in loaded["feature-map.json"]}
    return {"schema_version": 1,
            "purpose": "Exact reference gp17 mapping audit of completed whole-proteome features",
            "analysis_status": "post_hoc_after_feature_outputs; not a preregistered endpoint",
            "script_sha256": digest(Path(__file__).read_bytes()),
            "inputs": [{"path": MANIFEST, "sha256": MANIFEST_SHA256},
                       {"path": FEATURE_REPORT, "sha256": FEATURE_REPORT_SHA256}, *bindings],
            "new_searches_executed": False, "parameters_or_vocabulary_changed": False,
            "phenotype_labels_loaded": False,
            "results": [summarize_reference(name, manifest["reference_proteins"][name], loaded["protein-map.json"],
                                             loaded["best-query.json"], family_features, counts) for name in ["T7", "T3"]],
            "interpretation": "Whole-genome feature differences do not establish recognition-head representation. No returned gp17 hit under these frozen settings is a feature-mapping observation, not a conclusion about binding or infection."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/processed/genophi-gp17-feature-mapping.json"))
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Preserve existing outputs")
    result = run()
    protected = [Path(i["path"]).resolve() for i in result["inputs"]] + [Path(__file__).resolve()]
    if args.output.resolve() in protected:
        parser.error("Output overlaps input")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result["results"]))
