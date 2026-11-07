"""Group every authenticated phage representation without phenotype labels.

Connected components retain all source and target bridge nodes. These are
conservative relatedness/exclusion groups, not inferred taxonomy or host range.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path


FIELDS = ["query_id", "subject_id", "query_source", "subject_source",
          "similarity_percent_unrounded", "query_nident", "subject_nident", "query_length", "subject_length"]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def checked(path, expected):
    data = Path(path).read_bytes()
    if digest(data) != expected:
        raise ValueError(f"Checksum differs: {path}")
    return data


def _name(value):
    return (isinstance(value, str) and bool(value) and value == value.strip()
            and not any(ord(c) < 32 or ord(c) == 127 for c in value))


def parse_pairs(matrix_text, genomes):
    by_id = {g["id"]: g for g in genomes}
    if (not genomes or len(by_id) != len(genomes) or len({g["source_id"] for g in genomes}) != len(genomes)
            or any(not _name(g["id"]) or not _name(g["source_id"])
                   or g["role"] not in {"external_collection", "target_reference"}
                   or type(g["retained_bases"]) is not int or g["retained_bases"] < 1 for g in genomes)):
        raise ValueError("Invalid or duplicate phage manifest nodes")
    reader = csv.DictReader(io.StringIO(matrix_text), delimiter="\t")
    if reader.fieldnames != FIELDS:
        raise ValueError("Expected full unrounded phage-matrix schema")
    pairs = {}
    for row in reader:
        if None in row or any(v is None for v in row.values()):
            raise ValueError("Incomplete or extra matrix fields")
        a, b = row["query_id"], row["subject_id"]
        if a not in by_id or b not in by_id or a == b:
            raise ValueError("Unexpected matrix node or self pair")
        key = tuple(sorted((a, b)))
        if key in pairs:
            raise ValueError("Duplicate unordered matrix pair")
        try:
            value = float(row["similarity_percent_unrounded"])
            ai, bi, alen, blen = [int(row[k]) for k in ["query_nident", "subject_nident", "query_length", "subject_length"]]
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid numeric matrix value") from exc
        if (not math.isfinite(value) or not 0 <= value <= 100 or min(ai, bi) < 0
                or row["query_source"] != by_id[a]["source_id"] or row["subject_source"] != by_id[b]["source_id"]
                or alen != by_id[a]["retained_bases"] or blen != by_id[b]["retained_bases"]):
            raise ValueError("Matrix values/source names/lengths differ from manifest")
        if value != 100 * (ai + bi) / (alen + blen):
            raise ValueError("Unrounded similarity differs from directional identity formula")
        if a > b:
            a, b, ai, bi, alen, blen = b, a, bi, ai, blen, alen
        pairs[key] = {"query": a, "subject": b, "similarity_percent_unrounded": value,
                      "query_nident": ai, "subject_nident": bi, "query_length": alen, "subject_length": blen}
    if len(pairs) != len(genomes) * (len(genomes) - 1) // 2:
        raise ValueError("Incomplete all-node phage matrix")
    return dict(sorted(pairs.items()))


def group_phages(matrix_text, genomes, model_ids, near_gp17_seeds, threshold=70.0, buffer=0.01):
    if (not math.isfinite(threshold) or not math.isfinite(buffer) or threshold < 0 or threshold > 100
            or buffer < 0 or buffer > threshold):
        raise ValueError("Invalid grouping threshold/buffer")
    pairs = parse_pairs(matrix_text, genomes)
    by_id = {g["id"]: g for g in genomes}
    sources = {g["source_id"] for g in genomes if g["role"] == "external_collection"}
    if (len(set(model_ids)) != len(model_ids) or not set(model_ids) <= sources
            or len(set(near_gp17_seeds)) != len(near_gp17_seeds) or not set(near_gp17_seeds) <= sources):
        raise ValueError("Unmapped or duplicate model-table/near-gp17 source identifier")
    parent = {node: node for node in by_id}

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    edges, boundary = [], []
    for pair in pairs.values():
        a, b, value = pair["query"], pair["subject"], pair["similarity_percent_unrounded"]
        if value >= threshold - buffer:
            parent[find(a)] = find(b)
            names = sorted((by_id[a]["source_id"], by_id[b]["source_id"]))
            edge = {"source_a": names[0], "source_b": names[1], "similarity_percent_unrounded": value}
            edges.append(edge)
            if value < threshold:
                boundary.append(edge)
    members = {}
    for node in sorted(by_id):
        members.setdefault(find(node), []).append(node)
    groups, node_groups = {}, {}
    for nodes in members.values():
        names = sorted(by_id[node]["source_id"] for node in nodes)
        membership_hash = digest(json.dumps(names, ensure_ascii=True, separators=(",", ":")).encode())
        group_id = "phage-group-" + membership_hash[:16]
        if group_id in groups:
            raise ValueError("Group identifier hash collision")
        targets = sorted(by_id[n]["source_id"] for n in nodes if by_id[n]["role"] == "target_reference")
        seeds = sorted(set(names) & set(near_gp17_seeds))
        reasons = ([{"kind": "target_containing_component", "seed_source_id": t} for t in targets]
                   + [{"kind": "near_gp17_seed_component", "seed_source_id": s} for s in seeds])
        groups[group_id] = {"source_members": names, "membership_sha256": membership_hash,
                            "member_genome_ids": sorted(nodes),
                            "external_source_members": sorted(set(names) & sources), "target_members": targets,
                            "model_table_members": sorted(set(names) & set(model_ids)),
                            "near_gp17_seed_members": seeds, "excluded": bool(reasons), "exclusion_reasons": reasons}
        node_groups.update({n: group_id for n in nodes})
    source_groups = {g["source_id"]: node_groups[g["id"]] for g in genomes if g["role"] == "external_collection"}
    target_groups = {g["source_id"]: node_groups[g["id"]] for g in genomes if g["role"] == "target_reference"}
    excluded = sorted(s for s, group in source_groups.items() if groups[group]["excluded"])
    model_groups = {s: source_groups[s] for s in sorted(model_ids)}
    return {"source_genomes": len(sources), "target_genomes": len(target_groups), "all_nodes": len(genomes),
            "pair_count": len(pairs), "group_count": len(groups),
            "model_table_phages": len(model_ids), "model_table_group_count": len(set(model_groups.values())),
            "retained_model_table_group_count": len({group for group in model_groups.values() if not groups[group]["excluded"]}),
            "excluded_group_count": sum(g["excluded"] for g in groups.values()),
            "similarity_percent": threshold, "boundary_buffer_percentage_points": buffer,
            "edge_rule": f"similarity >= {threshold - buffer}; unrounded undirected connected components",
            "all_nodes_retained_as_bridges": True,
            "edges_at_or_above_buffered_threshold": sorted(edges, key=lambda x: (x["source_a"], x["source_b"])),
            "boundary_edges": sorted(boundary, key=lambda x: (x["source_a"], x["source_b"])),
            "groups": dict(sorted(groups.items())), "source_phage_groups": dict(sorted(source_groups.items())),
            "target_phage_groups": dict(sorted(target_groups.items())), "model_phage_groups": model_groups,
            "near_gp17_seed_source_ids": sorted(near_gp17_seeds), "excluded_source_ids": excluded,
            "excluded_model_table_phages": sorted(set(excluded) & set(model_ids)),
            "retained_model_table_phages": sorted(set(model_ids) - set(excluded))}


def run(policy_path):
    policy_bytes = Path(policy_path).read_bytes()
    policy = json.loads(policy_bytes)
    # All policy evidence files are hash-bound, but phenotype tables are never parsed.
    raw = {}
    for item in policy["inputs"]:
        path = Path(item["path"])
        if path.name in raw:
            raise ValueError("Ambiguous policy input basename")
        raw[path.name] = checked(path, item["sha256"])
    manifest = json.loads(raw["genophi-phage-inputs.json"])
    source_map = json.loads(raw["genophi-training-source-map.json"])
    overlap = json.loads(raw["genophi-phage-overlap.json"])
    validation = json.loads(raw["viridic-original-r-comparison-v2.json"])
    phage_policy = json.loads(raw["genophi-phage-policy-v2.json"])
    manifest_hash = digest(raw["genophi-phage-inputs.json"])
    if (manifest["external_genomes"] != 97 or manifest["target_references"] != 2 or len(manifest["genomes"]) != 99
            or len(source_map["phages"]) != 94 or overlap["normalized_manifest_sha256"] != manifest_hash
            or overlap["policy_sha256"] != digest(raw["genophi-phage-policy-v2.json"])
            or phage_policy["normalized_input_manifest_sha256"] != manifest_hash):
        raise ValueError("Source roster/report/policy/manifest bindings differ")
    for name, item in source_map["phages"].items():
        if item["status"] != "exact_identifier_match" or item["file"]["id"] != name:
            raise ValueError("Model phage lacks exact authenticated source identifier")
    comparison = validation["comparison"]
    if (comparison["status"] != "passed" or not comparison["original_R_function_executed"]
            or not comparison["frozen_dataset_directional_equivalence_verified"]
            or comparison["unordered_pairs_compared"] != 4851
            or comparison["ordered_nonself_directions_compared"] != 9702 or comparison["discrepancies"]):
        raise ValueError("Original-R dataset validation incomplete")
    calculator_spec = validation["inputs"]["adapter_output"]
    calculator_bytes = checked(calculator_spec["path"], calculator_spec["sha256"])
    if overlap["calculator_output_sha256"] != calculator_spec["sha256"]:
        raise ValueError("Overlap and original-R validation refer to different calculators")
    pairs = parse_pairs(raw["genophi-phage-similarities.tsv"].decode(), manifest["genomes"])
    calculator_pairs = {}
    for pair in json.loads(calculator_bytes)["pairs"]:
        a, b = pair["query"], pair["subject"]
        if a >= b or (a, b) in calculator_pairs:
            raise ValueError("Unexpected/duplicate validated calculator pair")
        calculator_pairs[(a, b)] = {k: pair[k] for k in next(iter(pairs.values()))}
    if pairs != calculator_pairs:
        raise ValueError("Numeric matrix differs from original-R-validated calculator output")
    grouping = policy["phage_grouping"]
    whole = phage_policy["whole_genome_exclusion"]
    if any(grouping[k] != whole[k] for k in ["similarity_percent", "boundary_buffer_percentage_points"]):
        raise ValueError("Grouping thresholds differ from frozen source audit")
    result = group_phages(raw["genophi-phage-similarities.tsv"].decode(), manifest["genomes"], sorted(source_map["phages"]),
                          overlap["gp17"]["near_match_excluded_source_ids"], grouping["similarity_percent"],
                          grouping["boundary_buffer_percentage_points"])
    return {"schema_version": 1, "purpose": "Complete source phage relatedness groups for frozen source-roster admission",
            "policy_path": str(policy_path), "policy_sha256": digest(policy_bytes),
            "script_sha256": digest(Path(__file__).read_bytes()), "inputs": policy["inputs"],
            "validated_calculator_artifact": calculator_spec,
            "original_R_validation_scope": "Directional calculator agrees on this frozen HSP dataset; not universal pipeline equivalence",
            "original_R_validation_sha256": digest(raw["viridic-original-r-comparison-v2.json"]),
            "matrix_matches_original_R_validated_calculator": True,
            "group_id_rule": "phage-group- plus first16 SHA256 hex of compact JSON sorted complete source_members (targets included)",
            **result, "phenotype_tables_parsed": False, "training_or_scoring_executed": False,
            "source_representation_caveat": "Deposited representations retain the documented +60/+62bp source-length discrepancy; no trimming/substitution or stock-identity inference",
            "interpretation": "Relatedness groups and whole-component exclusions are conservative project policies; no taxonomy, binding, host-range or assembly inference. Final roster admission and production run specification remain separate."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    policy = json.loads(args.policy.read_bytes())
    protected = {args.policy.resolve(), Path(__file__).resolve()} | {Path(i["path"]).resolve() for i in policy["inputs"]}
    if args.output.exists() or args.output.resolve() in protected:
        parser.error("Output must be new and separate from inputs")
    report = run(args.policy)
    if args.output.resolve() == Path(report["validated_calculator_artifact"]["path"]).resolve():
        parser.error("Output cannot overwrite validated calculator input")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
