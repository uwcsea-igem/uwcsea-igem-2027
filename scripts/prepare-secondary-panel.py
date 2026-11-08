"""Freeze source-reference memberships and grouped blocks; never fit a model."""
import argparse
from collections import Counter
import csv
from dataclasses import asdict
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path

from tipcraft.learning import Sample
from tipcraft.partition import assign_group_folds, build_blocked_folds

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/secondary-panel-policy.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked(binding):
    path = ROOT / binding["path"]
    if sha(path) != binding["sha256"]:
        raise ValueError("Bound source checksum differs: " + str(path))
    return path


def unique(rows, key):
    values = {r[key]: r for r in rows}
    if len(values) != len(rows):
        raise ValueError("Duplicate source identifier: " + key)
    return values


def host_mapping(report, authority, authentication):
    part, summary = report["partition"], report["comparison_summary"]
    if (part["status"] != "complete_host_reference_grouping" or summary["status"] != "complete_reference_triage"
            or summary["failed_self_controls"] or summary["missing_directed_pairs"]
            or report["prior25_control_exact"] is not True):
        raise ValueError("Complete checked host grouping required")
    sources = authority["host_id_to_wgs_master_accession"]
    auth = unique(authentication["entries"], "host")
    if set(auth) != set(sources) or authentication["unresolved_hosts"]:
        raise ValueError("Host source roster incomplete")
    groups, seen, mapping = part["source_groups"], set(), {}
    for gid, group in groups.items():
        members = group["members"]
        if not members or len(set(members)) != len(members) or seen & set(members):
            raise ValueError("Duplicate/empty host group members")
        seen.update(members)
        target = any(m.startswith("original:") for m in members)
        source = sorted(m for m in members if m.startswith("ecor:"))
        if group["target_containing"] is not target or group["source_members"] != source:
            raise ValueError("Host target/source group flags inconsistent")
        for node in source:
            name = node.removeprefix("ecor:")
            if name not in sources or name in mapping:
                raise ValueError("Unknown or duplicate host node")
            row = auth[name]
            if row["source_wgs_accession"] != sources[name] or row["status"] != "authenticated_official_wgs_contig_set":
                raise ValueError("Host accession correspondence differs")
            mapping[name] = {"genotype_key": node, "group": gid, "excluded": target,
                             "source_accession": sources[name], "wgs_master": row["master"]["wgs_master"],
                             "wgs_build": row["master"]["wgs_build"], "fasta": row["fasta"]}
    if set(mapping) != set(sources) or summary["directed_pair_count"] != len(seen)**2:
        raise ValueError("Incomplete host group coverage or numeric comparison inventory")
    excluded = sorted(v["genotype_key"] for v in mapping.values() if v["excluded"])
    candidates = sorted(v["genotype_key"] for v in mapping.values() if not v["excluded"])
    if (part["target_related_source_hosts"] != excluded or part["candidate_source_hosts"] != candidates
            or part["candidate_source_group_count"] != len({v["group"] for v in mapping.values() if not v["excluded"]})):
        raise ValueError("Host exclusion summaries disagree")
    return mapping


def phage_mapping(report, authority, replay):
    if (report["status"] != "completed_original_R_checked_official_variant_groups"
            or report["numeric_matrix_for_grouping_verified"] is not True
            or report["numeric_validation"]["dataset_specific_directional_equivalence_verified"] is not True):
        raise ValueError("Completed original-R-checked official variant groups required")
    roster = authority["phage_id_to_accession"]
    observations = unique(replay["observations"], "phage")
    if (set(observations) != set(roster) or replay["official_record_format_pairs_verified"] != len(roster)):
        raise ValueError("Official accession format inventory incomplete")
    expected, mapping = {}, {}
    for name, row in observations.items():
        if (row["requested_accession"] != roster[name] or row["official_formats_agree"] is not True
                or not row["official_versioned_accession"].startswith(roster[name] + ".")):
            raise ValueError("Official accession correspondence differs")
        expected[name] = ("source:" + name if row["exact_author_DNA_match"] else
                          "official:" + name + ":" + row["official_versioned_accession"])
    if report["official_panel_node_map"] != expected:
        raise ValueError("Exact official panel node map required; no name substitution")
    containers = unique(report["containers"], "container_key")
    grouping = report["grouping"]
    if grouping["all_nodes_retained_as_bridges"] is not True:
        raise ValueError("All nonmodel bridge nodes must remain in grouping")
    members = {}
    for gid, group in grouping["groups"].items():
        keys = group["members"]
        if not keys or len(keys) != len(set(keys)):
            raise ValueError("Duplicate/empty phage group")
        if set(keys) - set(containers) or set(keys) & set(members):
            raise ValueError("Unknown or multiply grouped phage node")
        target = sorted(k for k in keys if containers[k]["role"] == "target_reference")
        near = group["direct_gp17_seeds"]
        if (set(near) - set(keys) or group["target_seeds"] != target
                or group["target_only_excluded"] is not bool(target)
                or group["primary_excluded"] is not bool(target or near)):
            raise ValueError("Phage exclusion flags disagree with seeds")
        members.update({key: gid for key in keys})
    if set(members) != set(containers) or grouping["container_groups"] != members:
        raise ValueError("Complete consistent phage group mapping required")
    for name, key in sorted(expected.items()):
        row = observations[name]
        if key not in members or containers[key]["source_id"] != name or containers[key]["length"] != row["official_length"]:
            raise ValueError("Official genotype node missing or length differs")
        gid = members[key]
        mapping[name] = {"genotype_key": key, "group": gid, "excluded": grouping["groups"][gid]["primary_excluded"],
                         "source_accession": roster[name], "versioned_accession": row["official_versioned_accession"],
                         "sequence_sha256": row["official_sequence_sha256"],
                         "exact_author_DNA_match": row["exact_author_DNA_match"]}
    return mapping


def parse_observations(rows, authority):
    hosts, phages = authority["host_id_to_wgs_master_accession"], authority["phage_id_to_accession"]
    pairs = {}
    for row in rows:
        if set(row) != {"strain", "phage", "interaction"}:
            raise ValueError("Exact source observation columns required")
        host, phage = row["strain"], row["phage"]
        if host not in hosts or phage not in phages or (host, phage) in pairs:
            raise ValueError("Unknown/duplicate source observation")
        try:
            value = Decimal(row["interaction"])
        except (InvalidOperation, TypeError) as exc:
            raise ValueError("Invalid source interaction") from exc
        if not value.is_finite() or value not in (0, 1):
            raise ValueError("Require finite binary observed labels; missing pairs remain unknown")
        pairs[host, phage] = int(value)
    if (len(pairs) != authority["panel_pairs"] or {h for h, p in pairs} != set(hosts)
            or {p for h, p in pairs} != set(phages)):
        raise ValueError("Entire published observation roster required")
    return pairs


def build_panel(rows, authority, host_report, phage_report, host_auth, phage_replay, design):
    hosts = host_mapping(host_report, authority, host_auth)
    phages = phage_mapping(phage_report, authority, phage_replay)
    pairs = parse_observations(rows, authority)
    samples, memberships, labels = [], [], {}
    for (host, phage), label in sorted(pairs.items()):
        row_id = "basel-ecor:" + host + ":" + phage
        reasons = (["target_related_host_component"] if hosts[host]["excluded"] else []) + (
            ["target_or_direct_gp17_seed_phage_component"] if phages[phage]["excluded"] else [])
        sample = Sample(row_id, host, phage, hosts[host]["group"], phages[phage]["group"], authority["assay_source"])
        memberships.append({**asdict(sample), "eligible": not reasons, "exclusion_reasons": reasons,
                            "host_genotype_key": hosts[host]["genotype_key"], "phage_genotype_key": phages[phage]["genotype_key"]})
        labels[row_id] = label
        if not reasons:
            samples.append(sample)
    retained_hosts = sorted(h for h, v in hosts.items() if not v["excluded"])
    retained_phages = sorted(p for p, v in phages.items() if not v["excluded"])
    if ({s.host_id for s in samples} != set(retained_hosts) or {s.phage_id for s in samples} != set(retained_phages)):
        raise ValueError("Eligible source genome lacks observations; no silent roster reduction")
    settings = design["partition"]
    if (settings["host_folds"], settings["phage_folds"], settings["seed"], settings["regime"], settings["mixed_quadrants"]) != (3, 3, 2027, "both_unseen", "embargo"):
        raise ValueError("Require frozen 3x3 both-unseen partition design")
    hf = assign_group_folds([s.host_group for s in samples], axis="host", folds=3, seed=2027)
    pf = assign_group_folds([s.phage_group for s in samples], axis="phage", folds=3, seed=2027)
    blocks = build_blocked_folds(samples, hf, pf)  # Labels never enter assignment/partition APIs.
    diagnostics, invalid = [], []
    by_id = {s.row_id: s for s in samples}
    for block in blocks:
        train = Counter(labels[i] for i in block["train_ids"])
        by_phage = {}
        for row_id in block["test_ids"]:
            by_phage.setdefault(by_id[row_id].phage_id, Counter())[labels[row_id]] += 1
        informative = sorted(p for p, counts in by_phage.items() if counts[0] and counts[1])
        diagnostics.append({"fold_id": block["fold_id"], "training_label_counts": {str(k): train[k] for k in [0, 1]},
                            "informative_test_phages": informative, "informative_test_phage_count": len(informative),
                            "positive_negative_test_comparisons": sum(c[0]*c[1] for c in by_phage.values())})
        if set(train) != {0, 1}:
            invalid.append(block["fold_id"])
    unknown = [{"host_id": h, "phage_id": p, "status": "unknown_not_negative",
                "eligible_genomes": not hosts[h]["excluded"] and not phages[p]["excluded"]}
               for h in sorted(hosts) for p in sorted(phages) if (h, p) not in pairs]
    return {"status": "frozen_secondary_reference_panel" if not invalid else "HOLD_single_class_training_blocks",
            "host_genomes": hosts, "phage_genomes": phages, "memberships": memberships,
            "samples": [asdict(s) for s in samples], "labels": labels, "unknown_pairs": unknown,
            "host_group_folds": hf, "phage_group_folds": pf, "blocks": blocks,
            "post_partition_label_diagnostics": diagnostics, "invalid_training_blocks": invalid,
            "counts": {"source_hosts": len(hosts), "source_phages": len(phages), "source_observed_rows": len(pairs),
                       "source_unknown_pairs": len(unknown), "eligible_hosts": len(retained_hosts),
                       "eligible_phages": len(retained_phages), "eligible_observed_rows": len(samples),
                       "excluded_observed_rows": len(pairs)-len(samples), "host_groups": len(hf), "phage_groups": len(pf)},
            "label_independent_partitioning": True, "diagnostics_used_for_selection": False,
            "experimental_stock_identity_verified": False, "training_authorized": False,
            "model_fitted_or_scored": False, "new_original_benchmark_coverage": 0,
            "interpretation": "Published source-ID labels associated with authenticated official reference genotypes; physical experimental stocks unverified. Spot-clearance endpoint; no productive-infection or original T7/T3 performance claim."}


def join_annotations(result, manifest):
    entries = {}
    for row in manifest["entries"]:
        item = row["item"]
        if item["role"] == "technical_control":
            continue
        key = (item["kind"], item["source_id"])
        if item["role"] != "secondary_reference_not_admitted" or key in entries:
            raise ValueError("Unknown/duplicate annotation source role")
        entries[key] = row
    expected = {(kind, name) for kind, field in [("host", "host_genomes"), ("phage", "phage_genomes")]
                for name in result[field]}
    if set(entries) != expected:
        raise ValueError("Exact 77-reference annotation join required; controls are not source members")
    for kind, field in [("host", "host_genomes"), ("phage", "phage_genomes")]:
        for name, genome in result[field].items():
            row = entries[kind, name]
            item = row["item"]
            if row["status"] != "annotation_validated" or row["FAA_GFF_inventory_and_coordinates_agree"] is not True or row["proteins"] < 1:
                raise ValueError("Complete validated source annotation required")
            if kind == "host":
                if item["wgs_accession"] != genome["source_accession"] or item["fasta"]["sha256"] != genome["fasta"]["sha256"]:
                    raise ValueError("Host annotation genotype correspondence differs")
            elif (item["versioned_accession"] != genome["versioned_accession"]
                  or item["expected_sequence_sha256"] != genome["sequence_sha256"]):
                raise ValueError("Phage annotation genotype correspondence differs")
            for binding in [row["input"], row["contig_map"], row["protein_coordinate_inventory"], row["outputs"]["FAA"], row["outputs"]["GFF"]]:
                checked(binding)
            genome["annotation"] = {"id": item["id"], "proteins": row["proteins"],
                                     "FAA": row["outputs"]["FAA"], "GFF": row["outputs"]["GFF"],
                                     "normalized_input": row["input"], "contig_map": row["contig_map"],
                                     "protein_coordinate_inventory": row["protein_coordinate_inventory"],
                                     "proteome_multiset_sha256": row["proteome_multiset_sha256"]}
    result["annotation_join_complete"] = True
    result["ready_for_feature_preparation"] = not result["invalid_training_blocks"]
    result["feature_vocabularies_or_models_prepared"] = False


def prepare(policy_path, expected_phage_sha=None, check_only=False):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Panel helper differs from frozen policy")
    paths = {name: checked(binding) for name, binding in policy["inputs"].items()}
    design = json.loads(paths["design"].read_text())
    for binding in design["bound_design_inputs"]:
        checked(binding)
    future = ROOT / policy["future_phage_groups"]["path"]
    if not future.exists() or not expected_phage_sha:
        return {"status": "HOLD_missing_prerequisite", "missing": [
            name for name, absent in [(str(future.relative_to(ROOT)), not future.exists()),
                                      ("externally_frozen_phage_group_report_sha256", not expected_phage_sha)] if absent],
                "membership_generated": False, "training_authorized": False}
    checked({"path": str(future), "sha256": expected_phage_sha})
    phage = json.loads(future.read_text())
    if (phage["wrapper"]["sha256"] != policy["inputs"]["phage_group_helper"]["sha256"]
            or phage["policy_sha256"] not in [policy["inputs"][k]["sha256"] for k in policy["accepted_group_policies"]]):
        raise ValueError("Future phage grouping producer/policy differs")
    if (len(phage["containers"]) != 266 or phage["plan"]["expected_unordered_pairs"] != 35245
            or sorted(c["container_key"] for c in phage["containers"] if c["role"] == "target_reference") != ["reference:T3", "reference:T7"]):
        raise ValueError("Complete common-context official-variant grouping required")
    for binding in phage["inputs"].values():
        checked(binding)
    checked(phage["numeric_matrix"])
    host = json.loads(paths["host_groups"].read_text())
    replay = json.loads(paths["host_groups_replay"].read_text())
    if (replay["status"] != "exact_complete_matrix_and_partition_replay"
            or replay["report"]["sha256"] != policy["inputs"]["host_groups"]["sha256"]):
        raise ValueError("Host grouping replay differs")
    checked(replay["matrix"])
    authority = json.loads(paths["authority"].read_text())
    if (authority["panel_pairs"], len(authority["host_id_to_wgs_master_accession"]),
            len(authority["phage_id_to_accession"])) != (1240, 25, 52):
        raise ValueError("Full 1240-row/25-host/52-phage authority required")
    with paths["observations"].open(newline="") as stream:
        result = build_panel(list(csv.DictReader(stream)), authority, host, phage,
                             json.loads(paths["host_authentication"].read_text()),
                             json.loads(paths["phage_accession_replay"].read_text()), design)
    if result["counts"]["source_unknown_pairs"] != 60:
        raise ValueError("Expected 60 unknown source pairs")
    annotation = json.loads(paths["annotation_report"].read_text())
    annotation_replay = json.loads(paths["annotation_replay"].read_text())
    if (annotation["annotation_manifest"] != policy["inputs"]["annotation_manifest"]
            or annotation_replay["input_report"] != policy["inputs"]["annotation_report"]
            or annotation_replay["status"] != "all77_secondary_reference_annotations_replayed_offline"):
        raise ValueError("Annotation replay/manifest bindings differ")
    join_annotations(result, json.loads(paths["annotation_manifest"].read_text()))
    for binding in policy["inputs"].values():
        checked(binding)
    checked({"path": str(future), "sha256": expected_phage_sha})
    if policy_path.read_bytes() != policy_bytes:
        raise ValueError("Panel policy changed during preparation")
    result.update({"schema_version": 1, "helper_sha256": policy["helper_sha256"], "policy_sha256": hashlib.sha256(policy_bytes).hexdigest(),
                   "membership_generated": True,
                   "inputs": {**policy["inputs"], "phage_groups": {"path": str(future.relative_to(ROOT)), "sha256": expected_phage_sha}}})
    if check_only:
        return {"status": result["status"], "counts": result["counts"], "membership_generated": False,
                "training_authorized": False, "check_only": True}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--expected-phage-groups-sha256")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.output is None and not args.check_only:
        parser.error("--output required outside --check-only")
    if args.output and args.output.exists():
        parser.error("Use a fresh output; existing artifacts are preserved")
    policy = json.loads(args.policy.read_text())
    if args.output and args.output.resolve() == (ROOT / policy["future_phage_groups"]["path"]).resolve():
        parser.error("Panel output cannot replace the future grouping prerequisite")
    result = prepare(args.policy, args.expected_phage_groups_sha256, args.check_only)
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    print(json.dumps({k: result[k] for k in ["status", "counts", "missing", "membership_generated"] if k in result}))
    return 0 if result["status"] == "frozen_secondary_reference_panel" else 2


if __name__ == "__main__":
    raise SystemExit(main())
