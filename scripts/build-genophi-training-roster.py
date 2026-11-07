"""Admit authenticated external source rows, or write an explicit HOLD report.

No target phenotypes, model fitting, feature selection or prediction are inputs
to this stage. A ready source roster is not a complete production run spec.
"""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODEL_TABLE = ROOT / "data/raw/genophi-audit/data/interaction_matrices/ecoli_interaction_matrix.csv"
TARGET_MANIFEST = ROOT / "software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json"
SKANI_SHA256 = "1f1d08abe20d592522faf3f2454634e815269f61e37bb0adda4608599e6a4efb"


def module(filename):
    spec = importlib.util.spec_from_file_location(filename, Path(__file__).with_name(filename))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def pinned_inputs(policy):
    result = {}
    for item in policy["inputs"]:
        path = ROOT / item["path"]
        if sha(path) != item["sha256"]:
            raise ValueError(f"Frozen input hash mismatch: {item['path']}")
        if path.name in result:
            raise ValueError("Duplicate frozen input basename")
        result[path.name] = path
    split = policy["group_splits"]
    if ((split["regime"], split["host_folds"], split["phage_folds"], split["seed"])
            != ("both_unseen", 3, 3, 2027) or policy["production_run_spec_complete"] is not False):
        raise ValueError("Unsupported frozen roster/split specification")
    return result


def crosswalk(rows, audit, source_map, inventory):
    by_source = {item["source_id"]: item for item in inventory if item["kind"] == "host"}
    if len(by_source) != sum(item["kind"] == "host" for item in inventory):
        raise ValueError("Duplicate deposited host source identifier")
    names = audit["crosschecked_host_correspondence"]
    model_hosts, model_phages = {r["strain"] for r in rows}, {r["phage"] for r in rows}
    if set(names) != model_hosts or set(source_map["hosts"]) != model_hosts or set(source_map["phages"]) != model_phages:
        raise ValueError("Source-map/label-audit identifiers disagree with the model table")
    if len(set(names.values())) != len(names) or set(names.values()) - set(by_source):
        raise ValueError("Original-source ID must map exactly and uniquely to a deposited source ID")
    return {name: {"original_source_id": source, "deposit_file_id": by_source[source]["deposit_file_id"],
                   "path": by_source[source]["path"], "stock_identity_verified": False}
            for name, source in sorted(names.items())}


def propagate_host_exclusions(core, seeds):
    source_groups = {}
    for group, info in core["groups"].items():
        for source in info["source_members"]:
            if source in source_groups:
                raise ValueError("A source host belongs to multiple core groups")
            source_groups[source] = group
    if set(seeds) - set(source_groups):
        raise ValueError("Host exclusion seed is absent from full core-distance groups")
    excluded = {}
    for source, reasons in sorted(seeds.items()):
        excluded.setdefault(source_groups[source], {})[source] = sorted(set(reasons))
    return {"excluded_groups": excluded,
            "excluded_source_ids": sorted(source for source, group in source_groups.items() if group in excluded),
            "source_groups": source_groups}


def curate_rows(rows, label_audit, name_crosswalk, core, host_exclusions, phage_groups, excluded_phage_groups):
    blanks = {(r["strain"], r["phage"]) for r in label_audit["source_blank_pairs"]}
    conflicts = {(r["strain"], r["phage"]) for r in label_audit["conflicting_nonblank_pairs"]}
    seen, result = set(), []
    for row in rows:
        pair = (row["strain"], row["phage"])
        if pair in seen or row["interaction"] not in ("0", "1"):
            raise ValueError("Duplicate or nonbinary source interaction row")
        seen.add(pair)
        source = name_crosswalk.get(pair[0], {}).get("original_source_id")
        hg, pg = core["training_host_groups"].get(pair[0]), phage_groups.get(pair[1])
        if source and hg is not None and host_exclusions["source_groups"].get(source) != hg:
            raise ValueError("Training host group disagrees with its original-source membership")
        reasons = []
        if pair in blanks:
            reasons.append("original_source_empty")
        if pair in conflicts:
            reasons.append("nonblank_label_conflict")
        if not source or hg is None:
            reasons.append("unmapped_host_input_or_group")
        if pg is None:
            reasons.append("unmapped_phage_input_or_group")
        if hg in host_exclusions["excluded_groups"]:
            reasons.append("excluded_host_core_group")
        if pg in excluded_phage_groups:
            reasons.append("excluded_phage_component")
        identifier = hashlib.sha256(json.dumps(pair, separators=(",", ":")).encode()).hexdigest()[:24]
        result.append({"row_id": "source-row-" + identifier, "host_id": pair[0], "phage_id": pair[1],
                       "original_source_host_id": source, "label": int(row["interaction"]),
                       "host_group": hg, "phage_group": pg,
                       "publication": "Gaborieau_collection_doi_10.1038_s41564-024-01832-5",
                       "exclusion_reasons": reasons, "retained": not reasons})
    if (blanks | conflicts) - seen:
        raise ValueError("Label exclusions contain pairs absent from the source table")
    if len({r["row_id"] for r in result}) != len(result):
        raise ValueError("Stable row identifier collision")
    return sorted(result, key=lambda r: r["row_id"])


def label_summary(rows):
    retained = [r for r in rows if r["retained"]]
    return {"candidate_rows": len(rows), "retained_candidate_rows": len(retained),
            "excluded_candidate_rows": len(rows) - len(retained),
            "reason_counts_nonexclusive": dict(sorted(Counter(reason for r in rows for reason in r["exclusion_reasons"]).items())),
            "retained_candidate_label_counts": dict(sorted(Counter(str(r["label"]) for r in retained).items()))}


def verify_hosts(inventory, genomes, prepare):
    hashes, missing = {}, []
    for item in inventory:
        if item["kind"] != "host":
            continue
        path = genomes / item["path"]
        if not path.exists():
            if path.is_symlink():
                raise ValueError(f"Dangling host input symlink: {item['source_id']}")
            missing.append(item["path"])
        else:
            hashes[item["source_id"]] = prepare.verify(path, item)["file_sha256"]
    return hashes, sorted(missing)


def verify_phages(inventory, genomes, target_genomes, manifest, provenance, target_manifest):
    normalizer = module("normalize-genophi-phages.py")
    if (manifest["source_provenance"] != provenance or manifest["target_manifest_sha256"] != sha(target_manifest)
            or manifest["normalizer_sha256"] != sha(Path(normalizer.__file__))):
        raise ValueError("Normalized phage provenance/target/normalizer bindings differ")
    by_source = {item["source_id"]: item for item in manifest["genomes"]}
    if len(by_source) != len(manifest["genomes"]):
        raise ValueError("Duplicate phage manifest source identifiers")
    hashes, missing = {}, []
    for item in inventory:
        if item["kind"] != "phage":
            continue
        path = genomes / item["path"]
        if not path.exists():
            missing.append(str(path))
            continue
        data = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if len(data) != item["expected_bytes"] or blob != item["expected_git_blob_sha1"]:
            raise ValueError(f"Deposited phage bytes differ: {item['source_id']}")
        hashes[item["source_id"]] = sha(path)
        record = by_source[item["source_id"]]
        sequence, changes = normalizer.parse_complete_genome(data, normalizer.REVIEWED_EMPTY_RECORDS.get(sha(path)))
        if (record["source_file_sha256"] != sha(path) or record["git_blob_sha1"] != blob
                or record["sequence_sha256"] != hashlib.sha256(sequence.encode()).hexdigest()
                or any(record[key] != value for key, value in changes.items())):
            raise ValueError("Actual deposited phage differs from normalized audit input")
    for record in manifest["genomes"]:
        if record["role"] == "target_reference":
            path = target_genomes / "phages" / f"{record['source_id']}.fna"
            if not path.exists():
                missing.append(str(path))
            elif sha(path) != record["source_file_sha256"]:
                raise ValueError("Target phage file differs from frozen relatedness input")
    return hashes, sorted(missing)


def verify_phage_groups(path, policy_path):
    report = read_json(path)
    replay = module("group-genophi-phages.py").run(policy_path)
    if (Path(report["policy_path"]).resolve() != policy_path.resolve()
            or {k: v for k, v in report.items() if k != "policy_path"}
            != {k: v for k, v in replay.items() if k != "policy_path"}):
        raise ValueError("Phage component report differs from replayed authenticated evidence")
    return report


def verify_host_audit(directory, inventory, provenance, hashes, genomes, target_genomes, target_manifest, policy_path):
    audit = module("audit-genophi-host-overlap.py")
    expected_hosts = {item["source_id"] for item in inventory if item["kind"] == "host"}
    if set(hashes) != expected_hosts:
        raise ValueError("Host audit admission requires every authenticated source-host hash")
    report_path, matrix = directory / "report.json", directory / "skani.tsv"
    report, policy = read_json(report_path), read_json(policy_path)
    manifest = read_json(target_manifest)
    target_hashes = {}
    for name in manifest["hosts"]:
        path = target_genomes / "hosts" / f"{name}.fna"
        stats = module("prepare-genophi-genomes.py").fasta_stats(path.read_bytes())
        if sha(path) != manifest["genome_sources"][name]["file_sha256"] or stats["total_bases"] != manifest["genome_sources"][name]["length"]:
            raise ValueError(f"Frozen target reference differs: {name}")
        target_hashes[name] = sha(path)
    bindings = {"method": "skani 0.4.0", "wrapper_sha256": sha(Path(audit.__file__)),
                "policy_sha256": sha(policy_path), "target_manifest_sha256": sha(target_manifest),
                "source_provenance": provenance, "input_file_sha256": hashes,
                "target_file_sha256": target_hashes, "matrix_sha256": sha(matrix),
                "executable_sha256": SKANI_SHA256, "no_phenotype_labels_loaded": True,
                "infection_or_assembly_scores_produced": False}
    for key, expected in bindings.items():
        if report.get(key) != expected:
            raise ValueError(f"Host audit binding mismatch: {key}")
    hosts = {str((genomes / item["path"]).resolve()): item["source_id"] for item in inventory if item["kind"] == "host"}
    targets = {str((target_genomes / "hosts" / f"{name}.fna").resolve()): name for name in manifest["hosts"]}
    def references(path, expected):
        names = path.read_text().splitlines()
        if len(names) != len(expected) or len(set(names)) != len(names) or {str(Path(n).resolve()) for n in names} != set(expected):
            raise ValueError("Host audit input list does not cover exactly the authenticated genomes")
        return {name: expected[str(Path(name).resolve())] for name in names}
    refs = references(directory / "references.txt", hosts)
    queries = references(directory / "queries.txt", targets)
    command = report["command"]
    params = policy["host_parameters"]
    expected_tail = ["-c", str(params["compression"]), "-m", str(params["marker_compression"]),
                     "-s", str(params["screen_percent"]), "--min-af", str(params["minimum_reported_af_percent"]),
                     "-t", str(params["threads"])]
    if (len(command) != 18 or command[1:3] != ["dist", "--rl"] or command[4] != "--ql"
            or command[6:16] != expected_tail or command[16] != "-o"
            or Path(command[3]).resolve() != (directory / "references.txt").resolve()
            or Path(command[5]).resolve() != (directory / "queries.txt").resolve()
            or Path(command[17]).resolve() != matrix.resolve() or sha(Path(command[0])) != SKANI_SHA256):
        raise ValueError("Host audit command/executable differs from fixed parameters or actual inputs")
    with matrix.open() as stream:
        recomputed = audit.summarize(list(csv.DictReader(stream, delimiter="\t")), refs, queries, policy["host_exclusion"])
    if any(report.get(key) != value for key, value in recomputed.items()):
        raise ValueError("Host audit flags disagree with recomputed raw matrix summary")
    seeds = {}
    for host, info in recomputed["hosts"].items():
        reasons = (["named_target_identifier"] if info["named_target_flag"] else [])
        reasons.extend(target + ":" + pair["status"] for target, pair in info["comparisons"].items()
                       if pair["status"] != "outside_declared_near_target_rule")
        if reasons:
            seeds[host] = reasons
    return seeds, {"report_sha256": sha(report_path), "matrix_sha256": sha(matrix), **bindings,
                   "reported_pairs": recomputed["reported_pairs"], "expected_pairs": recomputed["expected_pairs"]}


def partitions(rows, policy):
    from tipcraft.learning import Sample
    from tipcraft.partition import assign_group_folds, build_blocked_folds
    from tipcraft import partition
    retained = [r for r in rows if r["retained"]]
    samples = [Sample(**{key: r[key] for key in ("row_id", "host_id", "phage_id", "host_group", "phage_group", "publication")}) for r in retained]
    split = policy["group_splits"]
    hf = assign_group_folds([s.host_group for s in samples], axis="host", folds=split["host_folds"], seed=split["seed"])
    pf = assign_group_folds([s.phage_group for s in samples], axis="phage", folds=split["phage_folds"], seed=split["seed"])
    blocks = build_blocked_folds(samples, hf, pf)
    labels = {r["row_id"]: r["label"] for r in retained}
    for block in blocks:
        for role in ("train", "test"):
            counts = Counter(str(labels[row_id]) for row_id in block[role + "_ids"])
            block[role + "_label_counts"] = dict(sorted(counts.items()))
            block["single_class_" + ("training" if role == "train" else "test")] = len(counts) == 1
        block["single_class_fit_must_abstain"] = block["single_class_training"]
    return {"host_group_folds": hf, "phage_group_folds": pf, "folds": blocks,
            "partition_source_sha256": sha(Path(partition.__file__)), "random_fallback_used": False,
            "publication_overlap": "Shared source publication; these are not leave-publication-out folds"}


def run(args):
    from tipcraft import partition
    policy = read_json(args.policy)
    report = {"schema_version": 1, "status": "HOLD", "policy_sha256": sha(args.policy),
              "builder_sha256": sha(Path(__file__)), "production_run_spec_complete": False,
              "partition_source_sha256": sha(Path(partition.__file__)),
              "training_or_scoring_executed": False, "target_scores_generated": False,
              "final_roster_written": False, "hold_reasons": [],
              "interpretation": "Source ID correspondence does not establish stock identity. Admission is separate from feature preparation and fitting."}
    try:
        pinned = pinned_inputs(policy)
        report["authenticated_policy_inputs"] = policy["inputs"]
        prepare = module("prepare-genophi-genomes.py")
        inventory, provenance = prepare.load_inventory(args.metadata_cache)
        label = read_json(pinned["genophi-label-audit.json"])
        if module("audit-genophi-labels.py").run(args.metadata_cache / "original-interaction-matrix.csv", MODEL_TABLE) != label:
            raise ValueError("Replayed source-label audit differs from frozen evidence")
        core = read_json(pinned["genophi-host-core-groups.json"])
        if module("group-genophi-hosts.py").run(args.metadata_cache / "original-host-core-distances.tsv", MODEL_TABLE) != core:
            raise ValueError("Replayed host core grouping differs from frozen evidence")
        source_map = read_json(pinned["genophi-training-source-map.json"])
        if any(source_map[key] != provenance[key] for key in ("source_commit", "source_tree_sha1", "reviewed_host_file_inventory_sha256")):
            raise ValueError("Source map is bound to different deposited source metadata")
        if source_map["genophi_table_sha256"] != label["genophi_table_sha256"] or sha(MODEL_TABLE) != label["genophi_table_sha256"]:
            raise ValueError("Source-map/label-audit model-table bindings differ")
        with MODEL_TABLE.open() as stream:
            rows = list(csv.DictReader(stream))
        names = crosswalk(rows, label, source_map, inventory)
        report.update({"source_provenance": provenance, "model_table_sha256": sha(MODEL_TABLE),
                       "name_crosswalk": names, "exact_original_source_id_matches": len(names)})
        hashes, missing = verify_hosts(inventory, args.genomes, prepare)
        report["host_inputs"] = {"expected": sum(x["kind"] == "host" for x in inventory), "verified": len(hashes),
                                 "missing_paths": missing, "verified_file_sha256": hashes}
        if missing:
            report["hold_reasons"].append("Missing authenticated host genome files")
        report["candidate_source_label_counts"] = {"all_rows": len(rows), "source_empty_excluded": len(label["source_blank_pairs"]),
                                                   "nonblank_conflicts": len(label["conflicting_nonblank_pairs"]),
                                                   "consistent_nonblank_rows": label["concordant_nonblank_pairs"]}
        phages = verify_phage_groups(args.phage_groups, args.policy)
        report["phage_groups_sha256"] = sha(args.phage_groups)
        report["phage_group_builder_sha256"] = phages["script_sha256"]
        phage_hashes, missing_phages = verify_phages(inventory, args.genomes, args.target_genomes,
                                                   read_json(pinned["genophi-phage-inputs.json"]), provenance, TARGET_MANIFEST)
        report["phage_inputs"] = {"verified_source_files": len(phage_hashes), "missing_paths": missing_phages,
                                  "verified_source_file_sha256": phage_hashes}
        if missing_phages:
            report["hold_reasons"].append("Missing authenticated phage/reference genome files")
        overlap_policy = pinned["genophi-overlap-policy.json"]
        named_targets = read_json(overlap_policy)["host_exclusion"]["named_target_identifiers"]
        all_core_sources = {source for group in core["groups"].values() for source in group["source_members"]}
        seeds = {name: ["named_target_identifier"] for name in named_targets if name in all_core_sources}
        if args.host_audit_dir is None:
            report["hold_reasons"].append("Complete target-host audit has not been supplied")
        elif missing:
            report["hold_reasons"].append("Supplied host audit cannot be admitted until every host file is authenticated")
        else:
            audit_seeds, bindings = verify_host_audit(args.host_audit_dir, inventory, provenance, hashes,
                                                     args.genomes, args.target_genomes, TARGET_MANIFEST, overlap_policy)
            seeds.update(audit_seeds)
            report["host_audit_bindings"] = bindings
        host_exclusions = propagate_host_exclusions(core, seeds)
        excluded_pg = {group for group, info in phages["groups"].items() if info["excluded"]}
        candidates = curate_rows(rows, label, names, core, host_exclusions, phages["model_phage_groups"], excluded_pg)
        report["candidate_row_summary"] = {"provisional": bool(report["hold_reasons"]), **label_summary(candidates)}
        report["host_group_exclusions"] = host_exclusions
        report["phage_group_exclusions"] = {group: phages["groups"][group]["exclusion_reasons"] for group in sorted(excluded_pg)}
        if report["hold_reasons"]:
            return report, None, None
        try:
            folds = partitions(candidates, policy)
        except ValueError as exc:
            report["hold_reasons"].append("Frozen deterministic split is infeasible: " + str(exc))
            report["candidate_row_summary"]["provisional"] = True
            return report, None, None
        report.update({"status": "ready_for_feature_preparation", "source_admission_checks_passed": True,
                       "final_roster_written": True, "retained_rows": sum(r["retained"] for r in candidates),
                       "excluded_rows": sum(not r["retained"] for r in candidates),
                       "partition_source_sha256": folds["partition_source_sha256"]})
        return report, candidates, folds
    except FileNotFoundError as exc:
        report["hold_reasons"].append("Required input missing: " + str(exc.filename))
        return report, None, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("policy", "metadata-cache", "genomes", "target-genomes", "phage-groups", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--host-audit-dir", type=Path)
    args = parser.parse_args()
    inputs = [args.metadata_cache, args.genomes, args.target_genomes, args.policy, args.phage_groups]
    if args.host_audit_dir:
        inputs.append(args.host_audit_dir)
    if args.output_dir.exists() or any(args.output_dir.resolve().is_relative_to(path.resolve())
                                       or path.resolve().is_relative_to(args.output_dir.resolve()) for path in inputs):
        parser.error("Output must be a fresh directory separate from every input")
    try:
        report, roster, folds = run(args)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    args.output_dir.mkdir(parents=True, exist_ok=False)
    if roster is not None:
        (args.output_dir / "rows.json").write_text(json.dumps(roster, indent=2) + "\n")
        (args.output_dir / "partitions.json").write_text(json.dumps(folds, indent=2) + "\n")
        report["rows_sha256"] = sha(args.output_dir / "rows.json")
        report["partitions_sha256"] = sha(args.output_dir / "partitions.json")
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["status"] == "ready_for_feature_preparation" else 1


if __name__ == "__main__":
    raise SystemExit(main())
