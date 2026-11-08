"""Direct gp17 screening of an immutable, augmented nucleotide database.

Reuse its completed database construction read-only; do not run nucleotide
alignments, propagate graph exclusions, load phenotypes or fit a model.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

POLICY = Path("software/tipcraft/benchmarks/receptor-resource/official-variant-gp17-policy.json")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def corpus(mapping, sequences, old_mapping, old_sequences):
    records = mapping["records"]
    if (len(records) != 273 or len(sequences) != 273 or len(old_sequences) != 270
            or records[:270] != old_mapping["records"] or sequences[:270] != old_sequences
            or len({r["id"] for r in records}) != 273
            or len({r["container_key"] for r in records}) != 266):
        raise ValueError("Require all unchanged 270 records plus three distinct official records")
    official = records[270:]
    if ({r["source_id"] for r in official} != {"Bas06", "Bas21", "Bas68"}
            or any(r["role"] != "official_reference_version"
                   or r["container_key"] != f"official:{r['source_id']}:{r['source_record_id']}" for r in official)):
        raise ValueError("Official variant identity/namespace differs")
    for row, (identifier, sequence) in zip(records, sequences):
        if identifier != row["id"] or len(sequence) != row["retained_bases"]:
            raise ValueError("Subject record order/length differs from map")
        if "sequence_sha256" in row and hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
            raise ValueError("Mapped nucleotide sequence hash differs")
    targets = [r for r in records if r["role"] == "target_reference"]
    if {(r["id"], r["container_key"]) for r in targets} != {("t000", "reference:T7"), ("t001", "reference:T3")}:
        raise ValueError("Original controls differ")
    return records


def screen(rows, records, policy, lengths, logic):
    """Retain every namespace and record; validate with the frozen gp17 parser."""
    adapted = [dict(r, source_id=r["container_key"],
                    role="target_reference" if r["role"] == "target_reference" else "external_collection") for r in records]
    checked = logic.gp17_matches(rows, adapted, policy, lengths)
    by_id = {r["id"]: r for r in records}
    near, broad = [], []
    for row in rows:
        subject = by_id[row[1]]
        identity = float(row[2])
        coverage = 100 * (int(row[4]) - int(row[3]) + 1) / int(row[5])
        hit = {"query": row[0], "nucleotide_record_id": row[1], "container_key": subject["container_key"],
               "source_id": subject["source_id"], "role": subject["role"],
               "identity_percent": identity, "single_hsp_query_coverage_percent": coverage,
               "subject_start": int(row[6]), "subject_end": int(row[7]), "subject_frame": int(row[8])}
        if identity >= policy["gp17_near_match"]["identity_percent"] and coverage >= policy["gp17_near_match"]["minimum_query_coverage_percent"]:
            near.append(hit)
        if identity >= policy["broader_gp17_family_diagnostic"]["identity_percent"] and coverage >= policy["broader_gp17_family_diagnostic"]["minimum_query_coverage_percent"]:
            broad.append(hit)
    for actual, key in ((near, "near_matches"), (broad, "broader_family_matches")):
        if Counter(h["container_key"] for h in actual if h["role"] != "target_reference") != Counter(h["source_id"] for h in checked[key]):
            raise ValueError("Container-level report disagrees with unchanged frozen hit logic")
    near_keys, broad_keys = {r["container_key"] for r in near}, {r["container_key"] for r in broad}
    containers = {}
    for record in records:
        key = record["container_key"]
        item = containers.setdefault(key, {"container_key": key, "source_id": record["source_id"], "role": record["role"],
                "nucleotide_record_ids": [], "direct_near_hit": key in near_keys,
                "broader_diagnostic_hit": key in broad_keys, "whole_genome_group_clearance": False})
        if item["source_id"] != record["source_id"] or item["role"] != record["role"]:
            raise ValueError("A namespace key merges distinct source roles")
        item["nucleotide_record_ids"].append(record["id"])
    return {"near_matches": near, "broader_family_matches": broad,
            "direct_near_hit_container_keys": sorted(near_keys), "broader_diagnostic_container_keys": sorted(broad_keys),
            "target_self_controls_passed": checked["target_self_controls_passed"],
            "containers": [containers[k] for k in sorted(containers)]}


def old_changes(old, result):
    near = set(old["direct_near_hit_excluded_primary_ids"])
    broad = {r["source_id"] for r in old["gp17"]["broader_family_matches"]}
    changes = []
    for row in result["containers"]:
        if row["role"] != "external_collection":
            continue
        name = row["source_id"]
        before = {"direct_near_hit": name in near, "broader_diagnostic_hit": name in broad}
        after = {k: row[k] for k in before}
        if before != after:
            changes.append({"container_key": row["container_key"], "prior": before, "augmented_database": after})
    return {"compared_source_containers": sum(r["role"] == "external_collection" for r in result["containers"]),
            "changed_source_classifications": changes}


def authenticate(policy_path):
    policy = json.loads(policy_path.read_text())
    bindings = {policy_path: sha(policy_path), Path(__file__): sha(__file__)}
    if policy["wrapper_sha256"] != bindings[Path(__file__)]:
        raise ValueError("Frozen wrapper differs")
    for item in policy["inputs"]:
        path = Path(item["path"])
        if sha(path) != item["sha256"]:
            raise ValueError(f"Frozen input changed: {path}")
        bindings[path] = item["sha256"]
    base = json.loads(Path(policy["original_gp17_policy"]).read_text())
    for key in ("query_file", "tblastn_parameters", "tblastn_format", "gp17_near_match", "broader_gp17_family_diagnostic"):
        if policy[key] != base[key]:
            raise ValueError("Search queries/settings/thresholds differ from original screen")
    if policy["tblastn"] != base["executables"]["tblastn"]:
        raise ValueError("Translated-search executable differs")
    prep = load(policy["preparation_helper"], "official_gp17_input_reader")
    logic = load(base["logic_helper"], "official_gp17_frozen_logic")
    plan = json.loads(Path(policy["alignment_plan"]).read_text())
    step = json.loads(Path(policy["database_receipt"]).read_text())
    alignment_policy = Path(policy["alignment_policy"])
    if (plan["policy_sha256"] != sha(alignment_policy) or plan["nucleotide_queries"] != 273
            or plan["container_nodes"] != 266 or step["exit_code"] != 0 or step["argv"] != plan["commands"][0]
            or plan["script_sha256"] != sha(policy["alignment_helper"])
            or sha(step["log"]) != step["log_sha256"]):
        raise ValueError("Database construction does not match frozen augmented plan")
    build = step["argv"]
    if (Path(build[0]).resolve() != Path(base["executables"]["makeblastdb"]["path"]).resolve()
            or Path(build[-1]).resolve() != Path(policy["database"]).resolve()
            or Path(build[build.index("-in") + 1]).resolve() != Path(plan["subjects_fasta"]["path"]).resolve()
            or build[build.index("-dbtype") + 1] != "nucl" or "-parse_seqids" not in build):
        raise ValueError("Database receipt names a different constructor/corpus/database")
    for item in plan["inputs"] + [plan["record_manifest"], plan["subjects_fasta"]]:
        path = Path(item["path"])
        if path.resolve() not in {p.resolve() for p in bindings} or sha(path) != item["sha256"]:
            raise ValueError("Augmented plan input absent from frozen screen bindings")
    actual_db = sorted(str(p) for p in Path(policy["database"]).parent.glob("db.*") if p.is_file())
    if actual_db != sorted(policy["database_artifacts"]):
        raise ValueError("Database artifact set changed")
    prep.database_count(Path(step["log"]).read_text(), 273)
    database = json.loads(Path(policy["database"] + ".njs").read_text())
    if (database["number-of-sequences"] != 273 or database["dbtype"] != "Nucleotide"
            or Path(database["description"]).resolve() != Path(plan["subjects_fasta"]["path"]).resolve()):
        raise ValueError("Nucleotide database count/type differs")
    sequences = prep.fasta(Path(plan["subjects_fasta"]["path"]).read_bytes(), prep.DNA)
    old_sequences = prep.fasta(Path(policy["old_subjects"]).read_bytes(), prep.DNA)
    records = corpus(json.loads(Path(plan["record_manifest"]["path"]).read_text()), sequences,
                     json.loads(Path(policy["old_record_map"]).read_text()), old_sequences)
    if database["number-of-letters"] != sum(len(s) for _, s in sequences):
        raise ValueError("Database base count differs from authenticated subjects")
    _, target_sequences, lengths = prep.reference_inputs(base,
        json.loads(Path(policy["target_manifest"]).read_text()), Path(policy["query_file"]).read_bytes(), bindings)
    if any(dict(sequences).get(k) != v for k, v in target_sequences):
        raise ValueError("Original gp17 control DNA differs in augmented corpus")
    executable = Path(policy["tblastn"]["path"])
    version = subprocess.check_output([str(executable), "-version"], text=True).splitlines()[0]
    if version != "tblastn: 2.17.0+":
        raise ValueError("Translated-search version differs")
    prep.verify(bindings)
    return policy, bindings, prep, logic, records, lengths, version


def execute(policy_path, workdir, output):
    policy, bindings, prep, logic, records, lengths, version = authenticate(policy_path)
    prep.guard(workdir, output, bindings)
    workdir.mkdir(parents=True, exist_ok=False)
    hits = workdir / "tblastn.tsv"
    argv = [policy["tblastn"]["path"], "-query", policy["query_file"], "-db", policy["database"],
            *policy["tblastn_parameters"], "-outfmt", policy["tblastn_format"], "-out", str(hits)]
    plan = workdir / "frozen-plan.json"
    prep.dump(plan, {"policy_sha256": bindings[policy_path], "wrapper_sha256": bindings[Path(__file__)],
                    "command": argv, "query_lengths": lengths, "database_reused_read_only": True,
                    "inputs": [{"path": str(p), "sha256": h} for p, h in bindings.items()]})
    bindings[plan] = sha(plan)
    result = {"schema_version": 1, "status": "failed_direct_gp17_screen", "policy_sha256": bindings[policy_path],
              "wrapper_sha256": bindings[Path(__file__)], "tool_version": version,
              "nucleotide_records": len(records), "container_nodes": len({r["container_key"] for r in records}),
              "database_records_verified": 273, "database_reused_read_only": True,
              "direct_near_hit_container_keys": [], "broader_diagnostic_container_keys": [],
              "phenotype_tables_loaded": False, "training_or_scoring_executed": False,
              "whole_genome_group_analysis_complete": False, "graph_exclusions_applied": False,
              "homology_free_training_claimed": False, "stock_identity_verified": False,
              "interpretation": "Direct translated gp17 similarity only, preserving namespace-qualified controls, sources and official variants. No qualifying hit means no reported single HSP meeting fixed thresholds, not absence of homology or group clearance."}
    try:
        command = prep.command(argv, workdir, 1, policy["timeout_seconds"], bindings)
        result["command"] = command
        if command["exit_code"]:
            raise ValueError("Translated search failed or timed out")
        rows, completion = logic.completed_blast(hits.read_text(), set(lengths), 11)
        matches = screen(rows, records, policy, lengths, logic)
        prep.verify(bindings)
        result.update(matches)
        result.update({"status": "completed_direct_gp17_screen", "query_completion": completion,
                       "old_source_comparison": old_changes(json.loads(Path(policy["old_gp17_report"]).read_text()), matches),
                       "input_hashes_verified_after_execution": True})
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        result["error"] = str(exc)
        result["unassessed_containers"] = sorted({r["container_key"] for r in records})
    result["authenticated_inputs"] = [{"path": str(p), "sha256": h} for p, h in bindings.items()]
    result["raw_artifacts"] = [{"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(workdir.iterdir()) if p.is_file()]
    output.parent.mkdir(parents=True, exist_ok=True)
    prep.dump(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = execute(args.policy, args.workdir, args.output)
    print(json.dumps({k: report[k] for k in ("status", "direct_near_hit_container_keys", "broader_diagnostic_container_keys")}))
    raise SystemExit(0 if report["status"] == "completed_direct_gp17_screen" else 1)
