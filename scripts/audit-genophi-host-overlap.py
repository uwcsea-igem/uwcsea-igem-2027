"""Run frozen host-sequence triage; no phenotype input or infection predictions."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def classify_pair(ani: float, af_ref: float, af_query: float, policy: dict) -> str:
    if any(not math.isfinite(x) or not 0 <= x <= 100 for x in (ani, af_ref, af_query)):
        raise ValueError("ANI and aligned fractions must be finite percentages")
    threshold = policy["ani_percent"]
    af_threshold = policy["minimum_af_in_each_genome_percent"]
    if ani >= threshold and min(af_ref, af_query) >= af_threshold:
        return "near_target_exclude"
    # A buffer is an engineering triage rule, not an estimator confidence bound.
    if ani >= threshold - policy["ani_boundary_buffer_percent"] and min(af_ref, af_query) >= af_threshold - policy["af_boundary_buffer_percent"]:
        return "boundary_exclude_pending_review"
    return "outside_declared_near_target_rule"


def summarize(rows: list[dict], refs: dict[str, str], queries: dict[str, str], policy: dict) -> dict:
    comparisons = {}
    for row in rows:
        try:
            host, target = refs[row["Ref_file"]], queries[row["Query_file"]]
            values = [float(row[k]) for k in ("ANI", "Align_fraction_ref", "Align_fraction_query")]
        except (KeyError, TypeError) as exc:
            raise ValueError("Unexpected skani matrix identifiers/schema") from exc
        key = (host, target)
        if key in comparisons:
            raise ValueError("Duplicate skani pair")
        comparisons[key] = {"ani_percent": values[0], "af_host_percent": values[1], "af_target_percent": values[2],
                            "status": classify_pair(*values, policy)}
    hosts = {}
    for host in sorted(refs.values()):
        pairs = {target: comparisons.get((host, target), {"status": "unreported_exclude_pending_review"})
                 for target in sorted(queries.values())}
        named = host in policy["named_target_identifiers"]
        excluded = named or any(p["status"] != "outside_declared_near_target_rule" for p in pairs.values())
        hosts[host] = {"named_target_flag": named, "exclude_from_training": excluded, "comparisons": pairs}
    return {"hosts": hosts, "reported_pairs": len(comparisons), "expected_pairs": len(refs) * len(queries),
            "excluded_hosts": [h for h, r in hosts.items() if r["exclude_from_training"]],
            "outside_rule_hosts": [h for h, r in hosts.items() if not r["exclude_from_training"]]}


def run(args) -> dict:
    spec = importlib.util.spec_from_file_location("prepare_genomes", Path(__file__).with_name("prepare-genophi-genomes.py"))
    prepare = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prepare)
    inventory, provenance = prepare.load_inventory(args.metadata_cache)
    hosts = [r for r in inventory if r["kind"] == "host"]
    policy = json.loads(args.policy.read_text())
    manifest = json.loads(args.targets.read_text())
    params = policy["host_parameters"]
    version = subprocess.check_output([str(args.executable), "--version"], text=True).strip()
    if version != "skani 0.4.0":
        raise ValueError(f"Policy requires skani 0.4.0; got {version}")
    # Validate every input before external analysis. No phenotype table is read.
    refs, input_hashes = {}, {}
    for item in hosts:
        path = args.genomes / item["path"]
        prepare.verify(path, item)
        refs[str(path)] = item["source_id"]
        input_hashes[item["source_id"]] = sha256(path)
    queries = {}
    target_hashes = {}
    for name in manifest["hosts"]:
        path = args.target_genomes / "hosts" / f"{name}.fna"
        data = path.read_bytes()
        seq = "".join(data.decode().splitlines()[1:]).upper()
        if hashlib.sha256(seq.encode()).hexdigest() != manifest["genome_sources"][name]["sequence_sha256"]:
            raise ValueError(f"Target reference hash mismatch: {name}")
        queries[str(path)] = name
        target_hashes[name] = sha256(path)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    ref_list, query_list = args.output_dir / "references.txt", args.output_dir / "queries.txt"
    ref_list.write_text("\n".join(refs) + "\n")
    query_list.write_text("\n".join(queries) + "\n")
    matrix = args.output_dir / "skani.tsv"
    command = [str(args.executable), "dist", "--rl", str(ref_list), "--ql", str(query_list),
               "-c", str(params["compression"]), "-m", str(params["marker_compression"]),
               "-s", str(params["screen_percent"]), "--min-af", str(params["minimum_reported_af_percent"]),
               "-t", str(params["threads"]), "-o", str(matrix)]
    completed = subprocess.run(command, text=True, capture_output=True)
    (args.output_dir / "run.log").write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise ValueError(f"skani exited {completed.returncode}; inspect run.log")
    result = summarize(list(csv.DictReader(matrix.open(), delimiter="\t")), refs, queries, policy["host_exclusion"])
    return {"schema_version": 1, "method": version, "executable_sha256": sha256(args.executable),
            "policy_sha256": sha256(args.policy), "target_manifest_sha256": sha256(args.targets),
            "source_provenance": provenance, "input_file_sha256": input_hashes, "target_file_sha256": target_hashes,
            "matrix_sha256": sha256(matrix), "command": command, **result,
            "no_phenotype_labels_loaded": True, "infection_or_assembly_scores_produced": False,
            "disposition": "Host triage only; target-related phage exclusions and the complete training specification remain required"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("metadata-cache", "genomes", "target-genomes", "targets", "policy", "executable", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    # A fresh directory keeps input evidence and historical analyses intact.
    for source in (args.metadata_cache, args.genomes, args.target_genomes):
        if args.output_dir.resolve().is_relative_to(source.resolve()):
            parser.error("Output directory must be separate from input directories")
    result = run(args)
    (args.output_dir / "report.json").write_text(json.dumps(result, indent=2) + "\n")
