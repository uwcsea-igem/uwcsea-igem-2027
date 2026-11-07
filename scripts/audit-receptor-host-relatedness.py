"""Compare five existing reference genomes without inferring stock equivalence."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(rows, inputs, rule, classify):
    names = {str(Path(x["path"]).resolve()): x["id"] for x in inputs}
    if len(names) != len(inputs) or len(set(names.values())) != len(inputs):
        raise ValueError("Duplicate reference path or identifier")
    pairs = {}
    for row in rows:
        a, b = (names[str(Path(row[k]).resolve())] for k in ("Ref_file", "Query_file"))
        if (a, b) in pairs:
            raise ValueError("Duplicate directed host comparison")
        values = [float(row[k]) for k in ("ANI", "Align_fraction_ref", "Align_fraction_query")]
        pairs[a, b] = {"reference": a, "query": b, "ani_percent": values[0],
                       "alignment_fraction_reference_percent": values[1],
                       "alignment_fraction_query_percent": values[2],
                       "status": classify(*values, rule)}
    ids = sorted(names.values())
    missing = [[a, b] for a in ids for b in ids if (a, b) not in pairs]
    failed_self = [a for a in ids if pairs.get((a, a), {}).get("status") != "near_target_exclude"]
    parent = {a: a for a in ids}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    edges = []
    if not missing and not failed_self:
        for i, a in enumerate(ids):
            for b in ids[i+1:]:
                if any(pairs[x]["status"] != "outside_declared_near_target_rule" for x in ((a, b), (b, a))):
                    parent[find(a)] = find(b)
                    edges.append([a, b])
    groups = {}
    for a in ids:
        groups.setdefault(find(a), []).append(a)
    return {"directed_pair_count": len(pairs), "missing_directed_pairs": missing,
            "failed_self_controls": failed_self, "comparisons": [pairs[k] for k in sorted(pairs)],
            "near_reference_edges": edges if not missing and not failed_self else None,
            "reference_components": sorted(groups.values()) if not missing and not failed_self else None,
            "status": "complete_reference_triage" if not missing and not failed_self else "HOLD_incomplete_reference_triage"}


def run(policy_path, workdir, output):
    policy_bytes, script_bytes = policy_path.read_bytes(), Path(__file__).read_bytes()
    policy = json.loads(policy_bytes)
    if len(policy["files"]) != 5 or policy["expected_directed_pairs_including_self"] != 25:
        raise ValueError("This policy requires five references and 25 directed pairs")
    if workdir.exists() or output.exists() or workdir.is_symlink() or output.is_symlink():
        raise ValueError("Use fresh workdir and report")
    bindings = policy["inputs"] + [policy["executable"]]
    for item in bindings:
        if sha(item["path"]) != item["sha256"]:
            raise ValueError("Frozen host-triage dependency changed")
    reader_path = Path(__file__).with_name("audit-receptor-host-payloads.py")
    spec = importlib.util.spec_from_file_location("host_source_reader", reader_path)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    reader_hash = sha(reader_path)
    for item in policy["files"]:
        if sha(item["path"]) != item["file_sha256"]:
            raise ValueError("Frozen host source changed")
        records = reader.fasta(Path(item["path"]).read_bytes(), set("ACGTRYSWKMBDHVN"))
        if len(records) != 1 or hashlib.sha256(records[0][1].upper().encode()).hexdigest() != item["sequence_sha256"]:
            raise ValueError("Expected exact frozen single-record host sequence")
    helper_spec = next(x for x in policy["inputs"] if x["path"].endswith("audit-genophi-host-overlap.py"))
    spec = importlib.util.spec_from_file_location("host_near_reference_rule", helper_spec["path"])
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    executable = policy["executable"]["path"]
    version = subprocess.check_output([executable, "--version"], text=True).strip()
    if version != policy["executable"]["version"]:
        raise ValueError("Unexpected skani version")
    workdir.mkdir(parents=True)
    files = workdir / "files.txt"
    files.write_text("".join(str(Path(x["path"]).resolve()) + "\n" for x in policy["files"]))
    matrix = workdir / "skani.tsv"
    p = policy["parameters"]
    command = [executable, "dist", "--rl", str(files), "--ql", str(files), "-c", str(p["compression"]),
               "-m", str(p["marker_compression"]), "-s", str(p["screen_percent"]), "--min-af",
               str(p["minimum_reported_af_percent"]), "-t", str(p["threads"]), "-o", str(matrix)]
    plan = {"policy_sha256": hashlib.sha256(policy_bytes).hexdigest(), "script_sha256": hashlib.sha256(script_bytes).hexdigest(),
            "files": policy["files"], "command": command, "reader_sha256": reader_hash}
    (workdir / "frozen-plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    (workdir / "run.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise ValueError("Host reference triage failed; inspect retained log")
    with matrix.open() as stream:
        summary = summarize(list(csv.DictReader(stream, delimiter="\t")), policy["files"],
                            policy["near_reference_rule"], helper.classify_pair)
    for item in bindings:
        if sha(item["path"]) != item["sha256"]:
            raise ValueError("Dependency changed during host reference triage")
    for item in policy["files"]:
        if sha(item["path"]) != item["file_sha256"]:
            raise ValueError("Host source changed during triage")
    if policy_path.read_bytes() != policy_bytes or Path(__file__).read_bytes() != script_bytes or sha(reader_path) != reader_hash:
        raise ValueError("Policy or helper changed during triage")
    report = {"schema_version": 1, **plan, **summary, "matrix": {"path": str(matrix), "sha256": sha(matrix)},
              "genomic_reference_relatedness_only": True, "experimental_stock_identity_verified": False,
              "phenotype_labels_loaded": False, "models_trained_or_scored": False,
              "interpretation": "Relatedness of five declared reference representations, not strain identity, experimental stock equivalence or receptor/functional infection prediction. Near-host overlap must be disclosed when selecting a known-host/cold-phage versus unseen-host evaluation regime. Missing comparisons never supply zero similarity."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: report[k] for k in ("status", "directed_pair_count", "missing_directed_pairs", "failed_self_controls", "reference_components")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/host-relatedness-policy.json"))
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.policy, args.workdir, args.output)
