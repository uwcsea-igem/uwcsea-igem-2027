"""Run frozen, label-free protein-calling checks on authenticated public inputs.

This creates no interaction scores or training features. Raw sequences stay in
ignored storage. Exact demo comparisons are technical correspondence checks.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import resource
import subprocess
import time


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def checked(path: Path, expected: str) -> bytes:
    data = path.read_bytes()
    if digest(data) != expected:
        raise ValueError(f"Input checksum differs: {path}")
    return data


def load_module(name: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def records(data: bytes) -> list[tuple[str, str]]:
    result, identifier, chunks = [], None, []
    for line in data.decode("ascii").splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if not line[1:].strip():
                raise ValueError("Empty FASTA header")
            if identifier is not None:
                result.append((identifier, "".join(chunks).upper()))
            identifier, chunks = line[1:].split()[0], []
        else:
            if identifier is None:
                raise ValueError("Sequence before FASTA header")
            chunks.append("".join(line.split()))
    if identifier is not None:
        result.append((identifier, "".join(chunks).upper()))
    if not result or any(not s for _, s in result):
        raise ValueError("Empty FASTA input or record")
    if len({i for i, _ in result}) != len(result):
        raise ValueError("Duplicate record identifiers; normalize input explicitly")
    return result


def protein_multiset(data: bytes) -> Counter:
    sequences = []
    for _, sequence in records(data):
        sequence = sequence.removesuffix("*")
        if not sequence or not set(sequence) <= set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError("Empty protein, internal stop or invalid amino acid")
        sequences.append(sequence)
    return Counter(sequences)


def multiset_digest(proteins: Counter) -> str:
    # JSON length framing prevents collisions between different record splits;
    # multiplicities are retained while headers/order are ignored.
    return digest(json.dumps(sorted(proteins.items()), separators=(",", ":")).encode())


def compare_multisets(predicted: Counter, deposited: Counter) -> dict:
    shared = sum((predicted & deposited).values())
    return {"predicted_proteins": sum(predicted.values()),
            "deposited_proteins": sum(deposited.values()),
            "exact_shared_protein_copies": shared,
            "predicted_only_copies": sum((predicted - deposited).values()),
            "deposited_only_copies": sum((deposited - predicted).values()),
            "exact_multiset_match": predicted == deposited,
            "predicted_multiset_sha256": multiset_digest(predicted),
            "deposited_multiset_sha256": multiset_digest(deposited)}


def prepare_inputs(policy: dict) -> tuple[list[dict], dict, list[dict]]:
    spec = policy["phage_manifest"]
    manifest = json.loads(checked(Path(spec["path"]), spec["sha256"]))
    phages = records(checked(Path(spec["normalized_fasta"]), manifest["genomes_fasta_sha256"]))
    if {i for i, _ in phages} != {g["id"] for g in manifest["genomes"]}:
        raise ValueError("Normalized phage inventory differs")
    sequences = dict(phages)
    inputs = []
    for g in manifest["genomes"]:
        sequence = sequences[g["id"]]
        if digest(sequence.encode()) != g["sequence_sha256"]:
            raise ValueError("Normalized phage sequence differs")
        inputs.append({"id": "p_" + g["id"], "kind": "phage", "role": g["role"],
                       "source_id": g["source_id"], "sequences": [sequence],
                       "source_file_sha256": g["source_file_sha256"],
                       "source_sequence_sha256": g["sequence_sha256"]})
    if len(inputs) != 99:
        raise ValueError("Require all 97 deposited phage representations and two targets")
    host_report = policy["host_integrity_report"]
    checked(Path(host_report["path"]), host_report["sha256"])
    prepare = load_module("prepare-genophi-genomes.py")
    for i, h in enumerate(policy["available_hosts"]):
        data = checked(Path(policy["host_root"]) / h["path"], h["file_sha256"])
        stats = prepare.fasta_stats(data)
        if stats["sequence_sha256"] != h["sequence_sha256"]:
            raise ValueError("Host contig-framed hash differs")
        # Replace duplicate source headers only; keep contig order and boundaries.
        seqs = [s for _, s in records_with_duplicate_headers(data)]
        inputs.append({"id": f"h_{i:03d}", "kind": "host", "role": "available_external_subset",
                       "source_id": h["source_id"], "sequences": seqs,
                       "source_file_sha256": h["file_sha256"],
                       "source_sequence_sha256": h["sequence_sha256"]})
    target = policy["target_manifest"]
    targets = json.loads(checked(Path(target["path"]), target["sha256"]))
    for i, name in enumerate(targets["hosts"]):
        path = Path(target["genomes_root"]) / "hosts" / f"{name}.fna"
        data = path.read_bytes()
        seqs = [s for _, s in records(data)]
        if len(seqs) != 1 or digest(seqs[0].encode()) != targets["genome_sources"][name]["sequence_sha256"]:
            raise ValueError(f"Frozen target host differs: {name}")
        inputs.append({"id": f"h_t{i:03d}", "kind": "host", "role": "target_reference",
                       "source_id": name, "sequences": seqs,
                       "source_file_sha256": digest(data),
                       "source_sequence_sha256": digest(seqs[0].encode())})
    upstream = policy["upstream_repository"]
    tree = json.loads(Path(upstream["tree_cache"]).read_bytes())
    entries = load_module("audit-genophi-sources.py").validate_tree(tree, upstream["root_tree_sha1"])
    demo = []
    for path, entry in sorted(entries.items()):
        if not path.startswith("data/test_data/phage_AAs/") or not path.endswith(".faa"):
            continue
        data = (Path(policy["demo_comparison"]["directory"]) / Path(path).name).read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if entry["type"] != "blob" or blob != entry["sha"]:
            raise ValueError(f"Demo source Git blob differs: {path}")
        demo.append({"source_id": Path(path).stem, "path": path, "git_blob_sha1": blob,
                     "file_sha256": digest(data), "proteins": protein_multiset(data)})
    if len(demo) != policy["demo_comparison"]["expected_files"]:
        raise ValueError("Incomplete demo-proteome comparison inventory")
    if not {d["source_id"] for d in demo} <= {i["source_id"] for i in inputs if i["kind"] == "phage"}:
        raise ValueError("Demo proteome lacks an authenticated genome counterpart")
    return inputs, manifest, demo


def records_with_duplicate_headers(data: bytes) -> list[tuple[str, str]]:
    # Validation remains the strict deposited-DNA parser; after it passes, assign
    # neutral unique headers, avoiding interpretation of source names as features.
    load_module("prepare-genophi-genomes.py").fasta_stats(data)
    lines, index = [], 0
    for line in data.decode("ascii").splitlines():
        if line.startswith(">"):
            lines.append(f">c{index:06d}")
            index += 1
        else:
            lines.append(line)
    return records(("\n".join(lines) + "\n").encode())


def run(policy_path: Path, workdir: Path) -> dict:
    if workdir.exists():
        raise ValueError("Use a new output directory; preserve previous runs")
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    tools = policy["tools"]
    for tool in tools.values():
        provenance = json.loads(checked(Path(tool["provenance_path"]), tool["provenance_sha256"]))
        checked(Path(tool["executable"]), tool["executable_sha256"])
        if provenance["commit"] != tool["commit"]:
            raise ValueError("Tool source commit differs")
    inputs, manifest, demo = prepare_inputs(policy)
    # All sources and executables are verified before starting any annotation.
    workdir.mkdir(parents=True)
    report = {"schema_version": 1, "purpose": policy["purpose"],
              "policy_path": str(policy_path), "policy_sha256": digest(policy_bytes),
              "script_sha256": digest(Path(__file__).read_bytes()),
              "tools": tools, "execution": policy["execution"],
              "phenotype_tables_loaded": False, "features_clustered": False,
              "interaction_scores_generated": False, "training_executed": False,
              "source_length_caveat": policy["source_length_caveat"],
              "files": [], "demo_comparisons": [], "target_gp17_exact_recovery": {}}
    start = time.monotonic()
    usage_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    proteomes = {}
    for item in inputs:
        directory = workdir / item["id"]
        directory.mkdir()
        fasta, faa, gff = [directory / s for s in ("input.fna", "proteins.faa", "genes.gff")]
        fasta.write_text("".join(f">c{i:06d}\n{sequence}\n" for i, sequence in enumerate(item["sequences"])))
        command = [tools[item["kind"]]["executable"], *policy["commands"][item["kind"]],
                   "-i", str(fasta), "-a", str(faa), "-o", str(gff)]
        before = time.monotonic()
        with (directory / "run.log").open("wb") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True,
                           timeout=policy["execution"]["timeout_seconds_per_genome"])
        proteins = protein_multiset(faa.read_bytes())
        proteomes[item["source_id"]] = proteins
        report["files"].append({k: v for k, v in item.items() if k != "sequences"} | {
            "input_contigs": len(item["sequences"]), "input_bases": sum(map(len, item["sequences"])),
            "normalized_input_sha256": digest(fasta.read_bytes()), "command": command,
            "wall_seconds": time.monotonic() - before, "proteins": sum(proteins.values()),
            "protein_residues": sum(len(s) * n for s, n in proteins.items()),
            "proteome_multiset_sha256": multiset_digest(proteins),
            "faa_sha256": digest(faa.read_bytes()), "gff_sha256": digest(gff.read_bytes()),
            "log_sha256": digest((directory / "run.log").read_bytes())})
    usage_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    report["resources"] = {"platform": platform.platform(), "architecture": platform.machine(),
                           "wall_seconds": time.monotonic() - start,
                           "child_user_cpu_seconds": usage_end.ru_utime - usage_start.ru_utime,
                           "child_system_cpu_seconds": usage_end.ru_stime - usage_start.ru_stime,
                           "peak_child_rss_bytes": usage_end.ru_maxrss * (1 if platform.system() == "Darwin" else 1024),
                           "rss_scope": "Peak of serial gene-caller child processes; not total workflow/training memory"}
    for d in demo:
        report["demo_comparisons"].append({k: v for k, v in d.items() if k != "proteins"} |
                                           compare_multisets(proteomes[d["source_id"]], d["proteins"]))
    for name, protein in manifest["reference_proteins"].items():
        expected = protein["protein_sequence_sha256"]
        matches = [s for s in proteomes[name] if digest(s.encode()) == expected]
        report["target_gp17_exact_recovery"][name] = {
            "reference_accession": protein["protein_id"], "expected_sha256": expected,
            "exact_sequence_recovered": bool(matches), "copies": sum(proteomes[name][s] for s in matches)}
    report["input_genomes"] = len(inputs)
    report["all_demo_multisets_match"] = all(c["exact_multiset_match"] for c in report["demo_comparisons"])
    report["interpretation"] = policy["interpretation"]
    (workdir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.policy, args.workdir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
