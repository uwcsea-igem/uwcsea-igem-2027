"""Run the frozen public-phage relatedness audit, never phenotype prediction."""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def component_exclusions(pairs: list[dict], genomes: list[dict], threshold: float, buffer: float) -> dict:
    by_id = {g["id"]: g for g in genomes}
    if len(by_id) != len(genomes):
        raise ValueError("Duplicate genome ID")
    targets = {g["id"] for g in genomes if g["role"] == "target_reference"}
    parent = {i: i for i in by_id}
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen, edges, boundary = set(), [], []
    for p in pairs:
        a, b, value = p["query"], p["subject"], p["similarity_percent_unrounded"]
        if a not in by_id or b not in by_id or a == b or not math.isfinite(value) or not 0 <= value <= 100:
            raise ValueError("Invalid intergenomic pair")
        key = tuple(sorted((a, b)))
        if key in seen:
            raise ValueError("Duplicate intergenomic pair")
        seen.add(key)
        if value >= threshold - buffer:
            parent[find(a)] = find(b)
            edge = {"query": a, "subject": b, "similarity_percent_unrounded": value}
            edges.append(edge)
            if value < threshold:
                boundary.append(edge)
    if len(seen) != len(genomes) * (len(genomes) - 1) // 2:
        raise ValueError("Incomplete all-genome similarity matrix")
    target_groups = {find(i) for i in targets}
    excluded = [g["source_id"] for g in genomes if g["role"] == "external_collection" and find(g["id"]) in target_groups]
    return {"excluded_source_ids": sorted(excluded), "edges_at_or_above_buffered_threshold": edges,
            "boundary_edges": boundary,
            "target_components": {t: sorted(by_id[i]["source_id"] for i in by_id if find(i) == find(t)) for t in sorted(targets)}}


def completed_blast(text: str, expected_queries: set[str], fields: int) -> tuple[list[list[str]], dict]:
    """Require outfmt-7 query blocks, declared row counts and final completion."""
    rows, counts, query = [], {}, None
    reported, batches, closed_batch, last_batch_end = {}, [], False, 0
    for line in text.splitlines():
        if line.startswith("# Query: "):
            closed_batch = False
            query = line[len("# Query: "):].split()[0]
            if query in counts or query not in expected_queries:
                raise ValueError("Unexpected or repeated completed-search query")
            counts[query] = 0
        elif re.fullmatch(r"# \d+ hits found", line):
            if query is None or query in reported:
                raise ValueError("Unexpected BLAST hit-count record")
            reported[query] = int(line.split()[1])
        elif line.startswith("# BLAST processed "):
            count = int(line.split()[3])
            if closed_batch or count != len(counts) - last_batch_end:
                raise ValueError("Invalid BLAST batch completion record")
            batches.append(count)
            last_batch_end, closed_batch, query = len(counts), True, None
        elif line and not line.startswith("#"):
            row = line.split("\t")
            if query is None or query not in reported or row[0] != query or len(row) != fields:
                raise ValueError("Unexpected BLAST data row")
            rows.append(row)
            counts[query] += 1
    if set(counts) != expected_queries or set(reported) != expected_queries or sum(batches) != len(expected_queries) or not closed_batch:
        raise ValueError("Incomplete BLAST query completion records")
    if counts != reported:
        raise ValueError("BLAST row counts differ from query completion records")
    return rows, {"completed_queries": sum(batches), "completed_batch_query_counts": batches, "reported_rows_by_query": reported}


def gp17_matches(rows: list[list[str]], genomes: list[dict], policy: dict, expected_query_lengths: dict[str, int]) -> dict:
    by_id = {g["id"]: g for g in genomes}
    queries = {g["id"] + "_gp17" for g in genomes if g["role"] == "target_reference"}
    if set(expected_query_lengths) != queries:
        raise ValueError("Reference protein lengths missing")
    near, family, self_controls = [], [], set()
    for row in rows:
        if len(row) != 11 or row[0] not in queries or row[1] not in by_id:
            raise ValueError("Unexpected translated-alignment schema/IDs")
        identity, start, end, length = float(row[2]), int(row[3]), int(row[4]), int(row[5])
        if not math.isfinite(identity) or not 0 <= identity <= 100 or not 1 <= start <= end <= length or length != expected_query_lengths[row[0]]:
            raise ValueError("Invalid translated-alignment identity/query coordinates")
        subject_start, subject_end, frame = int(row[6]), int(row[7]), int(row[8])
        bits, evalue = float(row[9]), float(row[10])
        if (frame not in {-3, -2, -1, 1, 2, 3} or not 1 <= min(subject_start, subject_end) <= max(subject_start, subject_end) <= by_id[row[1]]["retained_bases"]
                or (frame > 0) != (subject_end >= subject_start) or not math.isfinite(bits) or bits < 0 or not math.isfinite(evalue) or evalue < 0):
            raise ValueError("Invalid translated-alignment subject coordinates/frame/statistics")
        coverage = 100 * (end - start + 1) / length
        is_near = identity >= policy["gp17_near_match"]["identity_percent"] and coverage >= policy["gp17_near_match"]["minimum_query_coverage_percent"]
        if row[1] == row[0].removesuffix("_gp17") and is_near:
            self_controls.add(row[0])
        if by_id[row[1]]["role"] != "external_collection":
            continue
        record = {"query": row[0], "source_id": by_id[row[1]]["source_id"], "identity_percent": identity,
                  "query_coverage_percent": coverage, "subject_start": int(row[6]), "subject_end": int(row[7]),
                  "subject_frame": int(row[8])}
        if is_near:
            near.append(record)
        if identity >= policy["broader_gp17_family_diagnostic"]["identity_percent"] and coverage >= policy["broader_gp17_family_diagnostic"]["minimum_query_coverage_percent"]:
            family.append(record)
    if self_controls != queries:
        raise ValueError("gp17 target self controls incomplete")
    return {"near_matches": near, "broader_family_matches": family, "target_self_controls_passed": sorted(self_controls),
            "near_match_excluded_source_ids": sorted({x["source_id"] for x in near})}


def run(args) -> dict:
    policy = json.loads(args.policy.read_text())
    inputs = args.normalized
    manifest = json.loads((inputs / "manifest.json").read_text())
    if sha(inputs / "manifest.json") != policy["normalized_input_manifest_sha256"]:
        raise ValueError("Normalized manifest differs from frozen phage policy")
    for name, key in [("genomes.fna", "genomes_fasta_sha256"), ("lengths.json", "lengths_sha256"),
                      ("reference-gp17.faa", "reference_gp17_fasta_sha256")]:
        if sha(inputs / name) != manifest[key]:
            raise ValueError(f"Normalized input differs: {name}")
    if sha(args.adapter) != policy["intergenomic_calculator"]["standalone_adapter_sha256"]:
        raise ValueError("Adapter differs from frozen policy")
    versions, tool_hashes = {}, {}
    for name in ("blastn", "makeblastdb", "tblastn"):
        executable = args.blast_bin / name
        version = subprocess.check_output([str(executable), "-version"], text=True).splitlines()[0]
        if policy["blast_version"] not in version:
            raise ValueError(f"Unexpected BLAST version: {version}")
        versions[name], tool_hashes[name] = version, sha(executable)
    if not args.replay:
        args.workdir.mkdir(parents=True, exist_ok=False)
    elif not args.workdir.is_dir():
        raise ValueError("Replay requires the existing alignment directory")
    database, nucleotide, protein = args.workdir / "db", args.workdir / "blastn.tsv", args.workdir / "tblastn.tsv"
    numeric_nucleotide = args.workdir / "blastn-data.tsv"
    commands = [
        [str(args.blast_bin / "makeblastdb"), "-in", str(inputs / "genomes.fna"), "-dbtype", "nucl", "-parse_seqids", "-out", str(database)],
        [str(args.blast_bin / "blastn"), "-query", str(inputs / "genomes.fna"), "-db", str(database), *policy["blastn_parameters"], "-outfmt", policy["blastn_format"], "-out", str(nucleotide)],
        [str(args.blast_bin / "tblastn"), "-query", str(inputs / "reference-gp17.faa"), "-db", str(database), *policy["tblastn_parameters"], "-outfmt", policy["tblastn_format"], "-out", str(protein)],
    ]
    if not args.replay:
        for i, command in enumerate(commands):
            result = subprocess.run(command, capture_output=True, text=True)
            (args.workdir / f"step-{i}.log").write_text(result.stdout + result.stderr)
            if result.returncode:
                raise ValueError(f"Audit step {i} failed; inspect step-{i}.log")
    genomes = manifest["genomes"]
    nucleotide_rows, nucleotide_completion = completed_blast(nucleotide.read_text(), {g["id"] for g in genomes}, 14)
    protein_lengths = {g["id"] + "_gp17": manifest["reference_proteins"][g["source_id"]]["amino_acids"]
                       for g in genomes if g["role"] == "target_reference"}
    protein_rows, protein_completion = completed_blast(protein.read_text(), set(protein_lengths), 11)
    numeric_nucleotide.write_text("".join("\t".join(r) + "\n" for r in nucleotide_rows))
    command = [sys.executable, str(args.adapter), "--blast-tsv", str(numeric_nucleotide), "--genome-lengths", str(inputs / "lengths.json"), "--output", str(args.workdir / "similarities.json")]
    result = subprocess.run(command, capture_output=True, text=True)
    (args.workdir / "step-3.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise ValueError("Directional calculation failed; inspect step-3.log")
    similarity = json.loads((args.workdir / "similarities.json").read_text())
    # A complete successful search must at least recover every genome against itself.
    self_ids = {r[0] for r in nucleotide_rows if r[0] == r[1]}
    if self_ids != {g["id"] for g in genomes}:
        raise ValueError("BLAST self controls incomplete; do not interpret absent pair hits")
    whole = policy["whole_genome_exclusion"]
    groups = component_exclusions(similarity["pairs"], genomes, whole["similarity_percent"], whole["boundary_buffer_percentage_points"])
    matches = gp17_matches(protein_rows, genomes, policy, protein_lengths)
    excluded = sorted(set(groups["excluded_source_ids"]) | set(matches["near_match_excluded_source_ids"]))
    return {"schema_version": 1, "purpose": "Sequence-relatedness audit; not a host-range predictor",
            "policy_sha256": sha(args.policy), "normalized_manifest_sha256": sha(inputs / "manifest.json"),
            "adapter_sha256": sha(args.adapter), "tool_versions": versions, "executable_sha256": tool_hashes,
            "blastn_sha256": sha(nucleotide), "tblastn_sha256": sha(protein),
            "query_completion": {"blastn": nucleotide_completion, "tblastn": protein_completion},
            "wrapper_sha256": sha(Path(__file__)),
            "replayed_existing_alignments": args.replay,
            "calculator_output_sha256": sha(args.workdir / "similarities.json"),
            "source_genomes": manifest["external_genomes"], "target_genomes": manifest["target_references"],
            "pair_count": len(similarity["pairs"]), "whole_genome": groups, "gp17": matches,
            "provisional_excluded_source_ids": excluded,
            "broader_family_retained_source_ids": sorted({m["source_id"] for m in matches["broader_family_matches"]} - set(excluded)),
            "upstream_calculator_equivalence_verified": similarity["upstream_equivalence_verified"],
            "phenotype_labels_loaded": False, "infection_or_assembly_predictions_produced": False,
            "decision": "Provisional phage exclusions; calculator validation, host exclusions and a final training specification remain required"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("policy", "normalized", "blast-bin", "adapter", "workdir", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--replay", action="store_true", help="Re-evaluate existing completed alignment files without rerunning BLAST")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a new report output")
    if args.workdir.resolve().is_relative_to(args.normalized.resolve()) or args.output.resolve().is_relative_to(args.normalized.resolve()):
        parser.error("Outputs must be separate from normalized inputs")
    report = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
