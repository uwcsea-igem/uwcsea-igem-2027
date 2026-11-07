"""Screen authenticated receptor-cohort DNA for direct target-gp17 near hits.

This source-relatedness audit does not assess host range, binding, assembly,
graph independence or suitability for model fitting.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

POLICY = Path("software/tipcraft/benchmarks/receptor-resource/gp17-overlap-policy.json")
DNA = set("ACGTRYSWKMBDHVN")
PROTEIN = set("ACDEFGHIKLMNPQRSTVWY")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def file_sha(path):
    return sha(path.read_bytes())


def fasta(data, alphabet):
    records, seen, name, sequence = [], set(), None, []
    for line in data.decode("ascii").splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if name is not None:
                if not sequence:
                    raise ValueError("Empty FASTA record")
                records.append((name, "".join(sequence)))
            fields = line[1:].split()
            if not fields or fields[0] in seen:
                raise ValueError("Missing or duplicate FASTA identifier")
            name, sequence = fields[0], []
            seen.add(name)
        else:
            value = "".join(line.split()).upper()
            if name is None or not set(value) <= alphabet:
                raise ValueError("Sequence before header or unsupported FASTA alphabet")
            sequence.append(value)
    if name is None or not sequence:
        raise ValueError("Missing or empty FASTA record")
    records.append((name, "".join(sequence)))
    return records


def source_inputs(entries, bindings):
    """Authenticate each DNA record independently of protein annotation status."""
    candidates, genomes, sequences, seen = [], [], [], set()
    for entry in entries:
        source = entry["source"]
        candidate = {"primary_id": source["primary_id"], "container_id": source["container_id"],
                     "cohort_member": source["cohort_member"], "published_model_role": source["published_model_role"],
                     "container_status": entry["status"], "protein_status": entry["protein_status"],
                     "container_error": entry.get("error"), "protein_errors": entry.get("protein_errors", []),
                     "nucleotide_ids": [], "nucleotide_status": "unassessed_input_error", "admitted": False}
        if candidate["primary_id"] in seen:
            raise ValueError("Duplicate cohort primary identifier")
        seen.add(candidate["primary_id"])
        try:
            if entry["nucleotide_status"] != "available":
                raise ValueError("Container has no authenticated available nucleotide output")
            container = Path(entry["container_path"])
            data = container.read_bytes()
            blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
            if len(data) != source["declared_size_bytes"] or blob != source["git_blob_sha1"]:
                raise ValueError("Original container bytes differ from frozen source identity")
            output = entry["outputs"]["nucleotide_fasta"]
            path = Path(output["path"])
            content = path.read_bytes()
            if sha(content) != output["sha256"] or len(content) != output["bytes"]:
                raise ValueError("Normalized nucleotide file differs from container inventory")
            records = fasta(content, DNA)
            if len(records) != len(entry["records"]) or len(records) != output["sequences"]:
                raise ValueError("Nucleotide record count differs from inventory")
            pending = []
            for (identifier, sequence), record in zip(records, entry["records"]):
                if (identifier != record["neutral_nucleotide_id"] or len(sequence) != record["length"]
                        or sha(sequence.encode("ascii")) != record["sequence_sha256"]):
                    raise ValueError("Nucleotide identifiers/order/bases differ from inventory")
                pending.append({"id": identifier, "source_id": source["primary_id"], "role": "external_collection",
                                "retained_bases": len(sequence), "source_record_id": record["record_id"]})
            bindings[container], bindings[path] = sha(data), sha(content)
            genomes.extend(pending)
            sequences.extend(records)
            candidate["nucleotide_ids"] = [r[0] for r in records]
            candidate["nucleotide_status"] = "authenticated_available"
        except (OSError, ValueError, KeyError) as exc:
            candidate["nucleotide_error"] = str(exc)
        candidates.append(candidate)
    if len({g["id"] for g in genomes}) != len(genomes):
        raise ValueError("Duplicate normalized nucleotide identifiers")
    return candidates, genomes, sequences


def reference_inputs(policy, manifest, query_bytes, bindings):
    queries = dict(fasta(query_bytes, PROTEIN))
    target_records = {g["source_id"]: g for g in manifest["genomes"] if g["role"] == "target_reference"}
    if set(target_records) != {"T7", "T3"} or set(queries) != {"t000_gp17", "t001_gp17"}:
        raise ValueError("Require the two original full-length target references")
    genomes, sequences, lengths = [], [], {}
    for item in policy["target_nucleotide_references"]:
        name, path = item["source_id"], Path(item["path"])
        data = path.read_bytes()
        record, protein = target_records[name], manifest["reference_proteins"][name]
        dna = fasta(data, DNA)
        if len(dna) != 1 or sha(data) != item["sha256"] or sha(data) != record["source_file_sha256"]:
            raise ValueError("Original target nucleotide file differs")
        sequence = dna[0][1]
        if sha(sequence.encode()) != record["sequence_sha256"] or len(sequence) != record["retained_bases"]:
            raise ValueError("Original target nucleotide sequence differs")
        query = record["id"] + "_gp17"
        if (not protein["cds_translation_verified"] or protein["genome_sequence_sha256"] != record["sequence_sha256"]
                or len(queries[query]) != protein["amino_acids"]
                or sha(queries[query].encode()) != protein["protein_sequence_sha256"]):
            raise ValueError("Full-length target gp17 differs from the verified translation")
        bindings[path] = sha(data)
        genomes.append(record)
        sequences.append((record["id"], sequence))
        lengths[query] = len(queries[query])
    if len(genomes) != 2 or {g["source_id"] for g in genomes} != {"T7", "T3"}:
        raise ValueError("Duplicate or missing target nucleotide reference")
    return genomes, sequences, lengths


def verify(bindings):
    for path, expected in bindings.items():
        if file_sha(path) != expected:
            raise ValueError(f"Frozen input/implementation changed: {path}")


def authenticated(policy_path):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    bindings = {policy_path: sha(policy_bytes), Path(__file__): file_sha(Path(__file__))}
    if policy["wrapper_sha256"] != bindings[Path(__file__)]:
        raise ValueError("Wrapper differs from the pre-search frozen policy")
    raw = {}
    for item in policy["inputs"]:
        path = Path(item["path"])
        data = path.read_bytes()
        if path.name in raw or sha(data) != item["sha256"]:
            raise ValueError("Duplicate or changed frozen source input")
        raw[path.name], bindings[path] = data, item["sha256"]
    report = json.loads(raw["receptor-cohort-container-audit.json"])
    manifest_path = Path(report["sequence_manifest"]["path"])
    manifest_bytes = raw["sequence-manifest.json"]
    if (manifest_path not in bindings or sha(manifest_bytes) != report["sequence_manifest"]["sha256"]
            or json.loads(manifest_bytes)["policy_sha256"] != report["policy_sha256"]
            or report["policy_sha256"] != sha(raw["cohort-container-policy.json"])):
        raise ValueError("Container report/manifest/policy bindings differ")
    entries = json.loads(manifest_bytes)["entries"]
    if len(entries) != policy["expected_source_candidates"] or len(entries) != 261:
        raise ValueError("Require every one of the 261 frozen source candidates")
    if (policy["gp17_near_match"]["identity_percent"] != 90.0
            or policy["gp17_near_match"]["minimum_query_coverage_percent"] != 80.0
            or policy["broader_gp17_family_diagnostic"]["identity_percent"] != 40.0
            or policy["broader_gp17_family_diagnostic"]["minimum_query_coverage_percent"] != 80.0):
        raise ValueError("Frozen gp17 thresholds differ from the unchanged 90/80 and 40/80 rules")
    if (Path(policy["query_file"]) not in bindings
            or bindings[Path(policy["query_file"])] != json.loads(raw["manifest.json"])["reference_gp17_fasta_sha256"]
            or policy["tblastn_parameters"] != ["-evalue", "1e-5", "-max_target_seqs", "10000", "-num_threads", "4", "-seg", "yes"]
            or policy["tblastn_format"] != "7 qseqid sseqid pident qstart qend qlen sstart send sframe bitscore evalue"
            or set(policy["executables"]) != {"makeblastdb", "tblastn"}):
        raise ValueError("Query file/search settings differ from the frozen full-reference search")
    helper_path = Path(policy["logic_helper"])
    if helper_path not in bindings:
        raise ValueError("Completion/near-match helper is not authenticated")
    spec = importlib.util.spec_from_file_location("frozen_gp17_logic", helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    versions = {}
    for name, item in policy["executables"].items():
        path = Path(item["path"])
        bindings[path] = item["sha256"]
        verify({path: item["sha256"]})
        versions[name] = subprocess.check_output([str(path), "-version"], text=True).splitlines()[0]
        if "2.17.0+" not in versions[name]:
            raise ValueError("BLAST version differs from the frozen binary")
    verify(bindings)
    return policy, entries, raw, bindings, helper, versions


def guard(workdir, output, protected):
    if any(p.exists() or p.is_symlink() for p in (workdir, output)):
        raise ValueError("Workdir and output must both be new")
    a, b = workdir.resolve(), output.resolve()
    if a.is_relative_to(b) or b.is_relative_to(a):
        raise ValueError("Workdir and report must be separate")
    for destination in (a, b):
        if any(destination.is_relative_to(p.resolve()) or p.resolve().is_relative_to(destination) for p in protected):
            raise ValueError("Output overlaps an authenticated input")


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def command(argv, directory, index, timeout, bindings):
    verify(bindings)
    log = directory / f"step-{index}.log"
    start = time.monotonic()
    with log.open("xb") as stream:
        child = subprocess.Popen(argv, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
            code = -signal.SIGKILL
    return {"argv": argv, "exit_code": code, "wall_seconds": time.monotonic() - start,
            "log": str(log), "log_sha256": file_sha(log)}


def database_count(log, expected):
    counts = [int(n) for n in re.findall(r"Adding sequences from FASTA; added (\d+) sequences", log)]
    if counts != [expected]:
        raise ValueError("makeblastdb did not report every authenticated nucleotide record")
    return expected


def mark_candidates(candidates, matches):
    near = set(matches["near_match_excluded_source_ids"])
    broad = {m["source_id"] for m in matches["broader_family_matches"]}
    for row in candidates:
        if row["nucleotide_status"] != "authenticated_available":
            row["direct_gp17_disposition"] = "unassessed_input_error"
        else:
            row["direct_gp17_disposition"] = ("direct_near_match_exclude" if row["primary_id"] in near else
                                                "broader_similarity_diagnostic_only" if row["primary_id"] in broad else
                                                "no_reported_hit_meeting_frozen_thresholds")
        row["whole_genome_or_group_clearance"] = False
    return candidates


def execute(policy_path, workdir, output):
    policy, entries, raw, bindings, helper, versions = authenticated(policy_path)
    candidates, genomes, sequences = source_inputs(entries, bindings)
    targets, target_sequences, lengths = reference_inputs(policy, json.loads(raw["manifest.json"]), raw["reference-gp17.faa"], bindings)
    if set(g["id"] for g in genomes) & set(g["id"] for g in targets):
        raise ValueError("Cohort nucleotide IDs collide with target controls")
    if len(genomes) + len(targets) > 10000:
        raise ValueError("Frozen max_target_seqs would not cover every subject record")
    guard(workdir, output, bindings)
    report = {"schema_version": 1, "status": "HOLD_incomplete_nucleotide_inputs", "purpose": policy["purpose"],
              "policy_sha256": bindings[policy_path], "wrapper_sha256": bindings[Path(__file__)],
              "authenticated_inputs": [{"path": str(p), "sha256": h} for p, h in bindings.items()],
              "tool_versions": versions, "expected_source_candidates": len(entries),
              "authenticated_nucleotide_candidates": sum(c["nucleotide_status"] == "authenticated_available" for c in candidates),
              "source_nucleotide_records": len(genomes), "separate_original_reference_controls": [g["source_id"] for g in targets],
              "expected_source_nucleotide_records": policy["expected_source_nucleotide_records"],
              "candidate_rows": candidates, "commands": [], "query_completion_verified": False,
              "direct_near_hit_excluded_primary_ids": [], "whole_genome_group_analysis_complete": False,
              "homology_free_training_claimed": False, "stock_identity_verified": False,
              "training_or_scoring_executed": False, "phenotype_tables_loaded": False,
              "binding_infection_or_assembly_predictions_produced": False, "all_rows_admitted": False,
              "interpretation": "Direct translated gp17 similarity flags only; no graph propagation, whole-genome/group clearance or model-run authorization. No-hit means no reported qualifying single alignment under this search, not absence of homology."}
    workdir.mkdir(parents=True, exist_ok=False)
    if report["authenticated_nucleotide_candidates"] == len(entries) and len(genomes) == policy["expected_source_nucleotide_records"]:
        subject = workdir / "subjects.fna"
        subject.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in sequences + target_sequences))
        bindings[subject] = file_sha(subject)
        database, hits = workdir / "db", workdir / "tblastn.tsv"
        argv = [[policy["executables"]["makeblastdb"]["path"], "-in", str(subject), "-dbtype", "nucl", "-parse_seqids", "-out", str(database)],
                [policy["executables"]["tblastn"]["path"], "-query", policy["query_file"], "-db", str(database),
                 *policy["tblastn_parameters"], "-outfmt", policy["tblastn_format"], "-out", str(hits)]]
        plan = {"policy_sha256": bindings[policy_path], "wrapper_sha256": bindings[Path(__file__)],
                "subjects_sha256": bindings[subject], "subjects": genomes + targets,
                "query_lengths": lengths, "commands": argv, "phenotype_tables_loaded": False}
        dump(workdir / "frozen-plan.json", plan)
        bindings[workdir / "frozen-plan.json"] = file_sha(workdir / "frozen-plan.json")
        try:
            for index, args in enumerate(argv, 1):
                record = command(args, workdir, index, policy["timeout_seconds_per_step"], bindings)
                report["commands"].append(record)
                if record["exit_code"]:
                    raise ValueError(f"Search step {index} failed or timed out")
                if index == 1:
                    report["database_records_verified"] = database_count(Path(record["log"]).read_text(), len(genomes) + len(targets))
                    db_files = [p for p in workdir.glob("db.*") if p.is_file()]
                    if not db_files:
                        raise ValueError("No nucleotide database artifacts")
                    bindings.update({p: file_sha(p) for p in db_files})
            verify(bindings)
            rows, completion = helper.completed_blast(hits.read_text(), set(lengths), 11)
            matches = helper.gp17_matches(rows, genomes + targets, policy, lengths)
            report.update({"status": "completed_direct_gp17_screen", "query_completion_verified": True,
                           "query_completion": completion, "gp17": matches,
                           "direct_near_hit_excluded_primary_ids": matches["near_match_excluded_source_ids"],
                           "candidate_rows": mark_candidates(candidates, matches),
                           "tblastn_sha256": file_sha(hits), "frozen_plan_sha256": file_sha(workdir / "frozen-plan.json")})
            verify(bindings)
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            report["status"], report["error"] = "failed_direct_gp17_screen", str(exc)
            report["direct_near_hit_excluded_primary_ids"], report["query_completion_verified"] = [], False
            report.pop("gp17", None)
            for row in candidates:
                row["direct_gp17_disposition"] = "unassessed_search_failure"
    else:
        for row in candidates:
            row["direct_gp17_disposition"] = "unassessed_incomplete_input_roster"
    try:
        verify(bindings)
        report["input_hashes_verified_after_execution"] = True
    except (OSError, ValueError) as exc:
        report.update({"status": "failed_direct_gp17_screen", "error": str(exc),
                       "input_hashes_verified_after_execution": False, "query_completion_verified": False,
                       "direct_near_hit_excluded_primary_ids": []})
        report.pop("gp17", None)
        for row in candidates:
            row["direct_gp17_disposition"] = "unassessed_input_changed"
    report["disposition_counts"] = dict(Counter(c["direct_gp17_disposition"] for c in candidates))
    report["raw_artifacts"] = [{"path": str(p), "sha256": file_sha(p)} for p in sorted(workdir.iterdir()) if p.is_file()]
    output.parent.mkdir(parents=True, exist_ok=True)
    dump(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=POLICY)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = execute(args.policy, args.workdir, args.output)
    print(json.dumps({"status": result["status"], "source_candidates": result["expected_source_candidates"],
                      "authenticated_nucleotide_candidates": result["authenticated_nucleotide_candidates"],
                      "direct_near_hit_excluded_primary_ids": result["direct_near_hit_excluded_primary_ids"]}))
    raise SystemExit(0 if result["status"] == "completed_direct_gp17_screen" else 1)
