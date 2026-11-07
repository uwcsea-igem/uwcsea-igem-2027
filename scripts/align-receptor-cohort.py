"""Produce validated all-record nucleotide HSPs; no aggregation or grouping."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import resource
import subprocess
import sys
import time

POLICY = Path("software/tipcraft/benchmarks/receptor-resource/whole-genome-alignment-policy.json")
FIELDS = ["qseqid", "sseqid", "evalue", "bitscore", "qlen", "slen", "qstart", "qend",
          "sstart", "send", "qseq", "sseq", "nident", "gaps"]
PARAMETERS = ["-evalue", "1", "-max_target_seqs", "10000", "-num_threads", "4", "-word_size", "7",
              "-reward", "2", "-penalty", "-3", "-gapopen", "5", "-gapextend", "2"]


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def record_map(genomes):
    records = []
    for genome in genomes:
        row = dict(genome)
        row["container_key"] = ("source:" if row["role"] == "external_collection" else "reference:") + row["source_id"]
        if row["role"] == "target_reference":
            row["source_record_id"] = row["accession"]
        records.append(row)
    if len({r["id"] for r in records}) != len(records):
        raise ValueError("Duplicate nucleotide record identifiers")
    return records


def validate_hsp(fields, lengths):
    if len(fields) != len(FIELDS):
        raise ValueError("Expected exactly fourteen nucleotide alignment fields")
    query, subject = fields[:2]
    if query not in lengths or subject not in lengths:
        raise ValueError("Unexpected query/subject nucleotide identifier")
    try:
        evalue, bits = float(fields[2]), float(fields[3])
        qlen, slen, qs, qe, ss, se = map(int, fields[4:10])
        nident, gaps = map(int, fields[12:])
    except ValueError as exc:
        raise ValueError("Invalid numeric nucleotide HSP fields") from exc
    if (not math.isfinite(evalue) or evalue < 0 or not math.isfinite(bits) or bits < 0
            or qlen != lengths[query] or slen != lengths[subject]
            or not 1 <= qs <= qlen or not 1 <= qe <= qlen or not 1 <= ss <= slen or not 1 <= se <= slen):
        raise ValueError("Invalid HSP values/coordinates or lengths differ from authenticated records")
    qseq, sseq = fields[10:12]
    allowed = set("ACGTRYSWKMBDHVN-")
    if (not qseq or len(qseq) != len(sseq) or not set(qseq.upper()) <= allowed or not set(sseq.upper()) <= allowed
            or len(qseq.replace("-", "")) != abs(qe - qs) + 1
            or len(sseq.replace("-", "")) != abs(se - ss) + 1
            or not 0 <= nident <= min(abs(qe - qs) + 1, abs(se - ss) + 1)
            or gaps != qseq.count("-") + sseq.count("-")
            or any(a == b == "-" for a, b in zip(qseq, sseq))):
        raise ValueError("HSP alignment strings/counts disagree with coordinate spans")


def validate_blast(lines, lengths, numeric_stream):
    """Stream comments/data without materializing or reordering HSP strings."""
    if isinstance(lines, (str, bytes)):
        raise ValueError("Pass a stream of text lines")
    counts, declared, self_ids, current, terminal, rows = {}, {}, set(), None, None, 0
    for line in lines:
        text = line.rstrip("\r\n")
        if text.startswith("# Query: "):
            if terminal is not None:
                raise ValueError("Query appears after terminal completion")
            parts = text[len("# Query: "):].split()
            if not parts:
                raise ValueError("Empty query block identifier")
            current = parts[0]
            if current not in lengths or current in counts:
                raise ValueError("Unexpected or repeated query block")
            counts[current] = 0
        elif re.fullmatch(r"# \d+ hits found", text):
            if current is None or current in declared or terminal is not None:
                raise ValueError("Unexpected or repeated hit-count record")
            declared[current] = int(text.split()[1])
        elif text.startswith("# BLAST processed "):
            if terminal is not None or current is None:
                raise ValueError("Unexpected terminal query count")
            completion = re.fullmatch(r"# BLAST processed (\d+) queries", text)
            if completion is None:
                raise ValueError("Malformed terminal query count")
            terminal = int(completion[1])
            current = None
        elif text and not text.startswith("#"):
            fields = text.split("\t")
            if terminal is not None or current is None or current not in declared or fields[0] != current:
                raise ValueError("HSP occurs outside its declared query block")
            validate_hsp(fields, lengths)
            numeric_stream.write(text + "\n")
            counts[current] += 1
            rows += 1
            if fields[0] == fields[1]:
                self_ids.add(fields[0])
    if (set(counts) != set(lengths) or set(declared) != set(lengths) or counts != declared
            or terminal != len(lengths) or self_ids != set(lengths)):
        raise ValueError("Incomplete query counts/completion or missing per-record self control")
    return {"completed_queries": terminal, "reported_hsp_rows": rows,
            "reported_rows_by_query": counts, "all_record_self_controls_passed": sorted(self_ids),
            "raw_hsp_order_preserved": True, "coordinate_offsets_applied": False}


def validate_for_aggregation(alignment_report, record_mapping):
    """Read-only streaming proof that numeric rows exactly reproduce raw HSPs."""
    if alignment_report["status"] != "validated_all_record_nucleotide_hsps":
        raise ValueError("Alignment report has not passed")
    artifacts = [alignment_report[key] for key in ("record_manifest", "raw_hsps", "numeric_hsps")]
    for artifact in artifacts:
        if sha(Path(artifact["path"])) != artifact["sha256"]:
            raise ValueError("Alignment artifact hash differs")
    if json.loads(Path(artifacts[0]["path"]).read_text()) != record_mapping:
        raise ValueError("Provided record map differs from authenticated mapping")
    lengths = {r["id"]: r["retained_bases"] for r in record_mapping["records"]}
    if len(lengths) != len(record_mapping["records"]) or len(lengths) != alignment_report["nucleotide_queries"]:
        raise ValueError("Duplicate or incomplete nucleotide record map")
    with Path(artifacts[1]["path"]).open() as raw, Path(artifacts[2]["path"]).open() as numeric:
        class ExactSink:
            def write(self, line):
                if numeric.readline() != line:
                    raise ValueError("Numeric HSP differs from raw row/order")
        result = validate_blast(raw, lengths, ExactSink())
        if numeric.read(1):
            raise ValueError("Extra numeric HSP rows")
    if result != alignment_report["query_completion"]:
        raise ValueError("Recomputed query completion differs from report")
    for artifact in artifacts:
        if sha(Path(artifact["path"])) != artifact["sha256"]:
            raise ValueError("Alignment artifact changed during verification")
    return result


def authenticate(policy_path):
    policy = json.loads(policy_path.read_text())
    bindings = {policy_path: sha(policy_path), Path(__file__): sha(Path(__file__))}
    if (bindings[Path(__file__)] != policy["wrapper_sha256"] or policy["blastn_parameters"] != PARAMETERS
            or policy["blastn_format"] != "7 " + " ".join(FIELDS)):
        raise ValueError("Wrapper/search parameters differ from frozen policy")
    for item in policy["inputs"]:
        path = Path(item["path"])
        if sha(path) != item["sha256"]:
            raise ValueError(f"Frozen preparation input changed: {path}")
        bindings[path] = item["sha256"]
    helper_path = Path(policy["preparation_helper"])
    if helper_path not in bindings:
        raise ValueError("Input preparation helper is not authenticated")
    spec = importlib.util.spec_from_file_location("frozen_receptor_dna_inputs", helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    _, entries, raw, inherited, _, preparation_versions = helper.authenticated(Path(policy["preparation_policy"]))
    bindings.update(inherited)
    candidates, source, sequences = helper.source_inputs(entries, bindings)
    target, target_sequences, _ = helper.reference_inputs(json.loads(Path(policy["preparation_policy"]).read_text()),
        json.loads(raw["manifest.json"]), raw["reference-gp17.faa"], bindings)
    records = record_map(source + target)
    if (len(candidates) != 261 or len(source) != 268 or len(target) != 2 or len(records) != 270
            or any(c["nucleotide_status"] != "authenticated_available" for c in candidates)
            or len({r["container_key"] for r in records}) != 263):
        raise ValueError("Require all 261 source containers, 268 separate records and two original controls")
    versions = {}
    for name, item in policy["executables"].items():
        path = Path(item["path"])
        if sha(path) != item["sha256"]:
            raise ValueError("Frozen nucleotide executable differs")
        bindings[path] = item["sha256"]
        versions[name] = subprocess.check_output([str(path), "-version"], text=True).splitlines()[0]
        if versions[name] != item["version"]:
            raise ValueError("Nucleotide executable version differs")
    helper.verify(bindings)
    return policy, helper, bindings, records, sequences + target_sequences, versions, preparation_versions


def run_command(argv, workdir, index, helper, bindings):
    helper.verify(bindings)
    log = workdir / f"step-{index}.log"
    start, before = time.monotonic(), resource.getrusage(resource.RUSAGE_CHILDREN)
    print(f"START step {index}: {argv[0]}", flush=True)
    with log.open("xb") as stream:
        child = subprocess.Popen(argv, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        while True:
            try:
                code = child.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                print(f"RUNNING step {index}: pid={child.pid}, elapsed={time.monotonic()-start:.0f}s; same process", flush=True)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {"argv": argv, "exit_code": code, "log": str(log), "log_sha256": sha(log),
            "wall_seconds": time.monotonic() - start, "child_user_cpu_seconds": after.ru_utime - before.ru_utime,
            "child_system_cpu_seconds": after.ru_stime - before.ru_stime,
            "cumulative_max_direct_child_rss_bytes": after.ru_maxrss * (1 if sys.platform == "darwin" else 1024)}


def execute(policy_path, workdir, output):
    policy, helper, bindings, records, sequences, versions, preparation_versions = authenticate(policy_path)
    helper.guard(workdir, output, bindings)
    workdir.mkdir(parents=True, exist_ok=False)
    subjects, mapping = workdir / "subjects.fna", workdir / "record-manifest.json"
    subjects.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in sequences))
    if sha(subjects) != policy["subjects_fasta_sha256"]:
        raise ValueError("All-record corpus differs from the frozen, previously authenticated DNA corpus")
    dump(mapping, {"schema_version": 1, "records": records, "source_containers": 261, "source_records": 268,
                   "original_reference_controls": 2, "dna_records_concatenated": False, "bases_trimmed": False})
    bindings[subjects], bindings[mapping] = sha(subjects), sha(mapping)
    database, raw, numeric = workdir / "db", workdir / "blastn.tsv", workdir / "blastn-data.tsv"
    commands = [[policy["executables"]["makeblastdb"]["path"], "-in", str(subjects), "-dbtype", "nucl", "-parse_seqids", "-out", str(database)],
                [policy["executables"]["blastn"]["path"], "-query", str(subjects), "-db", str(database),
                 *PARAMETERS, "-outfmt", policy["blastn_format"], "-out", str(raw)]]
    if commands != [[v.format(workdir=str(workdir)) for v in command] for command in policy["command_templates"]]:
        raise ValueError("Actual commands differ from frozen command templates")
    frozen = workdir / "frozen-plan.json"
    dump(frozen, {"policy_sha256": bindings[policy_path], "wrapper_sha256": bindings[Path(__file__)],
                  "subjects_sha256": bindings[subjects], "record_manifest_sha256": bindings[mapping],
                  "queries": [r["id"] for r in records], "commands": commands, "inputs": {str(p): h for p, h in bindings.items()}})
    bindings[frozen] = sha(frozen)
    report = {"schema_version": 1, "status": "failed_or_incomplete_hsp_audit", "purpose": policy["purpose"],
              "policy_sha256": bindings[policy_path], "wrapper_sha256": bindings[Path(__file__)],
              "authenticated_inputs": [{"path": str(p), "sha256": h} for p, h in bindings.items()],
              "tool_versions": versions, "preparation_tool_versions": preparation_versions,
              "source_containers": 261, "source_nucleotide_records": 268, "original_reference_controls": 2,
              "nucleotide_queries": 270, "container_nodes": 263,
              "subjects_fasta": {"path": str(subjects), "sha256": sha(subjects)},
              "record_manifest": {"path": str(mapping), "sha256": sha(mapping)},
              "frozen_plan": {"path": str(frozen), "sha256": sha(frozen)}, "commands": [],
              "resource_scope": "Serial command CPU/time and cumulative direct-child peak RSS; not aggregate workflow memory",
              "intergenomic_similarities_calculated": False, "whole_genome_groups_created": False,
              "graph_exclusions_or_admission_applied": False, "phenotype_tables_loaded": False,
              "training_or_scoring_executed": False, "biological_predictions_produced": False,
              "interpretation": "Validated per-record HSP output only. Container-coordinate aggregation, original-R numeric comparison, grouping and admission remain separate. No records joined, trimmed or filtered."}
    try:
        for index, argv in enumerate(commands, 1):
            entry = run_command(argv, workdir, index, helper, bindings)
            report["commands"].append(entry)
            if entry["exit_code"]:
                raise ValueError(f"Nucleotide audit command {index} failed; no automatic restart")
            if index == 1:
                report["database_records_verified"] = helper.database_count(Path(entry["log"]).read_text(), 270)
                db_files = [p for p in workdir.glob("db.*") if p.is_file()]
                if not db_files:
                    raise ValueError("No nucleotide database artifacts")
                bindings.update({p: sha(p) for p in db_files})
        helper.verify(bindings)
        lengths = {r["id"]: r["retained_bases"] for r in records}
        with raw.open() as stream, numeric.open("x") as target:
            report["query_completion"] = validate_blast(stream, lengths, target)
        helper.verify(bindings)
        report.update({"status": "validated_all_record_nucleotide_hsps", "raw_hsps": {"path": str(raw), "sha256": sha(raw), "bytes": raw.stat().st_size},
                       "numeric_hsps": {"path": str(numeric), "sha256": sha(numeric), "bytes": numeric.stat().st_size},
                       "input_hashes_verified_after_execution": True})
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        report["error"] = str(exc)
    report["raw_artifacts"] = [{"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(workdir.iterdir()) if p.is_file()]
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
    print(json.dumps({"status": result["status"], "queries": result["nucleotide_queries"],
                      "reported_hsps": result.get("query_completion", {}).get("reported_hsp_rows")}))
    raise SystemExit(0 if result["status"] == "validated_all_record_nucleotide_hsps" else 1)
