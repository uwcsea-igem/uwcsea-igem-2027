"""Protocol-only streaming replay of an already completed nucleotide search.

Accept valid BLAST output batches while keeping the executed v1 HSP validator,
search commands, sequence inputs, coordinates and numeric literals unchanged.
No executable, alignment, grouping or model is run by this replay.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re

POLICY = Path("software/tipcraft/benchmarks/receptor-resource/whole-genome-alignment-replay-policy.json")


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_helper(spec):
    path = Path(spec["path"])
    if sha(path) != spec["sha256"]:
        raise ValueError("Executed v1 HSP validator bytes differ")
    module_spec = importlib.util.spec_from_file_location("executed_nucleotide_validator_v1", path)
    result = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(result)
    return result


def validate_batches(lines, lengths, numeric_stream, validate_hsp):
    """Require exact global coverage and each batch's declared query/HSP counts."""
    if isinstance(lines, (str, bytes)):
        raise ValueError("Pass a stream of text lines")
    counts, declared, self_ids, batch_queries, batches = {}, {}, set(), [], []
    current, closed, last_nonblank, rows = None, False, None, 0
    for line in lines:
        text = line.rstrip("\r\n")
        if text.strip():
            last_nonblank = text
        if text.startswith("# Query: "):
            names = text[len("# Query: "):].split()
            if not names or names[0] not in lengths or names[0] in counts:
                raise ValueError("Unexpected, empty or globally repeated query block")
            current, closed = names[0], False
            counts[current] = 0
            batch_queries.append(current)
        elif re.fullmatch(r"# \d+ hits found", text):
            if current is None or current in declared or closed:
                raise ValueError("Unexpected or repeated query hit count")
            declared[current] = int(text.split()[1])
        elif text.startswith("# BLAST processed "):
            match = re.fullmatch(r"# BLAST processed (\d+) queries", text)
            if match is None or closed or not batch_queries or current is None:
                raise ValueError("Malformed or repeated batch completion")
            count = int(match[1])
            if count != len(batch_queries) or any(q not in declared or counts[q] != declared[q] for q in batch_queries):
                raise ValueError("Batch completion differs from actual query/HSP counts")
            batches.append(count)
            current, closed, batch_queries = None, True, []
        elif text and not text.startswith("#"):
            fields = text.split("\t")
            if current is None or current not in declared or closed or fields[0] != current:
                raise ValueError("HSP occurs outside its declared query block")
            validate_hsp(fields, lengths)
            numeric_stream.write(text + "\n")
            counts[current] += 1
            rows += 1
            if fields[0] == fields[1]:
                self_ids.add(fields[0])
    if (set(counts) != set(lengths) or set(declared) != set(lengths) or counts != declared
            or sum(batches) != len(lengths) or not closed or batch_queries
            or last_nonblank is None or re.fullmatch(r"# BLAST processed \d+ queries", last_nonblank) is None
            or self_ids != set(lengths)):
        raise ValueError("Incomplete global coverage, terminal-at-EOF or per-record self controls")
    return {"completed_queries": len(counts), "completed_batch_query_counts": batches,
            "reported_hsp_rows": rows, "reported_rows_by_query": counts,
            "all_record_self_controls_passed": sorted(self_ids),
            "raw_hsp_order_preserved": True, "coordinate_offsets_applied": False}


def validate_for_aggregation(report, mapping):
    if report["status"] != "validated_all_record_nucleotide_hsps" or report["search_rerun"] is not False:
        raise ValueError("Expected a successful protocol-only replay receipt")
    helper = load_helper(report["original_execution"]["helper_snapshot"])
    artifacts = [report[key] for key in ("record_manifest", "raw_hsps", "numeric_hsps")]
    for artifact in artifacts:
        if sha(Path(artifact["path"])) != artifact["sha256"]:
            raise ValueError("Replay artifact hash differs")
    if json.loads(Path(artifacts[0]["path"]).read_text()) != mapping:
        raise ValueError("Record map differs from authenticated input")
    lengths = {r["id"]: r["retained_bases"] for r in mapping["records"]}
    if len(lengths) != len(mapping["records"]) or len(lengths) != report["nucleotide_queries"]:
        raise ValueError("Duplicate or incomplete nucleotide record map")
    with Path(artifacts[1]["path"]).open() as raw, Path(artifacts[2]["path"]).open() as numeric:
        class ExactSink:
            def write(self, line):
                if numeric.readline() != line:
                    raise ValueError("Replay numeric HSP differs from raw row/order")
        completion = validate_batches(raw, lengths, ExactSink(), helper.validate_hsp)
        if numeric.read(1):
            raise ValueError("Extra numeric rows after completed raw search")
    if completion != report["query_completion"]:
        raise ValueError("Replayed completion disagrees with accepted receipt")
    for artifact in artifacts:
        if sha(Path(artifact["path"])) != artifact["sha256"]:
            raise ValueError("Replay artifact changed during verification")
    return completion


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def execute(policy_path, workdir, output):
    policy = json.loads(policy_path.read_text())
    bindings = {policy_path: sha(policy_path), Path(__file__): sha(Path(__file__))}
    if policy["replay_helper_sha256"] != bindings[Path(__file__)] or policy["search_rerun"] is not False:
        raise ValueError("Replay policy/helper binding differs")
    for item in policy["inputs"]:
        path = Path(item["path"])
        if sha(path) != item["sha256"]:
            raise ValueError("Original execution input/receipt differs from replay policy")
        bindings[path] = item["sha256"]
    original = json.loads(Path(policy["original_report"]).read_text())
    plan = json.loads(Path(original["frozen_plan"]["path"]).read_text())
    original_helper = {"path": policy["original_helper"], "sha256": original["wrapper_sha256"]}
    helper = load_helper(original_helper)
    if (original["status"] != "failed_or_incomplete_hsp_audit" or original["database_records_verified"] != 270
            or len(original["commands"]) != 2 or any(c["exit_code"] != 0 for c in original["commands"])
            or plan["commands"] != [c["argv"] for c in original["commands"]]
            or plan["policy_sha256"] != original["policy_sha256"]
            or plan["wrapper_sha256"] != original["wrapper_sha256"]
            or sha(Path(policy["original_policy"])) != original["policy_sha256"]
            or sha(Path(original["frozen_plan"]["path"])) != original["frozen_plan"]["sha256"]):
        raise ValueError("Original search did not complete unchanged commands/database construction")
    for item in original["authenticated_inputs"] + original["raw_artifacts"]:
        path = Path(item["path"])
        if sha(path) != item["sha256"]:
            raise ValueError("Original search input/raw artifact changed")
        bindings[path] = item["sha256"]
    mapping_spec = original["record_manifest"]
    mapping = json.loads(Path(mapping_spec["path"]).read_text())
    if (sha(Path(mapping_spec["path"])) != mapping_spec["sha256"]
            or plan["record_manifest_sha256"] != mapping_spec["sha256"]
            or plan["subjects_sha256"] != original["subjects_fasta"]["sha256"]):
        raise ValueError("Original record/corpus bindings differ")
    raw_specs = [a for a in original["raw_artifacts"] if Path(a["path"]).name == "blastn.tsv"]
    if len(raw_specs) != 1 or sha(Path(raw_specs[0]["path"])) != raw_specs[0]["sha256"]:
        raise ValueError("Expected one authenticated completed raw BLAST output")
    helper_module = helper
    # Same fresh-output guard; this does not invoke any executable.
    preparation = importlib.util.spec_from_file_location("frozen_output_guard", "scripts/audit-receptor-gp17-overlap.py")
    guard_module = importlib.util.module_from_spec(preparation)
    preparation.loader.exec_module(guard_module)
    guard_module.guard(workdir, output, bindings)
    workdir.mkdir(parents=True, exist_ok=False)
    snapshot = workdir / "executed-producer-v1.py"
    with snapshot.open("xb") as stream:
        stream.write(Path(original_helper["path"]).read_bytes())
    if sha(snapshot) != original_helper["sha256"]:
        raise ValueError("Executed producer snapshot differs")
    bindings[snapshot] = sha(snapshot)
    numeric = workdir / "blastn-data.tsv"
    lengths = {r["id"]: r["retained_bases"] for r in mapping["records"]}
    if len(lengths) != 270 or len(mapping["records"]) != 270:
        raise ValueError("Expected the exact all270-record query corpus")
    report = {"schema_version": 2, "status": "failed_or_incomplete_protocol_replay",
              "purpose": policy["purpose"], "policy_path": str(policy_path),
              "policy_sha256": bindings[policy_path], "wrapper_sha256": bindings[Path(__file__)],
              "original_execution": {"report": {"path": policy["original_report"], "sha256": sha(Path(policy["original_report"]))},
                    "status": original["status"], "parser_error": original.get("error"),
                    "policy": {"path": policy["original_policy"], "sha256": original["policy_sha256"]},
                    "helper_source": original_helper,
                    "policy_sha256": original["policy_sha256"], "wrapper_sha256": original["wrapper_sha256"],
                    "helper_snapshot": {"path": str(snapshot), "sha256": sha(snapshot)},
                    "frozen_plan": original["frozen_plan"], "commands": original["commands"]},
              "record_manifest": mapping_spec, "subjects_fasta": original["subjects_fasta"], "raw_hsps": raw_specs[0],
              "frozen_plan": original["frozen_plan"], "database_records_verified": original["database_records_verified"],
              "nucleotide_queries": 270, "source_containers": 261, "source_nucleotide_records": 268,
              "original_reference_controls": 2, "container_nodes": 263,
              "search_rerun": False, "hsp_validation_changed": False, "coordinates_or_literals_changed": False,
              "intergenomic_similarities_calculated": False, "whole_genome_groups_created": False,
              "graph_exclusions_or_admission_applied": False, "phenotype_tables_loaded": False,
              "training_or_scoring_executed": False,
              "interpretation": "Completion-protocol repair only; valid batches accepted with global unique270-query coverage, per-batch counts and terminal-at-EOF. Original search execution and failed v1 receipt are preserved separately. No numeric similarities, graph or biological validation."}
    try:
        with Path(raw_specs[0]["path"]).open() as raw, numeric.open("x") as target:
            completion = validate_batches(raw, lengths, target, helper_module.validate_hsp)
        for path, expected in bindings.items():
            if sha(path) != expected:
                raise ValueError("Original input/raw/implementation changed during protocol replay")
        report.update({"status": "validated_all_record_nucleotide_hsps", "query_completion": completion,
                       "numeric_hsps": {"path": str(numeric), "sha256": sha(numeric), "bytes": numeric.stat().st_size},
                       "input_hashes_verified_after_execution": True})
    except (OSError, ValueError) as exc:
        report["error"] = str(exc)
    report["authenticated_inputs"] = [{"path": str(p), "sha256": h} for p, h in bindings.items()]
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
                      "batches": result.get("query_completion", {}).get("completed_batch_query_counts"), "search_rerun": False}))
    raise SystemExit(0 if result["status"] == "validated_all_record_nucleotide_hsps" else 1)
