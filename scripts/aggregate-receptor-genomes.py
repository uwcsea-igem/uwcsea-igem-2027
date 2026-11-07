"""Aggregate independently aligned records into conservative container groups.

Offsets are numeric coordinates only. No nucleotide strings are joined or
searched across record boundaries. The external GPL calculator and original R
function must agree on this frozen transformed HSP dataset before grouping.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys


BRIDGE = '''# SPDX-License-Identifier: GPL-3.0-only
# External interoperability driver for the pinned VIRIDIC Python adaptation.
# Original authors/source/license: scripts/third_party/README.md and VIRIDIC-LICENSE.
import dataclasses, importlib.util, json, pathlib, sys
spec = importlib.util.spec_from_file_location("viridic_external", sys.argv[1])
tool = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = tool
spec.loader.exec_module(tool)
lengths = json.loads(pathlib.Path(sys.argv[3]).read_text())
hsps = tool.parse_blast_tsv(pathlib.Path(sys.argv[2]).read_text(), lengths)
pairs = tool.calculate_similarities(hsps, lengths)
self_groups = {}
for h in hsps:
    if h.qseqid == h.sseqid:
        self_groups.setdefault(h.qseqid, []).append(h)
selves = [{"id": k, **dataclasses.asdict(tool.directional_nident(v))}
          for k,v in sorted(self_groups.items())]
pathlib.Path(sys.argv[4]).write_text(json.dumps({"pairs": pairs, "self_counts": selves,
    "upstream_commit": tool.UPSTREAM_COMMIT, "universal_upstream_equivalence_verified": False},
    sort_keys=True, separators=(",", ":")) + "\\n")
'''

FIELDS = ["query_id", "subject_id", "query_source", "subject_source", "similarity_percent_unrounded",
          "query_nident", "subject_nident", "query_length", "subject_length"]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                       allow_nan=False) + "\n").encode()


def bound(path):
    with Path(path).open("rb") as stream:
        value = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "sha256": value}


def checked(spec):
    data = Path(spec["path"]).read_bytes()
    if digest(data) != spec["sha256"]:
        raise ValueError("Frozen input changed: " + spec["path"])
    return data


def checked_file(spec):
    if bound(spec["path"])["sha256"] != spec["sha256"]:
        raise ValueError("Frozen artifact changed: " + spec["path"])


def module(path):
    spec = importlib.util.spec_from_file_location("receptor_external_"+digest(str(path).encode())[:16], path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = result
    spec.loader.exec_module(result)
    return result


def build_mapping(sources, targets, observed):
    """Derive offsets from original source-record order, never alignment order."""
    expected, containers = {}, {}
    for entry in sources:
        source = entry["source"]
        key = "source:" + source["primary_id"]
        if key in containers or entry["nucleotide_status"] != "available":
            raise ValueError("Duplicate or unavailable source container")
        offset = 0
        containers[key] = {"container_key": key, "source_id": source["primary_id"],
                           "role": "external_collection", "cohort_member": source["cohort_member"],
                           "published_model_role": source["published_model_role"],
                           "ambiguous_bases": sum(r["ambiguous_bases"] for r in entry["records"]), "records": []}
        for record in entry["records"]:
            identifier = record["neutral_nucleotide_id"]
            if identifier in expected or type(record["length"]) is not int or record["length"] <= 0:
                raise ValueError("Invalid or duplicate source record")
            expected[identifier] = {"id": identifier, "container_key": key, "source_id": source["primary_id"],
                                    "role": "external_collection", "retained_bases": record["length"],
                                    "source_record_id": record["record_id"], "offset": offset,
                                    "sequence_sha256": record["sequence_sha256"]}
            containers[key]["records"].append(identifier)
            offset += record["length"]
        containers[key]["length"] = offset
    for target in targets:
        key, identifier = "reference:" + target["source_id"], target["id"]
        if key in containers or identifier in expected or target["role"] != "target_reference":
            raise ValueError("Invalid or duplicate target reference")
        expected[identifier] = {"id": identifier, "container_key": key, "source_id": target["source_id"],
                                "role": "target_reference", "retained_bases": target["retained_bases"],
                                "source_record_id": target["accession"], "offset": 0,
                                "sequence_sha256": target["sequence_sha256"]}
        containers[key] = {"container_key": key, "source_id": target["source_id"], "role": "target_reference",
                           "cohort_member": False, "published_model_role": None, "ambiguous_bases": 0,
                           "records": [identifier], "length": target["retained_bases"]}
    actual = {x["id"]: x for x in observed}
    if len(actual) != len(observed) or set(actual) != set(expected):
        raise ValueError("Record manifest incomplete, duplicated or includes unexpected IDs")
    for identifier, want in expected.items():
        if any(actual[identifier].get(k) != want[k] for k in
               ["id", "container_key", "source_id", "role", "retained_bases", "source_record_id"]):
            raise ValueError("Alignment record differs from source provenance: " + identifier)
    return expected, dict(sorted(containers.items()))


def verify_fasta(path, records):
    observed, identifier, sequence = {}, None, []
    def finish():
        if identifier is not None:
            text = "".join(sequence).upper()
            if identifier in observed or not text or not set(text) <= set("ACGTRYSWKMBDHVN"):
                raise ValueError("Invalid or duplicate alignment FASTA record")
            observed[identifier] = (len(text), digest(text.encode("ascii")))
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith(">"):
                finish()
                identifier, sequence = line[1:].strip(), []
            elif line.strip():
                if identifier is None:
                    raise ValueError("Sequence before FASTA header")
                sequence.append(line.strip())
    finish()
    expected = {k: (v["retained_bases"], v["sequence_sha256"]) for k, v in records.items()}
    if observed != expected:
        raise ValueError("Combined alignment FASTA differs from exact source sequences")


def transform_rows(lines, records, containers):
    """Keep HSP order and strings; offsets cannot make disjoint records overlap."""
    for number, line in enumerate(lines, 1):
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != 14:
            raise ValueError(f"Row{number}: expected14 BLAST fields")
        q, s = fields[:2]
        if q not in records or s not in records:
            raise ValueError("Unknown BLAST record ID")
        for identifier, length_index, positions in [(q, 4, [6, 7]), (s, 5, [8, 9])]:
            record = records[identifier]
            if int(fields[length_index]) != record["retained_bases"]:
                raise ValueError("BLAST record length differs from source")
            for position in positions:
                coordinate = int(fields[position])
                if not 1 <= coordinate <= record["retained_bases"]:
                    raise ValueError("BLAST coordinate crosses source-record boundary")
                fields[position] = str(coordinate + record["offset"])
            fields[length_index] = str(containers[record["container_key"]]["length"])
        fields[0], fields[1] = records[q]["container_key"], records[s]["container_key"]
        yield "\t".join(fields) + "\n"


def validate_pairs(pairs, containers):
    expected = {(a, b) for i, a in enumerate(sorted(containers)) for b in sorted(containers)[i+1:]}
    seen = set()
    for row in pairs:
        key = row["query"], row["subject"]
        if key not in expected or key in seen:
            raise ValueError("Unexpected or duplicate adapter pair")
        seen.add(key)
        a, b = (containers[k]["length"] for k in key)
        n, m = row["query_nident"], row["subject_nident"]
        value = row["similarity_percent_unrounded"]
        if (type(n) is not int or type(m) is not int or not 0 <= n <= a or not 0 <= m <= b
                or row["query_length"] != a or row["subject_length"] != b
                or not math.isfinite(value) or not 0 <= value <= 100 or value != 100*(n+m)/(a+b)):
            raise ValueError("Invalid adapter identity/length/similarity")
    if seen != expected:
        raise ValueError("Missing whole-container pair; no absent values are imputed")


def components(pairs, containers, near_seeds, broad_seeds, threshold=69.99):
    validate_pairs(pairs, containers)
    parent = {k: k for k in containers}
    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k
    edges = []
    for pair in pairs:
        if pair["similarity_percent_unrounded"] >= threshold:
            a, b = pair["query"], pair["subject"]
            parent[find(a)] = find(b)
            edges.append({"a": a, "b": b, "similarity_percent_unrounded": pair["similarity_percent_unrounded"]})
    membership = defaultdict(list)
    for key in sorted(containers):
        membership[find(key)].append(key)
    groups, node_groups = {}, {}
    for members in membership.values():
        signature = digest(encoded(members))
        gid = "receptor-phage-group-" + signature[:16]
        if gid in groups:
            raise ValueError("Group-ID collision")
        targets = [k for k in members if containers[k]["role"] == "target_reference"]
        near = [k for k in members if k in near_seeds]
        groups[gid] = {"members": members, "membership_sha256": signature, "target_seeds": targets,
                       "direct_gp17_seeds": near, "target_only_excluded": bool(targets),
                       "primary_excluded": bool(targets or near),
                       "added_by_gp17_seed_component": bool(near and not targets),
                       "broader_gp17_diagnostic_members": [k for k in members if k in broad_seeds]}
        node_groups.update({k: gid for k in members})
    source = [k for k in containers if containers[k]["role"] == "external_collection"]
    return {"groups": dict(sorted(groups.items())), "container_groups": node_groups,
            "edges": sorted(edges, key=lambda e:(e["a"], e["b"])), "all_nodes_retained_as_bridges": True,
            "target_only_excluded_source_keys": [k for k in source if groups[node_groups[k]]["target_only_excluded"]],
            "primary_excluded_source_keys": [k for k in source if groups[node_groups[k]]["primary_excluded"]],
            "added_by_near_seed_component_source_keys": [k for k in source if groups[node_groups[k]]["added_by_gp17_seed_component"]],
            "broader_diagnostic_components": {k: node_groups[k] for k in broad_seeds}}


def self_checks(adapter, r_rows, containers):
    py = {r["id"]: r for r in adapter["self_counts"]}
    rself = {r["query"]: r for r in r_rows if r["query"] == r["subject"]}
    if (len(py) != len(adapter["self_counts"]) or len(rself) != sum(r["query"] == r["subject"] for r in r_rows)
            or set(py) != set(containers) or set(rself) != set(containers)):
        raise ValueError("Missing or duplicate self-control groups")
    results = []
    for key, source in containers.items():
        actual = rself[key]
        try:
            values = float(actual["nident"]), float(actual["original_aligned_query_bases"])
            finite = all(math.isfinite(x) and x.is_integer() and x > 0 for x in values)
            if not all(math.isfinite(x) for x in values):
                values = (None, None)
        except ValueError:
            values, finite = (None, None), False
        same = finite and values == (py[key]["nident"], py[key]["original_aligned_query_bases"])
        valid = same and not actual["error"] and values[0] <= source["length"]
        if source["ambiguous_bases"] == 0:
            valid = valid and values == (source["length"], source["length"])
        results.append({"container_key": key, "ambiguous_bases": source["ambiguous_bases"],
                        "length": source["length"], "original_R_nident": values[0], "original_R_span_statistic": values[1],
                        "Python_counts": py[key], "exact_R_Python_agreement": same, "passed": bool(valid),
                        "required_exact_length": source["ambiguous_bases"] == 0})
    return results


def command(argv, log, timeout=600, env=None):
    try:
        result = subprocess.run(argv, capture_output=True, text=True, timeout=timeout, env=env)
        Path(log).write_text(result.stdout + result.stderr)
        return {"argv": argv, "returncode": result.returncode, "log": bound(log)}
    except subprocess.TimeoutExpired:
        Path(log).write_text("Bounded calculation exceeded600seconds; no parameter change or retry.\n")
        return {"argv": argv, "timeout_seconds": timeout, "returncode": None, "log": bound(log)}


def run(policy_path, alignment_report_path, workdir, output, matrix):
    policy_bytes, script_bytes = Path(policy_path).read_bytes(), Path(__file__).read_bytes()
    policy = json.loads(policy_bytes)
    if policy["grouping"] != {"similarity_percent": 70.0, "boundary_buffer_percentage_points": 0.01}:
        raise ValueError("Only the frozen70-minus0.01 grouping rule is supported")
    if digest(BRIDGE.encode()) != policy["GPL_external_boundary"]["bridge_sha256"]:
        raise ValueError("External adapter bridge differs from frozen policy")
    if Path(alignment_report_path).resolve() != Path(policy["alignment_report_path"]).resolve():
        raise ValueError("Unexpected alignment report source")
    inputs = {name: checked(spec) for name, spec in policy["inputs"].items()}
    if workdir.exists() or output.exists() or matrix.exists():
        raise ValueError("Fresh workdir/report/matrix required")
    if output.resolve().is_relative_to(workdir.resolve()) or matrix.resolve().is_relative_to(workdir.resolve()):
        raise ValueError("Public artifacts must be outside raw workdir")
    alignment_bytes = Path(alignment_report_path).read_bytes()
    alignment = json.loads(alignment_bytes)
    if (alignment["status"] != "validated_all_record_nucleotide_hsps"
            or alignment["query_completion"]["completed_queries"] != 270
            or alignment["database_records_verified"] != 270
            or not alignment["input_hashes_verified_after_execution"]):
        raise ValueError("Complete270-query BLAST absence validation required")
    if alignment["policy_sha256"] != policy["inputs"]["alignment_policy"]["sha256"]:
        raise ValueError("Alignment was not run under frozen source policy")
    if alignment["wrapper_sha256"] != policy["inputs"]["alignment_helper"]["sha256"]:
        raise ValueError("Alignment helper differs from frozen implementation")
    record_bytes = checked(alignment["record_manifest"])
    if alignment["record_manifest"] != policy["prepared_record_manifest"]:
        raise ValueError("Alignment record map differs from frozen preparation")
    raw_spec, numeric_spec = alignment["raw_hsps"], alignment["numeric_hsps"]
    checked_file(raw_spec); checked_file(numeric_spec)
    manifest = json.loads(inputs["source_manifest"])
    if json.loads(inputs["source_report"])["sequence_manifest"] != policy["inputs"]["source_manifest"]:
        raise ValueError("Container report and raw manifest differ")
    targets = [g for g in json.loads(inputs["target_manifest"])["genomes"] if g["role"] == "target_reference"]
    record_manifest = json.loads(record_bytes)
    records, containers = build_mapping(manifest["entries"], targets, record_manifest["records"])
    if len(containers) != 263 or len(records) != 270 or sum(x["cohort_member"] for x in containers.values()) != 260:
        raise ValueError("Frozen263-container/270-record/260-cohort coverage differs")
    ambiguity_ids = sorted(x["source_id"] for x in containers.values() if x["ambiguous_bases"])
    if ambiguity_ids != policy["self_controls"]["known_ambiguous_source_ids"]:
        raise ValueError("Ambiguous-source inventory differs from pre-calculation policy")
    if alignment["subjects_fasta"] != policy["prepared_alignment_fasta"]:
        raise ValueError("Unexpected alignment FASTA source")
    checked_file(alignment["subjects_fasta"])
    verify_fasta(policy["alignment_fasta_path"], records)
    # Reuse the producer's authenticated outfmt7 parser for independent completion,
    # reported-row counts and byte-exact numeric-row extraction checks.
    producer = module(policy["inputs"]["alignment_helper"]["path"])
    producer.validate_for_aggregation(alignment, record_manifest)
    gp17 = json.loads(inputs["gp17_report"])
    near = sorted(gp17["direct_near_hit_excluded_primary_ids"])
    if near != policy["direct_gp17_source_seeds"]:
        raise ValueError("Direct gp17 seeds differ from frozen policy")
    broader = sorted({r["source_id"] for r in gp17["gp17"]["broader_family_matches"]} - set(near))
    if broader != policy["broader_gp17_diagnostic_source_ids"]:
        raise ValueError("Broader gp17 diagnostic source set differs from frozen policy")
    near_keys = ["source:"+x for x in near]
    broad_keys = ["source:"+x for x in policy["broader_gp17_diagnostic_source_ids"]]
    if not set(near_keys+broad_keys) <= set(containers):
        raise ValueError("Unmapped gp17 source seeds")
    workdir.mkdir(parents=True)
    transformed = workdir/"container-hsps.tsv"
    with Path(numeric_spec["path"]).open() as source, transformed.open("x") as dest:
        rows = 0
        for line in transform_rows(source, records, containers):
            dest.write(line)
            rows += 1
    if rows != alignment["query_completion"]["reported_hsp_rows"]:
        raise ValueError("Transformed HSP count differs from complete raw alignment")
    lengths_path = workdir/"container-lengths.json"
    lengths_path.write_bytes(encoded({k:v["length"] for k,v in containers.items()}))
    map_path = workdir/"container-record-map.json"
    map_path.write_bytes(encoded({"records": records, "containers": containers}))
    bridge_path = workdir/"external-calculator.py"
    bridge_path.write_text(BRIDGE)
    compare_tool = module(policy["inputs"]["R_comparison_helper"]["path"])
    if digest(compare_tool.R_DRIVER.encode()) != policy["GPL_external_boundary"]["R_driver_sha256"]:
        raise ValueError("Original R driver differs from frozen policy")
    runtime = json.loads(inputs["R_runtime_manifest"])
    verified_runtime = compare_tool.validate_runtime(Path(policy["R_home"]), Path(policy["R_library"]), runtime)
    rdir = workdir/"original-R"
    rdir.mkdir()
    driver = rdir/"run-original.R"
    driver.write_text(compare_tool.R_DRIVER)
    (rdir/"lengths.tsv").write_text("id\tlength\n"+"".join(f"{k}\t{v['length']}\n" for k,v in containers.items()))
    adapter_output = workdir/"adapter-similarities.json"
    adapter_argv = [sys.executable, str(bridge_path), policy["inputs"]["adapter_source"]["path"],
                    str(transformed), str(lengths_path), str(adapter_output)]
    rargv = [str(Path(policy["R_home"])/"bin/exec/R"), "--vanilla", "--slave", "-f", str(driver),
             "--args", policy["inputs"]["R_upstream_source"]["path"], str(transformed), str(rdir)]
    plan = {"policy_sha256": digest(policy_bytes), "wrapper_sha256": digest(script_bytes),
            "inputs": policy["inputs"], "alignment_report": bound(alignment_report_path),
            "raw_hsps": raw_spec, "numeric_hsps": numeric_spec, "record_manifest": alignment["record_manifest"],
            "alignment_fasta": bound(policy["alignment_fasta_path"]), "transformed_hsps": bound(transformed),
            "container_lengths": bound(lengths_path), "container_record_map": bound(map_path),
            "adapter_bridge": bound(bridge_path), "original_R_driver": bound(driver),
            "adapter_command": adapter_argv, "R_command": rargv, "HSP_rows": rows,
            "expected_unordered_pairs": 34453, "records": 270, "containers": 263}
    (workdir/"frozen-calculation-plan.json").write_bytes(encoded(plan))
    runs = [command(adapter_argv, workdir/"adapter.log")]
    status, numeric_validation, grouping = "HOLD_calculation_unresolved", {}, None
    adapter = None
    if runs[0]["returncode"] == 0:
        try:
            candidate = json.loads(adapter_output.read_bytes())
            validate_pairs(candidate["pairs"], containers)
            adapter = candidate
        except (ValueError, TypeError, KeyError, OSError) as exc:
            numeric_validation = {"status": "HOLD_adapter_output_unresolved", "error": f"{type(exc).__name__}: {exc}"}
    if adapter is not None:
        # Freeze generated numeric comparison inputs before original R execution.
        comparison_policy = {"inputs": {"genome_lengths": bound(lengths_path), "adapter_output": bound(adapter_output)},
                             "scope": {"distinct_unordered_pairs": 34453},
                             "agreement_requirements": {"unrounded_similarity_absolute_tolerance_percentage_points": 1e-12}}
        (workdir/"R-comparison-inputs.json").write_bytes(encoded(comparison_policy))
        environment = os.environ.copy()
        environment.update(R_HOME=str(Path(policy["R_home"]).resolve()), R_LIBS_USER=str(Path(policy["R_library"]).resolve()),
                           DYLD_LIBRARY_PATH=str((Path(policy["R_home"])/"lib").resolve()))
        runs.append(command(rargv, rdir/"execution.log", env=environment))
        if runs[-1]["returncode"] == 0:
            try:
                generic = compare_tool.compare(comparison_policy, rdir)
                selves = self_checks(adapter, compare_tool.read_tsv(rdir/"directional.tsv"), containers)
                passed = (not generic["discrepancies"] and not generic["original_R_errors"]
                          and all(s["passed"] for s in selves))
                numeric_validation = {"generic_original_helper_result": generic, "dataset_specific_self_controls": selves,
                                      "dataset_specific_directional_equivalence_verified": passed,
                                      "universal_upstream_equivalence_verified": False,
                                      "span_statistic_meaning": "Pinned upstream original-span sum diagnostic; not union coverage or forced100% for ambiguous DNA"}
                if passed:
                    status = "completed_container_grouping"
                    grouping = components(adapter["pairs"], containers, near_keys, broad_keys,
                                          policy["grouping"]["similarity_percent"]-policy["grouping"]["boundary_buffer_percentage_points"])
            except (ValueError, KeyError, TypeError, OSError) as exc:
                numeric_validation = {"dataset_specific_directional_equivalence_verified": False,
                                      "status": "HOLD_validation_unresolved", "error": f"{type(exc).__name__}: {exc}"}
        matrix.parent.mkdir(parents=True, exist_ok=True)
        with matrix.open("x") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(FIELDS)
            for row in adapter["pairs"]:
                writer.writerow([row["query"], row["subject"], containers[row["query"]]["source_id"],
                                 containers[row["subject"]]["source_id"], repr(row["similarity_percent_unrounded"]),
                                 row["query_nident"], row["subject_nident"], row["query_length"], row["subject_length"]])
    if Path(policy_path).read_bytes() != policy_bytes or Path(__file__).read_bytes() != script_bytes:
        raise ValueError("Policy or wrapper changed during calculation")
    for spec in policy["inputs"].values():
        checked(spec)
    if Path(alignment_report_path).read_bytes() != alignment_bytes:
        raise ValueError("Alignment report changed during calculation")
    checked_file(raw_spec); checked_file(numeric_spec); checked(alignment["record_manifest"])
    report = {"schema_version": 1, "status": status, "policy_sha256": digest(policy_bytes), "wrapper_sha256": digest(script_bytes),
              "inputs": plan, "frozen_calculation_plan": bound(workdir/"frozen-calculation-plan.json"),
              "commands": runs, "runtime_content_verification": verified_runtime,
              "containers": list(containers.values()), "expected_unordered_pairs": 34453,
              "numeric_matrix": bound(matrix) if matrix.exists() else None, "numeric_validation": numeric_validation,
              "numeric_matrix_for_grouping_verified": grouping is not None,
              "unresolved_reason": None if grouping else "External calculation or dataset-specific original-R validation did not complete/pass; inspect hash-bound command logs and numeric_validation. No grouping or admission is asserted.",
              "grouping": grouping, "direct_gp17_source_seeds": near, "broader_gp17_diagnostic_source_ids": policy["broader_gp17_diagnostic_source_ids"],
              "record_searches_separate": True, "DNA_concatenated_or_designed": False,
              "coordinate_offset_scope": "Numeric overlap accounting only; no HSP crosses an original record boundary",
              "all_source_nodes_including_mutants_and_extras_retained": True, "final_source_admission_decided": False,
              "host_phenotypes_loaded_or_models_trained": False,
              "raw_artifacts": [bound(p) for p in sorted(workdir.rglob("*")) if p.is_file()]}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream:
        stream.write(encoded(report))
    print(json.dumps({"status": status, "report_sha256": digest(encoded(report)),
                      "groups": len(grouping["groups"]) if grouping else None}))
    return 0 if grouping else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/whole-genome-group-policy.json"))
    parser.add_argument("--alignment-report", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    args = parser.parse_args()
    raise SystemExit(run(args.policy, args.alignment_report, args.workdir, args.output, args.matrix))
