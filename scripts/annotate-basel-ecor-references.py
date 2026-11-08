"""Call genes independently for every authenticated secondary-panel reference.

This prepares annotations, not training membership, shared features or scores.
All source contigs and all77 references remain in the inventory.
"""
import argparse
from collections import Counter
import importlib.util
import json
import hashlib
from pathlib import Path
import platform
import re
import resource
import subprocess
import sys
import time


FLAGS = {"host": ["-p", "single", "-g", "11", "-q", "-f", "gff"],
         "phage": ["-p", "meta", "-q", "-f", "gff"]}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                       allow_nan=False)+"\n").encode()


def binding(path):
    path = Path(path)
    return {"path": str(path), "sha256": sha(path.read_bytes())}


def checked(spec):
    path = Path(spec["path"])
    if path.is_symlink():
        raise ValueError("Symlink input is not accepted")
    data = path.read_bytes()
    if sha(data) != spec["sha256"]:
        raise ValueError("Frozen input changed: "+str(path))
    return data


def load_module(spec):
    checked(spec)
    name = "secondary_annotation_"+sha(spec["path"].encode())[:16]
    definition = importlib.util.spec_from_file_location(name, spec["path"])
    result = importlib.util.module_from_spec(definition)
    sys.modules[name] = result
    definition.loader.exec_module(result)
    return result


def attributes(text):
    result = {}
    for item in text.rstrip(";").split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError("Malformed caller attribute")
        key, value = item.split("=", 1)
        if not key or key in result:
            raise ValueError("Duplicate or empty caller attribute")
        result[key] = value
    return result


def validate_outputs(faa, gff, contig_lengths, parser):
    """Match every FAA header to one bounded GFF CDS; preserve zero-CDS contigs."""
    proteins = parser.protein_multiset(faa)
    sequences = dict(parser.records(faa))
    headers = [line[1:] for line in faa.decode("ascii").splitlines() if line.startswith(">")]
    by_gene, inventory = {}, []
    for header in headers:
        parts = header.split(" # ")
        if len(parts) != 5:
            raise ValueError("Expected five-field Prodigal FAA header")
        identifier = parts[0]
        contig, ordinal = identifier.rsplit("_", 1)
        start, end, strand = map(int, parts[1:4])
        attrs = attributes(parts[4])
        gene_id, partial = attrs.get("ID"), attrs.get("partial")
        if (contig not in contig_lengths or not ordinal.isdigit() or int(ordinal) < 1
                or strand not in {-1, 1} or not 1 <= start <= end <= contig_lengths[contig]
                or not gene_id or gene_id in by_gene or partial not in {"00", "01", "10", "11"}):
            raise ValueError("Invalid or duplicate FAA CDS coordinates/ID")
        aa = sequences[identifier].removesuffix("*")
        row = {"protein_id": identifier, "caller_gene_id": gene_id, "contig_id": contig,
               "start_1based": start, "end_1based": end, "strand": strand, "partial": partial,
               "length": len(aa), "sequence_sha256": sha(aa.encode("ascii")),
               "terminal_stop_removed_for_multiset_only": sequences[identifier].endswith("*"),
               "caller_genetic_code": attrs.get("genetic_code")}
        by_gene[gene_id] = row
        inventory.append(row)
    seen = set()
    for line in gff.decode("ascii").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] != "CDS":
            raise ValueError("Expected nine-field caller CDS GFF")
        contig, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
        attrs = attributes(fields[8])
        gene_id = attrs.get("ID")
        if gene_id in seen or gene_id not in by_gene or strand not in {"+", "-"} or fields[7] not in {"0", "1", "2"}:
            raise ValueError("Duplicate/unmapped GFF CDS or invalid strand/phase")
        expected = by_gene[gene_id]
        if (contig, start, end, 1 if strand == "+" else -1, attrs.get("partial")) != (
                expected["contig_id"], expected["start_1based"], expected["end_1based"], expected["strand"], expected["partial"]):
            raise ValueError("FAA/GFF coordinate, strand or partial annotation disagreement")
        seen.add(gene_id)
    if seen != set(by_gene) or len(headers) != sum(proteins.values()):
        raise ValueError("FAA/GFF protein inventory is incomplete")
    counts = Counter(row["contig_id"] for row in inventory)
    return {"proteins": len(inventory), "protein_residues": sum(p["length"] for p in inventory),
            "proteome_multiset_sha256": parser.multiset_digest(proteins),
            "partial_cds_counts": dict(sorted(Counter(p["partial"] for p in inventory).items())),
            "contigs_with_no_CDS": sorted(set(contig_lengths)-set(counts)),
            "CDS_counts_by_contig": {key: counts[key] for key in contig_lengths},
            "FAA_GFF_inventory_and_coordinates_agree": True}, inventory


def normalize(records, identifier, parser):
    parsed = parser.records(records)
    if any(not set(seq) <= set("ACGTRYSWKMBDHVN") for _, seq in parsed):
        raise ValueError("Invalid IUPAC DNA")
    mapping = [{"source_record_id": old, "neutral_id": f"{identifier}c{i:06d}", "length": len(seq),
                "sequence_sha256": sha(seq.encode("ascii"))} for i, (old, seq) in enumerate(parsed)]
    data = "".join(f">{row['neutral_id']}\n{seq}\n" for row, (_, seq) in zip(mapping, parsed)).encode("ascii")
    return data, mapping


def validate_roster(genomes, controls, authority):
    if (len(genomes) != 77 or Counter(i["kind"] for i in genomes) != {"host": 25, "phage": 52}
            or len({i["id"] for i in genomes+controls}) != 79 or len(controls) != 2):
        raise ValueError("Require all25+52 genomes and exactly2 separate controls")
    for kind, field, prefix in [("host", "host_id_to_wgs_master_accession", "h"), ("phage", "phage_id_to_accession", "p")]:
        selected = [i for i in genomes if i["kind"] == kind]
        if (len({i["source_id"] for i in selected}) != len(selected)
                or {i["source_id"] for i in selected} != set(authority[field])
                or any(not re.fullmatch(prefix+r"\d{3}", i["id"]) or i["role"] != "secondary_reference_not_admitted" for i in selected)):
            raise ValueError("Source cohort incomplete/duplicated or identifiers unsafe")
    if ({c["id"] for c in controls} != {"ch000", "cp000"}
            or {c["kind"] for c in controls} != {"host", "phage"}
            or any(c["role"] != "technical_control" for c in controls)):
        raise ValueError("Unexpected technical control inventory")


def prepare(item, directory, helpers):
    """Authenticate existing per-genome representations; never select by grouping."""
    parser, ecor, phage = helpers
    data = checked(item["fasta"])
    if item["role"] == "technical_control":
        expected = parser.protein_multiset(checked(item["expected_faa"]))
        if parser.multiset_digest(expected) != item["expected_proteome_multiset_sha256"]:
            raise ValueError("Original control proteome differs")
    elif item["kind"] == "host":
        master = ecor.parse_master(checked(item["master_file"]), item["source_id"], item["wgs_accession"], item["project"])
        observed = ecor.parse_fasta(data, master)
        stored = json.loads(checked(item["contig_manifest"]))
        if stored["records"] != observed or stored["master"] != master or stored["DNA_concatenated_trimmed_or_designed"]:
            raise ValueError("Full ECOR contig inventory differs")
    else:
        official = phage.sequence_record(data, "fasta", item["requested_accession"])
        annotation = phage.sequence_record(checked(item["annotation"]), item["annotation_format"], item["requested_accession"])
        if (official["sequence"], official["versioned_accession"]) != (annotation["sequence"], annotation["versioned_accession"]):
            raise ValueError("Official phage formats differ")
        if (official["versioned_accession"], official["length"], official["sequence_sha256"]) != (
                item["versioned_accession"], item["expected_length"], item["expected_sequence_sha256"]):
            raise ValueError("Official phage version or DNA differs from frozen replay")
    normalized, mapping = normalize(data, item["id"], parser)
    if item["role"] != "technical_control" and item["kind"] == "host":
        if [(r["length"], r["sequence_sha256"]) for r in mapping] != [(r["length"], r["sequence_sha256"]) for r in observed]:
            raise ValueError("Host normalization changed source records")
    directory.mkdir(parents=True)
    (directory/"input.fna").write_bytes(normalized)
    (directory/"contig-map.json").write_bytes(encoded(mapping))
    return {"item": item, "status": "input_prepared", "directory": str(directory),
            "input": binding(directory/"input.fna"), "contig_map": binding(directory/"contig-map.json"),
            "input_contigs": len(mapping), "input_bases": sum(r["length"] for r in mapping)}


def call(prepared, tool, flags, timeout, parser):
    item = prepared["item"]
    if prepared["status"] != "input_prepared":
        return prepared
    directory = Path(prepared["directory"])
    fasta, faa, gff, log = [directory/name for name in ["input.fna", "proteins.faa", "genes.gff", "caller.log"]]
    argv = [tool["executable"], *flags, "-i", str(fasta), "-a", str(faa), "-o", str(gff)]
    result = {**prepared, "command": argv, "status": "caller_failed"}
    start, before = time.monotonic(), resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        checked(prepared["input"]); checked(item["fasta"])
        checked({"path": tool["executable"], "sha256": tool["executable_sha256"]})
        with log.open("xb") as stream:
            completed = subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
        result["exit_code"] = completed.returncode
        if completed.returncode:
            raise ValueError("Caller exited nonzero; no retry or alternate flags")
        checked(prepared["input"]); checked(item["fasta"])
        maps = json.loads(checked(prepared["contig_map"]))
        summary, proteins = validate_outputs(faa.read_bytes(), gff.read_bytes(), {r["neutral_id"]: r["length"] for r in maps}, parser)
        result.update(summary)
        coordinate_path = directory/"protein-coordinate-inventory.json"
        coordinate_path.write_bytes(encoded(proteins))
        result["protein_coordinate_inventory"] = binding(coordinate_path)
        result["status"] = "annotation_validated"
        if item["role"] == "technical_control":
            result["control_exact_proteome_multiset_match"] = summary["proteome_multiset_sha256"] == item["expected_proteome_multiset_sha256"]
            if not result["control_exact_proteome_multiset_match"]:
                result["status"] = "control_mismatch"
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    result["resources"] = {"wall_seconds": time.monotonic()-start,
        "child_user_cpu_seconds": after.ru_utime-before.ru_utime, "child_system_cpu_seconds": after.ru_stime-before.ru_stime,
        "cumulative_direct_child_peak_rss_bytes": after.ru_maxrss*(1 if sys.platform == "darwin" else 1024)}
    result["outputs"] = {name: binding(path) for name, path in [("FAA", faa), ("GFF", gff), ("log", log)] if path.exists()}
    (directory/"receipt.json").write_bytes(encoded(result))
    return result


def run(policy_path, workdir, output):
    import Bio
    if Bio.__version__ != "1.86":
        raise ValueError("Frozen sequence-record parser requires Biopython1.86")
    policy_bytes, code_bytes = policy_path.read_bytes(), Path(__file__).read_bytes()
    policy = json.loads(policy_bytes)
    if sha(code_bytes) != policy["wrapper_sha256"] or policy["flags"] != FLAGS or policy["execution"]["workers"] != 1:
        raise ValueError("Wrapper/flags/serial-execution rule differs from frozen policy")
    if workdir.exists() or workdir.is_symlink() or output.exists() or output.is_symlink():
        raise ValueError("Fresh workdir and report required")
    if output.resolve().is_relative_to(workdir.resolve()):
        raise ValueError("Public report must be outside raw workdir")
    for spec in policy["inputs"].values():
        checked(spec)
    for spec in policy["source_files"]:
        checked(spec)
    annotation = json.loads(checked(policy["inputs"]["prior_annotation_report"]))
    if policy["tools"] != annotation["tools"]:
        raise ValueError("Gene callers differ from prior actual annotation run")
    for kind in FLAGS:
        actual_flags = {tuple(f["command"][1:f["command"].index("-i")]) for f in annotation["files"] if f["kind"] == kind}
        if actual_flags != {tuple(FLAGS[kind])}:
            raise ValueError("Flags differ from prior actual per-kind commands")
    for report_name, helper_name in [("prior_annotation_report", "annotation_helper"), ("ecor_report", "ecor_helper"), ("phage_report", "phage_helper")]:
        if json.loads(checked(policy["inputs"][report_name]))["script_sha256"] != policy["inputs"][helper_name]["sha256"]:
            raise ValueError("Authentication helper differs from source receipt")
    versions = {}
    for kind, tool in policy["tools"].items():
        provenance = json.loads(checked({"path": tool["provenance_path"], "sha256": tool["provenance_sha256"]}))
        checked({"path": tool["executable"], "sha256": tool["executable_sha256"]})
        if provenance["commit"] != tool["commit"]:
            raise ValueError("Caller source provenance differs")
        version = subprocess.run([tool["executable"], "-v"], capture_output=True, text=True, timeout=10)
        if version.returncode or version.stderr != tool["version_stderr"]:
            raise ValueError("Actual gene-caller version differs")
        versions[kind] = {"stdout": version.stdout, "stderr": version.stderr}
    helpers = [load_module(policy["inputs"][name]) for name in ["annotation_helper", "ecor_helper", "phage_helper"]]
    genomes, controls = policy["genomes"], policy["controls"]
    validate_roster(genomes, controls, json.loads(checked(policy["inputs"]["roster_authority"])))
    workdir.mkdir(parents=True)
    prepared = []
    for item in controls+genomes:
        try:
            prepared.append(prepare(item, workdir/item["id"], helpers))
        except Exception as exc:
            prepared.append({"item": item, "status": "input_error", "error": f"{type(exc).__name__}: {exc}"})
    (workdir/"frozen-call-plan.json").write_bytes(encoded({"policy_sha256": sha(policy_bytes), "wrapper_sha256": sha(code_bytes),
        "versions": versions, "prepared_inputs": prepared, "flags": FLAGS,
        "execution": policy["execution"], "source_files": policy["source_files"]}))
    initial_ok = all(r["status"] == "input_prepared" for r in prepared)
    results = []
    for i, entry in enumerate(prepared):
        control_ok = all(r["status"] == "annotation_validated" and r.get("control_exact_proteome_multiset_match") for r in results[:2])
        if not initial_ok or (i >= 2 and not control_ok):
            results.append({**entry, "status": entry["status"] if entry["status"] == "input_error" else "not_run_input_or_control_failure"})
            continue
        item = entry["item"]
        print("CALL", item["id"], item["kind"], item["source_id"], flush=True)
        results.append(call(entry, policy["tools"][item["kind"]], FLAGS[item["kind"]], policy["execution"]["timeout_seconds_per_genome"], helpers[0]))
    for spec in policy["inputs"].values(): checked(spec)
    for spec in policy["source_files"]: checked(spec)
    if policy_path.read_bytes() != policy_bytes or Path(__file__).read_bytes() != code_bytes:
        raise ValueError("Policy or helper changed during execution")
    full = {"schema_version": 1, "policy_sha256": sha(policy_bytes), "wrapper_sha256": sha(code_bytes), "entries": results}
    manifest = workdir/"annotation-manifest.json"
    manifest.write_bytes(encoded(full))
    public = [{k:v for k,v in r.items() if k != "CDS_counts_by_contig"} for r in results]
    panel = public[2:]
    good = [r for r in panel if r["status"] == "annotation_validated"]
    report = {"schema_version": 1, "status": "all77_references_annotated" if len(good) == 77 else "annotation_incomplete",
        "policy_path": str(policy_path), "policy_sha256": sha(policy_bytes), "wrapper_sha256": sha(code_bytes),
        "source_inputs": policy["inputs"], "tools": policy["tools"], "versions": versions,
        "parser": {"Biopython": Bio.__version__, "Bio_module_sha256": sha(Path(Bio.__file__).read_bytes())},
        "frozen_call_plan": binding(workdir/"frozen-call-plan.json"), "annotation_manifest": binding(manifest),
        "summary": {"requested_hosts": 25, "requested_phages": 52, "validated_genomes": len(good),
                    "status_counts": dict(Counter(r["status"] for r in panel)), "proteins": sum(r["proteins"] for r in good),
                    "contigs_retained": sum(r["input_contigs"] for r in good),
                    "contigs_without_CDS": sum(len(r["contigs_with_no_CDS"]) for r in good)},
        "controls": public[:2], "entries": panel,
        "execution": {**policy["execution"], "platform": platform.platform(), "architecture": platform.machine(),
                      "resource_scope": "Serial caller child timings and cumulative child peak RSS; not total workflow memory"},
        "all_target_related_and_pending_cases_retained": True, "shared_vocabulary_constructed": False,
        "interaction_models_fitted_or_scores_generated": False, "training_membership_admitted": False,
        "Picard_collection_retrieved": False, "source_DNA_concatenated_or_edited": False,
        "FAA_translation_independently_revalidated_against_DNA": False,
        "interpretation": "Per-genome computational gene predictions using prior per-kind callers. Inputs and FAA/GFF coordinate inventories are verified; this is neither biological validation nor source admission. Prodigal single-mode internal training remains within each individual genome."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream: stream.write(encoded(report))
    print(json.dumps({"status": report["status"], "summary": report["summary"], "report_sha256": sha(encoded(report))}))
    return 0 if len(good) == 77 else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/secondary-annotation-policy.json"))
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("data/processed/secondary-reference-annotation.json"))
    args = parser.parse_args()
    raise SystemExit(run(args.policy, args.workdir, args.output))
