"""Replay completed secondary annotations without executing callers or models."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


class Artifacts:
    """Hash every actual file against its externally frozen provenance."""
    def __init__(self):
        self.files = {}

    def read(self, spec):
        path = Path(spec["path"])
        if path.is_symlink():
            raise ValueError("Symlink artifact is not accepted")
        key = str(path.resolve())
        expected = spec["sha256"]
        if key in self.files and self.files[key]["sha256"] != expected:
            raise ValueError("Conflicting artifact bindings")
        data = path.read_bytes()
        if sha(data) != expected:
            raise ValueError("Artifact bytes differ: "+str(path))
        self.files[key] = {"path": str(path), "sha256": expected}
        return data

    def verify_again(self):
        for spec in list(self.files.values()):
            self.read(spec)


def verify_entry(entry, prepared, producer, helpers, artifacts):
    """Reconstruct DNA normalization and every FAA/GFF coordinate record read-only."""
    parser, ecor, phage = helpers
    item = entry["item"]
    if entry["status"] != "annotation_validated" or entry["exit_code"] != 0:
        raise ValueError("This completed-inventory replay cannot validate an unsuccessful call")
    for key, value in prepared.items():
        if key != "status" and entry.get(key) != value:
            raise ValueError("Call input differs from pre-call frozen plan")
    if prepared["status"] != "input_prepared":
        raise ValueError("Input was not prepared before the call")
    data = artifacts.read(item["fasta"])
    if item["role"] == "technical_control":
        expected_control = parser.protein_multiset(artifacts.read(item["expected_faa"]))
        if parser.multiset_digest(expected_control) != item["expected_proteome_multiset_sha256"]:
            raise ValueError("Original control proteome differs")
    elif item["kind"] == "host":
        master = ecor.parse_master(artifacts.read(item["master_file"]), item["source_id"], item["wgs_accession"], item["project"])
        records = ecor.parse_fasta(data, master)
        stored = json.loads(artifacts.read(item["contig_manifest"]))
        if stored["master"] != master or stored["records"] != records or stored["DNA_concatenated_trimmed_or_designed"]:
            raise ValueError("Source ECOR record inventory differs")
    else:
        official = phage.sequence_record(data, "fasta", item["requested_accession"])
        annotation = phage.sequence_record(artifacts.read(item["annotation"]), item["annotation_format"], item["requested_accession"])
        if ((official["sequence"], official["versioned_accession"]) != (annotation["sequence"], annotation["versioned_accession"])
                or (official["versioned_accession"], official["length"], official["sequence_sha256"]) !=
                (item["versioned_accession"], item["expected_length"], item["expected_sequence_sha256"])):
            raise ValueError("Official phage reference identity differs")
    normalized, mapping = producer.normalize(data, item["id"], parser)
    if normalized != artifacts.read(entry["input"]) or mapping != json.loads(artifacts.read(entry["contig_map"])):
        raise ValueError("Normalized DNA or complete neutral contig mapping differs")
    if (entry["input_contigs"], entry["input_bases"]) != (len(mapping), sum(r["length"] for r in mapping)):
        raise ValueError("Normalized input summary differs")
    faa, gff = artifacts.read(entry["outputs"]["FAA"]), artifacts.read(entry["outputs"]["GFF"])
    artifacts.read(entry["outputs"]["log"])
    summary, inventory = producer.validate_outputs(faa, gff, {r["neutral_id"]: r["length"] for r in mapping}, parser)
    if inventory != json.loads(artifacts.read(entry["protein_coordinate_inventory"])):
        raise ValueError("Protein-coordinate inventory differs from actual FAA/GFF")
    if any(entry.get(key) != value for key, value in summary.items()):
        raise ValueError("Protein/contig summary differs from actual outputs")
    if item["role"] == "technical_control":
        if parser.protein_multiset(faa) != expected_control or not entry["control_exact_proteome_multiset_match"]:
            raise ValueError("Technical control equality does not replay")
    receipt = Path(entry["directory"])/"receipt.json"
    # The per-call receipt has the producer's exact canonical serialization.
    artifacts.read({"path": str(receipt), "sha256": sha(producer.encoded(entry))})
    return {"id": item["id"], "kind": item["kind"], "source_id": item["source_id"], "role": item["role"],
            "input_contigs": len(mapping), "input_bases": entry["input_bases"], "proteins": summary["proteins"],
            "protein_residues": summary["protein_residues"], "proteome_multiset_sha256": summary["proteome_multiset_sha256"],
            "contigs_without_CDS": len(summary["contigs_with_no_CDS"]), "artifacts_replayed": True}, mapping, inventory


def timing_evidence(entries, receipt, spec):
    """Describe original timing evidence; never invent a pause-adjusted wall time."""
    if (receipt.get("action") != "resume_existing_annotation_child_then_existing_wrapper"
            or not receipt.get("no_restart_or_new_caller_submission")
            or not receipt.get("wall_time_in_suspended_call_includes_user_pause")):
        raise ValueError("Unexpected resumption receipt")
    # h013 was identified by PID/command when the user pause was applied.
    paused = [e for e in entries if e["item"]["id"] == "h013"]
    if len(paused) != 1 or paused[0]["item"]["kind"] != "host":
        raise ValueError("Paused-call identity absent or ambiguous")
    record = paused[0]
    return {"resumption_receipt": spec, "recorded_receipt": receipt,
            "affected_call": {"id": "h013", "source_id": record["item"]["source_id"], **record["resources"]},
            "wall_time_includes_user_suspension": True,
            "wall_time_is_compute_benchmark": False, "adjusted_wall_time_calculated": False,
            "pause_interval_independently_reconstructed": False,
            "interpretation": "Preserve exact recorded monotonic wall time and caller CPU counters separately. The receipt has a resume epoch but no complete stop/resume interval; no approximate narrative pause duration is substituted for the measured field."}


def run(report_path, expected_report, resumption_path, expected_resumption, output):
    import Bio
    if Bio.__version__ != "1.86":
        raise ValueError("Require the recorded Biopython1.86 parser")
    if output.exists() or output.is_symlink():
        raise ValueError("Use a new replay report")
    artifacts = Artifacts()
    report_spec = {"path": str(report_path), "sha256": expected_report}
    report = json.loads(artifacts.read(report_spec))
    policy_spec = {"path": report["policy_path"], "sha256": report["policy_sha256"]}
    policy = json.loads(artifacts.read(policy_spec))
    producer_spec = {"path": "scripts/annotate-basel-ecor-references.py", "sha256": report["wrapper_sha256"]}
    artifacts.read(producer_spec)
    if policy["wrapper_sha256"] != producer_spec["sha256"]:
        raise ValueError("Executed producer and frozen policy differ")
    # Import authenticated helpers without calling run/prepare/call or any executable.
    import importlib.util
    definition = importlib.util.spec_from_file_location("secondary_annotation_producer_replay", producer_spec["path"])
    producer = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(producer)
    for item in policy["inputs"].values(): artifacts.read(item)
    for item in policy["source_files"]: artifacts.read(item)
    if report["source_inputs"] != policy["inputs"] or report["tools"] != policy["tools"] or policy["flags"] != producer.FLAGS:
        raise ValueError("Source/helper/tool/flag bindings differ")
    for tool in policy["tools"].values():
        artifacts.read({"path": tool["executable"], "sha256": tool["executable_sha256"]})
        provenance = json.loads(artifacts.read({"path": tool["provenance_path"], "sha256": tool["provenance_sha256"]}))
        if provenance["commit"] != tool["commit"]:
            raise ValueError("Caller provenance changed")
    for kind, tool in policy["tools"].items():
        if report["versions"][kind]["stderr"] != tool["version_stderr"]:
            raise ValueError("Originally observed tool version differs")
    if (report["parser"]["Biopython"] != Bio.__version__
            or report["parser"]["Bio_module_sha256"] != sha(Path(Bio.__file__).read_bytes())):
        raise ValueError("Sequence-record parser bytes differ")
    plan = json.loads(artifacts.read(report["frozen_call_plan"]))
    if (plan["policy_sha256"] != report["policy_sha256"] or plan["wrapper_sha256"] != report["wrapper_sha256"]
            or plan["flags"] != policy["flags"] or plan["execution"] != policy["execution"]
            or plan["source_files"] != policy["source_files"] or plan["versions"] != report["versions"]):
        raise ValueError("Pre-call plan differs from executed provenance")
    full = json.loads(artifacts.read(report["annotation_manifest"]))
    entries = full["entries"]
    if full["policy_sha256"] != report["policy_sha256"] or full["wrapper_sha256"] != report["wrapper_sha256"]:
        raise ValueError("Raw manifest lineage differs")
    expected_items = policy["controls"]+policy["genomes"]
    public = report["controls"]+report["entries"]
    if (len(entries) != 79 or len(public) != 79 or len(plan["prepared_inputs"]) != 79
            or [e["item"] for e in entries] != expected_items
            or [{k:v for k,v in e.items() if k != "CDS_counts_by_contig"} for e in entries] != public):
        raise ValueError("Complete79-call/control inventory differs")
    producer.validate_roster(policy["genomes"], policy["controls"], json.loads(artifacts.read(policy["inputs"]["roster_authority"])))
    helpers = [producer.load_module(policy["inputs"][name]) for name in ["annotation_helper", "ecor_helper", "phage_helper"]]
    summaries, nucleotide_ids, protein_ids = [], set(), set()
    for entry, prepared in zip(entries, plan["prepared_inputs"]):
        kind = entry["item"]["kind"]
        expected_command = [policy["tools"][kind]["executable"], *policy["flags"][kind],
                            "-i", entry["input"]["path"], "-a", entry["outputs"]["FAA"]["path"],
                            "-o", entry["outputs"]["GFF"]["path"]]
        if entry["command"] != expected_command:
            raise ValueError("Recorded per-genome command differs from frozen caller/flags/paths")
        summary, contigs, proteins = verify_entry(entry, prepared, producer, helpers, artifacts)
        nids, pids = {r["neutral_id"] for r in contigs}, {r["protein_id"] for r in proteins}
        if (len(nids) != len(contigs) or len(pids) != len(proteins) or nucleotide_ids & nids or protein_ids & pids):
            raise ValueError("Duplicate neutral nucleotide/protein IDs within or across genomes")
        nucleotide_ids.update(nids); protein_ids.update(pids); summaries.append(summary)
    panel, controls = summaries[2:], summaries[:2]
    observed = {"requested_hosts": 25, "requested_phages": 52, "validated_genomes": len(panel),
                "status_counts": {"annotation_validated": len(panel)}, "proteins": sum(r["proteins"] for r in panel),
                "contigs_retained": sum(r["input_contigs"] for r in panel),
                "contigs_without_CDS": sum(r["contigs_without_CDS"] for r in panel)}
    if observed != report["summary"]:
        raise ValueError("Public panel totals differ from actual artifacts")
    resume_spec = {"path": str(resumption_path), "sha256": expected_resumption}
    timing = timing_evidence(entries, json.loads(artifacts.read(resume_spec)), resume_spec)
    artifacts.verify_again()
    result = {"schema_version": 1, "status": "all77_secondary_reference_annotations_replayed_offline",
              "input_report": report_spec, "policy": policy_spec, "executed_producer": producer_spec,
              "replay_helper_sha256": sha(Path(__file__).read_bytes()), "parser": "Biopython1.86",
              "summary": observed, "per_kind": {kind: {"genomes": sum(r["kind"] == kind for r in panel),
                    "proteins": sum(r["proteins"] for r in panel if r["kind"] == kind)} for kind in ["host", "phage"]},
              "controls": controls, "entries": panel, "timing_evidence": timing,
              "verified_unique_artifact_files": len(artifacts.files), "artifacts": sorted(artifacts.files.values(), key=lambda r:r["path"]),
              "nucleotide_ids_globally_unique": len(nucleotide_ids), "protein_ids_globally_unique": len(protein_ids),
              "gene_caller_model_or_network_execution": False, "source_outputs_modified": False,
              "training_membership_admitted": False, "FAA_DNA_translation_independently_revalidated": False,
              "interpretation": "Offline input-preservation and annotation-artifact consistency replay. Gene-caller compatibility and coordinate correspondence are not evidence of host specificity, binding, assembly or experimental-stock identity."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream: stream.write(producer.encoded(result))
    print(json.dumps({"status": result["status"], "summary": observed, "verified_unique_artifact_files": len(artifacts.files),
                      "replay_report_sha256": sha(producer.encoded(result))}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=Path("data/processed/secondary-reference-annotation.json"))
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument("--resumption-receipt", type=Path, default=Path("data/raw/genophi-external-validation/secondary-annotation-v1/resumption-receipt.json"))
    parser.add_argument("--expected-resumption-sha256", required=True)
    parser.add_argument("--output", type=Path, default=Path("data/processed/secondary-reference-annotation-replay.json"))
    args = parser.parse_args()
    run(args.report, args.expected_report_sha256, args.resumption_receipt, args.expected_resumption_sha256, args.output)
