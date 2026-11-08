"""Replay every saved copied-proteome control with actual identity mode.

Diagnostic only: no phenotype inputs, new clustering, vocabulary admission,
interaction fitting, native scoring or continuation of the failed batch.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

from tipcraft import protein_features

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            value.update(chunk)
    return value.hexdigest()


def bound(path):
    return {"path": str(path), "sha256": sha(path)}


def verify(policy):
    for value in policy["inputs"]:
        if sha(ROOT / value["path"]) != value["sha256"]:
            raise ValueError("Diagnostic input changed: " + value["path"])
    if sha(__file__) != policy["helper_sha256"] or sha(protein_features.__file__) != policy["kernel_sha256"]:
        raise ValueError("Executed probe/kernel differs from policy")


def run(policy_path, output):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    verify(policy)
    if output.exists():
        raise ValueError("Use a fresh diagnostic output")
    work = ROOT / policy["workdir"]
    ledger_path = work / "command-ledger.json"
    if ledger_path.exists():
        raise ValueError("Do not reuse an existing diagnostic attempt")
    parser_path = ROOT / policy["protein_parser"]["path"]
    spec = importlib.util.spec_from_file_location("pinned_probe_protein_parser", parser_path)
    parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parser)
    records = json.loads((ROOT / policy["protein_map"]["path"]).read_text())
    controls = [r for r in records if r["role"] == "mapping_control"]
    if len(controls) != policy["copied_protein_count"]:
        raise ValueError("Complete copied proteome required")
    training = dict(parser.proteins((ROOT / policy["training_faa"]["path"]).read_bytes()))
    queries = dict(parser.proteins((ROOT / policy["query_faa"]["path"]).read_bytes()))
    for record in controls:
        original, copy = record["copied_training_protein"], record["protein"]
        if (training[original] != queries[copy] or len(queries[copy]) != record["length"]
                or hashlib.sha256(queries[copy].encode()).hexdigest() != record["sequence_sha256"]):
            raise ValueError("Actual copied protein sequence differs")
    ledger = []
    for step in policy["commands"]:
        log = work / (step["id"] + ".log")
        print("RUN", step["id"], flush=True)
        with log.open("x") as stream:
            result = subprocess.run(step["argv"], stdout=stream, stderr=subprocess.STDOUT,
                                    timeout=policy["timeout_seconds_per_command"])
        ledger.append({"id": step["id"], "argv": step["argv"], "exit_code": result.returncode, "log": bound(log)})
        ledger_path.write_text(json.dumps(ledger, indent=2) + "\n")
        if result.returncode != 0:
            raise ValueError("Diagnostic command failed; no automatic retry")
    lengths = {name: len(sequence) for name, sequence in training.items()}
    families = json.loads((ROOT / policy["vocabulary"]["path"]).read_text())["training_protein_families"]
    hits, diagnostics = {}, {}
    for context in policy["contexts"]:
        key = context["context"]
        subset = {r["copied_training_protein"] if key == "self" else r["protein"]: r["length"] for r in controls}
        with (ROOT / context["raw_tsv"]).open() as stream:
            hits[key], diagnostics[key] = protein_features.best_hits(stream, subset, lengths, families)
    mismatches = []
    for record in controls:
        a = hits["self"].get(record["copied_training_protein"])
        b = hits["external"].get(record["protein"])
        if a is None or b is None or a["family"] != b["family"]:
            mismatches.append({"training_protein": record["copied_training_protein"], "copied_query": record["protein"],
                               "training_hit": a, "external_hit": b})
    verify(policy)
    if policy_path.read_bytes() != policy_bytes:
        raise ValueError("Probe policy changed during execution")
    report = {"schema_version": 1, "status": "completed_copied_proteome_identity_probe",
              "policy": bound(policy_path), "copied_proteins": len(controls), "actual_sequence_copies_verified": True,
              "alignment_mode": 3, "contexts": {k: {"queries_with_hits": len(hits[k]),
                  "raw_rows": diagnostics[k]["raw_rows"],
                  "cross_family_best_hit_tie_queries": len(diagnostics[k]["cross_family_best_hit_ties"])} for k in hits},
              "family_mismatches": mismatches, "family_mismatch_count": len(mismatches),
              "all_copied_families_equal": not mismatches, "commands": ledger,
              "artifacts": [bound(ROOT / c["raw_tsv"]) for c in policy["contexts"]],
              "existing_v1_clusters_reused": True, "new_vocabulary_fitted": False,
              "independent_new_proteome_tested": False, "full_feature_pipeline_validated": False,
              "original_failed_batch_resumed_or_control_waived": False,
              "phenotype_or_native_labels_used": False, "model_fitted_or_scored": False,
              "interpretation": "Technical consistency replay of the existing copied proteome against frozen candidates, not host-prediction performance or full workflow admission."}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.policy, args.output)
    print(json.dumps({k: result[k] for k in ["status", "copied_proteins", "family_mismatch_count", "all_copied_families_equal"]}))
