"""Compare a complete copied proteome under ordinary external projection.

The saved target families/vocabulary and copied-query search are reused.
This diagnostic does not admit a revised production procedure.
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
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def binding(path):
    return {"path": str(path), "sha256": sha(path)}


def compare_assignments(controls, source_hits, copy_hits):
    comparisons, assigned, missing = [], 0, 0
    originals = [r["copied_training_protein"] for r in controls]
    copies = [r["protein"] for r in controls]
    if (not controls or len(set(originals)) != len(controls) or len(set(copies)) != len(controls)
            or not set(source_hits) <= set(originals) or not set(copy_hits) <= set(copies)):
        raise ValueError("Complete unique copied roster and in-roster hits required")
    for row in controls:
        original, copy = row["copied_training_protein"], row["protein"]
        left, right = source_hits.get(original), copy_hits.get(copy)
        if left is None and right is None:
            state = "both_unassigned"
            missing += 1
        elif left is not None and right is not None and left["family"] == right["family"]:
            state = "same_assigned_family"
            assigned += 1
        else:
            state = "assignment_discrepancy"
        comparisons.append({"training_protein": original, "copied_protein": copy, "state": state,
                            "training_family": left["family"] if left else None,
                            "copy_family": right["family"] if right else None})
    return {"proteins": len(controls), "same_assigned_family": assigned,
            "both_unassigned": missing, "assignment_discrepancies": len(controls) - assigned - missing,
            "comparisons": comparisons}


def run(policy_path, output):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    def verify():
        for item in policy["inputs"]:
            if sha(ROOT / item["path"]) != item["sha256"]:
                raise ValueError("Diagnostic input differs: " + item["path"])
        if (sha(__file__) != policy["helper_sha256"]
                or sha(protein_features.__file__) != policy["kernel_sha256"]):
            raise ValueError("Diagnostic implementation differs")
    verify()
    if output.exists():
        raise ValueError("Use a fresh diagnostic output")
    work = ROOT / policy["workdir"]
    ledger_path = work / "command-ledger.json"
    if ledger_path.exists():
        raise ValueError("Preserve existing attempt; no automatic retry")
    parser_spec = importlib.util.spec_from_file_location("symmetric_probe_parser", ROOT / policy["protein_parser"])
    parser = importlib.util.module_from_spec(parser_spec)
    parser_spec.loader.exec_module(parser)
    training = dict(parser.proteins((ROOT / policy["training_faa"]).read_bytes()))
    queries = dict(parser.proteins((ROOT / policy["query_faa"]).read_bytes()))
    source = dict(parser.proteins((work / "source.faa").read_bytes()))
    controls = [r for r in json.loads((ROOT / policy["protein_map"]).read_text()) if r["role"] == "mapping_control"]
    if len(controls) != policy["copied_proteins"] or set(source) != {r["copied_training_protein"] for r in controls}:
        raise ValueError("Exact complete copied proteome required")
    for row in controls:
        original, copy = row["copied_training_protein"], row["protein"]
        if (source[original] != training[original] or training[original] != queries[copy]
                or len(source[original]) != row["length"]
                or hashlib.sha256(source[original].encode()).hexdigest() != row["sequence_sha256"]):
            raise ValueError("Actual source/copy sequence equality failed")
    ledger = []
    for step in policy["commands"]:
        log = work / (step["id"] + ".log")
        print("RUN", step["id"], flush=True)
        with log.open("x") as stream:
            result = subprocess.run(step["argv"], stdout=stream, stderr=subprocess.STDOUT,
                                    timeout=policy["timeout_seconds_per_command"])
        ledger.append({**step, "exit_code": result.returncode, "log": binding(log)})
        ledger_path.write_text(json.dumps(ledger, indent=2) + "\n")
        if result.returncode != 0:
            raise ValueError("Diagnostic command failed; attempt preserved")
    vocabulary = json.loads((ROOT / policy["vocabulary"]).read_text())
    families = vocabulary["training_protein_families"]
    lengths = {name: len(sequence) for name, sequence in training.items()}
    copy_lengths = {r["protein"]: r["length"] for r in controls}
    with (work / "source.tsv").open() as stream:
        source_hits, source_diagnostics = protein_features.best_hits(
            stream, {n: len(s) for n, s in source.items()}, lengths, families)
    with (ROOT / policy["saved_query_tsv"]).open() as stream:
        copy_hits, copy_diagnostics = protein_features.best_hits(
            (line for line in stream if line.split("\t", 1)[0] in copy_lengths), copy_lengths, lengths, families)
    result = compare_assignments(controls, source_hits, copy_hits)
    _, left = protein_features.project(vocabulary, ["source"], {n: "source" for n in source}, source_hits)
    _, right = protein_features.project(vocabulary, ["copy"], {n: "copy" for n in copy_lengths}, copy_hits)
    result.update({"schema_version": 1, "status": "completed_complete_proteome_external_projection_probe",
        "policy": binding(policy_path), "commands": ledger,
        "contexts": {"training_external": source_diagnostics, "saved_copy_external": copy_diagnostics},
        "existing_vocabulary_vector_equal": left["source"] == right["copy"],
        "existing_vocabulary_vectors": {"training_external": left["source"], "copy_external": right["copy"]},
        "actual_amino_acid_copies_verified": True,
        "source_search_newly_executed": True, "copy_search_reused_from_failed_v2_task": True,
        "new_vocabulary_fitted": False, "full_production_pipeline_validated": False,
        "original_failed_control_waived_or_batch_resumed": False,
        "phenotype_or_native_labels_used": False, "biological_model_fitted_or_scored": False,
        "artifacts": [binding(work / "source.faa"), binding(work / "source.tsv")],
        "interpretation": "Complete copied-proteome diagnostic under a fixed existing vocabulary. Equal absence is reported separately from successful family assignments. This is not a rebuilt production vocabulary, independent new proteome, biological accuracy or authorization to waive a failed v2 control."})
    verify()
    if policy_path.read_bytes() != policy_bytes:
        raise ValueError("Diagnostic policy changed during execution")
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.policy, args.output)
    print(json.dumps({k: result[k] for k in ["status", "proteins", "same_assigned_family", "both_unassigned",
                                           "assignment_discrepancies", "existing_vocabulary_vector_equal"]}))
