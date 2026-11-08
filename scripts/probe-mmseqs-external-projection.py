"""Prospectively fixed focal external-projection diagnostic; no model fitting."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def binding(path):
    return {"path": str(path), "sha256": sha(path)}


def verify(policy):
    for item in policy["inputs"]:
        if sha(ROOT / item["path"]) != item["sha256"]:
            raise ValueError("Diagnostic input differs: " + item["path"])
    if sha(__file__) != policy["helper_sha256"]:
        raise ValueError("Diagnostic helper differs")


def run(policy_path, output):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    verify(policy)
    if output.exists():
        raise ValueError("Use a fresh diagnostic output")
    work = ROOT / policy["workdir"]
    ledger_path = work / "command-ledger.json"
    if ledger_path.exists():
        raise ValueError("Preserve prior attempt; no automatic retry")
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
    contexts = {}
    artifacts = []
    for context in policy["contexts"]:
        path = ROOT / context["tsv"]
        rows = []
        with path.open() as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 11 or fields[0] not in context["query_ids"]:
                    raise ValueError("Unexpected diagnostic alignment row")
                rows.append(fields)
        contexts[context["name"]] = {"raw_rows": rows,
            "query_ids_without_hits": [q for q in context["query_ids"] if not any(r[0] == q for r in rows)]}
        artifacts.append(binding(path))
    verify(policy)
    if policy_path.read_bytes() != policy_bytes:
        raise ValueError("Diagnostic policy changed during execution")
    result = {"schema_version": 1, "status": "completed_focal_external_projection_probe",
        "policy": binding(policy_path), "commands": ledger, "contexts": contexts, "artifacts": artifacts,
        "both_contexts_have_no_focal_hits": all(not c["raw_rows"] for c in contexts.values()),
        "existing_target_database_reused": True, "training_self_search_special_case_used": False,
        "source_proteins_tested": 2, "full_proteome_validated": False,
        "new_vocabulary_or_model_fitted": False, "phenotype_or_native_labels_used": False,
        "failed_control_waived_or_batch_resumed": False,
        "interpretation": "Focal technical probe only. Equal absence is not recovered family assignment, full-proteome consistency, method admission or host-prediction accuracy."}
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
    print(json.dumps({k: result[k] for k in ["status", "both_contexts_have_no_focal_hits", "full_proteome_validated"]}))
