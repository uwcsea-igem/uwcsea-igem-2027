"""One fixed uniform-caller control for the label-free host encoder pilot."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def execute(policy_path, workdir, output):
    policy_bytes, script_hash = policy_path.read_bytes(), sha(__file__)
    policy = json.loads(policy_bytes)
    if workdir.exists() or output.exists() or workdir.is_symlink() or output.is_symlink():
        raise ValueError("Use fresh control destinations")
    if output.resolve().is_relative_to(workdir.resolve()):
        raise ValueError("Report must be outside raw workdir")
    bindings = {Path(x["path"]): x["sha256"] for x in policy["inputs"] + policy["source_genomes"] + [policy["caller"]]}
    def verify():
        for path, expected in bindings.items():
            if sha(path) != expected:
                raise ValueError("Frozen uniform-control input changed")
        if policy_path.read_bytes() != policy_bytes or sha(__file__) != script_hash:
            raise ValueError("Control policy or wrapper changed")
    verify()
    mixed_path = Path("software/tipcraft/benchmarks/receptor-resource/host-feature-distinguishability-policy.json")
    pilot_policy = json.loads(mixed_path.read_bytes())
    if pilot_policy["fit_parameters"] != policy["feature_parameters"]:
        raise ValueError("Uniform control cannot change feature parameters")
    historical = json.loads(Path("data/raw/genophi-annotation/run-1/report.json").read_bytes())
    if historical["tools"]["host"]["executable_sha256"] != policy["caller"]["sha256"]:
        raise ValueError("Caller differs from original target annotations")
    for x in historical["files"]:
        if x["kind"] == "host" and x["role"] == "target_reference":
            if x["command"][1:8] != policy["caller"]["parameters"]:
                raise ValueError("Caller parameters differ from original target annotations")
    reader = load("scripts/audit-receptor-host-payloads.py", "uniform_host_reader")
    pilot = load("scripts/check-receptor-host-features.py", "uniform_host_pilot")
    workdir.mkdir(parents=True)
    commands = []
    for index, source in enumerate(policy["source_genomes"]):
        folder = workdir / f"author-{index}"
        folder.mkdir()
        argv = [policy["caller"]["path"], *policy["caller"]["parameters"], "-i", source["path"],
                "-a", str(folder/"proteins.faa"), "-o", str(folder/"genes.gff")]
        commands.append({"sample_id": source["sample_id"], "argv": argv, "folder": str(folder)})
    dump(workdir/"frozen-plan.json", {"policy_sha256": hashlib.sha256(policy_bytes).hexdigest(),
                                     "wrapper_sha256": script_hash, "commands": commands,
                                     "feature_parameters": policy["feature_parameters"]})
    call_records = []
    for source, item in zip(policy["source_genomes"], commands):
        start = time.monotonic()
        result = subprocess.run(item["argv"], text=True, capture_output=True, timeout=120)
        folder = Path(item["folder"])
        (folder/"run.log").write_text(result.stdout+result.stderr)
        if result.returncode:
            raise ValueError("Uniform gene call failed; retained log/source must not be repaired")
        fasta = folder/"proteins.faa"
        records = reader.fasta(fasta.read_bytes(), set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"))
        call_records.append({"sample_id": source["sample_id"], "source_genome": source, "argv": item["argv"],
                             "path": str(fasta), "sha256": sha(fasta), "expected_proteins": len(records),
                             "gff_sha256": sha(folder/"genes.gff"), "log_sha256": sha(folder/"run.log"),
                             "wall_seconds": time.monotonic()-start})
    verify()
    call_report = workdir/"uniform-calls.json"
    dump(call_report, {"policy_sha256": hashlib.sha256(policy_bytes).hexdigest(), "wrapper_sha256": script_hash,
                       "calls": call_records, "interaction_models_trained_or_scored": False})
    pilot_policy["purpose"] = "Label-free host encoder capacity under a uniform host caller; not infection prediction"
    pilot_policy["timing"] = "Control parameters fixed before calls; generated protein hashes/counts bound before feature fitting"
    pilot_policy["annotation_limit"] = "All5hostreferences use Prodigal2.6.3 single-mode table11. Stock equivalence remains unverified."
    pilot_policy["files"] = [{k:x[k] for k in ("sample_id", "path", "sha256", "expected_proteins")} |
                             {"role": "vocabulary_fit_only"} for x in call_records] + [
        x for x in pilot_policy["files"] if x["role"] == "projection_only"]
    pilot_policy["inputs"] += [{"path": str(policy_path), "sha256": hashlib.sha256(policy_bytes).hexdigest()},
                               {"path": str(call_report), "sha256": sha(call_report)},
                               {"path": __file__, "sha256": script_hash}]
    generated_policy = workdir/"uniform-feature-policy.json"
    dump(generated_policy, pilot_policy)
    pilot.execute(generated_policy, workdir/"feature-pilot", output)
    verify()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/host-feature-uniform-control-policy.json"))
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    execute(args.policy, args.workdir, args.output)
