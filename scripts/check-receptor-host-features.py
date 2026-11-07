"""Measure a declared host encoder's input collisions without phenotype labels."""
import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time

from tipcraft.kmer_features import fit_vocabulary, project_proteomes
import tipcraft.kmer_features as kernel


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def collisions(projection):
    groups = defaultdict(list)
    for name, row in zip(projection["sample_ids"], projection["matrix"]):
        groups[tuple(row)].append(name)
    return [{"sample_ids": names, "vector": list(vector)} for vector, names in groups.items() if len(names) > 1]


def execute(policy_path, workdir, output):
    policy_bytes, script_bytes = policy_path.read_bytes(), Path(__file__).read_bytes()
    policy = json.loads(policy_bytes)
    if workdir.exists() or output.exists() or workdir.is_symlink() or output.is_symlink():
        raise ValueError("Use fresh workdir/report")
    if output.resolve().is_relative_to(workdir.resolve()):
        raise ValueError("Report must be outside raw workdir")
    bindings = {Path(x["path"]): x["sha256"] for x in policy["inputs"] + policy["files"]}
    for path, expected in bindings.items():
        if sha(path.read_bytes()) != expected:
            raise ValueError("Frozen diagnostic input changed")
    if Path(kernel.__file__).resolve() not in {p.resolve() for p in bindings}:
        raise ValueError("Imported feature kernel is not frozen")
    spec = importlib.util.spec_from_file_location("host_inventory_reader", "scripts/audit-receptor-host-payloads.py")
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    report = json.loads(Path("data/raw/genophi-annotation/run-1/report.json").read_bytes())
    original = json.loads(Path("software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json").read_bytes())
    annotated = {x["source_id"]: x for x in report["files"] if x["kind"] == "host" and x["role"] == "target_reference"}
    proteomes, inventory = {}, []
    for item in policy["files"]:
        records = reader.fasta(Path(item["path"]).read_bytes(), set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"))
        proteins, stops = {}, 0
        for header, source in records:
            sequence = source.upper()
            if sequence.endswith("*"):
                sequence = sequence[:-1]
                stops += 1
            if not sequence or "*" in sequence:
                raise ValueError("Missing protein or internal stop; no silent repair")
            proteins[header.split()[0]] = sequence
        if len(proteins) != item["expected_proteins"]:
            raise ValueError("Source protein count differs")
        name = item["sample_id"]
        if name in proteomes:
            raise ValueError("Duplicate sample")
        if item["role"] == "projection_only":
            target = name.split(":", 1)[1]
            if (annotated[target]["faa_sha256"] != item["sha256"]
                    or annotated[target]["source_sequence_sha256"] != item["original_genome_sequence_sha256"]
                    or original["genome_sources"][target]["sequence_sha256"] != item["original_genome_sequence_sha256"]):
                raise ValueError("Target protein/genome provenance differs")
        proteomes[name] = proteins
        inventory.append({"sample_id": name, "role": item["role"], "proteins": len(proteins),
                          "terminal_stop_markers_removed": stops, "normalized_proteome_sha256": sha(encoded(proteins))})
    training = {x["sample_id"]: proteomes[x["sample_id"]] for x in policy["files"] if x["role"] == "vocabulary_fit_only"}
    targets = {x["sample_id"]: proteomes[x["sample_id"]] for x in policy["files"] if x["role"] == "projection_only"}
    if set(training) != {"author:BL21", "author:Keio"} or set(targets) != {"original:BL21", "original:BW25113", "original:MG1655"}:
        raise ValueError("Require exactly two fit-only and three projection-only inputs")
    workdir.mkdir(parents=True)
    plan = {"policy_sha256": sha(policy_bytes), "script_sha256": sha(script_bytes), "inputs": policy["inputs"],
            "files": policy["files"], "fit_parameters": policy["fit_parameters"], "proteome_inventory": inventory}
    (workdir / "frozen-plan.json").write_bytes(encoded(plan))
    start = time.monotonic()
    vocabulary = fit_vocabulary(training, **policy["fit_parameters"])
    training_projection = project_proteomes(vocabulary, training)
    target_projection = project_proteomes(vocabulary, targets)
    features = vocabulary.features
    if len(features) > 2**len(training)-1:
        raise ValueError("Merged nonzero occurrence-pattern count exceeds its theoretical bound")
    artifact = workdir / "vocabulary.json"
    artifact.write_bytes(encoded({"training_ids": vocabulary.training_ids, "features": features,
                                  "sequence_policy": vocabulary.sequence_policy}))
    collisions_target = collisions(target_projection)
    all_identical = len({tuple(row) for row in target_projection["matrix"]}) == 1
    for path, expected in bindings.items():
        if sha(path.read_bytes()) != expected:
            raise ValueError("Diagnostic input changed during fitting/projection")
    if policy_path.read_bytes() != policy_bytes or Path(__file__).read_bytes() != script_bytes:
        raise ValueError("Policy or wrapper changed during diagnostic")
    result = {"schema_version": 1, **plan, "status": "completed_label_free_host_representation_diagnostic",
              "retained_kmer_count": vocabulary.retained_kmer_count, "merged_feature_count": len(features),
              "nonzero_pattern_upper_bound": 2**len(training)-1,
              "feature_member_counts": [len(group) for group in features],
              "training_projection": training_projection, "target_projection": target_projection,
              "target_collision_groups": collisions_target, "all_original_target_host_vectors_identical": all_identical,
              "vocabulary": {"path": str(artifact), "sha256": sha(artifact.read_bytes())},
              "elapsed_seconds": time.monotonic()-start, "peak_process_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024),
              "phenotype_labels_loaded": False, "interaction_model_trained_or_scored": False,
              "source_stock_identity_verified": False, "new_benchmark_rows_admitted": 0,
              "interpretation": "A label-free capacity diagnostic of one declared encoder on these mixed-provenance annotated reference inputs. Identical host vectors cannot yield different host scores for a fixed phage in a predictor using only these sequence features. Collisions do not establish biological equivalence or measure infection prediction accuracy."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("retained_kmer_count", "merged_feature_count", "target_projection", "all_original_target_host_vectors_identical", "elapsed_seconds", "peak_process_rss_bytes")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/host-feature-distinguishability-policy.json"))
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    execute(args.policy, args.workdir, args.output)
