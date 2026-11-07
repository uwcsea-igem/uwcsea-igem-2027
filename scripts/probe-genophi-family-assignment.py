#!/usr/bin/env python3
"""Execute a pinned upstream projection function on synthetic data only."""
from __future__ import annotations

import argparse
import ast
import contextlib
import hashlib
import io
import json
import logging
import os
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "software/tipcraft/benchmarks/receptor-resource/family-assignment-probe-policy.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "data/raw/receptor-resource/genophi-assignment-source-v1/kmer_assign_features_workflow.py")
    parser.add_argument("--output", type=Path, default=ROOT / "data/processed/genophi-family-assignment-probe.json")
    args = parser.parse_args()
    policy_bytes = POLICY.read_bytes()
    policy = json.loads(policy_bytes)
    source = args.source.read_bytes()
    digest = hashlib.sha256(source).hexdigest()
    blob = hashlib.sha1(b"blob " + str(len(source)).encode() + b"\0" + source).hexdigest()
    if digest != policy["source_sha256"] or blob != policy["source_git_blob"]:
        raise ValueError("Pinned upstream source mismatch")

    import Bio
    from Bio import SeqIO
    import pandas as pd
    import tqdm as tqdm_package
    from tqdm import tqdm
    from tqdm.auto import tqdm as tqdm_auto

    tree = ast.parse(source)
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in policy["functions"]]
    if [node.name for node in nodes] != policy["functions"]:
        raise ValueError("Expected upstream function nodes missing or reordered")
    namespace = {"pd": pd, "SeqIO": SeqIO, "tqdm": tqdm, "tqdm_auto": tqdm_auto,
                 "logging": logging, "os": os}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(args.source), "exec"), namespace)

    sequences = {"both_a": "AAAAA", "both_b": "CCCCC", "partial_a": "AAAAA",
                 "negative_a": "TTTTT", "unassigned": "DDDDD"}
    mapping = {"both_a": "both", "both_b": "both", "partial_a": "partial",
               "negative_a": "negative", "unassigned": "unassigned"}
    capture = io.StringIO()
    with tempfile.TemporaryDirectory(prefix="tipcraft-assignment-probe-") as tmp:
        tmp = Path(tmp)
        fasta = tmp / "synthetic.faa"
        fasta.write_text("".join(f">{key}\n{seq}\n" for key, seq in sequences.items()))
        best = tmp / "best.tsv"
        best.write_text("both_a\talpha\nboth_b\tbeta\npartial_a\talpha\nnegative_a\talpha\n")
        features = tmp / "features.csv"
        pd.DataFrame({"Feature": ["f0", "f0"], "Cluster_Label": ["alpha_AAAAA", "beta_CCCCC"]}).to_csv(features, index=False)
        kmers = tmp / "kmers.csv"
        pd.DataFrame({"protein_family": ["alpha", "beta"], "Feature": ["f0", "f0"],
                      "cluster": ["alpha_AAAAA", "beta_CCCCC"], "kmer": ["AAAAA", "CCCCC"]}).to_csv(kmers, index=False)
        with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(capture):
            result = namespace["map_features_with_kmers_and_sequences"](
                str(best), str(features), str(kmers), mapping, "strain", [str(fasta)],
                str(fasta), policy["threshold"], threads=policy["threads"])
        if result is None:
            raise RuntimeError("Upstream synthetic projection failed: " + capture.getvalue())
        values = {str(row["strain"]): int(row["f0"]) for row in result.to_dict(orient="records")}
    expected = {"both": 1, "partial": 1, "negative": 0}
    if values != expected:
        raise ValueError(f"Upstream behavior changed: {values}")
    report = {
        "schema_version": 1, "status": "executed_synthetic_contract_probe",
        "policy_sha256": hashlib.sha256(policy_bytes).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": digest, "source_git_blob": blob,
        "upstream_commit": policy["upstream_commit"], "upstream_path": policy["upstream_path"],
        "executed_functions": policy["functions"],
        "dependencies": {"pandas": pd.__version__, "biopython": Bio.__version__, "tqdm": tqdm_package.__version__},
        "synthetic_inputs": {"sequences": sequences, "protein_genome_map": mapping,
                             "threshold": policy["threshold"], "members": ["alpha_AAAAA", "beta_CCCCC"]},
        "upstream_output": values, "omitted_genomes": sorted(set(mapping.values()) - set(values)),
        "fixed_vocabulary_denominator_illustration": {"both": 1, "partial": 0, "negative": 0, "unassigned": 0},
        "interpretation": [
            "Only members reached through assigned protein families enter the upstream denominator.",
            "A missing family is absent from that denominator, rather than automatically counted as a failed vocabulary member.",
            "A genome with no matched family rows is omitted by this function, rather than explicitly returned as an all-zero vector.",
            "The fixed-vocabulary denominator is a different illustrative rule, not a certified correction or the paper's measured inference protocol.",
            "No biological prediction, accuracy, threshold optimization or whole-workflow equivalence is established."
        ],
        "biological_sequences_or_phenotypes_loaded": False, "model_fitted": False
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "output": str(args.output), "values": values,
                      "omitted": report["omitted_genomes"]}))


if __name__ == "__main__":
    main()
