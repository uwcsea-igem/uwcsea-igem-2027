"""Replay authenticated phage alignments through the reusable feature kernel.

This verifies a computational component. It neither executes MMseqs nor builds
outer-fold host vocabularies, fits an interaction model or creates target scores.
"""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import resource
import sys
import time

from tipcraft import protein_features


SOURCE_REPORT = Path("data/processed/genophi-phage-feature-pilot-verified.json")
SOURCE_REPORT_SHA256 = "b8a4ea98682d64fe1bc6cb84bf5976fc578170d41745be7a6859cacf9e0bd7e8"
PARSER = Path("scripts/check-genophi-phage-features.py")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def authenticated_inputs():
    if sha(SOURCE_REPORT) != SOURCE_REPORT_SHA256:
        raise ValueError("Original verified pilot report differs")
    report = json.loads(SOURCE_REPORT.read_text())
    if sha(PARSER) != report["script_sha256"]:
        raise ValueError("Historical strict FASTA parser differs from verified execution")
    policy_path = Path(report["policy_path"])
    if sha(policy_path) != report["policy_sha256"]:
        raise ValueError("Original pilot policy differs")
    inputs = {SOURCE_REPORT: SOURCE_REPORT_SHA256, PARSER: report["script_sha256"],
              policy_path: report["policy_sha256"]}
    files = {}
    for item in report["artifacts"]:
        path = Path(item["path"])
        if path.name in files or sha(path) != item["sha256"]:
            raise ValueError("Duplicate or changed original pilot artifact")
        files[path.name] = path
        inputs[path] = item["sha256"]
    if report["phenotype_tables_loaded"] or report["interaction_models_fitted"] or report["interaction_scores_generated"]:
        raise ValueError("Expected the no-interaction-model representation pilot")
    policy = json.loads(policy_path.read_text())
    return report, policy, files, inputs


def parser():
    spec = importlib.util.spec_from_file_location("historical_protein_parser", PARSER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.proteins


def checked_proteins(path, records, parse):
    seqs = dict(parse(path.read_bytes()))
    if set(seqs) != set(records):
        raise ValueError("FASTA IDs disagree with authenticated protein mapping")
    for name, sequence in seqs.items():
        if (len(sequence) != records[name]["length"]
                or hashlib.sha256(sequence.encode()).hexdigest() != records[name]["sequence_sha256"]):
            raise ValueError("Protein sequence differs from authenticated mapping")
    return {name: len(sequence) for name, sequence in seqs.items()}


def execute(workdir, output):
    start = time.monotonic()
    wrapper_hash, kernel_hash = sha(Path(__file__)), sha(Path(protein_features.__file__))
    source, policy, files, inputs = authenticated_inputs()
    for destination in (workdir, output):
        if destination.exists() or destination.is_symlink():
            raise ValueError("Output locations must be new")
        for path in files.values():
            if (destination.resolve().is_relative_to(path.parent.resolve())
                    or path.resolve().is_relative_to(destination.resolve())):
                raise ValueError("Output must be separate from original input artifact directories")
    if output.resolve().is_relative_to(workdir.resolve()) or workdir.resolve().is_relative_to(output.resolve()):
        raise ValueError("Summary output must be separate from raw replay output")
    records = json.loads(files["protein-map.json"].read_text())
    if len({r["protein"] for r in records}) != len(records):
        raise ValueError("Duplicate authenticated protein ID")
    train = {r["protein"]: r for r in records if r["role"] == "train"}
    query = {r["protein"]: r for r in records if r["role"] in {"target", "control"}}
    if len(train) + len(query) != len(records) or set(train) & set(query):
        raise ValueError("Invalid protein roles")
    train_genomes = policy["training_phages"]
    query_genomes = policy["target_phages"] + ["mapping_control"]
    train_map = {name: r["genome"] for name, r in train.items()}
    query_map = {name: r["genome"] for name, r in query.items()}
    if set(train_map.values()) != set(train_genomes) or set(query_map.values()) != set(query_genomes):
        raise ValueError("Actual genome roles disagree with the pilot policy")
    parse = parser()
    train_lengths = checked_proteins(files["training.faa"], train, parse)
    with files["clusters.tsv"].open() as stream:
        families = protein_features.read_families(stream, train_lengths)
    with files["self.tsv"].open() as stream:
        self_hits, self_diagnostic = protein_features.best_hits(stream, train_lengths, train_lengths, families)
    if families != json.loads(files["member-family-map.json"].read_text()):
        raise ValueError("Kernel cluster membership differs from the original pilot")
    if self_hits != json.loads(files["best-self.json"].read_text()) or set(self_hits) != set(train_lengths):
        raise ValueError("Kernel training hit assignments differ or lack full self coverage")
    vocabulary = protein_features.fit_vocabulary(train_genomes, train_map, families, self_hits)
    old_map = json.loads(files["feature-map.json"].read_text())
    if (vocabulary["family_features"] != {r["family"]: r["feature"] for r in old_map}
            or vocabulary["retained_family_patterns"] != {r["family"]: r["training_occurrence"] for r in old_map}):
        raise ValueError("Train-only occurrence features differ from the original frozen vocabulary")
    workdir.mkdir(parents=True, exist_ok=False)
    vocabulary_path = workdir / "vocabulary.json"
    dump(vocabulary_path, vocabulary)
    frozen_hash = sha(vocabulary_path)
    dump(workdir / "pre-query-freeze.json", {
        "vocabulary_sha256": frozen_hash, "kernel_sha256": kernel_hash,
        "query_alignment_parsing_started": False,
        "training_genomes": train_genomes,
        "training_inputs": [{"path": str(files[name]), "sha256": inputs[files[name]]}
                            for name in ["training.faa", "clusters.tsv", "self.tsv"]],
    })
    # Query occurrence and best-hit processing cannot affect the fitted vocabulary.
    query_lengths = checked_proteins(files["queries.faa"], query, parse)
    with files["query.tsv"].open() as stream:
        query_hits, query_diagnostic = protein_features.best_hits(stream, query_lengths, train_lengths, families)
    if query_hits != json.loads(files["best-query.json"].read_text()):
        raise ValueError("Kernel query hit assignments differ from the original pilot")
    names = train_genomes + query_genomes
    features, vectors = protein_features.project(vocabulary, names, train_map | query_map, self_hits | query_hits)
    if vocabulary != json.loads(vocabulary_path.read_text()) or sha(vocabulary_path) != frozen_hash:
        raise ValueError("Vocabulary changed during query projection")
    matrix = io.StringIO(newline="")
    writer = csv.writer(matrix, delimiter="\t", lineterminator="\n")
    writer.writerow(["phage", *features])
    for genome in names:
        writer.writerow([genome, *vectors[genome]])
    matrix_bytes = matrix.getvalue().encode()
    if matrix_bytes != files["feature-matrix.tsv"].read_bytes():
        raise ValueError("Reusable-kernel matrix differs byte-for-byte from the authenticated pilot")
    (workdir / "feature-matrix.tsv").write_bytes(matrix_bytes)
    for path, expected in inputs.items():
        if sha(path) != expected:
            raise ValueError("Original input changed during replay")
    if sha(Path(__file__)) != wrapper_hash or sha(Path(protein_features.__file__)) != kernel_hash:
        raise ValueError("Replay implementation changed during execution")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    report = {
        "schema_version": 1, "status": "passed", "purpose": "Reusable train-only feature-kernel verification on authenticated existing phage alignments",
        "wrapper_sha256": wrapper_hash, "kernel_sha256": kernel_hash,
        "source_report": {"path": str(SOURCE_REPORT), "sha256": SOURCE_REPORT_SHA256},
        "authenticated_inputs": [{"path": str(p), "sha256": h} for p, h in inputs.items()],
        "training_genomes": len(train_genomes), "training_proteins": len(train_lengths),
        "query_genomes": len(query_genomes), "query_proteins": len(query_lengths),
        "families": len(set(families.values())), "retained_families": len(vocabulary["family_features"]),
        "features": len(features), "training_assignments_match": True, "query_assignments_match": True,
        "pre_query_vocabulary_freeze_sha256": sha(workdir / "pre-query-freeze.json"),
        "vocabulary_sha256": frozen_hash, "vocabulary_unchanged_after_query_projection": True,
        "matrix_sha256": hashlib.sha256(matrix_bytes).hexdigest(), "original_matrix_byte_identical": True,
        "mapping_control_nonzero_equal": bool(any(vectors[train_genomes[0]]) and vectors[train_genomes[0]] == vectors["mapping_control"]),
        "self_alignment_diagnostics": self_diagnostic, "query_alignment_diagnostics": query_diagnostic,
        "raw_workdir": str(workdir),
        "resources": {"wall_seconds": time.monotonic() - start,
                      "peak_replay_process_rss_bytes": peak if sys.platform == "darwin" else peak * 1024},
        "mmseqs_executed": False, "phenotype_tables_loaded": False, "interaction_model_fitted": False,
        "outer_fold_feature_preparation_executed": False, "host_feature_preparation_executed": False,
        "new_method_case_coverage": 0, "production_run_spec_complete": False,
        "interpretation": "Component replay; no new biological prediction or grouped validation. Complete sources, per-fold external-tool execution and full training settings remain required.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    dump(output, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = execute(args.workdir, args.output)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": report["status"], "features": report["features"],
                      "matrix_sha256": report["matrix_sha256"], "new_method_case_coverage": 0}))


if __name__ == "__main__":
    main()
