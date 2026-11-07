"""Authenticate source-length and terminal-overlap diagnostics without trimming."""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path

SOURCE_COMMIT = "99656d8c020fcd1dbbc30a02ee10274a07655562"
SOURCE_TREE = "550f10819c6d565214bc26d1e3fea59af17ec398"
T4LD_SHA256 = "3cbcc8a322d3d80a17df4c24a90fd948936005b7bb9ff7effa640fde2ecceab5"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def authenticate(data, entry):
    blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if entry["type"] != "blob" or blob != entry["sha"] or len(data) != entry["size"]:
        raise ValueError("Source bytes differ from pinned Git inventory")
    return blob


def inspect_sequence(data: bytes, limit: int = 120) -> dict:
    records = []
    for line in data.decode("ascii").splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if not line[1:].strip():
                raise ValueError("Empty FASTA header")
            records.append("")
        else:
            if not records:
                raise ValueError("Sequence before FASTA header")
            sequence = "".join(line.split()).upper()
            if not set(sequence) <= set("ACGTRYSWKMBDHVN"):
                raise ValueError("Non-IUPAC DNA source")
            records[-1] += sequence
    empty = [i for i, sequence in enumerate(records) if not sequence]
    if empty != ([0] if sha(data) == T4LD_SHA256 else []):
        raise ValueError("Unreviewed empty-record pattern")
    retained = [s for s in records if s]
    if len(retained) != 1:
        raise ValueError("Require one nonempty source record; no contig concatenation")
    sequence = retained[0]
    return {"retained_bases": len(sequence), "source_records": len(records),
            "skipped_empty_record_indices_for_inspection": empty,
            "sequence_sha256": sha(sequence.encode()),
            "exact_prefix_suffix_match_lengths": [n for n in range(1, min(limit, len(sequence)) + 1)
                                                   if sequence[:n] == sequence[-n:]],
            "leading_n_bases": len(sequence) - len(sequence.lstrip("N")),
            "trailing_n_bases": len(sequence) - len(sequence.rstrip("N")),
            "total_n_bases": sequence.count("N"),
            "total_non_acgt_bases": sum(base not in "ACGT" for base in sequence)}


def run(policy_path: Path, tree_path: Path, metadata_path: Path, genomes: Path, readme_path: Path) -> dict:
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if policy["source_commit"] != SOURCE_COMMIT or policy["source_tree_sha1"] != SOURCE_TREE:
        raise ValueError("Policy source identity changed")
    if policy["maximum_terminal_overlap_bases"] != 120 or policy["expected_source_files"] != 97 or policy["expected_metadata_rows"] != 96:
        raise ValueError("Diagnostic scope differs from frozen review")
    spec = importlib.util.spec_from_file_location("source_audit", Path(__file__).with_name("audit-genophi-sources.py"))
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    entries = audit.validate_tree(json.loads(tree_path.read_text()), SOURCE_TREE)
    metadata = metadata_path.read_bytes()
    metadata_blob = authenticate(metadata, entries[policy["metadata_path"]])
    if metadata_blob != policy["metadata_git_blob_sha1"]:
        raise ValueError("Metadata identity differs from frozen policy")
    rows = list(csv.DictReader(io.StringIO(metadata.decode()), delimiter=";"))
    if len(rows) != 96 or len({row["phage"] for row in rows}) != 96:
        raise ValueError("Unexpected or duplicate metadata roster")
    sizes = {row["phage"]: int(row["Genome_size"]) for row in rows}
    if any(value <= 0 for value in sizes.values()):
        raise ValueError("Invalid metadata genome length")
    readme = readme_path.read_bytes()
    authenticate(readme, entries["README.md"])
    files = []
    prefix = "data/genomics/phages/FNA/"
    for path, entry in sorted(entries.items()):
        if not path.startswith(prefix) or not path.endswith(".fna"):
            continue
        name = Path(path).stem
        data = (genomes / "phages" / Path(path).name).read_bytes()
        blob = authenticate(data, entry)
        record = inspect_sequence(data, 120)
        expected = sizes.get(name)
        files.append({"source_id": name, "source_path": path, "git_blob_sha1": blob,
                      "source_file_sha256": sha(data), "metadata_genome_size": expected,
                      "length_minus_metadata": record["retained_bases"] - expected if expected is not None else None,
                      **record})
    if len(files) != 97 or len([r for r in files if r["metadata_genome_size"] is not None]) != 96:
        raise ValueError("Unexpected source/metadata correspondence counts")
    generation_files = [p for p in entries if p.startswith("data/genomics/phages/")
                        and Path(p).suffix.lower() in {".py", ".sh", ".r", ".ipynb", ".rmd"}]
    return {"schema_version": 1, "purpose": policy["purpose"], "policy_sha256": sha(policy_bytes),
            "helper_sha256": sha(Path(__file__).read_bytes()), "source_commit": SOURCE_COMMIT,
            "source_tree_sha1": SOURCE_TREE, "metadata_git_blob_sha1": metadata_blob,
            "metadata_sha256": sha(metadata), "readme_git_blob_sha1": entries["README.md"]["sha"],
            "metadata_url": f"https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/{SOURCE_COMMIT}/{policy['metadata_path']}",
            "source_files": len(files), "metadata_matched_files": 96,
            "files_without_metadata_row": [r["source_id"] for r in files if r["metadata_genome_size"] is None],
            "length_minus_metadata_counts": dict(sorted(Counter(str(r["length_minus_metadata"]) for r in files if r["metadata_genome_size"] is not None).items())),
            "source_files_with_identical_terminal_60_bases": sum(60 in r["exact_prefix_suffix_match_lengths"] for r in files),
            "sources_with_terminal_ns": [r["source_id"] for r in files if r["leading_n_bases"] or r["trailing_n_bases"]],
            "files": files, "source_generation_script_paths_under_phage_directory": generation_files,
            "generation_rationale_status": "Not established: pinned README describes FASTA genome inputs without a terminal-extension explanation; no generation-code file under the phage data directory. This bounded inspection does not establish that no rationale exists elsewhere.",
            "sequence_changes_performed": False, "phenotype_table_loaded": False,
            "training_or_prediction_executed": False,
            "interpretation": "Preserved source representations have terminal overlap and length differences from metadata; do not infer natural termini, circularity, padding provenance, stock identity or trim sequences. Normalization and relatedness results remain source-representation-specific and provisional."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("policy", "tree", "metadata", "genomes", "readme", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or any(args.output.resolve() == p.resolve() for p in (args.policy, args.tree, args.metadata, args.readme)):
        parser.error("Use a new report path separate from inputs")
    result = run(args.policy, args.tree, args.metadata, args.genomes, args.readme)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
