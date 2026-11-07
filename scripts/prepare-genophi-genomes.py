"""Verify deposited genomes for an input audit; never train, score or resolve aliases.

Metadata must already be cached by the source-mapping audit. Genome files are
published atomically only after their expected checksum and FASTA validation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile


def _module(filename: str):
    spec = importlib.util.spec_from_file_location(filename, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_inventory(metadata: Path) -> tuple[list[dict], dict]:
    mapping = _module("map-genophi-training-sources.py")
    audit = _module("audit-genophi-sources.py")
    deposit_bytes = (metadata / "figshare-v1.json").read_bytes()
    tree_bytes = (metadata / "coli-tree.json").read_bytes()
    files = mapping.validate_deposit(json.loads(deposit_bytes))
    entries = audit.validate_tree(json.loads(tree_bytes), mapping.COLI_TREE)
    inventory = []
    for item in files:
        if item["name"].endswith(".fasta"):
            inventory.append({
                "kind": "host", "source_id": Path(item["name"]).stem,
                "source_name": item["name"], "deposit_file_id": item["id"],
                "path": f"hosts/{item['id']}.fasta", "url": item["download_url"],
                "expected_bytes": item["size"], "expected_md5": item["computed_md5"],
            })
    for path, entry in entries.items():
        if path.startswith("data/genomics/phages/FNA/") and path.endswith(".fna"):
            if entry["type"] != "blob":
                raise ValueError(f"Genome inventory entry is not a file: {path}")
            inventory.append({
                "kind": "phage", "source_id": Path(path).stem,
                "source_name": path, "path": f"phages/{Path(path).name}",
                "url": f"https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/{mapping.COLI_COMMIT}/{path}",
                "expected_bytes": entry["size"], "expected_git_blob_sha1": entry["sha"],
            })
    inventory.sort(key=lambda item: item["path"])
    if (sum(x["kind"] == "host" for x in inventory), sum(x["kind"] == "phage" for x in inventory)) != (403, 97):
        raise ValueError("Reviewed inventory must contain 403 host and 97 phage files")
    if len({x["path"] for x in inventory}) != len(inventory):
        raise ValueError("Duplicate genome destination")
    provenance = {
        "source_commit": mapping.COLI_COMMIT, "source_tree_sha1": mapping.COLI_TREE,
        "host_deposit_doi": "10.6084/m9.figshare.25941691.v1",
        "reviewed_host_file_inventory_sha256": mapping.HOST_INVENTORY_SHA256,
        "metadata_files": [
            {"path": "coli-tree.json", "sha256": hashlib.sha256(tree_bytes).hexdigest()},
            {"path": "figshare-v1.json", "sha256": hashlib.sha256(deposit_bytes).hexdigest()},
        ],
    }
    return inventory, provenance


def fasta_stats(data: bytes) -> dict:
    """Hash uppercase sequences in source contig order, with length framing.

    Headers and wrapping are excluded from the sequence hash. Boundaries are
    retained: [AC, GT] and [A, CGT] have different sequence hashes.
    """
    try:
        lines = data.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise ValueError("FASTA must be ASCII") from exc
    allowed = set("ACGTRYSWKMBDHVN")
    headers, duplicates, sequences = set(), set(), []
    current = None
    for line in lines:
        if not line.strip():
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            if not header:
                raise ValueError("Empty FASTA header")
            identifier = header.split()[0]
            if identifier in headers:
                duplicates.add(identifier)
            headers.add(identifier)
            if current is not None:
                if not current:
                    raise ValueError("Empty FASTA contig")
                sequences.append("".join(current))
            current = []
        else:
            if current is None:
                raise ValueError("FASTA sequence occurs before a header")
            sequence = "".join(line.split()).upper()
            if not set(sequence) <= allowed:
                raise ValueError("FASTA contains non-IUPAC DNA characters")
            current.append(sequence)
    if current is None or not current:
        raise ValueError("Empty FASTA sequence")
    sequences.append("".join(current))
    digest = hashlib.sha256(b"TipCraft-contig-sequences-v1\0")
    for sequence in sequences:
        digest.update(len(sequence).to_bytes(8, "big"))
        digest.update(sequence.encode("ascii"))
    return {
        "contigs": len(sequences), "contig_lengths": [len(x) for x in sequences],
        "total_bases": sum(map(len, sequences)),
        "ambiguous_bases": sum(base not in "ACGT" for seq in sequences for base in seq),
        "n_bases": sum(seq.count("N") for seq in sequences),
        "duplicate_header_identifiers": sorted(duplicates),
        "requires_header_normalization": bool(duplicates),
        "sequence_sha256": digest.hexdigest(),
        "sequence_hash_scheme": "uppercase-IUPAC; source-contig-order; 8-byte-big-endian-length-framing; headers-excluded; v1",
    }


def verify(path: Path, item: dict) -> dict:
    if path.stat().st_size != item["expected_bytes"]:
        raise ValueError("File size differs from deposited metadata")
    data = path.read_bytes()
    if len(data) != item["expected_bytes"]:
        raise ValueError("File size differs from deposited metadata")
    if "expected_md5" in item and hashlib.md5(data).hexdigest() != item["expected_md5"]:
        raise ValueError("File MD5 differs from deposited metadata")
    if "expected_git_blob_sha1" in item:
        digest = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if digest != item["expected_git_blob_sha1"]:
            raise ValueError("Git blob SHA1 differs from pinned repository")
    return {"file_sha256": hashlib.sha256(data).hexdigest(), **fasta_stats(data)}


def prepare_one(item: dict, genomes: Path, offline: bool = False) -> dict:
    destination = genomes / item["path"]
    if destination.exists():
        # A changed cache is evidence to inspect, never a reason to overwrite it.
        return {**item, "status": "verified", **verify(destination, item)}
    if destination.is_symlink():
        raise ValueError("Dangling genome symlink")
    if offline:
        raise ValueError("Genome missing from offline cache")
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_name = tempfile.mkstemp(prefix=destination.name + ".", suffix=".part", dir=destination.parent)
    os.close(handle)
    temporary = Path(temporary_name)
    try:
        result = subprocess.run([
            "curl", "--location", "--fail", "--silent", "--show-error",
            "--connect-timeout", "15", "--max-time", "180", "--retry", "2",
            "--retry-delay", "1", "--max-filesize", str(item["expected_bytes"]),
            "--output", str(temporary), item["url"],
        ], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, timeout=600)
        if result.returncode:
            raise ValueError(f"curl failed with exit status {result.returncode}")
        stats = verify(temporary, item)
        if destination.exists() or destination.is_symlink():
            raise ValueError("Genome destination appeared during download; rerun to verify")
        temporary.replace(destination)
        return {**item, "status": "verified", **stats}
    finally:
        temporary.unlink(missing_ok=True)


def prepare_inventory(inventory: list[dict], provenance: dict, genomes: Path,
                      offline: bool = False, workers: int = 4) -> dict:
    if not 1 <= workers <= 4:
        raise ValueError("Workers must be between 1 and 4")
    def task(item):
        try:
            return prepare_one(item, genomes, offline)
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            # Keep errors independent of temporary absolute paths and timings.
            message = str(exc)
            if isinstance(exc, OSError):
                message = f"{type(exc).__name__}: {exc.strerror}"
            elif isinstance(exc, subprocess.TimeoutExpired):
                message = "curl subprocess exceeded its timeout"
            return {**item, "status": "error", "error": message}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        files = list(pool.map(task, inventory))
    files.sort(key=lambda item: item["path"])
    verified = sum(item["status"] == "verified" for item in files)
    return {
        "schema_version": 1, **provenance,
        "purpose": "Public genome input audit; source identifiers retained without accepting training aliases",
        "expected_files": len(files), "verified_files": verified,
        "error_files": len(files) - verified, "complete": verified == len(files),
        "files_requiring_header_normalization": [item["path"] for item in files if item.get("requires_header_normalization")],
        "files": files,
        "training_or_scoring_executed": False, "aliases_resolved": False,
        "interpretation": "Checksums authenticate deposited file contents. Repeated source headers are reported without merging contigs or changing sequences; normalize them before tools requiring unique IDs. Stock identity and sequence-overlap exclusions remain separate audits",
    }


def protect_output(output: Path, metadata: Path, genomes: Path, inventory: list[dict]):
    for root in (metadata, genomes):
        if output.resolve().is_relative_to(root.resolve()):
            raise ValueError("Output must be outside metadata and genome input directories")
    inputs = [metadata / name for name in ("figshare-v1.json", "coli-tree.json")]
    inputs.extend(genomes / item["path"] for item in inventory)
    if output.exists() and any(path.exists() and output.samefile(path) for path in inputs):
        raise ValueError("Output must not be a hard link to an input")
    for item in inventory:
        destination = (genomes / item["path"]).resolve()
        if not destination.is_relative_to(genomes.resolve()) or destination.is_relative_to(metadata.resolve()):
            raise ValueError("Genome destination resolves outside its directory or into metadata")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-cache", type=Path, required=True)
    parser.add_argument("--genomes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 5))
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    try:
        inventory, provenance = load_inventory(args.metadata_cache)
        protect_output(args.output, args.metadata_cache, args.genomes, inventory)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    report = prepare_inventory(inventory, provenance, args.genomes, args.offline, args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    handle, name = tempfile.mkstemp(prefix=args.output.name + ".", suffix=".part", dir=args.output.parent)
    try:
        with os.fdopen(handle, "w") as stream:
            stream.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
        Path(name).replace(args.output)
    finally:
        Path(name).unlink(missing_ok=True)
    return 0 if report["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
