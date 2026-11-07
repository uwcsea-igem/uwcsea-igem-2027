"""Restore only three target-reference proteomes from the public annotation receipt.

No Picard host collection or phage corpus is required. Existing artifacts are
verified, never overwritten; fresh calls must reproduce the recorded FAA hash.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


def sha(data):
    return hashlib.sha256(data).hexdigest()


def execute(cache, output=None):
    if output is not None and (output.exists() or output.is_symlink()):
        raise ValueError("Restoration check report must be new")
    receipt_path = Path("data/processed/genophi-annotation-check.json")
    receipt_bytes = receipt_path.read_bytes()
    receipt = json.loads(receipt_bytes)
    manifest_path = Path("software/tipcraft/benchmarks/retro-t7-t3/wish-inputs.json")
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    caller = receipt["tools"]["host"]
    executable = Path(caller["executable"])
    spec = importlib.util.spec_from_file_location("reference_host_fasta", "scripts/audit-receptor-host-payloads.py")
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    items = [x for x in receipt["files"] if x["kind"] == "host" and x["role"] == "target_reference"]
    if len(items) != 3 or {x["source_id"] for x in items} != {"BL21", "BW25113", "MG1655"}:
        raise ValueError("Public receipt must describe exactly the three target hosts")
    results, commands = [], []
    for item in items:
        name = item["source_id"]
        source = Path("data/raw/wish/hosts") / (name + ".fna")
        original = source.read_bytes()
        source_spec = manifest["genome_sources"][name]
        if sha(original) != source_spec["file_sha256"]:
            raise ValueError("Original target reference file differs")
        records = reader.fasta(original, set("ACGTRYSWKMBDHVN"))
        if len(records) != 1 or sha(records[0][1].upper().encode()) != item["source_sequence_sha256"]:
            raise ValueError("Target sequence differs from public annotation receipt")
        normalized = f">c000000\n{records[0][1].upper()}\n".encode()
        if sha(normalized) != item["normalized_input_sha256"]:
            raise ValueError("Reconstructed normalized input differs")
        folder = cache / item["id"]
        fasta, protein, gff, log = [folder / x for x in ("input.fna", "proteins.faa", "genes.gff", "run.log")]
        for p in (folder, fasta, protein, gff, log):
            if p.is_symlink():
                raise ValueError("Symlink restoration destinations are not accepted")
        if protein.exists():
            if sha(protein.read_bytes()) != item["faa_sha256"]:
                raise ValueError("Cached target protein differs; no automatic repair")
            origin = "existing_verified_proteome"
        else:
            if sha(executable.read_bytes()) != caller["executable_sha256"]:
                raise ValueError("Original pinned Prodigal binary differs")
            if folder.exists() and any(folder.iterdir()):
                raise ValueError("Incomplete target cache must remain intact; use a fresh cache")
            folder.mkdir(parents=True, exist_ok=True)
            fasta.write_bytes(normalized)
            argv = [str(executable), *item["command"][1:8], "-i", str(fasta), "-a", str(protein), "-o", str(gff)]
            commands.append(argv)
            with log.open("xb") as stream:
                subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=120)
            if sha(protein.read_bytes()) != item["faa_sha256"]:
                raise ValueError("Fresh target call differs from recorded expected FAA hash")
            origin = "fresh_hash_matching_call"
        results.append({"source_id": name, "protein_path": str(protein), "sha256": item["faa_sha256"], "origin": origin})
    cache.mkdir(parents=True, exist_ok=True)
    receipt_copy = cache / "report.json"
    if receipt_copy.exists():
        if receipt_copy.read_bytes() != receipt_bytes:
            raise ValueError("Existing annotation receipt differs; no overwrite")
    else:
        with receipt_copy.open("xb") as stream:
            stream.write(receipt_bytes)
    if receipt_path.read_bytes() != receipt_bytes or manifest_path.read_bytes() != manifest_bytes:
        raise ValueError("Public receipt or target manifest changed during restoration")
    result = {"schema_version": 1, "script_sha256": sha(Path(__file__).read_bytes()),
              "target_host_proteomes": results, "public_receipt_sha256": sha(receipt_bytes),
              "source_manifest_sha256": sha(manifest_bytes), "fresh_call_commands": commands,
              "copied_report_is_historical_provenance_not_a_new_full_collection_run": True,
              "picard_host_inputs_needed": False, "phenotype_labels_loaded": False}
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("data/raw/genophi-annotation/run-1"))
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    execute(args.cache, args.report)
