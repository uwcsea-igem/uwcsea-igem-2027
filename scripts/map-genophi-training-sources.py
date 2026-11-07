"""Map published training IDs to deposited file metadata before sequence audit.

This reads an external phenotype table for its identifiers. It does not fit,
score, infer stock identity, or accept punctuation-normalized names as aliases.
"""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import urllib.request

COLI_COMMIT = "99656d8c020fcd1dbbc30a02ee10274a07655562"
COLI_TREE = "550f10819c6d565214bc26d1e3fea59af17ec398"
FIGSHARE = "https://api.figshare.com/v2/articles/25941691/versions/1"
HOST_INVENTORY_SHA256 = "fc5019c8de4a88d2432b298624374f1b7006ae8102fddbbc997e4dea9d20c24e"


def validate_deposit(deposit: dict) -> list[dict]:
    if deposit["id"] != 25941691 or deposit["version"] != 1:
        raise ValueError("Unexpected Figshare version")
    files = sorted([{k: f[k] for k in ("id", "name", "size", "download_url", "computed_md5")}
                    for f in deposit["files"]], key=lambda f: f["id"])
    digest = hashlib.sha256(json.dumps(files, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if digest != HOST_INVENTORY_SHA256:
        raise ValueError("Figshare file inventory differs from reviewed metadata")
    return files


def match_ids(ids: list[str], records: list[dict]) -> dict:
    """Exact ID/file matches and unaccepted alias candidates; never pick a tie."""
    if len({r["id"] for r in records}) != len(records):
        raise ValueError("Duplicate deposited file ID")
    normalized = lambda value: re.sub(r"[^a-z0-9]", "", value.casefold())
    result = {}
    for identifier in sorted(set(ids)):
        exact = [r for r in records if r["id"] == identifier]
        candidates = [r["id"] for r in records if normalized(r["id"]) == normalized(identifier)]
        result[identifier] = {
            "status": "exact_identifier_match" if exact else "unresolved",
            "file": exact[0] if exact else None,
            "unaccepted_alias_candidates": [] if exact else sorted(candidates),
        }
    return result


def run(cache: Path, matrix: Path, offline: bool = False) -> dict:
    spec = importlib.util.spec_from_file_location("source_audit", Path(__file__).with_name("audit-genophi-sources.py"))
    source_audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source_audit)
    artifacts = []

    def get(name: str, url: str) -> bytes:
        path = cache / name
        if not path.exists():
            if offline:
                raise ValueError(f"Missing cached artifact: {path}")
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "TipCraft-source-mapping"}), timeout=45) as r:
                data = r.read(10_000_001)
            if len(data) > 10_000_000:
                raise ValueError("Input metadata exceeds size limit")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        data = path.read_bytes()
        artifacts.append({"url": url, "cache_file": name, "sha256": hashlib.sha256(data).hexdigest()})
        return data

    tree = json.loads(get("coli-tree.json", f"https://api.github.com/repos/mdmparis/coli_phage_interactions_2023/git/trees/{COLI_COMMIT}?recursive=1"))
    entries = source_audit.validate_tree(tree, COLI_TREE)
    deposit = json.loads(get("figshare-v1.json", FIGSHARE))
    files = validate_deposit(deposit)
    hosts = [{"id": Path(f["name"]).stem, "file_id": f["id"], "url": f["download_url"],
              "bytes": f["size"], "deposited_md5": f["computed_md5"]}
             for f in files if f["name"].endswith(".fasta")]
    phages = [{"id": Path(p).stem, "url": f"https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/{COLI_COMMIT}/{p}",
               "bytes": e["size"], "git_blob_sha1": e["sha"]}
              for p, e in entries.items() if p.startswith("data/genomics/phages/FNA/") and p.endswith(".fna")]
    data = matrix.read_bytes()
    # Bind the phenotype identifiers to the already pinned GenoPHI inventory.
    expected = "c4795d786a1e13ac998e263dbf4aaf95fc0f94c52b842d02b6f7fc6e99ccfef2"
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError("Training table differs from reviewed GenoPHI source")
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    host_map = match_ids([r["strain"] for r in rows], hosts)
    phage_map = match_ids([r["phage"] for r in rows], phages)
    summary = lambda mapping: {"total": len(mapping), "exact_identifier_matches": sum(v["file"] is not None for v in mapping.values()),
                               "unresolved": [k for k, v in mapping.items() if v["file"] is None]}
    return {
        "schema_version": 1, "reviewed_date": "2026-11-07", "source_commit": COLI_COMMIT,
        "source_tree_sha1": COLI_TREE, "genophi_table_sha256": expected,
        "source_artifacts": artifacts, "host_deposit_doi": "10.6084/m9.figshare.25941691.v1",
        "reviewed_host_file_inventory_sha256": HOST_INVENTORY_SHA256,
        "available_file_metadata": {"hosts": len(hosts), "phages": len(phages)},
        "host_summary": summary(host_map), "phage_summary": summary(phage_map),
        "hosts": host_map, "phages": phage_map,
        "phenotype_values_used_for_prediction": False, "training_or_scoring_executed": False,
        "sequence_overlap_audit": "pending; exact identifiers are not sequence/stock authentication",
        "disposition": "Metadata mapping advances input preparation; HOLD target inference until aliases, sequence exclusions and training specification are resolved",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    if args.output.resolve().is_relative_to(args.cache.resolve()) or args.output.resolve() == args.matrix.resolve():
        parser.error("Output must be separate from inputs")
    report = run(args.cache, args.matrix, args.offline)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
