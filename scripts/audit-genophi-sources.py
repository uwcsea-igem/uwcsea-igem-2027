"""Inventory pinned public inputs, without training or predicting interactions.

The upstream interaction tables contain labels. This source audit reads them;
it must never be described as a label-isolated scoring run. Name matches are
discovery aids, not sequence-identity or stock-equivalence tests.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile


COMMIT = "f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f"
REPO = "Noonanav/GenoPHI"
ZENODO = "https://zenodo.org/api/records/21696332/files/GenoPHI_processed_data.zip/content"
ZIP_MD5 = "157de25b9c075a93b94781db04182d7b"
TABLES = (
    "data/interaction_matrices/ecoli_interaction_matrix.csv",
    "data/interaction_matrices/ecoli_interaction_matrix_subset.csv",
    "data/test_data/ecoli_test_interaction_matrix.csv",
)


def validate_tree(tree: dict, expected_root: str) -> dict:
    """Authenticate every inventory path against a Git root tree."""
    if tree.get("truncated") or not isinstance(tree.get("tree"), list):
        raise ValueError("Incomplete repository inventory")
    entries = {entry["path"]: entry for entry in tree["tree"]}
    if len(entries) != len(tree["tree"]):
        raise ValueError("Duplicate repository inventory path")
    for path in entries:
        if not path or any(part in {"", ".", ".."} for part in path.split("/")):
            raise ValueError(f"Invalid repository inventory path: {path}")
        if "/" in path:
            parent = path.rsplit("/", 1)[0]
            if parent not in entries or entries[parent]["type"] != "tree":
                raise ValueError(f"Missing directory parent: {path}")
    for directory in sorted([p for p, e in entries.items() if e["type"] == "tree"] + [""],
                            key=lambda p: p.count("/") + bool(p), reverse=True):
        children = [(p.rsplit("/", 1)[-1], e) for p, e in entries.items()
                    if (p.rsplit("/", 1)[0] if "/" in p else "") == directory]
        children.sort(key=lambda item: item[0] + ("/" if item[1]["type"] == "tree" else ""))
        payload = b"".join((e["mode"].lstrip("0") + " " + name).encode() + b"\0" + bytes.fromhex(e["sha"])
                           for name, e in children)
        digest = hashlib.sha1(b"tree " + str(len(payload)).encode() + b"\0" + payload).hexdigest()
        if digest != (entries[directory]["sha"] if directory else expected_root):
            raise ValueError(f"Repository tree hash mismatch: {directory or '/'}")
    return entries


def audit(cache: Path, offline: bool = False) -> dict:
    artifacts = []

    def read(name: str, url: str) -> bytes:
        path = cache / name
        if not path.exists():
            if offline:
                raise ValueError(f"Missing cached artifact: {path}")
            request = urllib.request.Request(url, headers={"User-Agent": "TipCraft-source-audit"})
            with urllib.request.urlopen(request, timeout=45) as response:
                data = response.read(10_000_001)
            if len(data) > 10_000_000:
                raise ValueError("Source exceeds audit size limit")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        data = path.read_bytes()
        artifacts.append({"cache_file": name, "url": url, "bytes": len(data),
                          "sha256": hashlib.sha256(data).hexdigest()})
        return data

    tree = json.loads(read("tree.json", f"https://api.github.com/repos/{REPO}/git/trees/{COMMIT}?recursive=1"))
    # Reconstruct Git's tree hashes bottom-up, binding all paths to the pinned
    # commit's independently recorded root tree, including offline cache replay.
    expected_root = "e42e64caaaea9efb8aca2af00d2dd84541bfe946"
    entries = validate_tree(tree, expected_root)

    summaries = {}
    for path in TABLES:
        data = read(path, f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{path}")
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if blob != entries[path]["sha"]:
            raise ValueError(f"Pinned Git blob mismatch: {path}")
        rows = list(csv.DictReader(io.StringIO(data.decode())))
        if not rows or set(rows[0]) != {"strain", "phage", "interaction"}:
            raise ValueError(f"Unexpected table schema: {path}")
        pairs = {(r["strain"], r["phage"]) for r in rows}
        if len(pairs) != len(rows) or any(r["interaction"] not in {"0", "1"} for r in rows):
            raise ValueError(f"Invalid interaction table: {path}")
        strains = sorted({r["strain"] for r in rows})
        phages = sorted({r["phage"] for r in rows})
        summaries[path] = {
            "rows": len(rows), "strains": len(strains), "phages": len(phages),
            "target_name_candidates_only": {name: [x for x in strains + phages if name.casefold() in x.casefold()]
                                             for name in ("T7", "T3", "BL21", "BW25113", "MG1655")},
            "sequence_overlap_audit": "not_performed; names cannot establish absence of homologs or stock identity",
        }
    archive = read("GenoPHI_processed_data.zip", ZENODO)
    if hashlib.md5(archive).hexdigest() != ZIP_MD5:
        raise ValueError("Zenodo archive checksum differs from reviewed record")
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        archive_paths = sorted(x.filename for x in z.infolist() if not x.is_dir())
    model_suffixes = (".cbm", ".pkl", ".pickle", ".joblib", ".pt", ".pth", ".onnx")
    model_paths = [p for p in entries if p.lower().endswith(model_suffixes)]
    archive_models = [p for p in archive_paths if p.lower().endswith(model_suffixes)]
    return {
        "schema_version": 1, "reviewed_date": "2026-11-07", "upstream_commit": COMMIT,
        "repository_tree_sha1": expected_root, "source_artifacts": artifacts,
        "interaction_tables": summaries,
        "bundled_demo_proteomes": {
            kind: sum(p.startswith(f"data/test_data/{kind}_AAs/") and p.endswith(".faa") for p in entries)
            for kind in ("strain", "phage")},
        "repository_model_filename_matches": model_paths,
        "zenodo_model_filename_matches": archive_models,
        "zenodo_archive_paths": archive_paths,
        "scope": "Pinned repository tree and one processed-data archive; external deposits and unusual model formats not excluded",
        "training_or_scoring_executed": False,
        "new_method_covered_cases": 0,
        "disposition": "HOLD: no audited model/feature-reference bundle and no completed target/homology exclusion audit",
        "interpretation": "A reproducible input inventory, not a failed GenoPHI prediction or proof that required public data do not exist",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    if args.output.resolve().is_relative_to(args.cache.resolve()):
        parser.error("Keep the output outside the source cache")
    report = audit(args.cache, args.offline)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
