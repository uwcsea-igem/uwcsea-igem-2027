"""Import a downloaded Figshare host ZIP without extracting its paths.

The reviewed deposit inventory authenticates original filenames, byte sizes and
MD5s. Only verified, unchanged FASTA bytes enter the file-ID genome cache. This
does not establish experimental stock identity or execute training or scoring.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import stat
import tempfile
import zipfile


SPEC = importlib.util.spec_from_file_location(
    "prepare_genophi_for_import", Path(__file__).with_name("prepare-genophi-genomes.py"))
PREPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)


def host_inventory(inventory: list[dict]) -> dict[str, dict]:
    hosts = {}
    paths = set()
    for item in inventory:
        if item["kind"] != "host":
            continue
        name, path = item["source_name"], PurePosixPath(item["path"])
        if (name != PurePosixPath(name).name or "\\" in name
                or path.parts != ("hosts", f"{item['deposit_file_id']}.fasta")
                or not isinstance(item["expected_bytes"], int) or item["expected_bytes"] <= 0
                or not item.get("expected_md5")):
            raise ValueError("Invalid authenticated host inventory record")
        if name in hosts or item["path"] in paths:
            raise ValueError("Duplicate host inventory filename or destination")
        hosts[name] = item
        paths.add(item["path"])
    if not hosts:
        raise ValueError("Host inventory is empty")
    return dict(sorted(hosts.items()))


def member_basename(info: zipfile.ZipInfo) -> str:
    """Validate even ignored members; never pass an archive path to extraction."""
    name = info.orig_filename
    parts = (name[:-1] if info.is_dir() else name).split("/")
    if (name != info.filename or "\\" in name or name.startswith("/")
            or any(part in ("", ".", "..") for part in parts)
            or ":" in parts[0]):
        raise ValueError("Unsafe ZIP member path")
    mode = stat.S_IFMT(info.external_attr >> 16)
    if mode not in (0, stat.S_IFREG, stat.S_IFDIR):
        raise ValueError("ZIP links and special files are not accepted")
    if mode == stat.S_IFDIR and not info.is_dir():
        raise ValueError("ZIP directory attributes disagree with its name")
    return parts[-1]


def destination_for(item: dict, genomes: Path, archive: Path) -> Path:
    destination = genomes / item["path"]
    resolved = destination.resolve()
    if not resolved.is_relative_to(genomes.resolve()) or destination.is_symlink():
        raise ValueError("Genome destination is a symlink or escapes its cache")
    if resolved == archive.resolve() or (destination.exists() and destination.samefile(archive)):
        raise ValueError("Genome destination overlaps the source archive")
    return destination


def import_member(bundle: zipfile.ZipFile, info: zipfile.ZipInfo, item: dict,
                  destination: Path) -> dict:
    # Declared uncompressed size is checked before opening/decompressing a member.
    if info.file_size != item["expected_bytes"]:
        raise ValueError("ZIP member size differs from deposited metadata")
    existing = None
    if destination.exists():
        existing = PREPARE.verify(destination, item)  # A corrupt cache is never replaced.
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle, name = tempfile.mkstemp(prefix=destination.name + ".", suffix=".part", dir=destination.parent)
    temporary = Path(name)
    try:
        with os.fdopen(handle, "wb") as target, bundle.open(info) as source:
            remaining = item["expected_bytes"]
            while remaining:
                block = source.read(min(1024 * 1024, remaining))
                if not block:
                    raise ValueError("ZIP member ended before its expected size")
                target.write(block)
                remaining -= len(block)
            if source.read(1):
                raise ValueError("ZIP member exceeds its expected size")
        stats = PREPARE.verify(temporary, item)
        if existing is not None:
            return {"status": "preserved_existing", **existing}
        try:
            # A hard link publishes atomically and cannot overwrite a racing file.
            os.link(temporary, destination)
        except FileExistsError:
            if destination.is_symlink():
                raise ValueError("Genome symlink appeared during import")
            return {"status": "preserved_existing", **PREPARE.verify(destination, item)}
        return {"status": "imported", **stats}
    finally:
        temporary.unlink(missing_ok=True)


def import_archive(archive: Path, genomes: Path, inventory: list[dict], provenance: dict) -> dict:
    hosts = host_inventory(inventory)
    files, archive_errors, ignored, members = [], [], [], {}
    with archive.open("rb") as archive_stream:
        digest = hashlib.file_digest(archive_stream, "sha256").hexdigest()
        archive_bytes = os.fstat(archive_stream.fileno()).st_size
        archive_stream.seek(0)
        with zipfile.ZipFile(archive_stream) as bundle:
            for info in bundle.infolist():
                try:
                    basename = member_basename(info)
                    if info.is_dir():
                        continue
                    if basename in hosts:
                        members.setdefault(basename, []).append(info)
                    elif Path(basename).suffix.lower() in (".fasta", ".fa", ".fna", ".fas", ".fsa"):
                        raise ValueError("FASTA filename is absent from the authenticated host inventory")
                    else:
                        ignored.append({"member": info.filename, "uncompressed_bytes": info.file_size})
                except ValueError as exc:
                    archive_errors.append({"member": info.orig_filename, "error": str(exc)})
            for name, item in hosts.items():
                entries = members.get(name, [])
                try:
                    destination = destination_for(item, genomes, archive)
                    if len(entries) > 1:
                        raise ValueError("Duplicate expected filename in ZIP: " + ", ".join(x.filename for x in entries))
                    if entries:
                        result = import_member(bundle, entries[0], item, destination)
                    elif destination.exists():
                        result = {"status": "preserved_existing", **PREPARE.verify(destination, item)}
                    else:
                        result = {"status": "missing", "error": "Expected filename absent from ZIP and genome cache"}
                except (OSError, ValueError, RuntimeError, NotImplementedError, zipfile.BadZipFile) as exc:
                    result = {"status": "error", "error": str(exc)}
                files.append({**item, "archive_member": entries[0].filename if len(entries) == 1 else None, **result})
    verified = sum(x["status"] in ("imported", "preserved_existing") for x in files)
    missing_names = sorted(set(hosts) - set(members))
    archive_complete = not missing_names and not archive_errors and all(
        x["status"] in ("imported", "preserved_existing") for x in files)
    return {
        "schema_version": 1, **provenance,
        "purpose": "Authenticate and import a manually downloaded official host archive; no sequence modification",
        "archive": {"path": str(archive.resolve()), "bytes": archive_bytes, "sha256": digest},
        "expected_files": len(hosts), "verified_files": verified,
        "imported_files": sum(x["status"] == "imported" for x in files),
        "preserved_existing_files": sum(x["status"] == "preserved_existing" for x in files),
        "error_files": sum(x["status"] == "error" for x in files),
        "missing_files": [x["path"] for x in files if x["status"] == "missing"],
        "archive_missing_expected_names": missing_names,
        "archive_errors": archive_errors, "ignored_non_fasta_members": ignored,
        "archive_complete": archive_complete, "cache_complete": verified == len(hosts),
        "complete": archive_complete and verified == len(hosts), "files": files,
        "training_or_scoring_executed": False, "aliases_resolved": False,
        "interpretation": "Checksums authenticate deposited bytes, not experimental stock identity. A partial archive remains incomplete even when missing members already exist in the verified cache.",
    }


def protect_paths(output: Path, metadata: Path, genomes: Path, archive: Path, inventory: list[dict]):
    PREPARE.protect_output(output, metadata, genomes, inventory)
    for source in (metadata, genomes, archive):
        if source.resolve().is_relative_to(output.resolve()):
            raise ValueError("Report output overlaps an input path")
    if output.resolve() == archive.resolve() or (output.exists() and output.samefile(archive)):
        raise ValueError("Report output overlaps the source archive")
    for item in inventory:
        if item["kind"] == "host":
            destination_for(item, genomes, archive)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-cache", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--genomes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        inventory, provenance = PREPARE.load_inventory(args.metadata_cache)
        host_inventory(inventory)
        protect_paths(args.output, args.metadata_cache, args.genomes, args.archive, inventory)
        report = import_archive(args.archive, args.genomes, inventory, provenance)
    except (OSError, ValueError, KeyError, zipfile.BadZipFile) as exc:
        parser.error(str(exc))
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
