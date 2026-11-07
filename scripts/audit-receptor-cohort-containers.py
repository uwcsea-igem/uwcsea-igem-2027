"""Authenticate a frozen phage-container cohort and inventory deposited annotations."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tempfile
import warnings


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                       allow_nan=False) + "\n").encode()


def verify(data, item):
    blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if len(data) != item["declared_size_bytes"] or blob != item["git_blob_sha1"]:
        raise ValueError("Container size or Git blob SHA1 differs from frozen source")
    return {"size_bytes": len(data), "git_blob_sha1": blob, "sha256": sha(data)}


def publish(path, data):
    """Publish exact verified bytes atomically; never repair a changed cache."""
    if path.is_symlink():
        raise ValueError("Symlink cache/output is not accepted")
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError("Existing artifact differs from expected bytes")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".part", dir=path.parent)
    temp = Path(name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
        # The run owns this directory; refuse a late conflicting writer.
        if path.exists() or path.is_symlink():
            raise ValueError("Artifact appeared while publishing")
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def inventory(roster):
    items = []
    for row in roster["records"]:
        matches = row["exact_primary_id_file_matches"]
        if row["mapping_status"] != "unique_candidate" or len(matches) != 1:
            raise ValueError("Every frozen cohort row must have one exact container")
        items.append({"primary_id": row["primary_id"], "cohort_member": True,
                      "published_model_role": row["published_model_role"], **matches[0]})
    extras = [x for x in roster["unused_catalog_superset"] if x["filename"] == "T3.gbk"]
    if len(items) != 260 or len(extras) != 1:
        raise ValueError("Expected all 260 cohort containers and separate T3")
    if Counter(x["published_model_role"] for x in items) != {"training": 206, "testing": 49, "mutant": 5}:
        raise ValueError("Published cohort roles differ from frozen 206/49/5")
    items.append({"primary_id": "T3", "cohort_member": False,
                  "published_model_role": "catalogue_extra_target_candidate", **extras[0]})
    items.sort(key=lambda x: x["primary_id"])
    if len({x["primary_id"] for x in items}) != 261:
        raise ValueError("Duplicate primary ID")
    for index, item in enumerate(items):
        item["container_id"] = f"c{index:03d}"
        path = PurePosixPath(item["repository_path"])
        if path.parts != ("Ecoli_phages", "data", "phage_genomes", item["filename"]):
            raise ValueError("Container path escapes pinned phage directory")
        if path.suffix != ".gbk" or path.stem != item["primary_id"]:
            raise ValueError("Exact primary filename match required")
        if item["format"] != "GenBank" or item["declared_size_bytes"] <= 0:
            raise ValueError("Invalid frozen container metadata")
        # Remove obsolete metadata-only observation flags, preserving source identity.
        item.pop("sequence_payload_read", None)
        item.pop("payload_bytes_authenticated", None)
    return items


def retrieve(item, policy, workdir, offline=False):
    path = workdir / "containers" / item["filename"]
    source = Path(policy["existing_target_cache"]) / item["filename"]
    url = ("https://raw.githubusercontent.com/mjohnson11/PhageDataSheets/"
           + policy["source_commit"] + "/" + item["repository_path"])
    base = {"source": item, "container_path": str(path), "url": url}
    temp = None
    try:
        if path.is_symlink() or source.is_symlink():
            raise ValueError("Symlink input is not accepted")
        if path.exists():
            data = path.read_bytes()
            origin = "existing_cohort_cache"
        elif source.exists() and item["primary_id"] in {"T7", "T3"}:
            data = source.read_bytes()
            origin = "existing_authenticated_target_candidate_reverified"
        elif offline:
            raise ValueError("Container absent in offline cache")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix=item["filename"] + ".", suffix=".part", dir=path.parent)
            os.close(fd)
            temp = Path(name)
            command = ["curl", "--fail", "--location", "--silent", "--show-error",
                       "--max-time", "45", "--retry", "0", "--max-filesize",
                       str(item["declared_size_bytes"]), "--output", str(temp), url]
            result = subprocess.run(command, capture_output=True, text=True, timeout=50)
            if result.returncode:
                raise ValueError(f"curl exit {result.returncode}: {result.stderr.strip()}")
            data = temp.read_bytes()
            origin = "downloaded"
        checks = verify(data, item)
        publish(path, data)
        return {**base, "payload_status": "authenticated", "retrieval": origin, **checks}
    except Exception as exc:
        return {**base, "payload_status": "error", "error": f"{type(exc).__name__}: {exc}"}
    finally:
        if temp is not None:
            temp.unlink(missing_ok=True)


def fasta(entries):
    return "".join(f">{name}\n{seq}\n" for name, seq in entries).encode("ascii")


def framed_sequences(sequences):
    digest = hashlib.sha256(b"TipCraft-contig-sequences-v1\0")
    for seq in sequences:
        data = seq.encode("ascii")
        digest.update(len(data).to_bytes(8, "big"))
        digest.update(data)
    return digest.hexdigest()


def parse_container(data, container_id):
    """Only inspect sequence, feature types, and whitelisted CDS provenance fields."""
    from Bio import SeqIO
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        parsed = list(SeqIO.parse(io.StringIO(data.decode("utf-8")), "genbank"))
    if not parsed:
        raise ValueError("GenBank has no records")
    records, nucleotides, proteins, protein_errors = [], [], [], []
    for index, record in enumerate(parsed):
        nid = f"{container_id}r{index:03d}"
        seq = str(record.seq).upper()
        if not seq or not set(seq) <= set("ACGTRYSWKMBDHVN"):
            raise ValueError("Missing or non-IUPAC DNA sequence")
        nucleotides.append((nid, seq))
        cds = [feature for feature in record.features if feature.type == "CDS"]
        annotations = []
        for ci, feature in enumerate(cds):
            for ti, value in enumerate(feature.qualifiers.get("translation", [])):
                pid = f"{nid}p{ci:05d}v{ti:02d}"
                aa = value.upper()
                annotation = {"protein_id": pid, "cds_index": ci, "translation_index": ti,
                              "source_protein_ids": feature.qualifiers.get("protein_id", []),
                              "source_locus_tags": feature.qualifiers.get("locus_tag", []),
                              "source_location": str(feature.location),
                              "source_transl_table": feature.qualifiers.get("transl_table", []),
                              "source_codon_start": feature.qualifiers.get("codon_start", []),
                              "has_translation_exception": "transl_except" in feature.qualifiers,
                              "source_pseudo_annotation": "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers,
                              "translation_agreement_with_DNA": "not_verified"}
                if not aa or not set(aa) <= set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"):
                    annotation["error"] = "Empty or non-alphabetic deposited translation"
                    protein_errors.append({"protein_id": pid, "error": annotation["error"]})
                else:
                    proteins.append((pid, aa))
                    annotation.update({"length": len(aa), "sequence_sha256": sha(aa.encode("ascii")),
                                       "letters_outside_standard_20": dict(sorted(Counter(x for x in aa if x not in "ACDEFGHIKLMNPQRSTVWY").items())),
                                       "internal_stop_count": aa[:-1].count("*"),
                                       "terminal_stop_present": aa.endswith("*")})
                annotations.append(annotation)
        records.append({"record_id": record.id, "neutral_nucleotide_id": nid, "length": len(seq),
                        "sequence_sha256": sha(seq.encode("ascii")),
                        "nucleotide_alphabet_counts": dict(sorted(Counter(seq).items())),
                        "ambiguous_bases": sum(x not in "ACGT" for x in seq),
                        "feature_counts": dict(sorted(Counter(f.type for f in record.features).items())),
                        "cds_count": len(cds),
                        "cds_without_translation": sum(not f.qualifiers.get("translation") for f in cds),
                        "cds_with_multiple_translation_values": sum(len(f.qualifiers.get("translation", [])) > 1 for f in cds),
                        "deposited_translation_count": len(annotations), "proteins": annotations})
    return {"records": records, "parser_warnings": [str(w.message) for w in caught],
            "protein_errors": protein_errors,
            "nucleotide_sequence_sha256": framed_sequences([s for _, s in nucleotides])}, nucleotides, proteins


def parse_one(item, workdir):
    base = {**item, "nucleotide_status": "unavailable", "protein_status": "unavailable", "outputs": {}}
    if item["payload_status"] != "authenticated":
        return {**base, "status": "payload_error"}
    try:
        data = Path(item["container_path"]).read_bytes()
        verify(data, item["source"])
        summary, dna, aa = parse_container(data, item["source"]["container_id"])
        base.update(summary)
        for kind, entries, extension in [("nucleotide", dna, ".fna"), ("protein", aa, ".faa")]:
            if kind == "protein" and summary["protein_errors"]:
                base["protein_status"] = "annotation_error"
                continue
            path = workdir / ("nucleotides" if kind == "nucleotide" else "proteins") / (item["source"]["container_id"] + extension)
            content = fasta(entries)
            publish(path, content)
            base["outputs"][kind + "_fasta"] = {"path": str(path), "sha256": sha(content), "bytes": len(content), "sequences": len(entries)}
            base[kind + "_status"] = "available" if entries else "no_deposited_translations"
        return {**base, "status": "annotation_error" if summary["protein_errors"] else "inventoried"}
    except Exception as exc:
        return {**base, "status": "parse_or_artifact_error", "error": f"{type(exc).__name__}: {exc}"}


def compact_report(report):
    """Keep complete container coverage; bind protein details in the raw manifest."""
    compact = {**report, "entries": []}
    nonstandard, ambiguous = Counter(), Counter()
    multirecord = []
    for entry in report["entries"]:
        reduced = {**entry, "records": []}
        for record in entry.get("records", []):
            proteins = record.get("proteins", [])
            letters = sum((Counter(p.get("letters_outside_standard_20", {})) for p in proteins), Counter())
            nonstandard.update(letters)
            ambiguous.update({k: v for k, v in record["nucleotide_alphabet_counts"].items() if k not in "ACGT"})
            reduced["records"].append({**{k: v for k, v in record.items() if k != "proteins"},
                "deposited_protein_metadata_sha256": sha(encoded(proteins)),
                "deposited_protein_metadata_hash_scheme": "SHA256 of canonical sorted-key compact JSON list plus newline; full list retained in sequence manifest",
                "deposited_translation_residues": sum(p.get("length", 0) for p in proteins),
                "translation_letters_outside_standard_20": dict(sorted(letters.items())),
                "internal_stop_count": sum(p.get("internal_stop_count", 0) for p in proteins),
                "terminal_stop_count": sum(p.get("terminal_stop_present", False) for p in proteins)})
        compact["entries"].append(reduced)
        if len(entry.get("records", [])) > 1:
            multirecord.append({"primary_id": entry["source"]["primary_id"],
                                "record_lengths": [r["length"] for r in entry["records"]]})
    compact["summary"] = {**report["summary"], "nonstandard_amino_acid_characters": dict(sorted(nonstandard.items())),
                          "ambiguous_nucleotide_characters": dict(sorted(ambiguous.items())),
                          "multi_record_containers": multirecord,
                          "total_nucleotide_bases": sum(r["length"] for e in report["entries"] for r in e.get("records", [])),
                          "total_deposited_translation_residues": sum(r["deposited_translation_residues"] for e in compact["entries"] for r in e["records"])}
    return compact


def compact_existing(source, workdir, output):
    """Refine an archived full report without fetching or re-parsing containers."""
    original = source.read_bytes()
    report = json.loads(original)
    archive_script = workdir / "executed-audit-helper-v1.py"
    if sha(archive_script.read_bytes()) != report["script_sha256"]:
        raise ValueError("Archived executed helper differs from initial report")
    manifest = report["sequence_manifest"]
    if sha(Path(manifest["path"]).read_bytes()) != manifest["sha256"]:
        raise ValueError("Full raw manifest changed")
    if output.resolve() == source.resolve() or output.is_symlink():
        raise ValueError("Cannot replace archived report or symlink")
    if not output.exists() or output.read_bytes() != original:
        raise ValueError("Only the byte-identical initial versioned report may be compacted")
    result = compact_report(report)
    result["inventory_execution_script_sha256"] = report["script_sha256"]
    result["script_sha256"] = sha(Path(__file__).read_bytes())
    result["report_refinement"] = {
        "type": "post_execution_compaction_without_new_fetch_or_parse",
        "supersedes": {"path": str(source), "sha256": sha(original)},
        "initial_executed_helper": {"path": str(archive_script), "sha256": report["script_sha256"]},
        "rationale": "Versioned report retains every container and aggregate record evidence; 48083 detailed deposited protein annotations remain in the unchanged hash-bound ignored manifest.",
        "postprocessing_python": sys.version,
        "sequence_manifest_unchanged": True,
        "new_downloads_or_parser_execution": False}
    data = encoded(result)
    fd, name = tempfile.mkstemp(prefix=output.name + ".", suffix=".part", dir=output.parent)
    temp = Path(name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
        if output.read_bytes() != original:
            raise ValueError("Original report changed during compaction")
        os.replace(temp, output)
    finally:
        temp.unlink(missing_ok=True)
    print(json.dumps({"status": result["status"], "report_bytes": len(data), "sha256": sha(data), "summary": result["summary"]}))
    return 0


def run(policy_path, workdir, output, offline=False):
    import Bio
    from Bio.GenBank import Scanner
    from Bio.SeqIO import InsdcIO
    if Bio.__version__ != "1.86":
        raise ValueError("Frozen parser requires Biopython 1.86")
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    roster_path = Path(policy["roster_path"])
    roster_bytes = roster_path.read_bytes()
    if sha(roster_bytes) != policy["roster_sha256"]:
        raise ValueError("Frozen roster changed")
    roster = json.loads(roster_bytes)
    if policy["source_commit"] != roster["sources"]["commit"]:
        raise ValueError("Source commit differs from frozen roster")
    items = inventory(roster)
    if sha(encoded(items)) != policy["inventory_sha256"]:
        raise ValueError("Frozen inventory changed")
    if policy["download_workers"] != 4 or policy["parse_policy"]["biopython_version"] != Bio.__version__:
        raise ValueError("Unsupported frozen policy")
    if output.exists() or output.is_symlink() or output.resolve().is_relative_to(workdir.resolve()):
        raise ValueError("Report must be new and outside raw workdir")
    protected = [policy_path, roster_path, Path(__file__)]
    if any(p.resolve().is_relative_to(workdir.resolve()) for p in protected):
        raise ValueError("Workdir must not contain source inputs")
    script_bytes = Path(__file__).read_bytes()
    runtime = {"python": sys.version, "biopython": Bio.__version__, "parser_files_sha256": {
        name: sha(Path(module.__file__).read_bytes()) for name, module in [("Bio", Bio), ("Bio.GenBank.Scanner", Scanner), ("Bio.SeqIO.InsdcIO", InsdcIO)]}}
    plan = {"policy_sha256": sha(policy_bytes), "roster_sha256": sha(roster_bytes),
            "script_sha256": sha(script_bytes), "parser": runtime, "inventory": items,
            "download_command_template": ["curl", "--fail", "--location", "--silent", "--show-error", "--max-time", "45", "--retry", "0", "--max-filesize", "<declared_size_bytes>", "--output", "<temporary-path>", "<pinned-raw-url>"]}
    publish(workdir / "frozen-plan.json", encoded(plan))
    with ThreadPoolExecutor(max_workers=4) as pool:
        retrieved = list(pool.map(lambda item: retrieve(item, policy, workdir, offline), items))
    observations = [parse_one(item, workdir) for item in retrieved]
    if Path(__file__).read_bytes() != script_bytes or policy_path.read_bytes() != policy_bytes or roster_path.read_bytes() != roster_bytes:
        raise ValueError("Source inputs changed during execution")
    manifest = {"schema_version": 1, "policy_sha256": sha(policy_bytes), "plan_sha256": sha(encoded(plan)),
                "sequence_normalization": policy["parse_policy"],
                "entries": [{k: v for k, v in item.items() if k != "retrieval"} for item in observations]}
    manifest_bytes = encoded(manifest)
    manifest_path = workdir / "sequence-manifest.json"
    publish(manifest_path, manifest_bytes)
    errors = [x for x in observations if x["status"] != "inventoried"]
    report = {"schema_version": 1, "status": "complete_container_inventory" if not errors else "incomplete_container_inventory",
              "policy_path": str(policy_path), "policy_sha256": sha(policy_bytes), "script_sha256": sha(script_bytes),
              "roster_sha256": sha(roster_bytes), "inventory_sha256": sha(encoded(items)), "runtime": runtime,
              "frozen_plan": {"path": str(workdir / "frozen-plan.json"), "sha256": sha(encoded(plan))},
              "sequence_manifest": {"path": str(manifest_path), "sha256": sha(manifest_bytes)},
              "summary": {"requested_containers": len(items), "cohort_containers": 260, "separate_target_candidates": 1,
                          "authenticated_payloads": sum(x["payload_status"] == "authenticated" for x in observations),
                          "nucleotide_available_containers": sum(x["nucleotide_status"] == "available" for x in observations),
                          "records": sum(len(x.get("records", [])) for x in observations),
                          "deposited_translations": sum(r["deposited_translation_count"] for x in observations for r in x.get("records", [])),
                          "cds_without_translation": sum(r["cds_without_translation"] for x in observations for r in x.get("records", [])),
                          "containers_with_parser_warnings": sum(bool(x.get("parser_warnings")) for x in observations),
                          "errors": len(errors)},
              "entries": observations, "translation_agreement_with_DNA": "not_verified",
              "experimental_stock_identity_verified": False, "relatedness_exclusions_completed": False,
              "training_or_scoring_executed": False, "new_sequence_design_or_annotation_filling": False,
              "evidence_limit": "Authenticated public containers and their deposited annotations; source references and Ando experimental stocks remain distinct. Missing CDS translations remain missing."}
    report = compact_report(report)
    publish(output, encoded(report))
    print(json.dumps({"status": report["status"], "summary": report["summary"], "report_sha256": sha(encoded(report))}))
    return 1 if errors else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/cohort-container-policy.json"))
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--compact-existing-report", type=Path,
                        help="Compact an archived full report; do not download or parse")
    args = parser.parse_args()
    if args.compact_existing_report:
        raise SystemExit(compact_existing(args.compact_existing_report, args.workdir, args.output))
    raise SystemExit(run(args.policy, args.workdir, args.output, args.offline))
