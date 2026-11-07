"""Check deposited CDS translations against authenticated source DNA, without exports."""
import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import sys
import warnings


UNSUPPORTED = {"transl_except", "pseudo", "pseudogene", "ribosomal_slippage",
               "exception", "trans_splicing", "artificial_location", "partial"}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                       allow_nan=False) + "\n").encode()


def qualifiers(q):
    unsupported = sorted(UNSUPPORTED.intersection(q))
    if unsupported:
        return {"status": "unverified", "reason": "unsupported_qualifiers", "qualifier_keys": unsupported}
    table = q.get("transl_table", [])
    if len(table) != 1 or not str(table[0]).isdigit():
        return {"status": "unverified", "reason": "missing_multiple_or_invalid_transl_table"}
    codon = q.get("codon_start", [])
    if ("codon_start" in q and len(codon) != 1) or (codon and str(codon[0]) not in {"1", "2", "3"}):
        return {"status": "unverified", "reason": "multiple_or_invalid_codon_start"}
    return {"status": "eligible", "transl_table": int(table[0]),
            "codon_start": int(codon[0]) if codon else 1,
            "codon_start_source": "declared" if codon else "Biopython_1.86_absence_default_1"}


def check(feature, dna, deposited):
    """Return hashes/status only; no translated amino-acid string leaves this call."""
    from Bio.Data import CodonTable
    from Bio.SeqFeature import ExactPosition
    settings = qualifiers(feature.qualifiers)
    if settings["status"] != "eligible":
        return settings
    if len(feature.qualifiers.get("translation", [])) != 1:
        return {**settings, "status": "unverified", "reason": "multiple_or_missing_translation_values"}
    loc = feature.location
    if loc is None or getattr(loc, "operator", None) not in {None, "join"}:
        return {**settings, "status": "unverified", "reason": "missing_or_unsupported_location"}
    parts = list(loc.parts)
    if any(type(p.start) is not ExactPosition or type(p.end) is not ExactPosition for p in parts):
        return {**settings, "status": "unverified", "reason": "partial_or_fuzzy_location"}
    if any(p.ref or p.ref_db for p in parts):
        return {**settings, "status": "unverified", "reason": "remote_location"}
    if any(p.strand not in {-1, 1} or not 0 <= int(p.start) < int(p.end) <= len(dna) for p in parts):
        return {**settings, "status": "unverified", "reason": "invalid_bounds_or_strand"}
    if len({p.strand for p in parts}) != 1:
        return {**settings, "status": "unverified", "reason": "mixed_strand_join"}
    intervals = sorted((int(p.start), int(p.end)) for p in parts)
    if any(b[0] < a[1] for a, b in zip(intervals, intervals[1:])):
        return {**settings, "status": "unverified", "reason": "overlapping_join"}
    if settings["transl_table"] not in CodonTable.unambiguous_dna_by_id:
        return {**settings, "status": "unverified", "reason": "unknown_transl_table"}
    coding = feature.extract(dna)[settings["codon_start"] - 1:]
    if not set(str(coding).upper()) <= set("ACGT"):
        return {**settings, "status": "unverified", "reason": "ambiguous_coding_DNA"}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            translated = str(coding.translate(table=settings["transl_table"], cds=True))
        except Exception as exc:
            return {**settings, "status": "cds_validation_failed", "reason": type(exc).__name__,
                    "detail": str(exc), "coding_nucleotide_length": len(coding)}
    if caught:
        return {**settings, "status": "unverified", "reason": "translation_warning",
                "warnings": [str(w.message) for w in caught]}
    return {**settings, "status": "exact_match" if translated == deposited else "mismatch",
            "coding_nucleotide_length": len(coding), "computed_length": len(translated),
            "computed_sha256": sha(translated.encode("ascii"))}


def authenticated_inputs(policy):
    report_bytes = Path(policy["container_report"]["path"]).read_bytes()
    if sha(report_bytes) != policy["container_report"]["sha256"]:
        raise ValueError("Frozen container report changed")
    report = json.loads(report_bytes)
    spec = policy["sequence_manifest"]
    if report["sequence_manifest"] != spec:
        raise ValueError("Report/manifest cross-binding differs")
    manifest_bytes = Path(spec["path"]).read_bytes()
    if sha(manifest_bytes) != spec["sha256"]:
        raise ValueError("Frozen raw sequence manifest changed")
    manifest = json.loads(manifest_bytes)
    items = manifest["entries"]
    frozen = [{"container_path": x["container_path"], "sha256": x["sha256"],
               "git_blob_sha1": x["git_blob_sha1"], "size_bytes": x["size_bytes"]} for x in items]
    if sha(encoded(frozen)) != policy["container_inventory_sha256"] or len(items) != 261:
        raise ValueError("Frozen261-container inventory changed")
    count = sum(len(r["proteins"]) for x in items for r in x["records"])
    if count != policy["expected_translation_entries"]:
        raise ValueError("Frozen translation-entry count changed")
    payloads = []
    for item in items:
        data = Path(item["container_path"]).read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if len(data) != item["size_bytes"] or sha(data) != item["sha256"] or blob != item["git_blob_sha1"]:
            raise ValueError("Authenticated container bytes changed")
        payloads.append(data)
    return items, payloads, manifest_bytes


def run(policy_path, workdir, output):
    import Bio
    from Bio import SeqIO, SeqFeature, Seq
    from Bio.Data import CodonTable
    from Bio.GenBank import Scanner
    if Bio.__version__ != "1.86":
        raise ValueError("Frozen parser requires Biopython1.86")
    if output.exists() or output.is_symlink() or workdir.exists() or workdir.is_symlink():
        raise ValueError("New workdir and output required")
    if output.resolve().is_relative_to(workdir.resolve()):
        raise ValueError("Versioned report must be outside raw workdir")
    policy_bytes, script_bytes = policy_path.read_bytes(), Path(__file__).read_bytes()
    policy = json.loads(policy_bytes)
    items, payloads, original_manifest = authenticated_inputs(policy)
    runtime = {"python": sys.version, "biopython": Bio.__version__, "module_sha256": {
        name: sha(Path(module.__file__).read_bytes()) for name, module in
        [("Bio", Bio), ("Bio.SeqFeature", SeqFeature), ("Bio.Seq", Seq), ("Bio.Data.CodonTable", CodonTable), ("Bio.GenBank.Scanner", Scanner)]}}
    plan = {"policy_sha256": sha(policy_bytes), "script_sha256": sha(script_bytes),
            "input_manifest_sha256": sha(original_manifest), "runtime": runtime,
            "container_inventory_sha256": policy["container_inventory_sha256"]}
    workdir.mkdir(parents=True)
    (workdir / "frozen-plan.json").write_bytes(encoded(plan))
    details, summaries = [], []
    for item, data in zip(items, payloads):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            records = list(SeqIO.parse(io.StringIO(data.decode("utf-8")), "genbank"))
        if len(records) != len(item["records"]):
            raise ValueError("Parsed record count differs from frozen manifest")
        current = []
        for record, reference in zip(records, item["records"]):
            if record.id != reference["record_id"] or sha(str(record.seq).upper().encode("ascii")) != reference["sequence_sha256"]:
                raise ValueError("Parsed DNA record differs from frozen manifest")
            cds = [f for f in record.features if f.type == "CDS"]
            observed = [(ci, ti, f, value.upper()) for ci, f in enumerate(cds) for ti, value in enumerate(f.qualifiers.get("translation", []))]
            if len(observed) != len(reference["proteins"]):
                raise ValueError("Deposited translation count differs from frozen manifest")
            for (ci, ti, feature, deposited), protein in zip(observed, reference["proteins"]):
                if (ci, ti, str(feature.location), sha(deposited.encode("ascii"))) != (protein["cds_index"], protein["translation_index"], protein["source_location"], protein["sequence_sha256"]):
                    raise ValueError("Deposited translation mapping differs from frozen manifest")
                result = {"primary_id": item["source"]["primary_id"], "protein_id": protein["protein_id"],
                          "record_id": reference["record_id"], "neutral_nucleotide_id": reference["neutral_nucleotide_id"],
                          "cds_index": ci, "translation_index": ti, "source_location": str(feature.location),
                          "deposited_length": len(deposited), "deposited_sha256": protein["sequence_sha256"],
                          **check(feature, record.seq, deposited)}
                if caught:
                    result = {**result, "status": "unverified", "reason": "container_parser_warning"}
                current.append(result)
        details.extend(current)
        summaries.append({"primary_id": item["source"]["primary_id"], "cohort_member": item["source"]["cohort_member"],
                          "container_sha256": item["sha256"], "translation_entries": len(current),
                          "status_counts": dict(sorted(Counter(x["status"] for x in current).items())),
                          "reason_counts": dict(sorted(Counter(x["reason"] for x in current if "reason" in x).items())),
                          "parser_warnings": [str(w.message) for w in caught]})
    if len(details) != policy["expected_translation_entries"] or len({x["protein_id"] for x in details}) != len(details):
        raise ValueError("Translation-entry coverage is incomplete or duplicated")
    if policy_path.read_bytes() != policy_bytes or Path(__file__).read_bytes() != script_bytes:
        raise ValueError("Policy or helper changed during checking")
    # Re-authenticate all source bytes and the unchanged raw manifest after checks.
    authenticated_inputs(policy)
    detail_path = workdir / "translation-check-details.json"
    detail_bytes = encoded({"schema_version": 1, "plan_sha256": sha(encoded(plan)), "entries": details})
    detail_path.write_bytes(detail_bytes)
    counts = dict(sorted(Counter(x["status"] for x in details).items()))
    reasons = dict(sorted(Counter(x["reason"] for x in details if "reason" in x).items()))
    report = {"schema_version": 1, "status": "completed_translation_consistency_audit",
              "policy_path": str(policy_path), "policy_sha256": sha(policy_bytes), "script_sha256": sha(script_bytes),
              "inputs": {"container_report": policy["container_report"], "sequence_manifest": policy["sequence_manifest"],
                         "container_inventory_sha256": policy["container_inventory_sha256"]}, "runtime": runtime,
              "frozen_plan": {"path": str(workdir / "frozen-plan.json"), "sha256": sha(encoded(plan))},
              "details": {"path": str(detail_path), "sha256": sha(detail_bytes)},
              "summary": {"containers": len(items), "translation_entries": len(details), "status_counts": counts,
                          "reason_counts": reasons, "all_entries_exactly_verified": counts.get("exact_match", 0) == len(details)},
              "containers": summaries,
              "nonmatch_identifiers": [{k: x[k] for k in ["primary_id", "protein_id", "status", "reason"] if k in x}
                                       for x in details if x["status"] in {"mismatch", "cds_validation_failed"}],
              "computed_protein_sequences_exported": False, "missing_annotations_filled": False,
              "original_sequence_manifest_unchanged": True, "training_or_scoring_executed": False,
              "interpretation": "Annotation/DNA consistency only. Unsupported annotations remain unverified; exact translation agreement does not establish experimental-stock identity, receptor usage or infection."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream:
        stream.write(encoded(report))
    print(json.dumps({"summary": report["summary"], "report_sha256": sha(encoded(report))}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/cds-translation-policy.json"))
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.policy, args.workdir, args.output)
