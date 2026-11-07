"""Authenticate author host reference candidates without assigning stock identity."""
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fasta(data, alphabet):
    records, header, sequence = [], None, []
    for line in data.decode("ascii").splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        elif header is None:
            raise ValueError("Expected a FASTA header before sequence")
        else:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    if not records or len({h.split()[0] for h, _ in records if h.split()}) != len(records):
        raise ValueError("Missing or duplicate FASTA identifier")
    if any(not s or not set(s.upper()) <= alphabet for _, s in records):
        raise ValueError("Missing sequence or unsupported source alphabet")
    return records


def run(policy_path, cache, output):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if output.exists() or output.is_symlink():
        raise ValueError("Report must be new")
    for key in ("candidate_report", "target_manifest"):
        item = policy[key]
        if sha(Path(item["path"]).read_bytes()) != item["sha256"]:
            raise ValueError("Frozen source/target manifest changed")
    references = json.loads(Path(policy["target_manifest"]["path"]).read_bytes())["genome_sources"]
    cache.mkdir(parents=True, exist_ok=True)
    outputs = []
    for spec in policy["files"]:
        host, name = spec["source_path"].split("/")[-2:]
        path = cache / host / name
        url = "https://raw.githubusercontent.com/mjohnson11/PhageDataSheets/" + policy["source_commit"] + "/" + spec["source_path"]
        if path.exists():
            raw = path.read_bytes()
        else:
            with urllib.request.urlopen(url, timeout=45) as response:
                raw = response.read(spec["bytes"] + 1)
        blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        if blob != spec["git_blob_sha1"] or len(raw) != spec["bytes"]:
            raise ValueError("Host asset differs from pinned source tree")
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("xb") as stream:
                stream.write(raw)
        entry = {"source_path": spec["source_path"], "path": str(path), "url": url,
                 "bytes": len(raw), "git_blob_sha1": blob, "sha256": sha(raw)}
        try:
            if name == "genes.tab":
                reader = csv.reader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
                entry["header"] = next(reader)
                entry["data_rows"] = sum(1 for _ in reader)
            else:
                records = fasta(raw, set("ACGTRYSWKMBDHVN") if name == "genome.fna" else set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"))
                entry["record_count"] = len(records)
                entry["total_residues"] = sum(len(s) for _, s in records)
                entry["alphabet_counts"] = dict(sorted(Counter(c for _, s in records for c in s.upper()).items()))
                entry["records"] = [{"header": h, "length": len(s), "sequence_sha256": sha(s.upper().encode())}
                                    for h, s in records] if name == "genome.fna" else None
                if name == "genome.fna":
                    target = "BW25113" if host == "Keio" else "BL21"
                    entry["original_reference_comparison"] = {
                        "target_name": target, "original_accession": references[target]["accession"],
                        "original_sequence_sha256": references[target]["sequence_sha256"],
                        "exact_single_record_sequence_match": len(records) == 1 and sha(records[0][1].upper().encode()) == references[target]["sequence_sha256"],
                        "comparison_does_not_establish_experimental_stock_identity": True}
        except (ValueError, UnicodeError, StopIteration) as exc:
            entry["parse_status"] = "unresolved"
            entry["parse_error"] = str(exc)
        else:
            entry["parse_status"] = "inventoried"
        outputs.append(entry)
    if policy_path.read_bytes() != policy_bytes:
        raise ValueError("Policy changed during inspection")
    for item in outputs:
        if sha(Path(item["path"]).read_bytes()) != item["sha256"]:
            raise ValueError("Cached source changed during inspection")
    report = {"schema_version": 1, "policy_sha256": sha(policy_bytes), "script_sha256": sha(Path(__file__).read_bytes()),
              "authenticated_sources": outputs, "host_stock_correspondence_verified": False,
              "host_reference_substitution_performed": False, "phenotype_labels_created": False,
              "models_trained_or_scored": False, "new_benchmark_rows_admitted": 0,
              "interpretation": "Byte-authenticated author host-reference assets and deposited sequence inventory. Browser names and file headers do not authenticate experimental stocks; source sequence differences remain explicit and are not repaired."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"authenticated_files": len(outputs), "parse_status_counts": dict(Counter(x["parse_status"] for x in outputs)), "new_benchmark_rows_admitted": 0}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/host-payload-policy.json"))
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.policy, args.cache, args.output)
