"""Authenticate existing T7/T3 GenBank containers; inspect deposited annotations."""
import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import warnings

import Bio
from Bio import SeqIO


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main(cache, output):
    if Bio.__version__ != "1.86":
        raise ValueError("Use the frozen Biopython 1.86 parser")
    policy_path = Path("software/tipcraft/benchmarks/receptor-resource/t7-t3-container-policy.json")
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    roster = Path("data/processed/receptor-phage-roster-audit.json")
    if sha(roster.read_bytes()) != policy["roster_sha256"]:
        raise ValueError("Frozen roster changed")
    if output.exists() or output.is_symlink():
        raise ValueError("Report must be new")
    cache.mkdir(parents=True, exist_ok=True)
    observations = []
    for spec in policy["files"]:
        path = cache / spec["filename"]
        url = f'https://raw.githubusercontent.com/mjohnson11/PhageDataSheets/{policy["source_commit"]}/{spec["repository_path"]}'
        if path.exists():
            data = path.read_bytes()
        else:
            request = urllib.request.Request(url, headers={"User-Agent": "TipCraft-source-audit"})
            with urllib.request.urlopen(request, timeout=45) as response:
                data = response.read(spec["declared_size_bytes"] + 1)
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if len(data) != spec["declared_size_bytes"] or blob != spec["git_blob_sha1"]:
            raise ValueError("GenBank bytes differ from frozen source tree")
        if not path.exists():
            with path.open("xb") as stream:
                stream.write(data)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            records = list(SeqIO.parse(io.StringIO(data.decode("utf-8")), "genbank"))
        summaries = []
        for record in records:
            cds = [feature for feature in record.features if feature.type == "CDS"]
            translations = [value for feature in cds for value in feature.qualifiers.get("translation", [])]
            letters = Counter("".join(translations))
            sequence = str(record.seq)
            summaries.append({
                "record_id": record.id, "name": record.name, "description": record.description,
                "sequence_length": len(sequence), "sequence_sha256": sha(sequence.encode("ascii")),
                "sequence_alphabet_counts": dict(sorted(Counter(sequence).items())),
                "molecule_type": record.annotations.get("molecule_type"),
                "organism": record.annotations.get("organism"),
                "feature_counts": dict(sorted(Counter(f.type for f in record.features).items())),
                "cds_with_translation": sum(bool(f.qualifiers.get("translation")) for f in cds),
                "cds_without_translation": sum(not f.qualifiers.get("translation") for f in cds),
                "cds_with_multiple_translation_values": sum(len(f.qualifiers.get("translation", [])) > 1 for f in cds),
                "deposited_translation_count": len(translations),
                "deposited_translation_lengths": [len(value) for value in translations],
                "translation_letters_outside_standard_20": dict(sorted((k, v) for k, v in letters.items() if k not in "ACDEFGHIKLMNPQRSTVWY")),
                "translation_sha256_in_feature_order": [sha(value.encode("ascii")) for value in translations],
                "translation_tables_declared": dict(Counter(value for f in cds for value in f.qualifiers.get("transl_table", []))),
                "cds_translation_exception_count": sum("transl_except" in f.qualifiers for f in cds),
                "source_qualifiers": [f.qualifiers for f in record.features if f.type == "source"],
            })
        observations.append({"path": str(path), "url": url, "bytes": len(data), "git_blob_sha1": blob,
                             "sha256": sha(data), "record_count": len(records), "records": summaries,
                             "parser_warnings": [str(w.message) for w in caught]})
    if policy_path.read_bytes() != policy_bytes or sha(roster.read_bytes()) != policy["roster_sha256"]:
        raise ValueError("Source policy/roster changed during audit")
    report = {"schema_version": 1, "status": "two_existing_containers_authenticated",
              "policy_sha256": sha(policy_bytes), "script_sha256": sha(Path(__file__).read_bytes()),
              "parser": policy["parse_policy"]["parser"], "files": observations,
              "experimental_stock_identity_verified": False, "homology_exclusions_completed": False,
              "deposited_translations_validated_against_nucleotide_sequence": False,
              "training_or_scoring_executed": False, "new_sequences_designed": False,
              "admitted_benchmark_rows": 0,
              "interpretation": "Byte-authenticated source containers and deposited annotation inventory only. Names and record accessions do not establish equivalence to Ando experimental stocks. Missing translations remain missing; none were generated."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": report["status"], "files": len(observations),
                      "records": sum(item["record_count"] for item in observations)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main(args.cache, args.output)
