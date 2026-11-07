"""Post hoc mismatch description; fixed existing translation parameters only."""
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import sys
import warnings

import Bio
from Bio import SeqIO, SeqFeature, Seq
from Bio.Data import CodonTable

ROOT = Path(__file__).parent
sha = lambda data: hashlib.sha256(data).hexdigest()
plan_bytes = (ROOT / "frozen-plan.json").read_bytes()
plan = json.loads(plan_bytes)
def read(spec):
    data = Path(spec["path"]).read_bytes()
    assert sha(data) == spec["sha256"]
    return json.loads(data)
report = read(plan["baseline_report"])
details = read(plan["baseline_details"])["entries"]
manifest = read(plan["sequence_manifest"])
assert Bio.__version__ == "1.86"
for name, module in [("Bio", Bio), ("Bio.SeqFeature", SeqFeature), ("Bio.Seq", Seq), ("Bio.Data.CodonTable", CodonTable)]:
    assert sha(Path(module.__file__).read_bytes()) == report["runtime"]["module_sha256"][name]
selected = {r["protein_id"]: r for r in details if r["status"] == "mismatch"}
assert len(selected) == 3666
rows, first_deposited, first_computed = [], Counter(), Counter()
for item in manifest["entries"]:
    wanted = {p["protein_id"] for r in item["records"] for p in r["proteins"]} & selected.keys()
    if not wanted:
        continue
    data = Path(item["container_path"]).read_bytes()
    assert sha(data) == item["sha256"]
    records = list(SeqIO.parse(io.StringIO(data.decode()), "genbank"))
    for record, source in zip(records, item["records"]):
        cds = [f for f in record.features if f.type == "CDS"]
        for protein in source["proteins"]:
            identifier = protein["protein_id"]
            if identifier not in selected:
                continue
            baseline = selected[identifier]
            feature = cds[protein["cds_index"]]
            deposited = feature.qualifiers["translation"][protein["translation_index"]].upper()
            assert sha(deposited.encode()) == baseline["deposited_sha256"]
            assert int(feature.qualifiers["transl_table"][0]) == baseline["transl_table"]
            coding = feature.extract(record.seq)[baseline["codon_start"]-1:]
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                computed = str(coding.translate(table=baseline["transl_table"], cds=True))
            assert not caught and sha(computed.encode()) == baseline["computed_sha256"] and computed != deposited
            terminal_only = ((computed.endswith("*") and computed[:-1] == deposited) or
                             (deposited.endswith("*") and deposited[:-1] == computed))
            first_only = len(computed) == len(deposited) and bool(computed) and computed[1:] == deposited[1:]
            category = ("terminal_stop_only" if terminal_only else "first_residue_only" if first_only else
                        "other_same_length" if len(computed) == len(deposited) else "length_difference_other")
            if first_only:
                first_deposited[deposited[0]] += 1
                first_computed[computed[0]] += 1
            rows.append({"primary_id": item["source"]["primary_id"], "protein_id": identifier,
                         "category": category, "deposited_length": len(deposited), "computed_length": len(computed),
                         "deposited_sha256": baseline["deposited_sha256"], "computed_sha256": baseline["computed_sha256"],
                         "same_length_differing_positions": sum(a != b for a, b in zip(computed, deposited)) if len(computed) == len(deposited) else None})
assert len(rows) == len(selected) and len({r["protein_id"] for r in rows}) == len(rows)
categories = Counter(r["category"] for r in rows)
errors = Counter()
for row in details:
    if row["status"] != "cds_validation_failed":
        continue
    message = row["detail"]
    category = ("first_codon_not_start" if message.startswith("First codon") else
                "final_codon_not_stop" if message.startswith("Final codon") else
                "length_not_multiple3" if "multiple of three" in message else
                "internal_stop" if "Extra in frame stop" in message else "other")
    errors[category] += 1
result = {"schema_version": 1, "analysis_type": "post_hoc_descriptive_followup_not_reclassification",
          "plan_sha256": sha(plan_bytes), "helper_sha256": sha(Path(__file__).read_bytes()),
          "baseline_report": plan["baseline_report"], "baseline_details": plan["baseline_details"],
          "sequence_manifest": plan["sequence_manifest"], "mismatch_entries": len(rows),
          "category_counts": dict(sorted(categories.items())),
          "first_residue_only_deposited_first_letter_counts": dict(sorted(first_deposited.items())),
          "first_residue_only_computed_first_letter_counts": dict(sorted(first_computed.items())),
          "existing_cds_validation_failure_counts": dict(sorted(errors.items())),
          "representatives": {c: sorted([r for r in rows if r["category"] == c], key=lambda r:r["protein_id"])[:3] for c in sorted(categories)},
          "entries": rows, "original_statuses_changed": False,
          "computed_AA_sequences_exported": False, "annotations_corrected": False}
assert (ROOT / "frozen-plan.json").read_bytes() == plan_bytes
read(plan["baseline_report"]); read(plan["baseline_details"]); read(plan["sequence_manifest"])
out = ROOT / "diagnostic.json"
with out.open("x") as stream:
    json.dump(result, stream, sort_keys=True, separators=(",", ":")); stream.write("\n")
print(json.dumps({k:result[k] for k in ["mismatch_entries", "category_counts", "first_residue_only_deposited_first_letter_counts", "first_residue_only_computed_first_letter_counts", "existing_cds_validation_failure_counts", "representatives"]}))
print("diagnostic_sha256", sha(out.read_bytes()))
