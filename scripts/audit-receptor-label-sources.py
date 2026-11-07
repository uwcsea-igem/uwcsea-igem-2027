"""Compare all published label-source literals; do not infer negative classes."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main(workdir, output):
    policy_path = Path("software/tipcraft/benchmarks/receptor-resource/label-source-policy.json")
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if any(p.exists() or p.is_symlink() for p in (workdir, output)):
        raise ValueError("Use fresh audit destinations")
    data = []
    for item in policy["inputs"]:
        raw = Path(item["path"]).read_bytes()
        if sha(raw) != item["sha256"]:
            raise ValueError("Frozen source input differs")
        data.append(raw)
    parser_path = Path(__file__).with_name("compare-receptor-measurements.py")
    parser_hash = sha(parser_path.read_bytes())
    spec = importlib.util.spec_from_file_location("source_ooxml", parser_path)
    parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parser)
    rows, details = parser.parse_workbook(data[0], policy["workbook"])
    fields = policy["fields"]
    source = [{"source_row": r["row"], **{name: r["cells"].get(col, "") for col, name in policy["workbook"]["columns"].items()}}
              for r in rows]
    browser = [{"source_row": i, **{name: r[name] for name in fields}}
               for i, r in enumerate(csv.DictReader(io.StringIO(data[1].decode("utf-8-sig")), delimiter="\t"), 2)]
    indexes = []
    for records in (source, browser):
        index = defaultdict(list)
        for r in records:
            index[r["Phage"]].append(r)
        indexes.append(index)
    differences, unresolved, paired = [], [], []
    for key in sorted(set(indexes[0]) | set(indexes[1])):
        left, right = indexes[0].get(key, []), indexes[1].get(key, [])
        if not key or len(left) != 1 or len(right) != 1:
            unresolved.append({"primary_id": key, "workbook_rows": left, "browser_rows": right})
            continue
        paired.append(key)
        for field in fields:
            if left[0][field] != right[0][field]:
                differences.append({"primary_id": key, "field": field, "workbook_raw": left[0][field], "browser_raw": right[0][field]})
    vocab = {name: {field: dict(sorted(Counter(r[field] for r in records).items())) for field in fields if field != "Phage"}
             for name, records in (("workbook", source), ("browser", browser))}
    for item in policy["inputs"]:
        if sha(Path(item["path"]).read_bytes()) != item["sha256"]:
            raise ValueError("Source changed during audit")
    if policy_path.read_bytes() != policy_bytes or sha(parser_path.read_bytes()) != parser_hash:
        raise ValueError("Policy/parser changed during audit")
    workdir.mkdir(parents=True, exist_ok=False)
    artifact = workdir / "literal-label-source-details.json"
    artifact.write_text(json.dumps({"source_rows": source, "browser_rows": browser, "differences": differences,
                                    "unresolved_keys": unresolved, "literal_vocabularies": vocab}, indent=2, ensure_ascii=False) + "\n")
    report = {"schema_version": 1, "purpose": policy["purpose"], "policy_sha256": sha(policy_bytes),
              "script_sha256": sha(Path(__file__).read_bytes()), "ooxml_parser_sha256": parser_hash,
              "authenticated_inputs": policy["inputs"], "workbook_parse": details,
              "workbook_rows": len(source), "browser_rows": len(browser), "unique_paired_ids": len(paired),
              "unresolved_keys": len(unresolved), "literal_difference_count": len(differences),
              "differences_by_field": dict(Counter(d["field"] for d in differences)),
              "literal_vocabularies": vocab,
              "detail_artifact": {"path": str(artifact), "sha256": sha(artifact.read_bytes())},
              "missing_or_unlisted_receptors_converted_to_negative": False, "training_or_scoring_executed": False,
              "new_benchmark_rows_admitted": 0,
              "interpretation": "Literal source comparison and vocabulary inventory only. Matching source strings do not establish negative-label completeness, direct binding, stock identity or independent test eligibility."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: report[k] for k in ("workbook_rows", "browser_rows", "unique_paired_ids", "unresolved_keys", "literal_difference_count", "differences_by_field")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.workdir, args.output)
