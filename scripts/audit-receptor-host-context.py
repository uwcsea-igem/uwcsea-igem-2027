"""Inventory host contexts without promoting propagation metadata to labels."""
from collections import Counter
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def audit(policy_path, output):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if output.exists() or output.is_symlink():
        raise ValueError("Use a fresh report destination")
    inputs = {}
    for item in policy["inputs"]:
        raw = Path(item["path"]).read_bytes()
        if sha(raw) != item["sha256"]:
            raise ValueError("Frozen host-context source differs")
        inputs[item["name"]] = raw
    reader = csv.DictReader(io.StringIO(inputs["metadata"].decode("utf-8-sig")), delimiter="\t")
    header = reader.fieldnames
    if not header or not set(policy["phenotype_columns"].values()) <= set(header):
        raise ValueError("Declared phenotype fields are missing")
    rows = list(reader)
    primary_ids = [r["Phage"] for r in rows]
    if len(rows) != 260 or len(set(primary_ids)) != len(rows):
        raise ValueError("Require all 260 unique metadata records")
    parser_path = Path(__file__).with_name("compare-receptor-measurements.py")
    parser_hash = sha(parser_path.read_bytes())
    spec = importlib.util.spec_from_file_location("host_context_ooxml", parser_path)
    parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parser)
    strains, _ = parser.parse_workbook(inputs["workbook"], policy["strain_table"])
    mg_rows = [r for r in strains if "MG1655" in r["cells"].get("A", "")]
    report = {
        "schema_version": 1, "purpose": policy["purpose"], "policy_sha256": sha(policy_bytes),
        "script_sha256": sha(Path(__file__).read_bytes()), "ooxml_parser_sha256": parser_hash,
        "authenticated_inputs": policy["inputs"], "metadata_records": len(rows),
        "declared_susceptibility_fields": policy["phenotype_columns"],
        "susceptibility_state_counts": {host: dict(sorted(Counter(r[column] for r in rows).items()))
                                        for host, column in policy["phenotype_columns"].items()},
        "target_hosts_without_a_declared_susceptibility_field": sorted(set(policy["original_target_hosts"]) - set(policy["phenotype_columns"])),
        "production_host_literal_counts": dict(sorted(Counter(r["Production host strain"] for r in rows).items())),
        "exact_original_target_production_host_counts": {host: sum(r["Production host strain"] == host for r in rows)
                                                         for host in policy["original_target_hosts"]},
        "supplementary_strain_rows_containing_literal_MG1655": mg_rows,
        "production_hosts_used_as_phenotype_labels_or_predictive_features": False,
        "MG1655_derivative_names_equated_to_original_target": False,
        "unknown_or_not_assayed_states_converted_to_negatives": False,
        "experimental_stock_identity_verified": False, "training_or_scoring_executed": False,
        "new_benchmark_rows_admitted": 0,
        "interpretation": "The metadata schema declares susceptibility only for BW25113 and BL21. Production-host strings are provenance, not independent target labels or predictor inputs. Derivative names do not supply missing MG1655 observations or establish Ando stock identity. State counts are not infection classifications."
    }
    for item in policy["inputs"]:
        if sha(Path(item["path"]).read_bytes()) != item["sha256"]:
            raise ValueError("Source changed during host-context audit")
    if policy_path.read_bytes() != policy_bytes or sha(parser_path.read_bytes()) != parser_hash:
        raise ValueError("Policy/parser changed during host-context audit")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/host-context-policy.json"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.policy, args.output)
    print(json.dumps({"metadata_records": result["metadata_records"],
                      "target_hosts_without_field": result["target_hosts_without_a_declared_susceptibility_field"],
                      "new_benchmark_rows_admitted": 0}))
