"""Audit published source roles against verified phage groups, without labels."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def partition(containers, membership, excluded):
    keys = [x["container_key"] for x in containers]
    if len(keys) != len(set(keys)) or not set(keys) <= set(membership):
        raise ValueError("Duplicate or ungrouped container")
    cohort = {x["container_key"]: x for x in containers if x["cohort_member"]}
    if any(x["published_model_role"] not in {"training", "testing", "mutant"} for x in cohort.values()):
        raise ValueError("Unknown published source role")
    eligible = {k: x for k, x in cohort.items() if k not in excluded}
    training = {k for k, x in eligible.items() if x["published_model_role"] == "training"}
    testing = {k for k, x in eligible.items() if x["published_model_role"] == "testing"}
    test_groups = {membership[k] for k in testing}
    train_groups = {membership[k] for k in training}
    overlap = test_groups & train_groups
    quarantined = {k for k in training if membership[k] in overlap}
    candidate_training = training - quarantined
    mutants = {k for k, x in eligible.items() if x["published_model_role"] == "mutant"}
    buckets = {"target_related_cohort_exclusions": set(cohort) & excluded,
               "candidate_training": candidate_training, "candidate_testing": testing,
               "training_group_quarantine": quarantined, "mutant_nontraining": mutants}
    union = set().union(*buckets.values())
    if union != set(cohort) or sum(map(len, buckets.values())) != len(cohort):
        raise ValueError("Cohort partition is incomplete or overlapping")
    if {membership[k] for k in candidate_training} & test_groups:
        raise ValueError("Candidate training/test groups overlap")
    return {"cohort_records": len(cohort), "published_eligible_training": len(training),
            "published_eligible_testing": len(testing), "published_training_groups": len(train_groups),
            "published_testing_groups": len(test_groups), "shared_published_role_groups": sorted(overlap),
            "testing_members_in_shared_groups": sorted(k for k in testing if membership[k] in overlap),
            "counts": {name: len(values) for name, values in buckets.items()},
            "buckets": {name: sorted(values) for name, values in buckets.items()},
            "candidate_training_group_count": len({membership[k] for k in candidate_training}),
            "candidate_testing_group_count": len(test_groups),
            "all_cohort_records_accounted_once": True, "candidate_train_test_group_overlap": []}


def run(policy_path, output):
    policy_bytes = policy_path.read_bytes()
    policy = json.loads(policy_bytes)
    if output.exists() or output.is_symlink():
        raise ValueError("Use a fresh role-audit report")
    for name in ("group_report", "numeric_matrix"):
        if sha(policy[name]["path"]) != policy[name]["sha256"]:
            raise ValueError("Frozen verified grouping changed")
    group = json.loads(Path(policy["group_report"]["path"]).read_bytes())
    if (group["status"] != "completed_container_grouping" or not group["numeric_matrix_for_grouping_verified"]
            or not group["numeric_validation"]["dataset_specific_directional_equivalence_verified"]
            or group["numeric_matrix"] != policy["numeric_matrix"]):
        raise ValueError("Verified complete container grouping is required")
    summary = partition(group["containers"], group["grouping"]["container_groups"],
                        set(group["grouping"]["primary_excluded_source_keys"]))
    if summary["cohort_records"] != 260:
        raise ValueError("Expected all 260 source cohort records")
    for name in ("group_report", "numeric_matrix"):
        if sha(policy[name]["path"]) != policy[name]["sha256"]:
            raise ValueError("Group/matrix changed during role audit")
    if policy_path.read_bytes() != policy_bytes:
        raise ValueError("Policy changed during role audit")
    report = {"schema_version": 1, "status": "sequence_partition_candidate_not_model_admission",
              "policy_sha256": hashlib.sha256(policy_bytes).hexdigest(), "script_sha256": sha(__file__),
              "group_report": policy["group_report"], "numeric_matrix": policy["numeric_matrix"], **summary,
              "phenotype_labels_read": False, "feature_vocabulary_fitted": False,
              "models_trained_or_scored": False, "final_training_admission": False,
              "interpretation": "Published training/testing roles are preserved as metadata. Project-relatedness overlap is disclosed, and a deterministic candidate split keeps all eligible testing phages while quarantining related training phages. This is not a blinded test, final label/feature admission or a claim against the source study's stated prospective-isolate validation."}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"counts": report["counts"], "testing_members_in_shared_groups": len(report["testing_members_in_shared_groups"]),
                      "shared_groups": len(report["shared_published_role_groups"]), "final_training_admission": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("software/tipcraft/benchmarks/receptor-resource/source-role-overlap-policy.json"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.policy, args.output)
