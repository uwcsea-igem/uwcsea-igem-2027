"""Post-score evaluation. This module may read labels; the WIsH runner does not."""
import hashlib
import json
import math
from pathlib import Path
from .benchmark import load_cases
from .rank import classification_metrics
from .wish import parse_matrix, sha256


def evaluate_scores(cases: list[dict], inference: dict) -> dict:
    scores = inference["scores"]
    if set(scores) & set(inference.get("excluded_constructs", {})):
        raise ValueError("Excluded constructs cannot have scores")
    if not scores or any(not hosts for hosts in scores.values()):
        raise ValueError("Scores cannot be empty")
    if any(not isinstance(v, (float, int)) or not math.isfinite(v)
           for hosts in scores.values() for v in hosts.values()):
        raise ValueError("Scores must be finite numbers")
    covered = [c for c in cases if c["phage"] in scores and c["host"] in scores[c["phage"]]]
    if not covered:
        raise ValueError("No covered cases")
    winners, per_phage = {}, {}
    for phage, by_host in sorted(scores.items()):
        if len(by_host) < 2:
            raise ValueError("Need at least two candidate hosts")
        best = max(by_host.values())
        winners[phage] = [h for h, score in by_host.items() if score == best]
        rows = [c for c in covered if c["phage"] == phage]
        positives = [c for c in rows if c["label"] == "1"]
        negatives = [c for c in rows if c["label"] == "0"]
        pairs = [(p, n) for p in positives for n in negatives]
        concordant = sum(1.0 if by_host[p["host"]] > by_host[n["host"]] else
                         0.5 if by_host[p["host"]] == by_host[n["host"]] else 0.0
                         for p, n in pairs)
        labels = {c["host"]: int(c["label"]) for c in rows}
        top_success = (sum(labels[h] for h in winners[phage]) / len(winners[phage])
                       if all(h in labels for h in winners[phage]) else None)
        per_phage[phage] = {
            "ranked_hosts": sorted(by_host, key=lambda h: (-by_host[h], h)),
            "top_hosts": winners[phage], "top1_documented_positive": top_success,
            "positive_negative_pairs": len(pairs),
            "pair_concordance": concordant / len(pairs) if pairs else None,
            "labels": labels,
        }
    labels = [int(c["label"]) for c in covered]
    # Explicit sensitivity analysis: WIsH itself does not predict a host-range matrix.
    forced = [float(c["host"] in winners[c["phage"]]) for c in covered]
    successes = [x["top1_documented_positive"] for x in per_phage.values()
                 if x["top1_documented_positive"] is not None]
    details = [{"case_id": c["case_id"], "label": int(c["label"]),
                "mean_log_likelihood": scores[c["phage"]][c["host"]],
                "selected_top1": c["host"] in winners[c["phage"]]} for c in covered]
    return {
        "schema_version": 1,
        "method": inference["method"], "total_cases": len(cases),
        "covered_cases": len(covered), "coverage": len(covered) / len(cases),
        "abstained_case_ids": [c["case_id"] for c in cases if c not in covered],
        "primary": {"top1_documented_positive_rate": sum(successes) / len(successes) if successes else None,
                    "phage_queries": len(successes)},
        "per_phage": per_phage,
        "forced_choice_sensitivity": classification_metrics(labels, forced),
        "always_positive_on_same_subset": classification_metrics(labels, [1.0] * len(labels)),
        "cases": details,
        "gate": "HOLD: this genome-composition method has not established strain-level discrimination or assembly prediction",
    }


def leakage_audit(cases: list[dict], inference: dict, anonymous: dict) -> dict:
    if inference.get("no_labels_loaded") is not True or anonymous.get("no_labels_loaded") is not True:
        raise ValueError("Scoring runs must explicitly declare that no labels were loaded")
    if set(anonymous["scores"]) & set(anonymous.get("excluded_constructs", {})):
        raise ValueError("Excluded constructs cannot have anonymous scores")
    if inference.get("anonymized") is not False or anonymous.get("anonymized") is not True:
        raise ValueError("Leakage audit requires distinct named and anonymized scoring runs")
    for key in ("method", "order", "alpha", "upstream_commit", "executable_sha256", "manifest_sha256"):
        if not inference.get(key) or inference.get(key) != anonymous.get(key):
            raise ValueError(f"Scoring runs have incompatible provenance: {key}")
    original = inference["scores"]
    name_invariance = original == anonymous["scores"]
    original_json = json.dumps(original, sort_keys=True)
    frozen_hash = hashlib.sha256(original_json.encode()).hexdigest()
    # Exhaustively permute the one positive among T3's three candidate hosts.
    t3 = [c for c in cases if c["phage"] == "T3"]
    if len(t3) != 3 or sum(c["label"] == "1" for c in t3) != 1:
        raise ValueError("This permutation audit requires three T3 cases with exactly one positive")
    placements = []
    for positive in t3:
        altered = [dict(c, label=str(int(c["case_id"] == positive["case_id"])))
                   if c["phage"] == "T3" else dict(c) for c in cases]
        evaluated = evaluate_scores(altered, inference)
        placements.append({"positive_host": positive["host"],
                           "pair_concordance": evaluated["per_phage"]["T3"]["pair_concordance"]})
    if json.dumps(inference["scores"], sort_keys=True) != original_json:
        raise ValueError("Evaluator modified scores")
    observed = evaluate_scores(cases, inference)["per_phage"]["T3"]["pair_concordance"]
    return {
        "run_provenance_matches": True,
        "name_and_fasta_header_invariance": name_invariance,
        "status": "passed_with_limitations" if name_invariance else "failed",
        "frozen_scores_sha256": frozen_hash,
        "label_permutation_changes_only_evaluation": True,
        "t3_all_single_positive_placements": placements,
        "t3_exact_permutation_upper_tail_fraction": sum(p["pair_concordance"] >= observed for p in placements) / len(placements),
        "selection_after_labels_visible": True,
        "prospective_blind_validation": False,
        "training_uses_interaction_labels": not inference["no_labels_loaded"],
        "conclusion": ("No label/name input path found in scoring; post-hoc method selection and stock mismatch remain limitations"
                       if name_invariance else "FAILED name invariance: do not claim label/name independence"),
    }


def validate_artifacts(inference: dict, anonymous: dict, manifest_path: Path,
                       named_matrix: Path, anonymous_matrix: Path) -> dict:
    """Bind JSON scores to the frozen manifest and upstream numeric matrices.

    This establishes artifact consistency, not authenticity of experimental stocks.
    """
    manifest = json.loads(manifest_path.read_text())
    expected_hash = sha256(manifest_path)
    for run in (inference, anonymous):
        if run.get("manifest_sha256") != expected_hash:
            raise ValueError("Run does not match the supplied manifest hash")
        for key in ("order", "alpha", "upstream_commit"):
            if run.get(key) != manifest[key]:
                raise ValueError(f"Run differs from frozen manifest: {key}")
        if run.get("excluded_constructs") != manifest["excluded_constructs"]:
            raise ValueError("Excluded construct declarations differ from manifest")
        if set(run["scores"]) != set(manifest["phages"]):
            raise ValueError("Scored phages differ from manifest")
        for hosts in run["scores"].values():
            if set(hosts) != set(manifest["hosts"]):
                raise ValueError("Scored candidate hosts differ from manifest")
    expected_named = parse_matrix(named_matrix)
    raw_anonymous = parse_matrix(anonymous_matrix)
    inverse_phages = {f"p{i}": name for i, name in enumerate(manifest["phages"])}
    inverse_hosts = {f"h{i}": name for i, name in enumerate(manifest["hosts"])}
    try:
        expected_anonymous = {
            inverse_phages[p]: {inverse_hosts[h]: v for h, v in hosts.items()}
            for p, hosts in raw_anonymous.items()}
    except KeyError as error:
        raise ValueError("Anonymous matrix identifiers differ from manifest") from error
    if inference["scores"] != expected_named or anonymous["scores"] != expected_anonymous:
        raise ValueError("JSON scores differ from the upstream likelihood matrices")
    return {"status": "verified", "manifest_sha256": expected_hash,
            "named_matrix_sha256": sha256(named_matrix),
            "anonymous_matrix_sha256": sha256(anonymous_matrix),
            "candidate_sets_match": True, "exclusions_match": True,
            "scope": "artifact consistency only; no stock identity or model validity claim"}
