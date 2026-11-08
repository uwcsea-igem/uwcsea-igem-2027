"""Label-free descriptions of ordering within one fitted model at a time.

Callers authenticate score receipts and selected-feature manifests. This module
neither verifies that provenance nor supplies an accuracy test or admission gate.
"""
from itertools import combinations
import math
from numbers import Real

from .learning import _feature_names, _samples, _text


def diagnose_conditional_ordering(samples, paired_predictions, selected_feature_names_by_block):
    """Describe matched host-pair differences across phages within each block.

    Every source row must occur in exactly one nonempty prediction block. A host
    pair is eligible only when both hosts have observed scores under at least two
    common phages in that same block. Exact stored-score ties are retained; no
    tolerance, missing-cell imputation, labels or cross-model pooling is used.
    """
    rows = _samples(samples)
    by_id = {row.row_id: row for row in rows}
    if not isinstance(paired_predictions, dict) or not paired_predictions:
        raise ValueError("Nonempty explicit paired-prediction blocks required")
    if (not isinstance(selected_feature_names_by_block, dict)
            or set(selected_feature_names_by_block) != set(paired_predictions)):
        raise ValueError("Selected-feature manifests must cover exactly the prediction blocks")
    seen, selected = set(), {}
    for block, record in paired_predictions.items():
        if (not _text(block) or not isinstance(record, dict) or set(record) != {"model_id", "scores"}
                or not _text(record["model_id"])):
            raise ValueError("Each valid block requires exactly one model identity and score map")
        scores = record["scores"]
        if not isinstance(scores, dict) or not scores or not set(scores) <= set(by_id) or seen & set(scores):
            raise ValueError("Unknown, duplicate or empty scored-row coverage")
        if any(isinstance(score, bool) or not isinstance(score, Real) or not math.isfinite(score)
               for score in scores.values()):
            raise ValueError("Scores must be finite numeric values")
        seen.update(scores)
        selected[block] = _feature_names(selected_feature_names_by_block[block])
    if seen != set(by_id):
        raise ValueError("Every source row must be scored exactly once")

    blocks = {}
    for block, record in sorted(paired_predictions.items()):
        by_phage = {}
        for row_id, score in record["scores"].items():
            row = by_id[row_id]
            by_phage.setdefault(row.phage_id, {})[row.host_id] = (row_id, float(score))
        hosts = sorted({by_id[row_id].host_id for row_id in record["scores"]})
        eligible, unassessed = [], []
        for first, second in combinations(hosts, 2):
            common = sorted(p for p, values in by_phage.items() if first in values and second in values)
            if len(common) < 2:
                unassessed.append({"host_pair": [first, second], "jointly_observed_phages": common,
                                   "reason": "fewer_than_two_jointly_observed_phages"})
                continue
            differences = {}
            for phage in common:
                left, right = by_phage[phage][first], by_phage[phage][second]
                difference = left[1] - right[1]
                if not math.isfinite(difference):
                    raise ValueError("Derived score difference is not finite")
                differences[phage] = {"row_ids": [left[0], right[0]], "difference": difference,
                                      "order_state": (difference > 0) - (difference < 0)}
            values = [v["difference"] for v in differences.values()]
            lower, upper = min(values), max(values)
            spread = upper - lower
            if not math.isfinite(spread):
                raise ValueError("Derived score-difference range is not finite")
            states = sorted({v["order_state"] for v in differences.values()})
            eligible.append({"host_pair": [first, second], "phage_differences": differences,
                             "minimum_difference": lower, "maximum_difference": upper, "difference_range": spread,
                             "order_states": states, "state_change": len(states) > 1,
                             "strict_reversal": -1 in states and 1 in states})
        possible = len(eligible) + len(unassessed)
        blocks[block] = {
            "model_id": record["model_id"], "scored_rows": len(record["scores"]),
            "observed_hosts": len(hosts), "observed_phages": len(by_phage),
            "selected_host_columns": sum(n.startswith("host:") for n in selected[block]),
            "selected_phage_columns": sum(n.startswith("phage:") for n in selected[block]),
            "possible_host_pairs": possible, "eligible_host_pairs": len(eligible),
            "host_pair_comparison_coverage": len(eligible) / possible if possible else None,
            "host_pairs_with_no_jointly_observed_phage": sum(not p["jointly_observed_phages"] for p in unassessed),
            "host_pairs_with_one_jointly_observed_phage": sum(len(p["jointly_observed_phages"]) == 1 for p in unassessed),
            "eligible_phage_pair_comparisons": sum(math.comb(len(p["phage_differences"]), 2) for p in eligible),
            "host_pairs_with_state_change": sum(p["state_change"] for p in eligible),
            "host_pairs_with_strict_reversal": sum(p["strict_reversal"] for p in eligible),
            "eligible_pairs": eligible, "unassessed_pairs": unassessed}
    return {
        "schema_version": 1, "status": "completed_label_free_conditional_ordering_diagnostic",
        "scored_rows": len(rows), "row_coverage": 1.0, "blocks": blocks,
        "order_state_definition": {"-1": "first canonical host scores lower", "0": "exact score tie",
                                   "1": "first canonical host scores higher"},
        "labels_accepted": False, "cross_model_comparisons": False, "missing_cells_imputed": False,
        "accuracy_or_admission_gate_defined": False,
        "interpretation": "Descriptive behavior of authenticated predictions only. Feature counts do not establish conditional ordering. Exact state changes, reversals and raw difference ranges can reflect tiny floating-point effects; no numerical or biological significance threshold is applied. They do not establish prediction accuracy, a biological mechanism, global score non-additivity or independent replication. No observed change cannot establish global additivity. Repeated host pairs remain separate across blocks/models."}
