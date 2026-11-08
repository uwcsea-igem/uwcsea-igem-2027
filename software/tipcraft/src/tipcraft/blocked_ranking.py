"""Matched ranking diagnostics within one fitted model's test block.

Raw margins from different outer-fold models have no common ranking scale.
This evaluator deliberately never compares hosts across test blocks. Model
identities are provenance declarations; consumers must authenticate bundles.
"""
import math
from numbers import Real

from .learning import _samples, _text


def evaluate_blocked_rankings(samples, labels, test_blocks, methods):
    """Require complete, identical row coverage for every supplied method.

    ``methods[name][block]`` contains ``model_id`` and a row-ID-to-score map.
    The primary summary averages informative blocks per phage, then phages per
    relatedness group, then groups equally. Ties receive 0.5. No uncertainty or
    independent-replicate interpretation is attached to these descriptive means.
    """
    rows = _samples(samples)
    by_id = {row.row_id: row for row in rows}
    expected = set(by_id)
    if not isinstance(labels, dict) or set(labels) != expected:
        raise ValueError("Labels must cover exactly the evaluated rows")
    if any(isinstance(v, bool) or not isinstance(v, Real) or not math.isfinite(v)
           or v not in (0, 1) for v in labels.values()):
        raise ValueError("Labels must be finite binary numbers")
    if not isinstance(test_blocks, dict) or not test_blocks:
        raise ValueError("Explicit nonempty test blocks required")
    seen = set()
    for block, ids in test_blocks.items():
        if not _text(block) or not isinstance(ids, (list, tuple)) or not ids:
            raise ValueError("Each block needs a valid ID and nonempty row list")
        if any(not _text(i) for i in ids) or len(set(ids)) != len(ids):
            raise ValueError("Invalid or duplicate block rows")
        if not set(ids) <= expected or seen & set(ids):
            raise ValueError("Unknown or multiply tested rows")
        seen.update(ids)
    if seen != expected:
        raise ValueError("Every evaluated row must be tested exactly once")
    if not isinstance(methods, dict) or not methods or "tied" in methods:
        raise ValueError("Named methods required; tied is reserved")
    for name, predictions in methods.items():
        if not _text(name) or not isinstance(predictions, dict) or set(predictions) != set(test_blocks):
            raise ValueError("Each method must cover exactly the test blocks")
        for block, ids in test_blocks.items():
            record = predictions[block]
            if not isinstance(record, dict) or set(record) != {"model_id", "scores"} or not _text(record["model_id"]):
                raise ValueError("Each block requires one declared model identity and scores")
            scores = record["scores"]
            if not isinstance(scores, dict) or set(scores) != set(ids):
                raise ValueError("Every method must score exactly the same rows in each block")
            if any(isinstance(s, bool) or not isinstance(s, Real) or not math.isfinite(s) for s in scores.values()):
                raise ValueError("Scores must be finite numbers")
    units = []
    for block, ids in sorted(test_blocks.items()):
        for phage in sorted({by_id[i].phage_id for i in ids}):
            query = sorted(i for i in ids if by_id[i].phage_id == phage)
            positive = [i for i in query if labels[i] == 1]
            negative = [i for i in query if labels[i] == 0]
            pairs = [(p, n) for p in positive for n in negative]
            concordance = {"tied": 0.5 if pairs else None}
            for name, predictions in sorted(methods.items()):
                scores = predictions[block]["scores"]
                concordance[name] = (sum(1 if scores[p] > scores[n] else
                                          0.5 if scores[p] == scores[n] else 0
                                          for p, n in pairs) / len(pairs)) if pairs else None
            units.append({"block_id": block, "phage_id": phage,
                          "phage_group": by_id[query[0]].phage_group,
                          "row_ids": query, "positive_negative_pairs": len(pairs),
                          "concordance": concordance})
    per_phage, per_group, macro = {}, {}, {}
    for name in ["tied", *sorted(methods)]:
        queries = {}
        for unit in units:
            if unit["positive_negative_pairs"]:
                queries.setdefault(unit["phage_id"], []).append(unit["concordance"][name])
        per_phage[name] = {p: sum(v) / len(v) for p, v in sorted(queries.items())}
        groups = {}
        mapping = {row.phage_id: row.phage_group for row in rows}
        for phage, value in per_phage[name].items():
            groups.setdefault(mapping[phage], []).append(value)
        per_group[name] = {g: sum(v) / len(v) for g, v in sorted(groups.items())}
        values = list(per_group[name].values())
        macro[name] = sum(values) / len(values) if values else None
    return {"schema_version": 1, "covered_rows": len(rows), "coverage": 1.0,
            "cross_model_host_comparisons": False,
            "model_ids": {name: {block: rec["model_id"] for block, rec in sorted(pred.items())}
                          for name, pred in sorted(methods.items())},
            "units": units, "informative_units": sum(bool(u["positive_negative_pairs"]) for u in units),
            "informative_phages": len(per_phage["tied"]),
            "informative_phage_groups": len(per_group["tied"]),
            "per_phage": per_phage, "per_group": per_group, "group_macro_concordance": macro}
