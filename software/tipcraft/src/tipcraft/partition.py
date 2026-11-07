"""Deterministic relatedness-group blocks without reading phenotype labels."""
from collections import Counter
import hashlib
import json

from .learning import _samples, _text, validate_partition


def assign_group_folds(group_ids, *, axis: str, folds: int, seed: int) -> dict[str, int]:
    """Hash-order unique groups, then assign cyclically; never optimize labels."""
    if axis not in {"host", "phage"}:
        raise ValueError("axis must be host or phage")
    if isinstance(folds, bool) or not isinstance(folds, int) or folds < 2:
        raise ValueError("At least two folds are required")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    values = list(group_ids)
    if not values or any(not _text(value) for value in values):
        raise ValueError("Valid nonempty group IDs required")
    groups = sorted(set(values))
    if len(groups) < folds:
        raise ValueError("Fewer relatedness groups than requested folds; no fold-count fallback")
    def key(group):
        data = json.dumps([seed, axis, group], ensure_ascii=True, separators=(",", ":")).encode()
        return hashlib.sha256(data).hexdigest(), group
    ordered = sorted(groups, key=key)
    return dict(sorted((group, i % folds) for i, group in enumerate(ordered)))


def _fold_count(assignments, expected, axis):
    if set(assignments) != expected:
        raise ValueError(f"{axis} fold mapping must cover exactly the retained groups")
    values = list(assignments.values())
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in values):
        raise ValueError("Fold indices must be nonnegative integers")
    count = len(set(values))
    if count < 2 or set(values) != set(range(count)):
        raise ValueError("Fold indices must cover a contiguous range of at least two folds")
    return count


def build_blocked_folds(samples, host_folds: dict[str, int], phage_folds: dict[str, int]) -> list[dict]:
    """Cross all host/phage blocks; mixed quadrants are explicitly embargoed.

    Every retained pair is tested once. In each block training uses neither
    held-out host groups nor held-out phage groups. Labels are not an input.
    """
    rows = sorted(_samples(samples), key=lambda sample: sample.row_id)
    nh = _fold_count(host_folds, {s.host_group for s in rows}, "Host")
    np = _fold_count(phage_folds, {s.phage_group for s in rows}, "Phage")
    result, test_counts = [], Counter()
    all_ids = {s.row_id for s in rows}
    for host_fold in range(nh):
        for phage_fold in range(np):
            train, test, embargo = [], [], []
            for sample in rows:
                h = host_folds[sample.host_group] == host_fold
                p = phage_folds[sample.phage_group] == phage_fold
                (test if h and p else train if not h and not p else embargo).append(sample)
            identifier = f"host{host_fold}-phage{phage_fold}"
            if not train or not test:
                raise ValueError(f"Empty training/test block {identifier}; no random fallback")
            validation = validate_partition(train, test, regime="both_unseen")
            roles = [{s.row_id for s in subset} for subset in (train, test, embargo)]
            if set.union(*roles) != all_ids or any(roles[i] & roles[j] for i in range(3) for j in range(i)):
                raise ValueError("Block roles do not partition the retained rows")
            test_counts.update(roles[1])
            result.append({"fold_id": identifier, "host_fold": host_fold, "phage_fold": phage_fold,
                           "train_ids": [s.row_id for s in train], "test_ids": [s.row_id for s in test],
                           "embargo_ids": [s.row_id for s in embargo], "validation": validation})
    if set(test_counts) != all_ids or any(value != 1 for value in test_counts.values()):
        raise ValueError("Every retained row must occur in exactly one outer test block")
    return result
