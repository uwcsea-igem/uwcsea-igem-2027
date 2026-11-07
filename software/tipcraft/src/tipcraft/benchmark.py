"""Evaluate published binary observations against baselines or frozen predictions."""
import csv
import math
from pathlib import Path
from .assembly import AssemblyEvidence, summarize
from .rank import classification_metrics


def load_cases(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Case file is empty")
    seen = set()
    for row in rows:
        for field in ("case_id", "phage", "host", "label", "source", "locator", "interpretation"):
            if not (row.get(field) or "").strip():
                raise ValueError(f"Missing case field: {field}")
        if row["case_id"] in seen or row["label"] not in ("0", "1"):
            raise ValueError("Duplicate case ID or invalid label")
        seen.add(row["case_id"])
    return rows


def load_predictions(path: Path) -> tuple[dict[str, float], dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Prediction file is empty")
    scores, methods = {}, {}
    for row in rows:
        case_id = (row.get("case_id") or "").strip()
        method = (row.get("method") or "").strip()
        if not case_id or not method or case_id in scores:
            raise ValueError("Predictions require unique case_id and nonempty method")
        score = float(row.get("score", ""))
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Prediction score must be finite and in [0, 1]")
        scores[case_id], methods[case_id] = score, method
    return scores, methods


def evaluate(cases: list[dict], scores: dict[str, float] | None = None) -> dict:
    labels = [int(row["label"]) for row in cases]
    # A predeclared constant comparator, not a learned or structural predictor.
    baseline = classification_metrics(labels, [1.0] * len(cases))
    result = {
        "schema_version": 1,
        "dataset": "retro-t7-t3-v0.1",
        "label_type": "qualitative published plaque observation",
        "cases": len(cases),
        "baseline": {"method": "always-positive", **baseline},
        "assembly": summarize([AssemblyEvidence() for _ in cases]),
        "model": {"status": "not_evaluated", "prediction_coverage": 0.0},
        "gate": "HOLD: no independent model predictions; baseline is not validation",
        "limitations": [
            "Small correlated cases from one publication; no generalization claim",
            "Negative labels mean no plaques detected under the reported conditions",
            "No quantitative EOP values or direct assembly assays included",
            "Prediction independence and leakage require human review",
        ],
    }
    if scores is None:
        return result
    known = {row["case_id"] for row in cases}
    if set(scores) - known:
        raise ValueError("Unknown prediction case IDs")
    selected = [row for row in cases if row["case_id"] in scores]
    if not selected:
        raise ValueError("No predictions match the cases")
    metrics = classification_metrics([int(row["label"]) for row in selected],
                                     [scores[row["case_id"]] for row in selected])
    result["model"] = {"status": "evaluated_unreviewed", **metrics,
                       "prediction_coverage": len(selected) / len(cases)}
    result["gate"] = "HOLD: predictions evaluated; independence and gate criteria need review"
    return result
