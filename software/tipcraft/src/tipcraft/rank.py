"""Generic classification evaluation; no biological design ranking is implemented."""

def classification_metrics(labels: list[int], scores: list[float]) -> dict:
    if not labels or len(labels) != len(scores):
        raise ValueError("Need equally sized nonempty labels and scores")
    if any(type(y) is not int or y not in (0, 1) for y in labels):
        raise ValueError("Labels must be integer 0 or 1")
    if any(not 0 <= p <= 1 for p in scores):
        raise ValueError("Scores must be finite probabilities in [0, 1]")
    predictions = [int(p >= 0.5) for p in scores]
    recalls = [sum(y == c and p == c for y, p in zip(labels, predictions)) / labels.count(c)
               for c in (0, 1) if c in labels]
    return {"n": len(labels),
            "accuracy": sum(y == p for y, p in zip(labels, predictions)) / len(labels),
            "balanced_accuracy": sum(recalls) / 2 if len(recalls) == 2 else None,
            "threshold": 0.5}
