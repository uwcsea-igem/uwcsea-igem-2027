"""Command-line entry point for evidence evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
from .benchmark import evaluate, load_cases, load_predictions


def main() -> int:
    parser = argparse.ArgumentParser(prog="tipcraft")
    sub = parser.add_subparsers(dest="command", required=True)
    benchmark = sub.add_parser("benchmark", help="Evaluate a published-observation dataset")
    benchmark.add_argument("--cases", type=Path, required=True)
    benchmark.add_argument("--predictions", type=Path)
    benchmark.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        scores, methods = load_predictions(args.predictions) if args.predictions else (None, {})
        result = evaluate(load_cases(args.cases), scores)
        result["case_sha256"] = hashlib.sha256(args.cases.read_bytes()).hexdigest()
        if args.predictions:
            result["prediction_sha256"] = hashlib.sha256(args.predictions.read_bytes()).hexdigest()
            result["prediction_methods"] = methods
        payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.output:
            if args.output.resolve() in {p.resolve() for p in (args.cases, args.predictions) if p}:
                raise ValueError("Output cannot overwrite input data")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(payload, encoding="utf-8")
        print(payload, end="")
    except (OSError, ValueError) as error:
        parser.exit(2, f"tipcraft: {error}\n")
    return 0
