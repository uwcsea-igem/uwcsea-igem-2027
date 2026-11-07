"""Command-line entry point for evidence evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
from .benchmark import evaluate, load_cases, load_predictions
from .wish import run_scores, sha256
from .retrospective import evaluate_scores, leakage_audit


def main() -> int:
    parser = argparse.ArgumentParser(prog="tipcraft")
    sub = parser.add_subparsers(dest="command", required=True)
    benchmark = sub.add_parser("benchmark", help="Evaluate a published-observation dataset")
    benchmark.add_argument("--cases", type=Path, required=True)
    benchmark.add_argument("--predictions", type=Path)
    benchmark.add_argument("--output", type=Path)
    score = sub.add_parser("wish-score", help="Run external WIsH without phenotype labels")
    score.add_argument("--manifest", type=Path, required=True)
    score.add_argument("--genomes", type=Path, required=True)
    score.add_argument("--executable", type=Path, required=True)
    score.add_argument("--output-dir", type=Path, required=True)
    score.add_argument("--anonymize", action="store_true")
    report = sub.add_parser("wish-report", help="Evaluate already frozen WIsH scores")
    report.add_argument("--cases", type=Path, required=True)
    report.add_argument("--scores", type=Path, required=True)
    report.add_argument("--anonymous-scores", type=Path, required=True)
    report.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "wish-score":
            result = run_scores(args.manifest, args.genomes, args.executable,
                                args.output_dir, args.anonymize)
            print(json.dumps(result, indent=2, allow_nan=False))
            return 0
        if args.command == "wish-report":
            inputs = (args.cases, args.scores, args.anonymous_scores)
            if args.output.resolve() in {p.resolve() for p in inputs}:
                raise ValueError("Output cannot overwrite input data")
            cases = load_cases(args.cases)
            inference = json.loads(args.scores.read_text())
            anonymous = json.loads(args.anonymous_scores.read_text())
            result = evaluate_scores(cases, inference)
            result["leakage_audit"] = leakage_audit(cases, inference, anonymous)
            result["input_sha256"] = {"cases": sha256(args.cases), "scores": sha256(args.scores),
                                      "anonymous_scores": sha256(args.anonymous_scores)}
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
            print(json.dumps(result, indent=2, allow_nan=False))
            return 0
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
