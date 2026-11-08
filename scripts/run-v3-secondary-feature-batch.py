"""Execute the 17 remaining frozen v3 feature jobs after initial validation.

No search retry, scientific parameter selection, model fitting or scoring.
Already-running tasks finish after a failure; queued tasks stay unstarted.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import threading

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def checked(binding):
    path = ROOT / binding["path"]
    if sha(path) != binding["sha256"]:
        raise ValueError("Frozen artifact differs: " + str(path))
    return path


def now():
    return datetime.now(timezone.utc).isoformat()


def validate_report(report, job, preplan, plan):
    if (report["status"] != "completed_secondary_axis_features"
            or report["block_id"] != job["block_id"] or report["axis"] != job["axis"]
            or report["panel"] != plan["panel"]
            or report["helper_sha256"] != plan["producer"]["sha256"]
            or report["policy_sha256"] != plan["feature_policy"]["sha256"]
            or report["mapping_control_passed"] is not True or report["ready_for_fitting"] is not True
            or report["control_contract"] != "symmetric_mapping_outcomes_and_vector_equality_v3"
            or report["training_projection_mode"] != "ordinary_external_query_database_v3"):
        raise ValueError("Completed exact v3 block/axis with passing mapping control required")
    if len(report["commands"]) != 9 or len(preplan["commands"]) != 9:
        raise ValueError("Nine actual commands required")
    for actual, planned in zip(report["commands"], preplan["commands"]):
        if (actual["command"] != planned["argv"] or actual.get("stdin") != planned.get("stdin")
                or type(actual["exit_code"]) is not int or actual["exit_code"] != 0):
            raise ValueError("Actual command or observed exit differs from preplan")
    normalized = {Path(b["path"]).name: b["sha256"] for b in report["normalized_input_bindings"]}
    if (len(report["normalized_input_bindings"]) != 3
            or normalized != job["normalized_sha256"] or normalized != preplan["normalized_sha256"]):
        raise ValueError("Normalized inputs differ from pre-execution plan")


def run(plan_path, expected_sha):
    plan_binding = {"path": str(plan_path), "sha256": expected_sha}
    plan = json.loads(checked(plan_binding).read_text())
    if (plan["job_count"] != 18 or len(plan["jobs"]) != 18
            or len({j["id"] for j in plan["jobs"]}) != 18
            or plan["execution_stages"][1]["maximum_workers"] != 2):
        raise ValueError("Exact frozen 18-job/two-worker plan required")
    for key in ["panel", "producer", "feature_policy"]:
        checked(plan[key])
    preplans = {j["id"]: json.loads(checked(j["pre_execution_plan"]).read_text()) for j in plan["jobs"]}
    directory = checked(plan["jobs"][0]["pre_execution_plan"]).parent.parent
    ledger_path = directory / "execution-ledger.json"
    ledger = json.loads(ledger_path.read_text())
    if (ledger["plan"] != plan_binding or ledger["status"] != "completed_initial_validation"
            or ledger["remaining_tasks_started"] is not False
            or [j["id"] for j in ledger["jobs"]] != [j["id"] for j in plan["jobs"]]
            or any(j["status"] != "not_started" for j in ledger["jobs"][1:])):
        raise ValueError("Completed initial validation and entirely unstarted remainder required")
    first = ledger["jobs"][0]
    if type(first["exit_code"]) is not int or first["exit_code"] != 0:
        raise ValueError("Initial validation must have an observed zero exit")
    initial_report_path = checked(first["report"])
    if initial_report_path.resolve() != (ROOT / plan["jobs"][0]["output"]).resolve():
        raise ValueError("Initial report destination differs from frozen plan")
    validate_report(json.loads(initial_report_path.read_text()), plan["jobs"][0],
                    preplans[first["id"]], plan)
    snapshot = directory / "initial-validation-ledger.json"
    with snapshot.open("x") as stream:
        stream.write(ledger_path.read_text())
    ledger.update({"status": "running_remaining_v3_features", "remaining_tasks_started": True,
                   "remaining_stage_started_utc": now(), "maximum_workers": 2,
                   "initial_validation_snapshot": {"path": str(snapshot), "sha256": sha(snapshot)},
                   "remaining_runner": {"path": str(Path(__file__)), "sha256": sha(__file__)},
                   "execution_scope": "initial_validation_plus_remaining_fixed_v3_jobs"})
    lock, stop = threading.Lock(), threading.Event()
    records = {r["id"]: r for r in ledger["jobs"]}
    def save():
        ledger_path.write_text(json.dumps(ledger, indent=2, allow_nan=False) + "\n")
    save()
    def execute(job):
        record = records[job["id"]]
        with lock:
            if stop.is_set():
                return
            record.update({"status": "starting", "argv": job["argv"], "started_utc": now(), "log": job["log"]})
            save()
        try:
            checked(plan_binding)
            for key in ["producer", "feature_policy"]:
                checked(plan[key])
            checked(job["pre_execution_plan"])
            if (ROOT / job["output"]).exists():
                raise ValueError("Do not reuse an existing task report")
            with (ROOT / job["log"]).open("x") as log:
                with lock:
                    if stop.is_set():
                        record["status"] = "not_started"
                        save()
                        return
                    process = subprocess.Popen(job["argv"], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                    record.update({"status": "running", "pid": process.pid})
                    save()
                    print("START_V3", job["id"], "pid", process.pid, flush=True)
                code = process.wait()
            with lock:
                record.update({"exit_code": code, "finished_utc": now(), "log_sha256": sha(ROOT / job["log"])})
            report_path = ROOT / job["output"]
            if report_path.exists():
                report = json.loads(report_path.read_text())
                with lock:
                    record.update({"report": {"path": job["output"], "sha256": sha(report_path)}, "status": report["status"]})
            else:
                raise ValueError("Task ended without a report")
            if code != 0:
                raise ValueError("Task returned nonzero exit")
            validate_report(report, job, preplans[job["id"]], plan)
            for binding in report["normalized_input_bindings"]:
                checked(binding)
            with lock:
                record.update({"commands_match_preplan": True, "normalized_inputs_match_preplan": True,
                               "ready_for_fitting": True, "features": len(report["feature_names"])})
                save()
                print("DONE_V3", job["id"], "features", record["features"], flush=True)
        except Exception as error:
            with lock:
                stop.set()
                record.update({"admission_failed": True, "exception_type": type(error).__name__,
                               "error": str(error), "finished_utc": now()})
                if record["status"] in {"starting", "running"}:
                    record["status"] = "failed_before_completed_report"
                save()
                print("HOLD_V3", job["id"], type(error).__name__, str(error), flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(execute, plan["jobs"][1:]))
    try:
        checked(plan_binding)
        for key in ["producer", "feature_policy"]:
            checked(plan[key])
        if sha(__file__) != ledger["remaining_runner"]["sha256"]:
            raise ValueError("Remaining runner source changed during execution")
    except (OSError, ValueError) as error:
        stop.set()
        ledger["final_integrity_error"] = str(error)
    complete = not stop.is_set() and all(r.get("ready_for_fitting") is True for r in ledger["jobs"])
    ledger.update({"status": "completed_v3_feature_batch" if complete else "HOLD_feature_execution_failure",
                   "finished_utc": now(), "biological_model_fitted_or_scored": False,
                   "complete_fitting_experiment_manifest_still_required": True})
    save()
    print(ledger["status"], flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--expected-plan-sha256", required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.plan, args.expected_plan_sha256))
