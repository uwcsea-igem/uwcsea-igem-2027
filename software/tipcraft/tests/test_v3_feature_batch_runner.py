"""Reject partial, mismatched or unexecuted evidence before batch admission."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("v3_feature_batch", ROOT / "scripts/run-v3-secondary-feature-batch.py")
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


def fixture():
    panel = {"path": "panel.json", "sha256": "panel"}
    plan = {"panel": panel, "producer": {"sha256": "producer"}, "feature_policy": {"sha256": "policy"}}
    commands = [{"argv": ["tool", str(i)], "stdin": None} for i in range(9)]
    normalized = {"training.faa": "train", "queries.faa": "query", "protein-map.json": "map"}
    job = {"block_id": "host1-phage0", "axis": "host", "normalized_sha256": normalized}
    preplan = {"commands": commands, "normalized_sha256": normalized}
    report = {"status": "completed_secondary_axis_features", "block_id": job["block_id"], "axis": "host",
              "panel": panel, "helper_sha256": "producer", "policy_sha256": "policy",
              "mapping_control_passed": True, "ready_for_fitting": True,
              "control_contract": "symmetric_mapping_outcomes_and_vector_equality_v3",
              "training_projection_mode": "ordinary_external_query_database_v3",
              "commands": [{"command": c["argv"], "stdin": None, "exit_code": 0} for c in commands],
              "normalized_input_bindings": [{"path": "raw/" + k, "sha256": v} for k, v in normalized.items()]}
    return report, job, preplan, plan


class V3FeatureBatchTests(unittest.TestCase):
    def test_complete_actual_command_receipt_is_admitted(self):
        RUN.validate_report(*fixture())

    def test_incomplete_false_exit_or_changed_scope_evidence_is_rejected(self):
        mutations = [lambda r: r["commands"].pop(), lambda r: r["commands"][0].update(exit_code=False),
                     lambda r: r["commands"][0].update(command=["changed"]),
                     lambda r: r.update(mapping_control_passed=False), lambda r: r.update(ready_for_fitting=False),
                     lambda r: r.update(axis="phage"), lambda r: r.update(helper_sha256="changed"),
                     lambda r: r.update(control_contract="v2"),
                     lambda r: r["normalized_input_bindings"][0].update(sha256="changed"),
                     lambda r: r["normalized_input_bindings"].append(r["normalized_input_bindings"][0])]
        for mutate in mutations:
            args = deepcopy(fixture())
            mutate(args[0])
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                RUN.validate_report(*args)

    def test_running_initial_validation_never_starts_remaining_tools(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "plans").mkdir()
            receipt = directory / "plans/receipt.json"
            receipt.write_text("{}")
            bound = {"path": str(receipt), "sha256": RUN.sha(receipt)}
            plan = {"job_count": 18, "jobs": [{"id": str(i), "pre_execution_plan": bound} for i in range(18)],
                    "execution_stages": [{}, {"maximum_workers": 2}],
                    "panel": bound, "producer": bound, "feature_policy": bound}
            path = directory / "plan.json"
            path.write_text(json.dumps(plan))
            ledger = directory / "execution-ledger.json"
            original = json.dumps({"plan": {"path": str(path), "sha256": RUN.sha(path)},
                                   "status": "running_initial_validation"})
            ledger.write_text(original)
            with patch.object(RUN.subprocess, "Popen") as popen, self.assertRaisesRegex(ValueError, "initial validation"):
                RUN.run(path, RUN.sha(path))
            popen.assert_not_called()
            self.assertEqual(ledger.read_text(), original)
            self.assertFalse((directory / "initial-validation-ledger.json").exists())


if __name__ == "__main__":
    unittest.main()
