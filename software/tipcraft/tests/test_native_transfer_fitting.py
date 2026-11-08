"""Final-source/native assembly checks on synthetic metadata and vectors."""
from copy import deepcopy
from dataclasses import asdict
import importlib.util
from pathlib import Path
import unittest

from tipcraft.learning import Sample

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("native_transfer_fit", ROOT / "scripts/fit-native-transfer.py")
NATIVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NATIVE)


def fixture():
    samples = [Sample(f"s{i}{j}", f"sh{i}", f"sp{j}", f"shg{i}", f"spg{j}", "synthetic-source")
               for i in range(2) for j in range(2)]
    panel = {"status": "frozen_secondary_reference_panel", "membership_generated": True,
             "annotation_join_complete": True, "ready_for_feature_preparation": True,
             "samples": [asdict(s) for s in samples], "labels": {s.row_id: i % 2 for i, s in enumerate(samples)}}
    design = {"native_candidates": {"hosts": NATIVE.HOSTS, "phages": NATIVE.PHAGES}}
    host_groups = {"partition": {"source_groups": {
        "native-bl": {"members": ["original:BL21"]},
        "native-k": {"members": ["original:BW25113", "original:MG1655"]}}}}
    phage_groups = {"grouping": {"container_groups": {"reference:T7": "native-p", "reference:T3": "native-p"}}}
    reports = []
    for axis, prefix, queries in [("host", "sh", NATIVE.HOSTS), ("phage", "sp", NATIVE.PHAGES)]:
        training = [prefix + str(i) for i in range(2)]
        vectors = {n: [i % 2, 1-i % 2] for i, n in enumerate(training + queries)}
        reports.append({"status": "completed_native_transfer_axis_features", "axis": axis,
                        "training_genomes": training, "query_genomes": queries, "vectors": vectors,
                        "feature_names": ["f0", "f1"], "mapping_control_passed": True})
    return panel, *reports, design, host_groups, phage_groups


class NativeTransferAssemblyTests(unittest.TestCase):
    def test_all_source_rows_and_exact_six_native_pairs_without_native_labels(self):
        args = fixture()
        data = NATIVE.assemble_native(*args)
        self.assertEqual(len(data["source"]), 4)
        self.assertEqual(len(data["target"]), 6)
        self.assertEqual({(s.phage_id, s.host_id) for s in data["target"]},
                         {(p, h) for p in NATIVE.PHAGES for h in NATIVE.HOSTS})
        self.assertFalse(data["partition"]["overlaps"]["host_group"])
        self.assertFalse(data["partition"]["overlaps"]["phage_group"])
        self.assertTrue(data["partition"]["publication_holdout"])
        self.assertNotIn("y_query", data)
        args[0]["labels"]["native:T3:BL21"] = "unused native annotation"
        self.assertEqual(data, NATIVE.assemble_native(*args))

    def test_changed_candidates_groups_or_source_feature_rosters_fail(self):
        for mutate in [lambda a: a[0].update(status="HOLD"),
                       lambda a: a[3].update(native_candidates={"hosts": NATIVE.HOSTS, "phages": ["T7"]}),
                       lambda a: a[1].update(query_genomes=list(reversed(NATIVE.HOSTS))),
                       lambda a: a[1]["training_genomes"].append("BL21"),
                       lambda a: a[2]["vectors"].update(extra=[0, 1]),
                       lambda a: a[1].update(mapping_control_passed=False),
                       lambda a: a[5]["grouping"]["container_groups"].update({"reference:T7": "spg0"})]:
            args = deepcopy(fixture())
            mutate(args)
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                NATIVE.assemble_native(*args)


if __name__ == "__main__":
    unittest.main()
