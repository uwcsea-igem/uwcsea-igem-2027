"""Host-relatedness grouping cannot hide missing or duplicate comparisons."""
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "host_relatedness", Path(__file__).resolve().parents[3] / "scripts/audit-receptor-host-relatedness.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
SOURCE = Path(__file__).resolve().parents[3] / "scripts/audit-genophi-host-overlap.py"
SPEC2 = importlib.util.spec_from_file_location("host_near_rule", SOURCE)
RULE = importlib.util.module_from_spec(SPEC2)
SPEC2.loader.exec_module(RULE)
POLICY = {"ani_percent": 99.5, "minimum_af_in_each_genome_percent": 90,
          "ani_boundary_buffer_percent": 0.05, "af_boundary_buffer_percent": 0.5}


def row(a, b, ani=100):
    return {"Ref_file": a, "Query_file": b, "ANI": str(ani),
            "Align_fraction_ref": "99", "Align_fraction_query": "99"}


class HostRelatednessTests(unittest.TestCase):
    def test_complete_graph_preserves_transitive_reference_bridges(self):
        files = [{"id": x, "path": x} for x in ("a", "b", "c")]
        rows = [row(a, b, 100 if a == b else 99.6 if {a, b} != {"a", "c"} else 99)
                for a in ("a", "b", "c") for b in ("a", "b", "c")]
        result = AUDIT.summarize(rows, files, POLICY, RULE.classify_pair)
        self.assertEqual(result["reference_components"], [["a", "b", "c"]])

    def test_missing_direction_remains_hold_without_a_zero_or_group(self):
        files = [{"id": x, "path": x} for x in ("a", "b")]
        result = AUDIT.summarize([row("a", "a"), row("b", "b"), row("a", "b")], files, POLICY, RULE.classify_pair)
        self.assertEqual(result["missing_directed_pairs"], [["b", "a"]])
        self.assertIsNone(result["reference_components"])

    def test_duplicate_paths_ids_and_rows_are_rejected(self):
        for files in ([{"id": "a", "path": "same"}, {"id": "b", "path": "same"}],
                      [{"id": "a", "path": "x"}, {"id": "a", "path": "y"}]):
            with self.subTest(files=files), self.assertRaises(ValueError):
                AUDIT.summarize([], files, POLICY, RULE.classify_pair)
        with self.assertRaises(ValueError):
            AUDIT.summarize([row("a", "a"), row("a", "a")], [{"id": "a", "path": "a"}], POLICY, RULE.classify_pair)


if __name__ == "__main__":
    unittest.main()
