import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("phage_overlap", Path(__file__).resolve().parents[3] / "scripts/audit-genophi-phage-overlap.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
GENOMES = [{"id": "t0", "source_id": "T7", "role": "target_reference", "retained_bases": 1500},
           {"id": "a", "source_id": "first", "role": "external_collection", "retained_bases": 1500},
           {"id": "b", "source_id": "second", "role": "external_collection", "retained_bases": 1500}]


class PhageOverlapTests(unittest.TestCase):
    def test_completed_no_hit_differs_from_truncated_search(self):
        text = "# BLASTN\n# Query: a\n# 1 hits found\na\ta\n# Query: b\n# 0 hits found\n# BLAST processed 2 queries\n"
        rows, completion = AUDIT.completed_blast(text, {"a", "b"}, 2)
        self.assertEqual(rows, [["a", "a"]])
        self.assertEqual(completion["reported_rows_by_query"], {"a": 1, "b": 0})
        batched = text.replace("# Query: b", "# BLAST processed 1 queries\n# Query: b").replace("# BLAST processed 2 queries", "# BLAST processed 1 queries")
        self.assertEqual(AUDIT.completed_blast(batched, {"a", "b"}, 2)[1]["completed_batch_query_counts"], [1, 1])
        for malformed in [text.split("# Query: b")[0], text.replace("# 1 hits found", "# 2 hits found"), text + "# Query: x\n"]:
            with self.subTest(malformed=malformed), self.assertRaises(ValueError):
                AUDIT.completed_blast(malformed, {"a", "b"}, 2)
    def test_target_group_excludes_transitive_related_source(self):
        pairs = [{"query": a, "subject": b, "similarity_percent_unrounded": value}
                 for a, b, value in [("t0", "a", 71), ("a", "b", 72), ("t0", "b", 30)]]
        result = AUDIT.component_exclusions(pairs, GENOMES, 70, 0.01)
        self.assertEqual(result["excluded_source_ids"], ["first", "second"])

    def test_boundary_is_excluded_and_incomplete_or_duplicate_matrix_rejected(self):
        pair = {"query": "t0", "subject": "a", "similarity_percent_unrounded": 69.995}
        result = AUDIT.component_exclusions([pair], GENOMES[:2], 70, 0.01)
        self.assertEqual(result["excluded_source_ids"], ["first"])
        self.assertEqual(len(result["boundary_edges"]), 1)
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            AUDIT.component_exclusions([pair], GENOMES, 70, 0.01)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            AUDIT.component_exclusions([pair, pair], GENOMES[:2], 70, 0.01)

    def test_gp17_identity_and_query_coverage_are_separate_requirements(self):
        policy = {"gp17_near_match": {"identity_percent": 90, "minimum_query_coverage_percent": 80},
                  "broader_gp17_family_diagnostic": {"identity_percent": 40, "minimum_query_coverage_percent": 80}}
        row = ["t0_gp17", "a", "92", "1", "400", "500", "10", "1209", "1", "500", "1e-30"]
        control = ["t0_gp17", "t0", "100", "1", "500", "500", "1", "1500", "1", "600", "0"]
        result = AUDIT.gp17_matches([row, control], GENOMES, policy, {"t0_gp17": 500})
        self.assertEqual(result["near_match_excluded_source_ids"], ["first"])
        short = row.copy(); short[4] = "399"
        self.assertEqual(AUDIT.gp17_matches([short, control], GENOMES, policy, {"t0_gp17": 500})["near_match_excluded_source_ids"], [])
        distant = row.copy(); distant[2] = "45"
        result = AUDIT.gp17_matches([distant, control], GENOMES, policy, {"t0_gp17": 500})
        self.assertEqual(len(result["broader_family_matches"]), 1)
        self.assertEqual(result["near_match_excluded_source_ids"], [])
        with self.assertRaisesRegex(ValueError, "self controls"):
            AUDIT.gp17_matches([], GENOMES, policy, {"t0_gp17": 500})
        bad = row.copy(); bad[8] = "0"
        with self.assertRaisesRegex(ValueError, "frame"):
            AUDIT.gp17_matches([bad, control], GENOMES, policy, {"t0_gp17": 500})
