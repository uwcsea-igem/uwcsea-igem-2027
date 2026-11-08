import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


AUDIT = module("scripts/audit-official-variant-gp17.py", "official_variant_gp17_tests")
LOGIC = module("scripts/audit-genophi-phage-overlap.py", "official_variant_gp17_frozen_logic_tests")
POLICY = json.loads((ROOT / "software/tipcraft/benchmarks/receptor-resource/gp17-overlap-policy.json").read_text())


def hit(subject, identity=90, start=1, end=800, query="t000_gp17"):
    return [query, subject, str(identity), str(start), str(end), "1000", "1", str(3*(end-start+1)), "1", "100", "1e-20"]


def controls():
    return [hit("t000", 100, end=1000), hit("t001", 100, end=1000, query="t001_gp17")]


def screen_records():
    return [{"id": "t000", "source_id": "T7", "container_key": "reference:T7", "role": "target_reference", "retained_bases": 4000},
            {"id": "t001", "source_id": "T3", "container_key": "reference:T3", "role": "target_reference", "retained_bases": 4000},
            {"id": "source_a", "source_id": "Bas06", "container_key": "source:Bas06", "role": "external_collection", "retained_bases": 4000},
            {"id": "source_b", "source_id": "Bas06", "container_key": "source:Bas06", "role": "external_collection", "retained_bases": 4000},
            {"id": "official_a", "source_id": "Bas06", "container_key": "official:Bas06:ABC.1", "role": "official_reference_version", "retained_bases": 4000}]


def corpus_fixture():
    old = [{"id": f"c{i}", "source_id": f"S{i}", "container_key": f"source:S{i}",
            "role": "external_collection", "retained_bases": 3} for i in range(261)]
    old += [dict(old[0], id=f"extra{i}") for i in range(7)]
    old += [dict(row, retained_bases=3) for row in screen_records()[:2]]
    official = [{"id": f"official_{name}", "source_id": name, "source_record_id": accession,
                 "role": "official_reference_version", "retained_bases": 3,
                 "container_key": f"official:{name}:{accession}"}
                for name, accession in [("Bas06", "A.1"), ("Bas21", "B.1"), ("Bas68", "C.1")]]
    sequences = [(r["id"], "ACG") for r in old + official]
    return {"records": old + official}, sequences, {"records": old}, sequences[:270]


class OfficialVariantGp17Tests(unittest.TestCase):
    def screen(self, rows, records=None):
        return AUDIT.screen(rows, records or screen_records(), POLICY,
                            {"t000_gp17": 1000, "t001_gp17": 1000}, LOGIC)

    def test_primary_boundary_and_distinct_same_name_namespaces(self):
        result = self.screen(controls() + [hit("official_a")])
        self.assertEqual(result["direct_near_hit_container_keys"],
                         ["official:Bas06:ABC.1", "reference:T3", "reference:T7"])
        self.assertNotIn("source:Bas06", result["direct_near_hit_container_keys"])
        self.assertEqual(result["target_self_controls_passed"], ["t000_gp17", "t001_gp17"])

    def test_broad_boundary_does_not_become_primary(self):
        result = self.screen(controls() + [hit("source_a", identity=40), hit("official_a", identity=89.999)])
        self.assertIn("source:Bas06", result["broader_diagnostic_container_keys"])
        self.assertIn("official:Bas06:ABC.1", result["broader_diagnostic_container_keys"])
        self.assertNotIn("source:Bas06", result["direct_near_hit_container_keys"])
        self.assertNotIn("official:Bas06:ABC.1", result["direct_near_hit_container_keys"])

    def test_does_not_sum_separate_short_hsps(self):
        result = self.screen(controls() + [hit("source_a", end=400), hit("source_a", start=401, end=800)])
        self.assertEqual(result["direct_near_hit_container_keys"], ["reference:T3", "reference:T7"])
        self.assertNotIn("source:Bas06", result["broader_diagnostic_container_keys"])

    def test_multirecord_hit_flags_container_retains_record_id(self):
        result = self.screen(controls() + [hit("source_b")])
        row = next(r for r in result["containers"] if r["container_key"] == "source:Bas06")
        self.assertEqual(row["nucleotide_record_ids"], ["source_a", "source_b"])
        self.assertTrue(row["direct_near_hit"])
        self.assertEqual(result["near_matches"][-1]["nucleotide_record_id"], "source_b")

    def test_missing_control_unknown_subject_and_invalid_statistics_fail(self):
        for rows in [controls()[:1], controls()+[hit("unknown")], controls()+[hit("source_a", identity="nan")]]:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                self.screen(rows)

    def test_container_cannot_merge_roles(self):
        records = screen_records()
        records[-1]["container_key"] = "source:Bas06"
        with self.assertRaises(ValueError):
            self.screen(controls(), records)

    def test_old_classification_changes_are_explicit(self):
        result = self.screen(controls() + [hit("official_a")])
        old = {"direct_near_hit_excluded_primary_ids": ["Bas06"],
               "gp17": {"broader_family_matches": [{"source_id": "Bas06"}]}}
        changes = AUDIT.old_changes(old, result)
        self.assertEqual(changes["compared_source_containers"], 1)
        self.assertEqual(changes["changed_source_classifications"][0]["container_key"], "source:Bas06")
        self.assertFalse(changes["changed_source_classifications"][0]["augmented_database"]["direct_near_hit"])

    def test_corpus_retains_all_prior_records_and_three_official_nodes(self):
        self.assertEqual(len(AUDIT.corpus(*corpus_fixture())), 273)

    def test_corpus_rejects_reordering_changed_sequence_and_role_alias(self):
        for mutation in ("order", "sequence", "role", "duplicate"):
            args = copy.deepcopy(corpus_fixture())
            if mutation == "order":
                args[0]["records"][0], args[0]["records"][1] = args[0]["records"][1], args[0]["records"][0]
            elif mutation == "sequence":
                args[1][0] = (args[1][0][0], "TGC")
            elif mutation == "role":
                args[0]["records"][-1]["role"] = "external_collection"
            else:
                args[0]["records"][-1]["id"] = args[0]["records"][-2]["id"]
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                AUDIT.corpus(*args)


if __name__ == "__main__":
    unittest.main()
