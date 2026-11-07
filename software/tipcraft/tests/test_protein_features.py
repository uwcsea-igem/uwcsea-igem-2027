import copy
import json
import unittest

from tipcraft.protein_features import best_hits, fit_vocabulary, project, read_families


def hit(target, family=None):
    return {"target": target, "family": family or target, "identity_float32": 1.0,
            "evalue_float32": 0.0, "raw_row": 1}


class ProteinFeatureTests(unittest.TestCase):
    def setUp(self):
        self.owners = {"a1": "g1", "a2": "g2", "b1": "g1", "b2": "g2", "rare": "g1"}
        self.families = {"a1": "a1", "a2": "a1", "b1": "b1", "b2": "b1", "rare": "rare"}
        self.hits = {p: hit(p, self.families[p]) for p in self.owners}

    def test_train_only_consolidation_heldout_occurrences_do_not_rescue_rare_family(self):
        vocabulary = fit_vocabulary(["g1", "g2", "g3"], self.owners, self.families, self.hits)
        self.assertEqual(vocabulary["family_features"], {"a1": "pc_0", "b1": "pc_0"})
        self.assertEqual(vocabulary["retained_family_patterns"], {"a1": [1, 1, 0], "b1": [1, 1, 0]})
        frozen = copy.deepcopy(vocabulary)
        features, vectors = project(vocabulary, ["heldout1", "heldout2", "zero"],
                                    {"q1": "heldout1", "q2": "heldout2", "q3": "heldout1"},
                                    {"q1": hit("a1"), "q2": hit("b1"), "q3": hit("rare")})
        self.assertEqual(features, ["pc_0"])
        self.assertEqual(vectors, {"heldout1": [1], "heldout2": [1], "zero": [0]})
        self.assertEqual(vocabulary, frozen)
        # Rare family is still excluded although a new heldout genome carries it;
        # a1/b1 remain consolidated although heldout presence differs.
        _, rare_only = project(vocabulary, ["new"], {"q": "new"}, {"q": hit("rare")})
        self.assertEqual(rare_only, {"new": [0]})

    def test_zero_hits_and_zero_feature_vocabulary_are_explicit(self):
        best, diagnostics = best_hits(iter([]), {"q": 3}, {"a": 3}, {"a": "a"})
        self.assertEqual(best, {})
        self.assertEqual(diagnostics["queries_without_hits"], 1)
        vocabulary = fit_vocabulary(["g1", "g2"], self.owners, self.families, {})
        self.assertEqual(vocabulary["training_proteins_without_hits"], 5)
        self.assertEqual(vocabulary["features"], [])
        self.assertEqual(project(vocabulary, ["x", "y"], {}, {}), ([], {"x": [], "y": []}))
        json.dumps(vocabulary, allow_nan=False)

    def test_streaming_hit_order_preserves_float32_raw_first_cross_family_ties(self):
        rows = iter(["q\tb\t40\t0.800000001\t1e-60\t0\t2\t3\t0\t2\t3\n",
                     "q\ta\t999\t0.800000002\t2e-60\t0\t2\t3\t0\t2\t3\n"])
        best, diagnostics = best_hits(rows, {"q": 3}, {"a": 3, "b": 3}, {"a": "a", "b": "b"})
        self.assertEqual(set(best["q"]), {"target", "family", "identity_float32", "evalue_float32", "raw_row"})
        self.assertEqual(best["q"]["target"], "b")
        self.assertEqual(best["q"]["raw_row"], 1)
        self.assertEqual(diagnostics["cross_family_best_hit_ties"],
                         [{"query": "q", "families": ["a", "b"], "chosen": best["q"]}])
        self.assertEqual(diagnostics["raw_rows"], 2)

    def test_best_hits_identity_then_evalue_not_bit_score(self):
        rows = ["q\ta\t999\t0.7\t1e-30\t0\t2\t3\t0\t2\t3",
                "q\tb\t10\t0.8\t1e-5\t0\t2\t3\t0\t2\t3",
                "q\tc\t1\t0.8\t1e-10\t0\t2\t3\t0\t2\t3"]
        best, diagnostics = best_hits(rows, {"q": 3}, {"a": 3, "b": 3, "c": 3}, {p: p for p in "abc"})
        self.assertEqual(best["q"]["target"], "c")
        self.assertEqual(diagnostics["cross_family_best_hit_ties"], [])

    def test_cluster_membership_complete_exact_and_training_only(self):
        self.assertEqual(read_families(iter(["a\ta\n", "a\tb\n"]), ["a", "b"]), {"a": "a", "b": "a"})
        for rows in [["a\ta"], ["a\ta", "a\ta"], ["a\ta", "heldout\tb"], ["a\tb", "b\ta"]]:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                read_families(iter(rows), ["a", "b"])
        with self.assertRaises(ValueError):
            read_families(iter([]), ["a", "a"])

    def test_malformed_coordinates_nonfinite_fields_unknown_targets_and_queries_fail(self):
        base = ["q", "a", "10", "0.8", "1e-5", "0", "2", "3", "0", "2", "3"]
        cases = [(0, "unknown"), (1, "heldout"), (2, "inf"), (3, "nan"), (3, "2"),
                 (4, "-1"), (5, "-1"), (6, "3"), (7, "4"), (8, "3"), (9, "-1"), (10, "4")]
        for index, value in cases:
            row = list(base)
            row[index] = value
            with self.subTest(index=index, value=value), self.assertRaises(ValueError):
                best_hits(iter(["\t".join(row)]), {"q": 3}, {"a": 3}, {"a": "a"})
        with self.assertRaises(ValueError):
            best_hits(iter(["q\ta"]), {"q": 3}, {"a": 3}, {"a": "a"})
        with self.assertRaises(ValueError):
            best_hits([], {"q": 0}, {"a": 3}, {"a": "a"})

    def test_fit_and_projection_validate_training_boundary_and_genome_membership(self):
        for owners, families, hits in [({"a1": "heldout"}, {"a1": "a1"}, {"a1": hit("a1")}),
                                       (self.owners, self.families, {"heldout_query": hit("a1")}),
                                       (self.owners, self.families, {"a1": hit("heldout_target", "a1")}),
                                       (self.owners, self.families, {"a1": hit("a1", "b1")})]:
            with self.assertRaises(ValueError):
                fit_vocabulary(["g1", "g2"], owners, families, hits)
        vocabulary = fit_vocabulary(["g1", "g2"], self.owners, self.families, self.hits)
        for genomes, owners, hits in [(["new"], {"q": "unlisted"}, {}),
                                       (["new", "new"], {}, {}),
                                       (["new"], {}, {"unknown": hit("a1")}),
                                       (["new"], {"q": "new"}, {"q": hit("heldout", "a1")})]:
            with self.assertRaises(ValueError):
                project(vocabulary, genomes, owners, hits)

    def test_feature_names_sort_lexically_and_vocabulary_drift_is_rejected(self):
        owners, families, hits = {}, {}, {}
        for i in range(12):
            family, member = f"f{i:02d}", f"m{i:02d}"
            owners[family], owners[member] = "anchor", f"g{i:02d}"
            families[family], families[member] = family, family
            hits[family], hits[member] = hit(family), hit(member, family)
        vocabulary = fit_vocabulary(["anchor"] + [f"g{i:02d}" for i in range(12)], owners, families, hits)
        self.assertEqual(vocabulary["features"], sorted(f"pc_{i}" for i in range(12)))
        self.assertLess(vocabulary["features"].index("pc_10"), vocabulary["features"].index("pc_2"))
        damaged = copy.deepcopy(vocabulary)
        damaged["family_features"]["f00"] = "pc_999"
        with self.assertRaises(ValueError):
            project(damaged, ["new"], {}, {})


if __name__ == "__main__":
    unittest.main()
