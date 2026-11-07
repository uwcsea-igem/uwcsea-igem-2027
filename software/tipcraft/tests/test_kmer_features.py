"""Test scientific boundaries of train-only global protein k-mer features."""
import unittest

from tipcraft.kmer_features import fit_vocabulary, project_proteomes


class KmerFeatureTests(unittest.TestCase):
    def test_never_join_protein_boundaries(self):
        model = fit_vocabulary({"A": {"p": "ACD", "q": "EFG"}}, minimum_protein_occurrences=1)
        self.assertEqual(model.features, ())
        self.assertEqual(project_proteomes(model, {"B": {"r": "ACDEF"}})["matrix"], [[]])

    def test_occurrence_counts_distinct_proteins_not_repetitions_or_genomes(self):
        repeated = fit_vocabulary({"A": {"p": "AAAAAAAAAA"}})
        two_proteins = fit_vocabulary({"A": {"p": "AAAAA", "q": "AAAAA"}})
        self.assertEqual(repeated.features, ())
        self.assertEqual(two_proteins.features, (("AAAAA",),))

    def test_heldout_only_and_heldout_rescued_singletons_never_enter_vocabulary(self):
        model = fit_vocabulary({"A": {"p": "AAAAA"}, "B": {"p": "CCCCC", "q": "CCCCC"}})
        before = model.features
        result = project_proteomes(model, {"test": {"p": "AAAAA", "q": "DDDDD"}})
        self.assertEqual(model.features, before)
        self.assertEqual(before, (("CCCCC",),))
        self.assertEqual(result["matrix"], [[0]])
        self.assertEqual(result["zero_vector_sample_ids"], ["test"])

    def test_identical_training_patterns_merge_and_project_union_of_members(self):
        model = fit_vocabulary({"A": {"p": "AAAAA", "q": "CCCCC"},
                                "B": {"p": "AAAAA", "q": "CCCCC"},
                                "C": {"p": "DDDDD", "q": "DDDDD"}})
        self.assertEqual(model.retained_kmer_count, 3)
        self.assertEqual(model.features, (("AAAAA", "CCCCC"), ("DDDDD",)))
        result = project_proteomes(model, {"only_nonrepresentative": {"p": "CCCCC"}})
        self.assertEqual(result["matrix"], [[1, 0]])

    def test_training_sample_order_and_input_protein_order_do_not_change_features(self):
        a = {"B": {"p": "AAAAA", "q": "CCCCC"}, "A": {"p": "AAAAA", "q": "CCCCC"}}
        b = {"A": {"q": "CCCCC", "p": "AAAAA"}, "B": {"q": "CCCCC", "p": "AAAAA"}}
        self.assertEqual(fit_vocabulary(a), fit_vocabulary(b))
        self.assertEqual(project_proteomes(fit_vocabulary(a), a)["sample_ids"], ["A", "B"])

    def test_invalid_or_missing_proteins_are_not_silently_cleaned(self):
        for value in ({}, {"A": {}}, {"A": {"p": ""}}, {"A": {"p": "acdef"}},
                      {"A": {"p": "ACDXF"}}, {"A": {"p": "ACDE*"}}, {" A": {"p": "AAAAA"}}):
            with self.subTest(value=value), self.assertRaises(ValueError):
                fit_vocabulary(value)

    def test_invalid_parameters_and_nonfitted_projection_fail(self):
        for params in ({"k": True}, {"k": 0}, {"minimum_protein_occurrences": 0},
                       {"minimum_protein_occurrences": 1.5}):
            with self.subTest(params=params), self.assertRaises(ValueError):
                fit_vocabulary({"A": {"p": "AAAAA"}}, **params)
        with self.assertRaises(ValueError):
            project_proteomes({}, {"A": {"p": "AAAAA"}})

    def test_explicit_deposited_policy_retains_U_and_skips_only_X_windows(self):
        proteins = {"A": {"p": "AAAAAU", "q": "AAAAAU", "r": "AAAAAXAAAAA"}}
        model = fit_vocabulary(proteins, sequence_policy="deposited21_skip_X")
        members = {mer for group in model.features for mer in group}
        self.assertIn("AAAAU", members)
        self.assertNotIn("AAAAX", members)
        result = project_proteomes(model, proteins)
        self.assertEqual(result["window_counts"]["A"], {
            "total_windows": 11, "skipped_X_windows": 5,
            "deposited_X_residues": 1, "deposited_U_residues": 2})
        with self.assertRaises(ValueError):
            fit_vocabulary(proteins)

    def test_skipped_X_does_not_join_fragments_or_inflate_protein_frequency(self):
        model = fit_vocabulary({"A": {"p": "ACDXEFG"}}, minimum_protein_occurrences=1,
                               sequence_policy="deposited21_skip_X")
        self.assertEqual(model.features, ())
        repeated = fit_vocabulary({"A": {"p": "AAAAAXAAAAA"}}, sequence_policy="deposited21_skip_X")
        self.assertEqual(repeated.features, ())

    def test_sequence_policy_remains_fixed_during_heldout_projection(self):
        model = fit_vocabulary({"A": {"p": "AAAAA", "q": "AAAAA"}})
        with self.assertRaises(ValueError):
            project_proteomes(model, {"B": {"p": "AAAAU"}})
        with self.assertRaises(ValueError):
            fit_vocabulary({"A": {"p": "AAAAA"}}, sequence_policy="guess_or_clean")


if __name__ == "__main__":
    unittest.main()
