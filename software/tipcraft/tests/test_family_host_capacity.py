"""Family conditioning and missing-member arithmetic of the declared pilot."""
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    'family_capacity', Path(__file__).resolve().parents[3] / 'scripts/check-receptor-family-host-capacity.py')
PILOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PILOT)


class FamilyHostCapacityTests(unittest.TestCase):
    def test_repetition_in_one_protein_cannot_pass_rare_filter(self):
        members, expected = PILOT.fit(['g'], {'p': 'AAAAAAA'}, {'p': 'g'}, {'p': {'family': 'a'}})
        self.assertEqual(members, {})
        self.assertEqual(expected, {})

    def test_same_kmer_in_different_families_does_not_rescue_singletons(self):
        members, _ = PILOT.fit(['g'], {'p': 'AAAAA', 'q': 'AAAAA'},
                               {'p': 'g', 'q': 'g'}, {'p': {'family': 'a'}, 'q': {'family': 'b'}})
        self.assertEqual(members, {})

    def test_self_projection_preserves_fit_patterns(self):
        sequences = {'p': 'AAAAA', 'q': 'AAAAA', 'r': 'CCCCC', 's': 'CCCCC'}
        genomes = {'p': 'g1', 'q': 'g2', 'r': 'g1', 's': 'g1'}
        hits = {p: {'family': 'a' if p in ['p', 'q'] else 'b'} for p in sequences}
        members, expected = PILOT.fit(['g1', 'g2'], sequences, genomes, hits)
        actual = PILOT.project(['g1', 'g2'], sequences, genomes, hits, members, expected, .5)
        self.assertEqual(actual['matrix'], {'g1': [1, 1], 'g2': [1, 0]})

    def test_reached_denominator_and_zero_rows_match_declared_adaptation(self):
        members = {('a', 'AAAAA'): 'f0', ('b', 'CCCCC'): 'f0'}
        seqs = {'p': 'AAAAA', 'q': 'TTTTT', 'r': 'AAAAA'}
        genomes = {'p': 'partial', 'q': 'negative', 'r': 'unassigned'}
        hits = {'p': {'family': 'a'}, 'q': {'family': 'a'}}
        result = PILOT.project(['partial', 'negative', 'unassigned'], seqs, genomes, hits,
                               members, {'f0': (1, 1)}, .75)
        self.assertEqual(result['matrix'], {'partial': [1], 'negative': [0], 'unassigned': [0]})
        self.assertEqual(result['member_match_counts']['partial']['f0']['fraction'], 1)
        self.assertIsNone(result['member_match_counts']['unassigned']['f0']['fraction'])

    def test_duplicate_query_proteins_count_each_member_once(self):
        result = PILOT.project(['g'], {'p': 'AAAAA', 'q': 'AAAAA'}, {'p': 'g', 'q': 'g'},
                               {'p': {'family': 'a'}, 'q': {'family': 'a'}},
                               {('a', 'AAAAA'): 'f0', ('a', 'CCCCC'): 'f0'}, {'f0': (1,)}, .5)
        counts = result['member_match_counts']['g']['f0']
        self.assertEqual(counts, {'matched_members': 1, 'reached_members': 2, 'fraction': .5})
        self.assertEqual(result['matrix']['g'], [1])


if __name__ == '__main__':
    unittest.main()
