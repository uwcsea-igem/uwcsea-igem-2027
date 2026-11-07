import importlib.util
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location('ecor_contrasts',Path(__file__).resolve().parents[3]/'scripts/audit-basel-ecor-contrasts.py')
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def row(host,phage,value):return {'strain':host,'phage':phage,'interaction':value}


class BaselEcorContrastTests(unittest.TestCase):
    def test_opposite_orientations_require_phage_dependent_host_order(self):
        result=AUDIT.summarize([row('a','p','1'),row('b','p','0'),row('a','q','0.0'),row('b','q','1.0')])
        self.assertEqual(result['host_pairs_with_opposing_directions'],1)
        self.assertEqual(result['informative_query_count'],2)

    def test_missing_observation_never_becomes_negative(self):
        result=AUDIT.summarize([row('a','p','1'),row('b','p','0'),row('a','q','1')])
        self.assertEqual(result['missing_pairs'],1)
        self.assertEqual(result['queries']['q']['negative'],0)
        self.assertEqual(result['informative_query_count'],1)

    def test_all_positive_query_contributes_no_ordering_comparisons(self):
        result=AUDIT.summarize([row('a','p','1'),row('b','p','1')])
        self.assertEqual(result['total_positive_negative_comparisons'],0)

    def test_duplicate_unknown_and_nonfinite_values_fail(self):
        for rows in [[row('a','p','1')]*2,[row('a','p','')],[row('a','p','NaN')],[row('a','p','0.5')]]:
            with self.assertRaises(ValueError):AUDIT.summarize(rows)


if __name__=='__main__':unittest.main()
