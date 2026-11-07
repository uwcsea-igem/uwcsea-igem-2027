"""A collision is an endpoint limitation only for contrasts it erases."""
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location('collision_endpoint',Path(__file__).resolve().parents[3]/'scripts/audit-host-collisions-against-endpoint.py')
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class HostCollisionEndpointTests(unittest.TestCase):
    def test_identical_negative_hosts_do_not_force_a_positive_negative_tie(self):
        rows=[{'phage':'p','host':h,'label':label} for h,label in [('a',1),('b',0),('c',0)]]
        result=AUDIT.diagnose(rows,{'a':[1,0],'b':[0,1],'c':[0,1]})
        self.assertEqual(result['forced_tied_primary_pair_count'],0)
        self.assertEqual(result['collision_groups'],[['b','c']])

    def test_mixed_label_collisions_preserve_forced_ties(self):
        rows=[{'phage':'p','host':h,'label':label} for h,label in [('a',1),('b',0),('c',0)]]
        result=AUDIT.diagnose(rows,dict.fromkeys(['a','b','c'],[1]))
        self.assertEqual(result['forced_tied_primary_pair_count'],2)

    def test_all_positive_query_supplies_no_concordance_contrast(self):
        result=AUDIT.diagnose([{'phage':'p','host':'a','label':1}],{'a':[1]})
        self.assertEqual(result['informative_queries'],[])

    def test_duplicate_observations_or_missing_vectors_fail(self):
        row={'phage':'p','host':'a','label':1}
        for rows,vectors in [([row,row],{'a':[1]}),([row],{'b':[1]})]:
            with self.assertRaises(ValueError):AUDIT.diagnose(rows,vectors)

    def test_multiple_queries_can_require_every_host_pair_to_separate(self):
        rows=[{'phage':q,'host':h,'label':int(h==q)} for q in ['a','b','c'] for h in ['a','b','c']]
        result=AUDIT.diagnose(rows,{'a':[1,0,0],'b':[0,1,0],'c':[0,0,1]})
        self.assertTrue(result['all_hosts_distinct_required_by_this_endpoint'])
        self.assertEqual(len(result['required_host_separations']),3)


if __name__=='__main__':
    unittest.main()
