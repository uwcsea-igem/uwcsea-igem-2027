"""Target-related connected components cannot leak into candidate hosts."""
import importlib.util
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location('ecor_groups',Path(__file__).resolve().parents[3]/'scripts/group-ecor-host-references.py')
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class EcorHostGroupTests(unittest.TestCase):
    def test_target_component_excludes_source_even_with_author_reference_bridge(self):
        groups=[['original:MG1655','author:Keio','ecor:a'],['ecor:b','ecor:c']]
        result=AUDIT.partition_components(groups,[name for group in groups for name in group])
        self.assertEqual(result['target_related_source_hosts'],['ecor:a'])
        self.assertEqual(result['candidate_source_hosts'],['ecor:b','ecor:c'])
        self.assertEqual(result['candidate_source_group_count'],1)

    def test_missing_comparison_state_produces_no_candidate_partition(self):
        result=AUDIT.partition_components(None,['ecor:a','original:BL21'])
        self.assertIsNone(result['candidate_source_hosts'])
        self.assertIsNone(result['target_related_source_hosts'])

    def test_duplicate_or_missing_component_member_fails(self):
        for groups in [[['ecor:a','ecor:a']],[['ecor:a']]]:
            with self.assertRaises(ValueError):AUDIT.partition_components(groups,['ecor:a','ecor:b'])

    def test_group_identity_is_independent_of_member_order(self):
        first=AUDIT.partition_components([['ecor:a','ecor:b']],['ecor:a','ecor:b'])
        second=AUDIT.partition_components([['ecor:b','ecor:a']],['ecor:a','ecor:b'])
        self.assertEqual(first['source_groups'],second['source_groups'])


if __name__=='__main__':unittest.main()
