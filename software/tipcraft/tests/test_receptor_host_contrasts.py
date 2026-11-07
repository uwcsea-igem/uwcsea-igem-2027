"""Count source contrasts without treating unmeasured annotations as negatives."""
import importlib.util
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location('host_contrasts',Path(__file__).resolve().parents[3]/'scripts/audit-receptor-host-contrasts.py')
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def metadata(**pairs):
    return {k:dict(zip(AUDIT.COLUMNS,v)) for k,v in pairs.items()}


class ReceptorHostContrastTests(unittest.TestCase):
    def test_nonbinary_and_unmeasured_literals_are_retained(self):
        result=AUDIT.summarize(['source:a','source:b'],metadata(a=('Yes','Not assayed'),b=('No','Resistant')),
                               {'source:a':'g','source:b':'h'})
        self.assertEqual(result['complete_literal_pairs'],0)
        self.assertEqual(len(result['incomplete_or_nonbinary_pairs']),2)
        self.assertEqual(result['literal_state_counts'][AUDIT.COLUMNS[1]],{'Not assayed':1,'Resistant':1})
        self.assertFalse(result['has_any_mixed_state_comparisons'])

    def test_one_orientation_cannot_identify_phage_dependent_ordering(self):
        result=AUDIT.summarize(['source:a','source:b'],metadata(a=('Yes','No'),b=('Yes','Yes')),
                               {'source:a':'g','source:b':'h'})
        self.assertTrue(result['fixed_host_order_can_satisfy_every_mixed_state_comparison'])
        self.assertFalse(result['opposing_host_preferences_observed'])

    def test_opposing_orientations_require_different_orderings(self):
        result=AUDIT.summarize(['source:a','source:b'],metadata(a=('Yes','No'),b=('No','Yes')),
                               {'source:a':'g','source:b':'h'})
        self.assertFalse(result['fixed_host_order_can_satisfy_every_mixed_state_comparison'])
        self.assertTrue(result['opposing_host_preferences_observed'])

    def test_related_phages_are_not_counted_as_independent_groups(self):
        result=AUDIT.summarize(['source:a','source:b'],metadata(a=('Yes','No'),b=('Yes','No')),
                               {'source:a':'g','source:b':'g'})
        self.assertEqual(result['categories']['BW_Yes__BL_No']['count'],2)
        self.assertEqual(result['categories']['BW_Yes__BL_No']['group_count'],1)

    def test_duplicate_or_unmapped_sources_fail(self):
        for keys,groups in [(['source:a','source:a'],{'source:a':'g'}),(['source:a'],{})]:
            with self.assertRaises(ValueError):AUDIT.summarize(keys,metadata(a=('Yes','No')),groups)


if __name__=='__main__':unittest.main()
