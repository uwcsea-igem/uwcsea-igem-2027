"""A new official genotype must remain separate from its same-name author node."""
import importlib.util
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location('official_variant_groups',Path(__file__).resolve().parents[3]/'scripts/group-official-basel-variants.py')
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def source():
    return {'Bas06':{'official_versioned_accession':'PZ445410.1','official_length':6,
                     'official_sequence_sha256':'new-sequence','official_ambiguous_bases':0}}


def row():
    return {'id':'official_Bas06','source_id':'Bas06','role':'official_reference_version',
            'container_key':'official:Bas06:PZ445410.1','source_record_id':'PZ445410.1',
            'retained_bases':6,'sequence_sha256':'new-sequence'}


class OfficialVariantMappingTests(unittest.TestCase):
    def test_same_name_distinct_sequences_preserve_both_nodes(self):
        records={'old_record':{'container_key':'source:Bas06','retained_bases':5}}
        containers={'source:Bas06':{'length':5}}
        new_records,new_containers=AUDIT.add_official_records(records,containers,[row()],source())
        self.assertEqual(len(new_records),2)
        self.assertEqual(len(new_containers),2)
        self.assertEqual(new_records['official_Bas06']['offset'],0)
        self.assertEqual(containers,{'source:Bas06':{'length':5}})

    def test_different_accession_length_hash_or_namespace_fails(self):
        for field,value in [('container_key','source:Bas06'),('retained_bases',5),
                            ('sequence_sha256','old-sequence'),('source_record_id','PZ445410.2')]:
            candidate=row();candidate[field]=value
            with self.assertRaises(ValueError):AUDIT.add_official_records({}, {},[candidate],source())

    def test_duplicate_or_missing_official_records_fail(self):
        for rows in [[row(),row()],[]]:
            with self.assertRaises(ValueError):AUDIT.add_official_records({}, {},rows,source())

    def test_official_cannot_reuse_old_record_identifier(self):
        with self.assertRaises(ValueError):
            AUDIT.add_official_records({'official_Bas06':{}},{},[row()],source())

    def test_unknown_exit_is_never_silently_zero(self):
        frozen={'subjects_fasta':{},'record_manifest':{},'policy_sha256':'p','script_sha256':'s','commands':[['make'],['search']]}
        alignment={**frozen,'status':'validated_all_record_nucleotide_hsps','commands':[{'argv':['make'],'exit_code':0},{'argv':['search'],'exit_code':None}]}
        with self.assertRaises(ValueError):AUDIT.require_execution_evidence(alignment,frozen,{})
        alignment['commands'][1]['exit_code']=0
        self.assertEqual(AUDIT.require_execution_evidence(alignment,frozen,{}),'both_original_command_exit_codes_observed_zero')

    def test_recovered_data_requires_explicit_policy_and_preserves_unknown_exit(self):
        frozen={'subjects_fasta':{},'record_manifest':{},'policy_sha256':'p','script_sha256':'s','commands':[['make'],['search']]}
        alignment={**frozen,'status':'validated_recovered_all_record_nucleotide_hsps',
                   'commands':[{'argv':['make'],'exit_code':0},{'argv':['search'],'exit_code':None,'terminal_observed':True}],
                   'execution_status':'orphaned_controller_exit_code_unobservable','execution_exit_success_claimed':False,
                   'artifact_completion_verified':True,'search_rerun':False,
                   'recovery':{'helper':{'sha256':'h'},'policy':{'sha256':'r'}}}
        with self.assertRaises(ValueError):AUDIT.require_execution_evidence(alignment,frozen,{})
        policy={'allow_explicit_orphan_artifact_recovery':True,'inputs':{'recovery_helper':{'sha256':'h'},'recovery_policy':{'sha256':'r'}}}
        self.assertEqual(AUDIT.require_execution_evidence(alignment,frozen,policy),'artifact_validated_but_search_exit_code_unknown')
        alignment['commands'][1]['exit_code']=1
        with self.assertRaises(ValueError):AUDIT.require_execution_evidence(alignment,frozen,policy)


if __name__=='__main__':unittest.main()
