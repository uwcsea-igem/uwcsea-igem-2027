"""Official accession identity and DNA equality cannot be inferred from names."""
import importlib.util
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location('phage_accessions',Path(__file__).resolve().parents[3]/'scripts/authenticate-basel-phage-accessions.py')
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
HAS_BIO=importlib.util.find_spec('Bio') is not None


class BaselPhageAccessionTests(unittest.TestCase):
    def test_formats_and_author_must_match_exact_dna(self):
        annotation={'versioned_accession':'PZ123456.1','sequence':'ACGT'}
        self.assertTrue(AUDIT.compare_records(annotation,dict(annotation),{'sequence':'ACGT'}))
        with self.assertRaises(ValueError):AUDIT.compare_records(annotation,dict(annotation),{'sequence':'ACGA'})

    def test_official_format_version_disagreement_fails(self):
        a={'versioned_accession':'PZ123456.1','sequence':'ACGT'}
        b={'versioned_accession':'PZ123456.2','sequence':'ACGT'}
        with self.assertRaises(ValueError):AUDIT.compare_records(a,b,{'sequence':'ACGT'})

    @unittest.skipUnless(HAS_BIO,'Optional Biopython parser')
    def test_ena_wrapper_resolves_exact_requested_version(self):
        data=b'>ENA|PZ123456|PZ123456.2 example\nACGTN\n'
        result=AUDIT.sequence_record(data,'fasta','PZ123456')
        self.assertEqual(result['versioned_accession'],'PZ123456.2')
        self.assertEqual(result['ambiguous_bases'],1)

    @unittest.skipUnless(HAS_BIO,'Optional Biopython parser')
    def test_unversioned_wrong_accession_multi_record_and_invalid_dna_fail(self):
        for data in [b'>PZ123456\nACGT\n',b'>PZ999999.1\nACGT\n',
                     b'>PZ123456.1\nACGT\n>PZ123456.1\nACGT\n',b'>PZ123456.1\nACGX\n']:
            with self.assertRaises(ValueError):AUDIT.sequence_record(data,'fasta','PZ123456')


if __name__=='__main__':unittest.main()
