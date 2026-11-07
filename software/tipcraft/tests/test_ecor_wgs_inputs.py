"""Preserve complete versioned WGS sets rather than accepting partial genomes."""
import importlib.util
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location('ecor_inputs',Path(__file__).resolve().parents[3]/'scripts/authenticate-ecor-wgs-inputs.py')
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)

MASTER=b'''ID   QOWX01000000; SV 01; linear; genomic DNA; SET; PRO; 2 SQ.
PR   Project:PRJNA230969;
DR   BioSample; SAMN04158349.
FT                   /strain="MOD1-ECOR12"
WGS  QOWX01000001-QOWX01000002
//
'''
FASTA=b'''>ENA|QOWX01000001|QOWX01000001.1 Escherichia coli strain MOD1-ECOR12 contig1
ACGTN
>ENA|QOWX01000002|QOWX01000002.2 Escherichia coli strain MOD1-ECOR12 contig2
TTTTA
'''


def master(data=MASTER):return AUDIT.parse_master(data,'ECOR12','QOWX00000000','PRJNA230969')


class EcorWgsInputTests(unittest.TestCase):
    def test_master_build_is_distinct_from_individual_sequence_versions(self):
        m=master();records=AUDIT.parse_fasta(FASTA,m)
        self.assertEqual(m['wgs_build'],'01')
        self.assertEqual([r['sequence_version'] for r in records],[1,2])
        self.assertEqual(m['biosample'],'SAMN04158349')

    def test_every_contig_and_ambiguity_is_preserved(self):
        records=AUDIT.parse_fasta(FASTA,master())
        self.assertEqual(len(records),2)
        self.assertEqual(sum(r['length'] for r in records),10)
        self.assertEqual(sum(r['ambiguous_bases'] for r in records),1)

    def test_partial_set_is_rejected(self):
        with self.assertRaises(ValueError):AUDIT.parse_fasta(FASTA.split(b'>ENA|QOWX01000002')[0],master())

    def test_duplicate_or_unexpected_contigs_are_rejected(self):
        for data in [FASTA+FASTA,FASTA.replace(b'QOWX01000002',b'QOWX01000003')]:
            with self.assertRaises(ValueError):AUDIT.parse_fasta(data,master())

    def test_wrong_master_identity_strain_project_or_count_is_rejected(self):
        for data in [MASTER.replace(b'QOWX',b'QOWY'),MASTER.replace(b'ECOR12',b'ECOR13'),
                     MASTER.replace(b'PRJNA230969',b'PRJNA999999'),MASTER.replace(b'2 SQ.',b'3 SQ.')]:
            with self.assertRaises(ValueError):master(data)

    def test_unversioned_wrong_strain_or_invalid_dna_is_rejected(self):
        for data in [FASTA.replace(b'000001.1',b'000001'),FASTA.replace(b'ECOR12',b'ECOR13'),FASTA.replace(b'ACGTN',b'ACGT*')]:
            with self.assertRaises(ValueError):AUDIT.parse_fasta(data,master())


if __name__=='__main__':unittest.main()
