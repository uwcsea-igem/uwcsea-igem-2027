import importlib.util
import hashlib
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("phage_normalization", Path(__file__).resolve().parents[3] / "scripts/normalize-genophi-phages.py")
NORMALIZE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NORMALIZE)


class PhageNormalizationTests(unittest.TestCase):
    def test_reference_protein_requires_matching_accession_dna_and_annotation(self):
        data = (b'VERSION     REF.1\nFEATURES             Location/Qualifiers\n'
                b'     CDS             1..6\n                     /product="tail fiber protein"\n'
                b'                     /protein_id="P.1"\n                     /translation="M"\n'
                b'ORIGIN\n        1 atgtaa\n//\n')
        digest = hashlib.sha256(b"ATGTAA").hexdigest()
        protein, metadata = NORMALIZE.extract_gp17(data, "REF.1", "P.1", digest)
        self.assertEqual(protein, "M")
        self.assertEqual(metadata["cds_location"], [1, 6])
        with self.assertRaisesRegex(ValueError, "disagrees with CDS"):
            NORMALIZE.extract_gp17(data.replace(b'/translation="M"', b'/translation="C"'), "REF.1", "P.1", digest)
        for accession, protein_id, genome_hash in [("REF.2", "P.1", digest), ("REF.1", "P.2", digest), ("REF.1", "P.1", "wrong")]:
            with self.subTest(accession=accession, protein_id=protein_id), self.assertRaises(ValueError):
                NORMALIZE.extract_gp17(data, accession, protein_id, genome_hash)

    def test_exact_reviewed_empty_record_keeps_the_real_sequence(self):
        source = b">old\n>old\nacgt\nnn\n"
        sequence, changes = NORMALIZE.parse_complete_genome(source, [0])
        self.assertEqual(sequence, "ACGTNN")
        self.assertEqual(changes["removed_empty_record_indices"], [0])
        with self.assertRaises(ValueError):
            NORMALIZE.parse_complete_genome(source)

    def test_actual_contigs_are_never_concatenated(self):
        with self.assertRaisesRegex(ValueError, "never concatenate"):
            NORMALIZE.parse_complete_genome(b">first\nACGT\n>second\nACGT\n")

    def test_unexpected_empty_record_or_non_dna_is_rejected(self):
        for data, exceptions in [(b">first\nACGT\n>empty\n", [0]), (b">x\nACZ\n", []), (b"ACGT\n", [])]:
            with self.subTest(data=data), self.assertRaises(ValueError):
                NORMALIZE.parse_complete_genome(data, exceptions)
