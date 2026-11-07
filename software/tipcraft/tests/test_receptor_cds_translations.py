"""Conservative checking rules; optional parser is absent from core CI."""
import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("receptor_cds_check", ROOT / "scripts/check-receptor-cds-translations.py")
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)
HAS_BIO = importlib.util.find_spec("Bio") is not None


class QualifierTests(unittest.TestCase):
    def test_missing_table_is_not_defaulted(self):
        self.assertEqual(CHECK.qualifiers({})["status"], "unverified")
        self.assertEqual(CHECK.qualifiers({"transl_table": ["11", "11"]})["status"], "unverified")

    def test_absent_codon_start_default_is_disclosed(self):
        result = CHECK.qualifiers({"transl_table": ["11"]})
        self.assertEqual(result["codon_start"], 1)
        self.assertEqual(result["codon_start_source"], "Biopython_1.86_absence_default_1")
        self.assertEqual(CHECK.qualifiers({"transl_table": ["11"], "codon_start": ["3"]})["codon_start"], 3)

    def test_malformed_offsets_are_not_repaired(self):
        for offset in [[], ["0"], ["4"], ["1", "2"], ["1.0"]]:
            with self.subTest(offset=offset):
                self.assertEqual(CHECK.qualifiers({"transl_table": ["11"], "codon_start": offset})["status"], "unverified")

    def test_recoding_and_pseudo_qualifiers_are_unverified(self):
        for key in CHECK.UNSUPPORTED:
            with self.subTest(key=key):
                result = CHECK.qualifiers({"transl_table": ["11"], key: []})
                self.assertEqual(result["status"], "unverified")
                self.assertEqual(result["qualifier_keys"], [key])


@unittest.skipUnless(HAS_BIO, "optional Biopython1.86 absent")
class TranslationTests(unittest.TestCase):
    def feature(self, dna, translation="MK", strand=1, start=None, codon=None):
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        q = {"transl_table": ["11"], "translation": [translation]}
        if codon is not None:
            q["codon_start"] = [str(codon)]
        return SeqFeature(SimpleLocation(0 if start is None else start, len(dna), strand=strand), type="CDS", qualifiers=q), Seq(dna)

    def test_alternative_start_and_reverse_strand_match(self):
        for dna, strand in [("GTGAAATAA", 1), ("TTATTTCAC", -1)]:
            feature, seq = self.feature(dna, strand=strand)
            result = CHECK.check(feature, seq, "MK")
            self.assertEqual(result["status"], "exact_match")
            self.assertEqual(result["computed_sha256"], CHECK.sha(b"MK"))
            self.assertNotIn("computed_sequence", result)

    def test_mismatch_is_not_a_sequence_export(self):
        feature, seq = self.feature("ATGAAATAA", translation="MA")
        result = CHECK.check(feature, seq, "MA")
        self.assertEqual(result["status"], "mismatch")
        self.assertNotIn("MK", result.values())

    def test_cds_start_and_terminal_stop_checks_are_enforced(self):
        for dna in ["AAAAAATAA", "ATGAAA", "ATGTAATAA"]:
            feature, seq = self.feature(dna)
            self.assertEqual(CHECK.check(feature, seq, "MK")["status"], "cds_validation_failed")

    def test_partial_and_ambiguous_coding_dna_are_unverified(self):
        from Bio.SeqFeature import BeforePosition
        feature, seq = self.feature("ATGAAATAA", start=BeforePosition(0))
        self.assertEqual(CHECK.check(feature, seq, "MK")["reason"], "partial_or_fuzzy_location")
        feature, seq = self.feature("ATGNNNTAA")
        self.assertEqual(CHECK.check(feature, seq, "MX")["reason"], "ambiguous_coding_DNA")

    def test_declared_frame_offset_and_exact_join(self):
        from Bio.Seq import Seq
        from Bio.SeqFeature import CompoundLocation, SimpleLocation
        feature, seq = self.feature("AATGAAATAA", codon=2)
        self.assertEqual(CHECK.check(feature, seq, "MK")["status"], "exact_match")
        feature.location = CompoundLocation([SimpleLocation(0, 3, strand=1), SimpleLocation(6, 12, strand=1)], operator="join")
        feature.qualifiers["codon_start"] = ["1"]
        self.assertEqual(CHECK.check(feature, Seq("ATGCCCAAATAA"), "MK")["status"], "exact_match")


if __name__ == "__main__":
    unittest.main()
