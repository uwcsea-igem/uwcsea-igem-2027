import importlib.util
from pathlib import Path
import unittest


SPEC = importlib.util.spec_from_file_location("demo_cds", Path(__file__).resolve().parents[3] / "scripts/audit-genophi-demo-cds.py")
CDS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CDS)


def gene(protein, end=9, strand=1, partial="00", start_type="GTG"):
    return {"id": "example", "start": 1, "end": end, "strand": strand,
            "partial": partial, "start_type": start_type, "code": 11, "protein": protein}


class DemoCDSTests(unittest.TestCase):
    def test_alternative_initiator_and_reverse_complement(self):
        self.assertTrue(CDS.translate_check("GTGGGTTAA", gene("MG*"))["passed"])
        self.assertTrue(CDS.translate_check("TTAACCCAC", gene("MG*", strand=-1))["passed"])
        changed = CDS.translate_check("TTGGGTTAA", gene("MG*"))
        self.assertFalse(changed["passed"])
        self.assertFalse(changed["declared_start_consistent"])

    def test_partial_flags_follow_genomic_ends_on_both_strands(self):
        self.assertTrue(CDS.translate_check("GGTTAA", gene("G*", end=6, partial="10", start_type="Edge"))["passed"])
        self.assertTrue(CDS.translate_check("TTAACC", gene("G*", end=6, strand=-1, partial="01", start_type="Edge"))["passed"])
        self.assertTrue(CDS.translate_check("ATGGGT", gene("MG", end=6, partial="01", start_type="ATG"))["passed"])
        self.assertTrue(CDS.translate_check("ACCCAT", gene("MG", end=6, strand=-1, partial="10", start_type="ATG"))["passed"])

    def test_wrong_protein_stop_and_unsupported_coordinates_are_not_repaired(self):
        for example in (gene("MT*"), gene("MG"), gene("MG*", partial="01")):
            self.assertFalse(CDS.translate_check("GTGGGTTAA", example)["passed"])
        for example in (gene("MG*", end=12), gene("MG*", end=8), gene("MG*", start_type="AAA")):
            with self.assertRaises(ValueError):
                CDS.translate_check("GTGGGTTAA", example)
        with self.assertRaises(ValueError):
            CDS.translate_check("GTGNNNTAA", gene("MG*"))

    def test_coordinate_comparison_separates_changed_start_from_changed_stop_and_strand(self):
        original = [gene("MG*")]
        shifted = gene("MG*"); shifted["start"] = 4
        shared = CDS.coordinate_comparison(original, original)
        self.assertEqual(shared["exact_shared_coordinate_copies"], 1)
        self.assertEqual(shared["unshared_coordinates"], [])
        for result in CDS.coordinate_comparison(original, [shifted])["unshared_coordinates"]:
            self.assertTrue(result["same_strand_stop_in_other"])
        for result in CDS.coordinate_comparison(original, [gene("MG*", strand=-1)])["unshared_coordinates"]:
            self.assertFalse(result["same_strand_stop_in_other"])

    def test_coordinate_fasta_parser_rejects_empty_ambiguous_duplicate_and_unframed_entries(self):
        valid = b">p # 1 # 9 # 1 # partial=00;start_type=GTG;genetic_code=11\nMG*\n"
        self.assertEqual(CDS.parse_genes(valid)[0]["protein"], "MG*")
        for data in (valid + valid, b">p\nMG*\n", valid.replace(b"MG*", b"MX*"), valid.replace(b"MG*", b""), valid.replace(b"MG*", b"M*G")):
            with self.subTest(data=data), self.assertRaises(ValueError):
                CDS.parse_genes(data)
