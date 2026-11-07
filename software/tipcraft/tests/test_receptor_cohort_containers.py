"""Integrity and annotation boundaries for the public-container audit."""
import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("cohort_container_audit", ROOT / "scripts/audit-receptor-cohort-containers.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
HAS_BIO = importlib.util.find_spec("Bio") is not None


def item(data):
    return {"filename": "example.gbk", "primary_id": "example", "repository_path": "Ecoli_phages/data/phage_genomes/example.gbk",
            "declared_size_bytes": len(data), "git_blob_sha1": hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()}


def genbank_bytes():
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    first = SeqRecord(Seq("ACGTN"), id="first", name="first", description="unused phenotype text")
    first.annotations["molecule_type"] = "DNA"
    first.features = [SeqFeature(SimpleLocation(0, 3), type="CDS", qualifiers={"translation": ["MXX"], "transl_table": ["11"]}),
                      SeqFeature(SimpleLocation(2, 5), type="CDS")]
    second = SeqRecord(Seq("RY"), id="second", name="second", description="unused")
    second.annotations["molecule_type"] = "DNA"
    stream = io.StringIO()
    SeqIO.write([first, second], stream, "genbank")
    return stream.getvalue().encode()


class ContainerAuditTests(unittest.TestCase):
    def test_verify_rejects_same_size_corruption(self):
        with self.assertRaisesRegex(ValueError, "SHA1"):
            AUDIT.verify(b"abd", item(b"abc"))

    def test_framing_preserves_record_boundaries(self):
        self.assertNotEqual(AUDIT.framed_sequences(["AC", "GT"]), AUDIT.framed_sequences(["A", "CGT"]))
        self.assertNotEqual(AUDIT.framed_sequences(["AC", "GT"]), AUDIT.framed_sequences(["GT", "AC"]))

    def test_atomic_failure_leaves_no_partial_or_destination(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result"
            with mock.patch.object(AUDIT.os, "replace", side_effect=OSError("simulated rename failure")):
                with self.assertRaises(OSError):
                    AUDIT.publish(path, b"data")
            self.assertFalse(path.exists())
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_changed_cache_is_not_downloaded_or_repaired(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            path = tmp / "work/containers/example.gbk"
            path.parent.mkdir(parents=True)
            path.write_bytes(b"abd")
            with mock.patch.object(AUDIT.subprocess, "run", side_effect=AssertionError("must not download")) as run:
                result = AUDIT.retrieve(item(b"abc"), {"existing_target_cache": str(tmp / "other"), "source_commit": "fixed"}, tmp / "work")
            run.assert_not_called()
            self.assertEqual(result["payload_status"], "error")
            self.assertEqual(path.read_bytes(), b"abd")

    def test_download_failure_preserves_entry_and_cleans_partial(self):
        def fail(command, **kwargs):
            Path(command[command.index("--output") + 1]).write_bytes(b"partial")
            return SimpleNamespace(returncode=22, stderr="HTTP 404")
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            with mock.patch.object(AUDIT.subprocess, "run", side_effect=fail):
                result = AUDIT.retrieve(item(b"abc"), {"existing_target_cache": str(tmp / "other"), "source_commit": "fixed"}, tmp / "work")
            self.assertEqual(result["source"]["primary_id"], "example")
            self.assertEqual(result["payload_status"], "error")
            self.assertIn("404", result["error"])
            self.assertEqual(list((tmp / "work/containers").iterdir()), [])

    @unittest.skipUnless(HAS_BIO, "optional Biopython parser absent")
    def test_deposited_proteins_only_missing_annotations_remain_missing(self):
        summary, dna, aa = AUDIT.parse_container(genbank_bytes(), "c000")
        self.assertEqual(dna, [("c000r000", "ACGTN"), ("c000r001", "RY")])
        self.assertEqual(aa, [("c000r000p00000v00", "MXX")])
        self.assertEqual(summary["records"][0]["cds_without_translation"], 1)
        self.assertEqual(summary["records"][0]["proteins"][0]["translation_agreement_with_DNA"], "not_verified")
        self.assertEqual(summary["records"][0]["ambiguous_bases"], 1)
        self.assertNotIn("unused phenotype text", repr(summary))
        self.assertTrue(all(not name.startswith("t") for name, seq in dna))

    @unittest.skipUnless(HAS_BIO, "optional Biopython parser absent")
    def test_dna_available_despite_protein_annotation_errors(self):
        data = genbank_bytes().replace(b"MXX", b"M?X")
        source = {**item(data), "container_id": "c000"}
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            path = tmp / "container.gbk"
            path.write_bytes(data)
            result = AUDIT.parse_one({"source": source, "container_path": str(path), "payload_status": "authenticated"}, tmp)
            self.assertEqual(result["nucleotide_status"], "available")
            self.assertEqual(result["protein_status"], "annotation_error")
            self.assertIn("nucleotide_fasta", result["outputs"])
            self.assertNotIn("protein_fasta", result["outputs"])
            self.assertTrue(result["protein_errors"])

    def test_compact_report_keeps_failures_and_binds_full_proteins(self):
        protein = {"length": 3, "letters_outside_standard_20": {"X": 1}, "internal_stop_count": 0, "terminal_stop_present": False}
        record = {"length": 5, "nucleotide_alphabet_counts": {"A": 4, "N": 1}, "proteins": [protein]}
        original = {"summary": {"errors": 1}, "sequence_manifest": {"path": "full.json", "sha256": "bound"},
                    "entries": [{"source": {"primary_id": "ok"}, "records": [record]},
                                {"source": {"primary_id": "failed"}, "status": "payload_error", "error": "404"}]}
        result = AUDIT.compact_report(original)
        self.assertEqual(len(result["entries"]), 2)
        self.assertEqual(result["entries"][1]["error"], "404")
        self.assertEqual(result["entries"][0]["records"][0]["deposited_protein_metadata_sha256"], AUDIT.sha(AUDIT.encoded([protein])))
        self.assertNotIn("proteins", result["entries"][0]["records"][0])
        self.assertEqual(original["entries"][0]["records"][0]["proteins"], [protein])
        self.assertEqual(result["summary"]["nonstandard_amino_acid_characters"], {"X": 1})
        self.assertEqual(result["summary"]["ambiguous_nucleotide_characters"], {"N": 1})


if __name__ == "__main__":
    unittest.main()
