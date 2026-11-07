"""Exercise content authentication and atomic publication without network access."""
import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "prepare_genophi", Path(__file__).resolve().parents[3] / "scripts/prepare-genophi-genomes.py")
PREPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)
FASTA = b">one description\nACGTRN\n>two\nttac\n"


def host_record(data=FASTA):
    return {"kind": "host", "source_id": "published-id", "source_name": "published-id.fasta",
            "path": "hosts/1.fasta", "url": "https://example.invalid/published.fasta",
            "expected_bytes": len(data), "expected_md5": hashlib.md5(data).hexdigest()}


class PrepareGenoPHITests(unittest.TestCase):
    def test_sequence_hash_retains_boundaries_and_normalizes_wrapping(self):
        first = PREPARE.fasta_stats(b">a\nAC\n>b\nGT\n")
        changed_boundaries = PREPARE.fasta_stats(b">a\nA\n>b\nCGT\n")
        renamed_wrapped = PREPARE.fasta_stats(b">renamed\na\nc\n>other\ngt\n")
        self.assertNotEqual(first["sequence_sha256"], changed_boundaries["sequence_sha256"])
        self.assertEqual(first["sequence_sha256"], renamed_wrapped["sequence_sha256"])
        actual = PREPARE.fasta_stats(FASTA)
        self.assertEqual((actual["contigs"], actual["contig_lengths"], actual["ambiguous_bases"], actual["n_bases"]),
                         (2, [6, 4], 2, 1))

    def test_repeated_source_headers_preserve_both_contigs_and_report_issue(self):
        result = PREPARE.fasta_stats(b">a\nAC\n>a other-description\nGT\n")
        unique = PREPARE.fasta_stats(b">first\nAC\n>second\nGT\n")
        self.assertEqual(result["contigs"], 2)
        self.assertEqual(result["duplicate_header_identifiers"], ["a"])
        self.assertTrue(result["requires_header_normalization"])
        self.assertEqual(result["sequence_sha256"], unique["sequence_sha256"])

    def test_rejects_malformed_empty_and_non_dna_fasta(self):
        for malformed in (b"", b"ACGT\n", b">\nAC\n", b">a\n>b\nAC\n",
                          b">a\nACZ\n", b">a\nACU\n"):
            with self.subTest(malformed=malformed), self.assertRaises(ValueError):
                PREPARE.fasta_stats(malformed)

    def test_resumes_verified_files_and_rejects_changed_cache_without_download(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = host_record()
            target = root / item["path"]
            target.parent.mkdir()
            target.write_bytes(FASTA)
            with patch.object(PREPARE.subprocess, "run") as download:
                before = PREPARE.prepare_one(item, root)
                target.write_bytes(FASTA.replace(b"ACGTRN", b"ACGTNN"))
                with self.assertRaisesRegex(ValueError, "MD5"):
                    PREPARE.prepare_one(item, root)
                download.assert_not_called()
            self.assertEqual(before["status"], "verified")
            self.assertEqual(target.read_bytes(), FASTA.replace(b"ACGTRN", b"ACGTNN"))

    def test_git_blob_authentication_rejects_same_size_change(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "phage.fna"
            path.write_bytes(FASTA)
            digest = hashlib.sha1(b"blob " + str(len(FASTA)).encode() + b"\0" + FASTA).hexdigest()
            item = {"expected_bytes": len(FASTA), "expected_git_blob_sha1": digest}
            PREPARE.verify(path, item)
            path.write_bytes(FASTA.replace(b"ACGTRN", b"ACGTNN"))
            with self.assertRaisesRegex(ValueError, "Git blob"):
                PREPARE.verify(path, item)

    def test_download_is_verified_before_publication_and_offline_report_matches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = host_record()
            def download(argv, **kwargs):
                self.assertIn("--max-time", argv)
                self.assertIn("--retry", argv)
                self.assertFalse((root / item["path"]).exists())
                Path(argv[argv.index("--output") + 1]).write_bytes(FASTA)
                return subprocess.CompletedProcess(argv, 0, stderr="")
            with patch.object(PREPARE.subprocess, "run", side_effect=download):
                online = PREPARE.prepare_inventory([item], {"source": "pinned"}, root)
            with patch.object(PREPARE.subprocess, "run") as network:
                offline = PREPARE.prepare_inventory([item], {"source": "pinned"}, root, offline=True)
                network.assert_not_called()
            self.assertEqual(online, offline)
            self.assertTrue(offline["complete"])
            self.assertEqual(list((root / "hosts").glob("*.part")), [])

    def test_failed_or_corrupt_download_never_publishes_and_cleans_partial(self):
        for exit_status in (0, 22):
            with self.subTest(exit_status=exit_status), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                item = host_record()
                def download(argv, **kwargs):
                    Path(argv[argv.index("--output") + 1]).write_bytes(b"partial")
                    return subprocess.CompletedProcess(argv, exit_status, stderr="network failure")
                with patch.object(PREPARE.subprocess, "run", side_effect=download):
                    report = PREPARE.prepare_inventory([item], {}, root)
                self.assertFalse(report["complete"])
                self.assertEqual(report["error_files"], 1)
                self.assertFalse((root / item["path"]).exists())
                self.assertEqual(list((root / "hosts").iterdir()), [])

    def test_offline_missing_file_remains_error_and_parallelism_is_bounded(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(PREPARE.subprocess, "run") as network:
                report = PREPARE.prepare_inventory([host_record()], {}, Path(directory), offline=True)
                network.assert_not_called()
            self.assertEqual((report["verified_files"], report["error_files"], report["complete"]), (0, 1, False))
            with self.assertRaisesRegex(ValueError, "between 1 and 4"):
                PREPARE.prepare_inventory([], {}, Path(directory), workers=5)

    def test_output_cannot_overwrite_input_paths_or_hard_links(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, genomes = root / "metadata", root / "genomes"
            metadata.mkdir()
            source = metadata / "figshare-v1.json"
            source.write_text("source")
            with self.assertRaisesRegex(ValueError, "outside"):
                PREPARE.protect_output(source, metadata, genomes, [])
            hard_link = root / "report.json"
            os.link(source, hard_link)
            with self.assertRaisesRegex(ValueError, "hard link"):
                PREPARE.protect_output(hard_link, metadata, genomes, [])


if __name__ == "__main__":
    unittest.main()
