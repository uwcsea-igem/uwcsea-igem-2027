"""Exercise a manually downloaded archive importer using tiny synthetic ZIPs."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile


SPEC = importlib.util.spec_from_file_location(
    "host_import", Path(__file__).resolve().parents[3] / "scripts/import-genophi-hosts.py")
IMPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPORT)
FASTA = b">original source contig\nACGTRN\n>another\nttac\n"
PROVENANCE = {"reviewed_host_file_inventory_sha256": "reviewed-inventory-hash",
              "metadata_files": [{"path": "figshare-v1.json", "sha256": "metadata-hash"}]}


def host_record(name="original-host.fasta", file_id=123, data=FASTA):
    return {"kind": "host", "source_id": Path(name).stem, "source_name": name,
            "deposit_file_id": file_id, "path": f"hosts/{file_id}.fasta",
            "expected_bytes": len(data), "expected_md5": hashlib.md5(data).hexdigest()}


def archive_at(root, entries):
    archive = root / "download.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for name, data in entries:
            bundle.writestr(name, data)
    return archive


class HostImportTests(unittest.TestCase):
    def test_original_filename_maps_to_file_id_and_binds_archive_and_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = host_record()
            archive = archive_at(root, [("wrapper/nested/" + item["source_name"], FASTA),
                                        ("wrapper/README.txt", b"metadata")])
            report = IMPORT.import_archive(archive, root / "genomes", [item], PROVENANCE)
            self.assertTrue(report["complete"])
            self.assertEqual((root / "genomes/hosts/123.fasta").read_bytes(), FASTA)
            self.assertFalse((root / "genomes" / item["source_name"]).exists())
            self.assertEqual(report["archive"]["sha256"], hashlib.sha256(archive.read_bytes()).hexdigest())
            self.assertEqual(report["archive"]["bytes"], archive.stat().st_size)
            self.assertEqual(report["metadata_files"], PROVENANCE["metadata_files"])
            self.assertEqual(report["reviewed_host_file_inventory_sha256"], "reviewed-inventory-hash")
            self.assertEqual(report["ignored_non_fasta_members"],
                             [{"member": "wrapper/README.txt", "uncompressed_bytes": 8}])
            self.assertFalse(report["training_or_scoring_executed"])
            self.assertEqual(list((root / "genomes/hosts").glob("*.part")), [])

    def test_duplicate_expected_basename_is_not_selected_from_either_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = host_record()
            archive = archive_at(root, [("a/" + item["source_name"], FASTA),
                                        ("b/" + item["source_name"], FASTA)])
            report = IMPORT.import_archive(archive, root / "genomes", [item], {})
            self.assertFalse(report["complete"])
            self.assertEqual(report["error_files"], 1)
            self.assertIn("Duplicate expected filename", report["files"][0]["error"])
            self.assertFalse((root / "genomes/hosts/123.fasta").exists())

    def test_unexpected_fasta_rejected_but_safe_known_member_can_import(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = host_record()
            archive = archive_at(root, [(item["source_name"], FASTA), ("unknown.FNA", FASTA)])
            report = IMPORT.import_archive(archive, root / "genomes", [item], {})
            self.assertEqual(report["imported_files"], 1)
            self.assertFalse(report["complete"])
            self.assertEqual(report["archive_errors"][0]["member"], "unknown.FNA")
            self.assertEqual(list((root / "genomes/hosts").iterdir()), [root / "genomes/hosts/123.fasta"])

    def test_unsafe_paths_and_symlinks_are_rejected_even_for_metadata(self):
        unsafe = ("../original-host.fasta", "/original-host.fasta", "C:/original-host.fasta",
                  "wrapper\\original-host.fasta", "a//README.txt", "a/./README.txt", "a//")
        for name in unsafe:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                archive = archive_at(root, [(name, FASTA)])
                report = IMPORT.import_archive(archive, root / "genomes", [host_record()], {})
                self.assertFalse(report["complete"])
                self.assertEqual(report["verified_files"], 0)
                self.assertIn("Unsafe ZIP", report["archive_errors"][0]["error"])
                self.assertFalse((root / "original-host.fasta").exists())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            link = zipfile.ZipInfo("original-host.fasta")
            link.create_system = 3
            link.external_attr = (stat.S_IFLNK | 0o777) << 16
            archive = archive_at(root, [(link, b"outside.fasta")])
            report = IMPORT.import_archive(archive, root / "genomes", [host_record()], {})
            self.assertIn("links and special files", report["archive_errors"][0]["error"])
            self.assertFalse((root / "genomes/hosts/123.fasta").exists())

    def test_wrong_declared_size_is_rejected_before_member_is_opened(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = archive_at(root, [("original-host.fasta", FASTA * 100)])
            with patch.object(zipfile.ZipFile, "open", side_effect=AssertionError("Must not decompress")):
                report = IMPORT.import_archive(archive, root / "genomes", [host_record()], {})
            self.assertIn("size differs", report["files"][0]["error"])
            self.assertFalse((root / "genomes/hosts/123.fasta").exists())

    def test_bad_md5_and_bad_fasta_never_publish_and_clean_temporaries(self):
        malformed = b">bad\nACZ\n"
        cases = [(host_record(), FASTA.replace(b"ACGTRN", b"ACGTNN"), "MD5"),
                 (host_record(data=malformed), malformed, "non-IUPAC")]
        for item, data, error in cases:
            with self.subTest(error=error), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                archive = archive_at(root, [(item["source_name"], data)])
                report = IMPORT.import_archive(archive, root / "genomes", [item], {})
                self.assertFalse(report["complete"])
                self.assertIn(error, report["files"][0]["error"])
                self.assertEqual(list((root / "genomes/hosts").iterdir()), [])

    def test_preserves_verified_existing_inode_and_never_overwrites_corrupt_cache(self):
        for existing, success in [(FASTA, True), (FASTA.replace(b"ACGTRN", b"ACGTNN"), False)]:
            with self.subTest(success=success), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                archive = archive_at(root, [("original-host.fasta", FASTA)])
                destination = root / "genomes/hosts/123.fasta"
                destination.parent.mkdir(parents=True)
                destination.write_bytes(existing)
                inode = destination.stat().st_ino
                report = IMPORT.import_archive(archive, root / "genomes", [host_record()], {})
                self.assertEqual(report["complete"], success)
                self.assertEqual(destination.stat().st_ino, inode)
                self.assertEqual(destination.read_bytes(), existing)
                self.assertEqual(list(destination.parent.glob("*.part")), [])
                if success:
                    self.assertEqual((report["preserved_existing_files"], report["imported_files"]), (1, 0))
                else:
                    self.assertIn("MD5", report["files"][0]["error"])

    def test_partial_archive_has_explicit_missing_list_even_if_cache_is_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first, second = host_record(), host_record("second.fasta", 456)
            archive = archive_at(root, [(first["source_name"], FASTA)])
            report = IMPORT.import_archive(archive, root / "genomes", [first, second], {})
            self.assertFalse(report["complete"])
            self.assertEqual(report["missing_files"], ["hosts/456.fasta"])
            self.assertEqual(report["archive_missing_expected_names"], ["second.fasta"])
            (root / "genomes/hosts/456.fasta").write_bytes(FASTA)
            resumed = IMPORT.import_archive(archive, root / "genomes", [first, second], {})
            self.assertTrue(resumed["cache_complete"])
            self.assertFalse(resumed["complete"])
            self.assertEqual(resumed["missing_files"], [])
            self.assertEqual(resumed["archive_missing_expected_names"], ["second.fasta"])

    def test_destination_symlink_and_racing_existing_file_are_never_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = archive_at(root, [("original-host.fasta", FASTA)])
            destination = root / "genomes/hosts/123.fasta"
            destination.parent.mkdir(parents=True)
            outside = root / "outside.fasta"
            outside.write_bytes(FASTA)
            destination.symlink_to(outside)
            report = IMPORT.import_archive(archive, root / "genomes", [host_record()], {})
            self.assertIn("symlink", report["files"][0]["error"])
            self.assertEqual(outside.read_bytes(), FASTA)
            destination.unlink()
            def race(source, target):
                Path(target).write_bytes(b"racing unrelated file")
                raise FileExistsError()
            with patch.object(IMPORT.os, "link", side_effect=race):
                report = IMPORT.import_archive(archive, root / "genomes", [host_record()], {})
            self.assertFalse(report["complete"])
            self.assertEqual(destination.read_bytes(), b"racing unrelated file")
            self.assertEqual(list(destination.parent.glob("*.part")), [])

    def test_report_output_cannot_overlap_inputs_or_their_hard_links(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, genomes = root / "metadata", root / "genomes"
            metadata.mkdir()
            archive = archive_at(root, [])
            (metadata / "figshare-v1.json").write_text("cached input")
            for output in (metadata / "report.json", genomes / "report.json", archive, root):
                with self.subTest(output=output), self.assertRaises(ValueError):
                    IMPORT.protect_paths(output, metadata, genomes, archive, [host_record()])
            for source in (archive, metadata / "figshare-v1.json"):
                output = root / "report.json"
                os.link(source, output)
                with self.assertRaises(ValueError):
                    IMPORT.protect_paths(output, metadata, genomes, archive, [host_record()])
                output.unlink()

    def test_cli_authenticates_inventory_and_writes_explicit_incomplete_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = archive_at(root, [("original-host.fasta", FASTA)])
            output = root / "processed/report.json"
            argv = ["import-genophi-hosts.py", "--metadata-cache", str(root / "metadata"),
                    "--archive", str(archive), "--genomes", str(root / "genomes"), "--output", str(output)]
            with patch.object(sys, "argv", argv), patch.object(
                    IMPORT.PREPARE, "load_inventory", return_value=(
                        [host_record(), host_record("second.fasta", 456)], PROVENANCE)) as inventory:
                self.assertEqual(IMPORT.main(), 1)
            inventory.assert_called_once_with(root / "metadata")
            report = json.loads(output.read_text())
            self.assertEqual(report["imported_files"], 1)
            self.assertEqual(report["archive_missing_expected_names"], ["second.fasta"])
            self.assertEqual(list(output.parent.glob("*.part")), [])

    def test_inventory_ambiguity_is_rejected_before_archive_access(self):
        with self.assertRaisesRegex(ValueError, "Duplicate host inventory"):
            IMPORT.host_inventory([host_record(), host_record(file_id=456)])
        for change in ({"source_name": "../host.fasta"}, {"path": "../123.fasta"}, {"expected_bytes": 0}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                IMPORT.host_inventory([{**host_record(), **change}])


if __name__ == "__main__":
    unittest.main()
