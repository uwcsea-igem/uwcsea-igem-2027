"""Guard source-inventory authentication without network or raw datasets."""
import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("source_audit", ROOT / "scripts/audit-genophi-sources.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
EMPTY_TREE = "4b825dc642cb6eb9a060e54bf8d69288fbee4904"
EMPTY_BLOB = "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"


class SourceAuditTests(unittest.TestCase):
    def test_empty_git_tree_is_authenticated(self):
        self.assertEqual(AUDIT.validate_tree({"tree": []}, EMPTY_TREE), {})

    def test_orphan_file_cannot_bypass_root_hash(self):
        entry = {"path": "missing/fake.faa", "type": "blob", "mode": "100644", "sha": EMPTY_BLOB}
        with self.assertRaisesRegex(ValueError, "Missing directory parent"):
            AUDIT.validate_tree({"tree": [entry]}, EMPTY_TREE)

    def test_duplicate_paths_are_not_silently_collapsed(self):
        entry = {"path": "fake.faa", "type": "blob", "mode": "100644", "sha": EMPTY_BLOB}
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            AUDIT.validate_tree({"tree": [entry, entry]}, EMPTY_TREE)

    def test_incomplete_or_changed_tree_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            AUDIT.validate_tree({"tree": [], "truncated": True}, EMPTY_TREE)
        entry = {"path": "fake.faa", "type": "blob", "mode": "100644", "sha": EMPTY_BLOB}
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            AUDIT.validate_tree({"tree": [entry]}, EMPTY_TREE)
