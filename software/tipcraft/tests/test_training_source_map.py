import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("training_source_map", Path(__file__).resolve().parents[3] / "scripts/map-genophi-training-sources.py")
MAPPING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MAPPING)


class TrainingSourceMapTests(unittest.TestCase):
    def test_article_id_and_version_do_not_authenticate_changed_inventory(self):
        deposit = {"id": 25941691, "version": 1, "files": [
            {"id": 46741567, "name": "BL21.fasta", "size": 1,
             "download_url": "https://example.invalid/changed", "computed_md5": "changed"}]}
        with self.assertRaisesRegex(ValueError, "inventory differs"):
            MAPPING.validate_deposit(deposit)

    def test_exact_match_keeps_file_identity(self):
        record = {"id": "BL21", "url": "published-file"}
        result = MAPPING.match_ids(["BL21"], [record])
        self.assertEqual(result["BL21"]["file"], record)

    def test_unique_normalized_match_is_not_accepted_as_identity(self):
        result = MAPPING.match_ids(["025010"], [{"id": "025-010"}])
        self.assertIsNone(result["025010"]["file"])
        self.assertEqual(result["025010"]["unaccepted_alias_candidates"], ["025-010"])

    def test_ambiguous_aliases_and_duplicates_are_not_silently_resolved(self):
        records = [{"id": "A-1"}, {"id": "A_1"}]
        self.assertEqual(MAPPING.match_ids(["A1"], records)["A1"]["unaccepted_alias_candidates"], ["A-1", "A_1"])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            MAPPING.match_ids(["A-1"], [records[0], records[0]])
