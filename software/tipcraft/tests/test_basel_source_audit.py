import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("basel_sources", Path(__file__).resolve().parents[3] / "scripts/audit-basel-validation-sources.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class BaselSourceTests(unittest.TestCase):
    def setUp(self):
        self.matrix = [{"strain": "host", "phage": "phage"}]
        self.phages = [{"Phage ID": "phage", "Genome Accession": "ACC.1", "Interaction matrix": ""},
                       {"Phage ID": "extra", "Genome Accession": "ACC.2", "Interaction matrix": "x"}]
        self.hosts = [{"Strain": "host", "NCBI accession": "WGS000"}, {"Strain": "extra_host", "NCBI accession": "WGS001"}]

    def test_flags_are_not_panel_roster_authority(self):
        report = AUDIT.reconcile(self.matrix, self.phages, self.hosts, [{"strain": "host", "phage": "training_phage"}])
        self.assertEqual(report["phage_id_to_accession"], {"phage": "ACC.1"})
        self.assertEqual(report["flagged_phages_absent_from_matrix"], ["extra"])
        self.assertEqual(report["matrix_phages_without_metadata_flag"], ["phage"])
        self.assertEqual(report["metadata_hosts_absent_from_matrix"], ["extra_host"])
        self.assertEqual(report["training_host_id_overlap_count"], 1)
        self.assertEqual(report["training_phage_id_overlap_count"], 0)

    def test_duplicate_metadata_ids_and_missing_exact_mappings_rejected(self):
        for phages, hosts in [(self.phages + [self.phages[0]], self.hosts),
                              (self.phages, self.hosts + [self.hosts[0]]),
                              (self.phages[1:], self.hosts), (self.phages, self.hosts[1:])]:
            with self.subTest(phages=phages, hosts=hosts), self.assertRaises(ValueError):
                AUDIT.reconcile(self.matrix, phages, hosts, self.matrix)

    def test_outcomes_do_not_select_roster(self):
        base = AUDIT.reconcile(self.matrix, self.phages, self.hosts, self.matrix)
        altered = AUDIT.reconcile([{**row, "interaction": "arbitrary_uninspected_value"} for row in self.matrix], self.phages, self.hosts, self.matrix)
        self.assertEqual(base, altered)

    def test_duplicate_pairs_rejected(self):
        with self.assertRaises(ValueError):
            AUDIT.reconcile(self.matrix + self.matrix, self.phages, self.hosts, self.matrix)
