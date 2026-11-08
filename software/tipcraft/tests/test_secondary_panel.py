"""Membership and partition guards use synthetic source metadata only."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("secondary_panel", ROOT / "scripts/prepare-secondary-panel.py")
PANEL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PANEL)


def fixture():
    hosts = ["h0", "h1", "h2", "h3", "hx"]
    phages = ["p0", "p1", "p2", "p3", "px"]
    authority = {"host_id_to_wgs_master_accession": {h: "WGS" + h for h in hosts},
                 "phage_id_to_accession": {p: "ACC" + p for p in phages},
                 "panel_pairs": 23, "assay_source": "synthetic-source"}
    rows = [{"strain": h, "phage": p, "interaction": str((i+j) % 2)}
            for i, h in enumerate(hosts) for j, p in enumerate(phages)
            if (h, p) not in [("h0", "p0"), ("h1", "p1")]]
    hg = {"group-" + h: {"members": ["ecor:" + h] + (["original:target"] if h == "hx" else []),
                           "source_members": ["ecor:" + h], "target_containing": h == "hx"} for h in hosts}
    host = {"partition": {"status": "complete_host_reference_grouping", "source_groups": hg,
                          "target_related_source_hosts": ["ecor:hx"], "candidate_source_hosts": ["ecor:" + h for h in hosts[:-1]],
                          "candidate_source_group_count": 4},
            "comparison_summary": {"status": "complete_reference_triage", "failed_self_controls": [],
                                   "missing_directed_pairs": [], "directed_pair_count": 36}, "prior25_control_exact": True}
    auth = {"unresolved_hosts": [], "entries": [{"host": h, "source_wgs_accession": "WGS" + h,
            "status": "authenticated_official_wgs_contig_set", "master": {"wgs_master": "WGS" + h + "01", "wgs_build": "01"},
            "fasta": {"path": h + ".fna", "sha256": h}} for h in hosts]}
    observations = [{"phage": p, "requested_accession": "ACC" + p, "official_formats_agree": True,
                     "official_versioned_accession": "ACC" + p + ".1", "exact_author_DNA_match": p != "p3",
                     "official_length": 100, "official_sequence_sha256": p} for p in phages]
    node_map = {p: ("source:" + p if p != "p3" else "official:p3:ACCp3.1") for p in phages}
    containers = [{"container_key": n, "source_id": p, "role": "external_collection" if p != "p3" else "official_reference_version",
                   "length": 100} for p, n in node_map.items()]
    groups = {"group-" + p: {"members": [n], "target_seeds": [], "target_only_excluded": False,
                            "direct_gp17_seeds": [n] if p == "px" else [], "primary_excluded": p == "px"}
              for p, n in node_map.items()}
    phage = {"status": "completed_original_R_checked_official_variant_groups", "numeric_matrix_for_grouping_verified": True,
             "numeric_validation": {"dataset_specific_directional_equivalence_verified": True}, "official_panel_node_map": node_map,
             "containers": containers, "grouping": {"all_nodes_retained_as_bridges": True, "groups": groups,
             "container_groups": {n: "group-" + p for p, n in node_map.items()}}}
    replay = {"observations": observations, "official_record_format_pairs_verified": 5}
    design = {"partition": {"host_folds": 3, "phage_folds": 3, "seed": 2027, "regime": "both_unseen", "mixed_quadrants": "embargo"}}
    return rows, authority, host, phage, auth, replay, design


class SecondaryPanelTests(unittest.TestCase):
    def test_exact_official_variants_exclusions_unknowns_and_disjoint_blocks(self):
        result = PANEL.build_panel(*fixture())
        self.assertEqual(result["counts"]["source_observed_rows"], 23)
        self.assertEqual(result["counts"]["source_unknown_pairs"], 2)
        self.assertEqual(result["counts"]["eligible_observed_rows"], 14)
        self.assertEqual(result["phage_genomes"]["p3"]["genotype_key"], "official:p3:ACCp3.1")
        self.assertFalse(result["experimental_stock_identity_verified"])
        self.assertFalse(result["training_authorized"])
        self.assertEqual(len(result["blocks"]), 9)
        self.assertTrue(all("interaction" not in r and "label" not in r for r in result["samples"]))
        seen = []
        for block in result["blocks"]:
            roles = [set(block[k]) for k in ["train_ids", "test_ids", "embargo_ids"]]
            self.assertFalse(roles[0] & roles[1] or roles[1] & roles[2] or roles[0] & roles[2])
            self.assertFalse(block["validation"]["overlaps"]["host_group"])
            self.assertFalse(block["validation"]["overlaps"]["phage_group"])
            seen += block["test_ids"]
        self.assertEqual(sorted(seen), sorted(r["row_id"] for r in result["samples"]))
        both = next(r for r in result["memberships"] if r["host_id"] == "hx" and r["phage_id"] == "px")
        self.assertEqual(len(both["exclusion_reasons"]), 2)
        self.assertTrue(all(u["status"] == "unknown_not_negative" for u in result["unknown_pairs"]))

    def test_labels_cannot_change_membership_or_fold_assignment(self):
        args = fixture()
        original = PANEL.build_panel(*args)
        permuted_labels = list(reversed([row["interaction"] for row in args[0]]))
        for row, label in zip(args[0], permuted_labels):
            row["interaction"] = label
        permuted = PANEL.build_panel(*args)
        for key in ["samples", "memberships", "blocks", "host_group_folds", "phage_group_folds", "unknown_pairs"]:
            self.assertEqual(original[key], permuted[key], key)
        for row in args[0]:
            row["interaction"] = "0"
        changed = PANEL.build_panel(*args)
        for key in ["samples", "memberships", "blocks", "host_group_folds", "phage_group_folds", "unknown_pairs"]:
            self.assertEqual(original[key], changed[key], key)
        self.assertEqual(changed["status"], "HOLD_single_class_training_blocks")
        self.assertEqual(len(changed["invalid_training_blocks"]), 9)

    def test_incomplete_unchecked_or_name_substituted_groups_fail(self):
        for mutation in [lambda a: a[3].update(status="HOLD"),
                         lambda a: a[3].update(numeric_matrix_for_grouping_verified=False),
                         lambda a: a[3]["numeric_validation"].update(dataset_specific_directional_equivalence_verified=False),
                         lambda a: a[3]["official_panel_node_map"].update(p3="source:p3"),
                         lambda a: a[3]["grouping"]["groups"]["group-p0"].update(primary_excluded=True),
                         lambda a: a[2]["partition"].update(status="HOLD"),
                         lambda a: a[2]["comparison_summary"].update(directed_pair_count=35),
                         lambda a: a[2]["partition"]["source_groups"]["group-h0"].update(target_containing=True)]:
            args = fixture()
            mutation(args)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                PANEL.build_panel(*args)

    def test_duplicate_missing_nonbinary_and_nonfinite_observations_fail(self):
        for mutation in [lambda a: a[0].append(deepcopy(a[0][0])), lambda a: a[0].pop(),
                         lambda a: a[0][0].update(interaction="NaN"), lambda a: a[0][0].update(interaction="-1"),
                         lambda a: a[0][0].update(strain="unmapped")]:
            args = fixture()
            mutation(args)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                PANEL.build_panel(*args)

    def test_annotation_join_never_admits_technical_control_or_wrong_genotype(self):
        result = PANEL.build_panel(*fixture())
        entries = [{"item": {"role": "technical_control", "source_id": "T7", "kind": "phage"}}]
        for kind, field in [("host", "host_genomes"), ("phage", "phage_genomes")]:
            for name, genome in result[field].items():
                item = {"role": "secondary_reference_not_admitted", "kind": kind, "source_id": name, "id": kind+name,
                        "wgs_accession": genome["source_accession"], "fasta": genome.get("fasta"),
                        "versioned_accession": genome.get("versioned_accession"), "expected_sequence_sha256": genome.get("sequence_sha256")}
                entries.append({"item": item, "status": "annotation_validated", "FAA_GFF_inventory_and_coordinates_agree": True,
                                "proteins": 2, "input": {}, "contig_map": {}, "protein_coordinate_inventory": {},
                                "outputs": {"FAA": {}, "GFF": {}}, "proteome_multiset_sha256": name})
        with patch.object(PANEL, "checked"):
            PANEL.join_annotations(result, {"entries": entries})
            self.assertTrue(result["annotation_join_complete"])
            self.assertNotIn("T7", result["phage_genomes"])
            entries[-1]["item"]["versioned_accession"] = "wrong.1"
            with self.assertRaisesRegex(ValueError, "genotype"):
                PANEL.join_annotations(result, {"entries": entries})

    def test_missing_group_prerequisite_check_only_emits_no_memberships(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "design.json").write_text(json.dumps({"bound_design_inputs": []}))
            policy = {"helper_sha256": PANEL.sha(PANEL.__file__),
                      "inputs": {"design": {"path": "design.json", "sha256": PANEL.sha(root / "design.json")}},
                      "future_phage_groups": {"path": "not-yet-created.json"}}
            path = root / "policy.json"
            path.write_text(json.dumps(policy))
            with patch.object(PANEL, "ROOT", root):
                result = PANEL.prepare(path, check_only=True)
            self.assertEqual(result["status"], "HOLD_missing_prerequisite")
            self.assertFalse(result["membership_generated"])
            self.assertNotIn("samples", result)
            self.assertNotIn("labels", result)
            self.assertFalse((root / "not-yet-created.json").exists())


if __name__ == "__main__":
    unittest.main()
