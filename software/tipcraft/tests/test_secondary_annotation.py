"""DNA preservation, complete rosters, output correspondence and failed calls."""
import copy
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[3]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
CALL = load("secondary_annotation", ROOT/"scripts/annotate-basel-ecor-references.py")
PARSER = load("prior_annotation", ROOT/"scripts/check-genophi-annotation.py")
FAA = b">h000c000000_1 # 1 # 9 # 1 # ID=1_1;partial=00\nMK*\n"
GFF = b"##gff-version 3\nh000c000000\tProdigal\tCDS\t1\t9\t1.0\t+\t0\tID=1_1;partial=00;\n"


class SecondaryAnnotationTests(unittest.TestCase):
    def test_neutralization_preserves_boundaries_order_and_copies(self):
        data, rows = CALL.normalize(b">old\nACGT\n>other\nACGT\n", "h000", PARSER)
        self.assertEqual(PARSER.records(data), [("h000c000000", "ACGT"), ("h000c000001", "ACGT")])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["sequence_sha256"], rows[1]["sequence_sha256"])
        self.assertEqual([r["source_record_id"] for r in rows], ["old", "other"])
        with self.assertRaises(ValueError): CALL.normalize(b">bad\nAC?T\n", "h000", PARSER)

    def test_complete_pair_inventory_and_zero_CDS_contigs(self):
        summary, proteins = CALL.validate_outputs(FAA, GFF, {"h000c000000":9,"h000c000001":3}, PARSER)
        self.assertEqual(summary["proteins"], 1)
        self.assertEqual(summary["contigs_with_no_CDS"], ["h000c000001"])
        self.assertEqual(proteins[0]["length"], 2)
        self.assertTrue(proteins[0]["terminal_stop_removed_for_multiset_only"])

    def test_coordinate_strand_missing_extra_and_duplicate_GFF_fail(self):
        variants = [GFF.replace(b"\t1\t9\t",b"\t1\t8\t"), GFF.replace(b"\t+\t",b"\t-\t"),
                    b"##gff-version 3\n", GFF+GFF.splitlines(True)[1], GFF.replace(b"ID=1_1",b"ID=2_1")]
        for data in variants:
            with self.subTest(data=data), self.assertRaises(ValueError):
                CALL.validate_outputs(FAA,data,{"h000c000000":9},PARSER)

    def test_invalid_faa_and_out_of_bounds_are_rejected(self):
        for faa in [FAA.replace(b"MK*",b"M*K"), FAA+FAA, FAA.replace(b"# 9 #",b"# 10 #")]:
            with self.assertRaises(ValueError): CALL.validate_outputs(faa,GFF,{"h000c000000":9},PARSER)

    def test_all_authority_ids_retained_not_selected_by_status(self):
        host=[{"id":f"h{i:03d}","kind":"host","source_id":f"H{i}","role":"secondary_reference_not_admitted"} for i in range(25)]
        phage=[{"id":f"p{i:03d}","kind":"phage","source_id":f"P{i}","role":"secondary_reference_not_admitted"} for i in range(52)]
        controls=[{"id":"ch000","kind":"host","role":"technical_control"},{"id":"cp000","kind":"phage","role":"technical_control"}]
        authority={"host_id_to_wgs_master_accession":{r['source_id']:"meta" for r in host},"phage_id_to_accession":{r['source_id']:"meta" for r in phage}}
        CALL.validate_roster(host+phage,controls,authority)
        bad=copy.deepcopy(host+phage);bad[1]['source_id']=bad[0]['source_id']
        with self.assertRaises(ValueError): CALL.validate_roster(bad,controls,authority)
        with self.assertRaises(ValueError): CALL.validate_roster(host+phage[:-1],controls,authority)

    def test_failed_child_keeps_failure_no_retry_and_input_unmodified(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);fasta=p/'input.fna';fasta.write_bytes(b">h000c000000\nATGAAATAA\n")
            toolpath=p/'caller';toolpath.write_bytes(b"tool")
            mapping=p/'contig-map.json';mapping.write_bytes(CALL.encoded([{'neutral_id':'h000c000000','length':9}]))
            source=p/'source.fna';source.write_bytes(fasta.read_bytes()); before=fasta.read_bytes()
            entry={'status':'input_prepared','directory':str(p),'input':CALL.binding(fasta),'contig_map':CALL.binding(mapping),
                   'item':{'id':'h000','kind':'host','role':'secondary_reference_not_admitted','source_id':'test','fasta':CALL.binding(source)}}
            with mock.patch.object(CALL.subprocess,'run',return_value=SimpleNamespace(returncode=9)) as run:
                result=CALL.call(entry,{'executable':str(toolpath),'executable_sha256':CALL.sha(b'tool')},CALL.FLAGS['host'],180,PARSER)
            self.assertEqual(run.call_count,1)
            self.assertEqual(result['status'],'caller_failed')
            self.assertEqual(result['exit_code'],9)
            self.assertEqual(fasta.read_bytes(),before)
            self.assertTrue((p/'receipt.json').exists())


if __name__ == '__main__': unittest.main()
