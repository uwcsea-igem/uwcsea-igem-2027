"""Replay detects changed and internally inconsistent annotation artifacts."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[3]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
REPLAY = load("secondary_replay", ROOT/"scripts/replay-secondary-reference-annotations.py")
PRODUCER = load("secondary_original", ROOT/"scripts/annotate-basel-ecor-references.py")
PARSER = load("annotation_original", ROOT/"scripts/check-genophi-annotation.py")


def fixture(directory):
    p = Path(directory)
    source = p/"source.fna"; source.write_bytes(b">source\nATGAAATAA\n")
    normalized, mapping = PRODUCER.normalize(source.read_bytes(), "ch000", PARSER)
    dna = p/"input.fna"; dna.write_bytes(normalized)
    maps = p/"contig-map.json"; maps.write_bytes(PRODUCER.encoded(mapping))
    faa = p/"proteins.faa"; faa.write_bytes(b">ch000c000000_1 # 1 # 9 # 1 # ID=1_1;partial=00\nMK*\n")
    gff = p/"genes.gff"; gff.write_bytes(b"ch000c000000\tProdigal\tCDS\t1\t9\t1\t+\t0\tID=1_1;partial=00;\n")
    log = p/"caller.log"; log.write_bytes(b"")
    summary, inventory = PRODUCER.validate_outputs(faa.read_bytes(), gff.read_bytes(), {"ch000c000000":9}, PARSER)
    coordinates = p/"protein-coordinate-inventory.json"; coordinates.write_bytes(PRODUCER.encoded(inventory))
    item = {"id":"ch000","kind":"host","source_id":"control","role":"technical_control", "fasta":PRODUCER.binding(source),
            "expected_faa":PRODUCER.binding(faa),"expected_proteome_multiset_sha256":summary["proteome_multiset_sha256"]}
    prepared = {"item":item,"status":"input_prepared","directory":str(p),"input":PRODUCER.binding(dna),
                "contig_map":PRODUCER.binding(maps),"input_contigs":1,"input_bases":9}
    entry = {**prepared,**summary,"status":"annotation_validated","exit_code":0,"control_exact_proteome_multiset_match":True,
             "outputs":{"FAA":PRODUCER.binding(faa),"GFF":PRODUCER.binding(gff),"log":PRODUCER.binding(log)},
             "protein_coordinate_inventory":PRODUCER.binding(coordinates)}
    (p/"receipt.json").write_bytes(PRODUCER.encoded(entry))
    return entry, prepared


class SecondaryReplayTests(unittest.TestCase):
    def test_complete_read_only_replay_preserves_artifact_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            entry, prepared = fixture(directory)
            before = {p:p.read_bytes() for p in Path(directory).iterdir()}
            summary, contigs, proteins = REPLAY.verify_entry(entry,prepared,PRODUCER,[PARSER,None,None],REPLAY.Artifacts())
            self.assertEqual(summary["proteins"],1)
            self.assertEqual(len(contigs),1)
            self.assertEqual(len(proteins),1)
            self.assertEqual(before,{p:p.read_bytes() for p in Path(directory).iterdir()})

    def test_changed_source_bytes_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            entry, prepared = fixture(directory)
            Path(entry["item"]["fasta"]["path"]).write_bytes(b">source\nATGCCCTAA\n")
            with self.assertRaisesRegex(ValueError,"Artifact bytes differ"):
                REPLAY.verify_entry(entry,prepared,PRODUCER,[PARSER,None,None],REPLAY.Artifacts())

    def test_rehashed_coordinate_inventory_still_requires_actual_output_agreement(self):
        with tempfile.TemporaryDirectory() as directory:
            entry, prepared = fixture(directory)
            path = Path(entry["protein_coordinate_inventory"]["path"])
            rows = json.loads(path.read_bytes()); rows[0]["start_1based"] = 2
            path.write_bytes(PRODUCER.encoded(rows))
            entry["protein_coordinate_inventory"] = PRODUCER.binding(path)
            with self.assertRaisesRegex(ValueError,"Protein-coordinate inventory differs"):
                REPLAY.verify_entry(entry,prepared,PRODUCER,[PARSER,None,None],REPLAY.Artifacts())

    def test_rehashed_contig_map_does_not_hide_changed_boundaries(self):
        with tempfile.TemporaryDirectory() as directory:
            entry, prepared = fixture(directory)
            path = Path(entry["contig_map"]["path"]); rows=json.loads(path.read_bytes());rows[0]["length"] = 8
            path.write_bytes(PRODUCER.encoded(rows)); spec=PRODUCER.binding(path)
            entry["contig_map"]=spec; prepared["contig_map"]=spec
            with self.assertRaisesRegex(ValueError,"neutral contig mapping differs"):
                REPLAY.verify_entry(entry,prepared,PRODUCER,[PARSER,None,None],REPLAY.Artifacts())

    def test_pause_timing_is_preserved_without_adjustment(self):
        receipt={"action":"resume_existing_annotation_child_then_existing_wrapper", "no_restart_or_new_caller_submission":True,
                 "wall_time_in_suspended_call_includes_user_pause":True,"resumed_epoch":1}
        entry={"item":{"id":"h013","kind":"host","source_id":"ECOR25"},"resources":{"wall_seconds":43206.2,"child_user_cpu_seconds":5.4}}
        result=REPLAY.timing_evidence([entry],receipt,{"path":"receipt","sha256":"bound"})
        self.assertEqual(result["affected_call"]["wall_seconds"],43206.2)
        self.assertFalse(result["wall_time_is_compute_benchmark"])
        self.assertFalse(result["adjusted_wall_time_calculated"])


if __name__ == "__main__": unittest.main()
