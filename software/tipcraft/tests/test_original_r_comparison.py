"""Provenance and invalid-output regressions; no R runtime or network needed."""
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "original_r_v2", Path(__file__).resolve().parents[3]/"scripts/compare-viridic-original-r-v2.py")
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


def write_tsv(path, fields, rows):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream,fieldnames=fields,delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class OriginalRComparisonTests(unittest.TestCase):
    def runtime_fixture(self, root):
        home, library = root/"r-home",root/"library"
        (home/"bin/exec").mkdir(parents=True)
        (home/"lib").mkdir()
        (home/"bin/exec/R").write_bytes(b"authenticated-R")
        (home/"lib/libR.dylib").write_bytes(b"authenticated-libR")
        files = {"Package/DESCRIPTION":b"Package: Package\nVersion: 1.0\n",
                 "Package/R/function.R":b"function() 1\n"}
        archive=root/"packages/Package_1.0.tgz"
        archive.parent.mkdir()
        with tarfile.open(archive,"w:gz") as bundle:
            for name,data in files.items():
                item=tarfile.TarInfo(name);item.size=len(data)
                bundle.addfile(item,io.BytesIO(data))
                path=library/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        manifest={"runtime_root":str(root),"r_executable_sha256":COMPARE.sha(home/"bin/exec/R"),
                  "libR_sha256":COMPARE.sha(home/"lib/libR.dylib"),
                  "packages":[{"name":"Package","version":"1.0","md5_verified":True,
                               "sha256":COMPARE.sha(archive)}]}
        return home,library,manifest

    def test_supplied_runtime_bytes_cannot_hide_behind_valid_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            home,library,manifest=self.runtime_fixture(Path(directory))
            self.assertEqual(COMPARE.validate_runtime(home,library,manifest)["installed_regular_files_verified"],2)
            (home/"bin/exec/R").write_bytes(b"different-R")
            with self.assertRaisesRegex(ValueError,"Actual supplied R runtime"):
                COMPARE.validate_runtime(home,library,manifest)

    def test_installed_package_content_and_description_are_authenticated(self):
        for field in ("content","version"):
            with self.subTest(field=field),tempfile.TemporaryDirectory() as directory:
                home,library,manifest=self.runtime_fixture(Path(directory))
                if field=="content":
                    (library/"Package/R/function.R").write_bytes(b"function() 2\n")
                else:
                    (library/"Package/DESCRIPTION").write_bytes(b"Package: Package\nVersion: 2.0\n")
                with self.assertRaisesRegex(ValueError,"Installed R package"):
                    COMPARE.validate_runtime(home,library,manifest)

    def test_changed_package_archive_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);home,library,manifest=self.runtime_fixture(root)
            (root/"packages/Package_1.0.tgz").write_bytes(b"changed archive")
            with self.assertRaisesRegex(ValueError,"archive changed"):
                COMPARE.validate_runtime(home,library,manifest)

    def test_na_and_nonfinite_similarity_cannot_skip_numeric_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/"lengths.json").write_text(json.dumps({"A":10,"B":10}))
            (root/"adapter.json").write_text(json.dumps({"pairs":[{
                "query":"A","subject":"B","query_nident":3,"subject_nident":3,
                "query_original_aligned_bases":3,"subject_original_aligned_bases":3,
                "similarity_percent_unrounded":30,"similarity_percent":30}]}))
            fields=["query","subject","qlen","slen","nident","original_aligned_query_bases","error"]
            rows=[dict(zip(fields,[a,b,10,10,n,n,""])) for a,b,n in
                  [("A","A",10),("B","B",10),("A","B",3),("B","A",3)]]
            write_tsv(root/"directional.tsv",fields,rows)
            policy={"inputs":{"genome_lengths":{"path":str(root/"lengths.json")},
                              "adapter_output":{"path":str(root/"adapter.json")}},
                    "scope":{"distinct_unordered_pairs":1},
                    "agreement_requirements":{"unrounded_similarity_absolute_tolerance_percentage_points":1e-12}}
            fields=["query","subject","similarity_percent_unrounded","similarity_percent"]
            write_tsv(root/"pairs.tsv",fields,[dict(zip(fields,["A","B","30","30"]))])
            self.assertEqual(COMPARE.compare(policy,root)["status"],"passed")
            for bad in ("NA","NaN","Inf",""):
                with self.subTest(bad=bad):
                    write_tsv(root/"pairs.tsv",fields,[dict(zip(fields,["A","B",bad,"30"]))])
                    result=COMPARE.compare(policy,root)
                    self.assertEqual(result["status"],"discrepancy")
                    self.assertFalse(result["frozen_dataset_directional_equivalence_verified"])
                    self.assertTrue(result["discrepancies"])


if __name__=="__main__":
    unittest.main()
