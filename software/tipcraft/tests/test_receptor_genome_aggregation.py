"""Contig offsets preserve disjoint HSP domains and conservative grouping."""
import importlib.util
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("receptor_aggregate", ROOT/"scripts/aggregate-receptor-genomes.py")
AGG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AGG)


def mapping():
    records = {"a1": {"container_key": "source:A", "retained_bases": 3, "offset": 0},
               "a2": {"container_key": "source:A", "retained_bases": 3, "offset": 3},
               "b1": {"container_key": "source:B", "retained_bases": 4, "offset": 0},
               "b2": {"container_key": "source:B", "retained_bases": 4, "offset": 4}}
    containers = {"source:A": {"length": 6}, "source:B": {"length": 8}}
    return records, containers


def row(q, s, qs=1, qe=3, ss=4, se=1):
    return f"{q}\t{s}\t1e-9\t20\t3\t4\t{qs}\t{qe}\t{ss}\t{se}\tA-CG\tATCG\t3\t1\n"


def pairs(containers, links):
    result = []
    for a,b in itertools.combinations(sorted(containers),2):
        # Length10000 makes symmetric counts equal the similarity*100 exactly.
        n = int(round(links.get((a,b),0)*100))
        result.append({"query":a,"subject":b,"query_length":10000,"subject_length":10000,
                       "query_nident":n,"subject_nident":n,"similarity_percent_unrounded":100*(n+n)/20000})
    return result


class AggregationTests(unittest.TestCase):
    def test_record_map_is_exact_order_independent_and_namespaced(self):
        sources=[{'source':{'primary_id':'T','cohort_member':True,'published_model_role':'mutant'},
                  'nucleotide_status':'available','records':[{'neutral_nucleotide_id':'a1','record_id':'one','length':3,'sequence_sha256':'x','ambiguous_bases':0},
                  {'neutral_nucleotide_id':'a2','record_id':'two','length':4,'sequence_sha256':'y','ambiguous_bases':1}]}]
        targets=[{'id':'t000','source_id':'T','role':'target_reference','retained_bases':5,'accession':'ref','sequence_sha256':'z'}]
        observed=[{'id':identifier,'container_key':'source:T','source_id':'T','role':'external_collection',
                   'retained_bases':length,'source_record_id':original} for identifier,length,original in [('a2',4,'two'),('a1',3,'one')]]
        observed.append({'id':'t000','container_key':'reference:T','source_id':'T','role':'target_reference','retained_bases':5,'source_record_id':'ref'})
        records, containers=AGG.build_mapping(sources,targets,observed)
        self.assertEqual(set(containers),{'source:T','reference:T'})
        self.assertEqual(records['a2']['offset'],3)
        self.assertEqual(containers['source:T']['length'],7)
        for bad in [observed[:-1], observed+[observed[0]], [{**observed[0],'retained_bases':5},*observed[1:]]]:
            with self.assertRaises(ValueError): AGG.build_mapping(sources,targets,bad)

    def test_offsets_preserve_order_reverse_subject_and_gap_columns(self):
        records, containers = mapping()
        original = [row("a2","b2"),row("a1","b1")]
        transformed = list(AGG.transform_rows(original,records,containers))
        self.assertEqual(transformed[0].split("\t")[4:10],["6","8","4","6","8","5"])
        self.assertEqual(transformed[1].split("\t")[4:10],["6","8","1","3","4","1"])
        self.assertTrue(all(r.endswith("A-CG\tATCG\t3\t1\n") for r in transformed))

    def test_coordinate_cannot_cross_a_source_record(self):
        records, containers = mapping()
        for original in [row("a1","b1",qe=4),row("a1","b1",ss=5),row("unknown","b1")]:
            with self.assertRaises(ValueError):
                list(AGG.transform_rows([original],records,containers))

    def test_external_calculator_counts_disjoint_contigs_without_collision(self):
        records, containers = mapping()
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            (p/"bridge.py").write_text(AGG.BRIDGE)
            (p/"hsps.tsv").write_text("".join(AGG.transform_rows([row("a1","b1"),row("a2","b2")],records,containers)))
            (p/"lengths.json").write_text(json.dumps({k:v['length'] for k,v in containers.items()}))
            subprocess.run([sys.executable,str(p/"bridge.py"),str(ROOT/"scripts/third_party/viridic_directional.py"),
                            str(p/"hsps.tsv"),str(p/"lengths.json"),str(p/"out.json")],check=True,capture_output=True)
            result=json.loads((p/"out.json").read_text())["pairs"][0]
            self.assertEqual(result["query_nident"],6)
            self.assertEqual(result["similarity_percent_unrounded"],100*6/14)

    def test_components_retain_bridges_and_distinguish_seed_effects(self):
        keys=["reference:T","source:mutant","source:train","source:near","source:broad"]
        c={k:{"length":10000,"role":"target_reference" if k.startswith('reference') else "external_collection"} for k in keys}
        links={tuple(sorted(p)):v for p,v in [(('reference:T','source:mutant'),69.99),(('source:mutant','source:train'),70),
                                             (('source:near','source:broad'),69.99)]}
        matrix=pairs(c,links)
        a=AGG.components(matrix,c,['source:near'],['source:broad'])
        b=AGG.components(list(reversed(matrix)),dict(reversed(list(c.items()))),['source:near'],['source:broad'])
        self.assertEqual(a['groups'],b['groups'])
        self.assertEqual(set(a['target_only_excluded_source_keys']),{'source:mutant','source:train'})
        self.assertEqual(set(a['added_by_near_seed_component_source_keys']),{'source:near','source:broad'})
        self.assertEqual(len(a['primary_excluded_source_keys']),4)
        # A diagnostic broad-family node alone is not an exclusion seed.
        isolated=AGG.components(pairs(c,{}),c,['source:near'],['source:broad'])
        self.assertNotIn('source:broad',isolated['primary_excluded_source_keys'])

    def test_incomplete_duplicate_nonfinite_and_false_formula_rejected(self):
        c={k:{"length":10000} for k in ['a','b','c']}
        rows=pairs(c,{})
        invalid=[rows[:-1],rows+[rows[0]], [{**rows[0],'similarity_percent_unrounded':float('nan')},*rows[1:]],
                 [{**rows[0],'similarity_percent_unrounded':1},*rows[1:]]]
        for matrix in invalid:
            with self.assertRaises(ValueError): AGG.validate_pairs(matrix,c)

    def test_ambiguous_self_requires_agreement_bounds_not_false_full_identity(self):
        c={'canonical':{'length':10,'ambiguous_bases':0},'ambiguous':{'length':10,'ambiguous_bases':1}}
        py={'self_counts':[{'id':'canonical','nident':10,'original_aligned_query_bases':10},
                           {'id':'ambiguous','nident':9,'original_aligned_query_bases':12}]}
        r=[{'query':x['id'],'subject':x['id'],'nident':str(x['nident']),
            'original_aligned_query_bases':str(x['original_aligned_query_bases']),'error':''} for x in py['self_counts']]
        self.assertTrue(all(x['passed'] for x in AGG.self_checks(py,r,c)))
        r[1]['nident']='8'
        self.assertFalse({x['container_key']:x for x in AGG.self_checks(py,r,c)}['ambiguous']['passed'])
        r[1]['nident']='NA'
        self.assertFalse(all(x['passed'] for x in AGG.self_checks(py,r,c)))


if __name__ == '__main__':
    unittest.main()
