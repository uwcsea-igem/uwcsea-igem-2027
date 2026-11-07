"""Post hoc label-aware diagnostic of frozen representations, never a predictor."""
import argparse
import csv
import hashlib
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def diagnose(rows, vectors):
    """Count endpoint comparisons forced to tie, without selecting any scores."""
    seen, queries = set(), {}
    for row in rows:
        query, host = row['phage'], row['host']
        label = row['label']
        if label not in (0, 1) or (query, host) in seen or host not in vectors:
            raise ValueError('Invalid/duplicate observation or missing host vector')
        seen.add((query, host))
        queries.setdefault(query, []).append((host, label))
    groups = {}
    for host, vector in vectors.items():
        if not vector or any(type(x) is not int or x not in (0, 1) for x in vector):
            raise ValueError('Expected nonempty binary host vector')
        groups.setdefault(tuple(vector), []).append(host)
    if len({len(v) for v in vectors.values()}) != 1:
        raise ValueError('Host vector widths differ')
    diagnostics, required_separations = {}, set()
    for query, observations in sorted(queries.items()):
        positive = [h for h,label in observations if label == 1]
        negative = [h for h,label in observations if label == 0]
        required_separations.update(tuple(sorted((p,n))) for p in positive for n in negative)
        forced = [{'positive_host': p, 'negative_host': n} for p in positive for n in negative
                  if vectors[p] == vectors[n]]
        diagnostics[query] = {'observations': len(observations), 'positive_negative_pairs': len(positive)*len(negative),
                              'forced_tied_positive_negative_pairs': forced,
                              'informative_for_primary_concordance': bool(positive and negative)}
    informative = [q for q,d in diagnostics.items() if d['informative_for_primary_concordance']]
    return {'collision_groups': [sorted(gs) for gs in groups.values() if len(gs)>1],
            'queries': diagnostics, 'informative_queries': informative,
            'forced_tied_primary_pair_count': sum(len(diagnostics[q]['forced_tied_positive_negative_pairs']) for q in informative),
            'single_informative_query_cannot_establish_phage_dependent_ranking': len(informative)==1,
            'required_host_separations': [list(pair) for pair in sorted(required_separations)],
            'all_hosts_distinct_required_by_this_endpoint': len(vectors)>1 and all(
                pair in required_separations for pair in combinations(sorted(vectors),2))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'data/processed/host-collision-endpoint-audit.json')
    args = parser.parse_args()
    cases = ROOT/'software/tipcraft/benchmarks/retro-t7-t3/cases.csv'
    contract = ROOT/'software/tipcraft/benchmarks/retro-t7-t3/next-method-contract.json'
    sources = {
        'global_union_uniform': ROOT/'data/processed/receptor-host-feature-uniform-control.json',
        'family_conditioned_binary': ROOT/'data/processed/receptor-family-host-capacity.json'}
    with cases.open(newline='') as stream:
        all_rows = list(csv.DictReader(stream))
    rows = [{**r, 'label': int(r['label'])} for r in all_rows if r['phage'] in ('T7','T3')]
    if len(rows)!=6 or {r['host'] for r in rows}!={'BL21','BW25113','MG1655'}:
        raise ValueError('Original six-observation native panel changed')
    results = {}
    for name,path in sources.items():
        report = json.loads(path.read_text())
        projection = report['target_projection']
        if isinstance(projection['matrix'],dict):
            vectors = {k.removeprefix('original:'):v for k,v in projection['matrix'].items()}
        else:
            vectors = {k.removeprefix('original:'):v for k,v in zip(projection['sample_ids'],projection['matrix'],strict=True)}
        results[name] = diagnose(rows,vectors)
    bindings = [{'path': str(p.relative_to(ROOT)), 'sha256': sha(p)} for p in [cases,contract,*sources.values()]]
    result = {'schema_version':1,'status':'post_hoc_endpoint_collision_diagnostic',
              'script_sha256':sha(__file__),'inputs':bindings,'native_observations':6,
              'exchange_observations_excluded_from_diagnostic':len(all_rows)-len(rows),
              'target_labels_read':True,'target_label_use':'Interpret frozen collisions against the recorded endpoint only; not fit, select or calibrate any feature, threshold, model or score',
              'primary_endpoint':json.loads(contract.read_text())['primary_endpoint'],
              'representations':results,'scores_generated':False,'model_fitted':False,
              'existing_feature_reports_and_policies_modified':False,
              'new_method_covered_cases':0,
              'interpretation':'The family collision does not force a positive-negative tie on this panel; no accuracy is measured. Existing source, stock, overlap, split and training-admission gates remain unresolved. A host-only comparator is essential because only T3 contributes mixed-label native ranking comparisons.'}
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,sort_keys=True)
        stream.write('\n')
    print(json.dumps({name:r['forced_tied_primary_pair_count'] for name,r in results.items()}))


if __name__=='__main__':
    main()
