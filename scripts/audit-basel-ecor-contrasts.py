"""Measure secondary-panel contrast availability, never predictive performance."""
import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal, InvalidOperation
import hashlib
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/basel-ecor-contrast-policy.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(rows):
    seen,by_phage=set(),defaultdict(dict)
    literals=Counter()
    for row in rows:
        host,phage=row['strain'],row['phage']
        if not host or not phage or (host,phage) in seen:
            raise ValueError('Empty ID or duplicate observation')
        try:value=Decimal(row['interaction'])
        except InvalidOperation as exc:raise ValueError('Invalid interaction') from exc
        if not value.is_finite() or value not in (0,1):raise ValueError('Nonbinary or nonfinite interaction')
        literals[row['interaction']]+=1
        seen.add((host,phage));by_phage[phage][host]=int(value)
    if not seen:raise ValueError('Empty panel')
    hosts=sorted({host for host,_ in seen})
    directions=defaultdict(lambda:[[],[]])
    queries={}
    for phage,states in sorted(by_phage.items()):
        positive=sum(states.values());negative=len(states)-positive
        queries[phage]={'observed':len(states),'positive':positive,'negative':negative,
                        'missing':len(hosts)-len(states),'positive_negative_comparisons':positive*negative}
        for a,b in combinations(sorted(states),2):
            if states[a]!=states[b]:directions[(a,b)][0 if states[a]>states[b] else 1].append(phage)
    pairs=[{'hosts':[a,b],'first_host_positive_phages':forward,'second_host_positive_phages':reverse,
            'opposing_directions':bool(forward and reverse)} for (a,b),(forward,reverse) in sorted(directions.items())]
    return {'hosts':hosts,'phages':sorted(by_phage),'observed_pairs':len(seen),
            'cartesian_pairs':len(hosts)*len(by_phage),'missing_pairs':len(hosts)*len(by_phage)-len(seen),
            'raw_value_counts':dict(sorted(literals.items())),'queries':queries,
            'informative_query_count':sum(bool(q['positive'] and q['negative']) for q in queries.values()),
            'total_positive_negative_comparisons':sum(q['positive_negative_comparisons'] for q in queries.values()),
            'host_pairs_with_any_direction':len(pairs),
            'host_pairs_with_opposing_directions':sum(p['opposing_directions'] for p in pairs),
            'host_pair_contrasts':pairs}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/basel-ecor-contrast-audit.json')
    args=parser.parse_args()
    policy=json.loads(POLICY.read_text())
    for item in policy['inputs']:
        if sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Bound source differs')
    authority=json.loads((ROOT/policy['inputs'][1]['path']).read_text())
    with (ROOT/policy['inputs'][0]['path']).open(newline='') as stream:
        result=summarize(list(csv.DictReader(stream)))
    if (result['observed_pairs']!=authority['panel_pairs'] or
        set(result['hosts'])!=set(authority['host_id_to_wgs_master_accession']) or
        set(result['phages'])!=set(authority['phage_id_to_accession'])):
        raise ValueError('Entire source roster changed')
    result.update({'schema_version':1,'status':'secondary_source_contrast_audit_not_prediction',
                   'script_sha256':sha(__file__),'policy_sha256':sha(POLICY),'inputs':policy['inputs'],
                   'outcomes_inspected_for_source_feasibility':True,'outcome_selected_subset':False,
                   'target_exclusions_or_genomic_groups_applied':False,'model_fitted_or_scored':False,
                   'host_genomes_downloaded':False,'original_benchmark_replaced':False,
                   'new_original_benchmark_coverage':0,
                   'endpoint':'Published spot-clearance/no-clearance records; missing pairs unknown',
                   'interpretation':'Opposing host-pair directions establish contrast availability, not a learned score, genetic independence, host stock identity or a new blind holdout.'})
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['observed_pairs','missing_pairs','informative_query_count','host_pairs_with_opposing_directions']}))


if __name__=='__main__':main()
