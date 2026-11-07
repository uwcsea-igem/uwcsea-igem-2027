"""Literal host-state contrasts across frozen, disjoint source-phage buckets."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/source-host-contrast-policy.json'
COLUMNS=('BW25113 susceptibility','BL21 susceptibility')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(keys,metadata,groups):
    if len(set(keys))!=len(keys):
        raise ValueError('Duplicate source member')
    categories=defaultdict(list)
    exclusions=[]
    states={column:Counter() for column in COLUMNS}
    for key in keys:
        if not key.startswith('source:') or key[7:] not in metadata or key not in groups:
            raise ValueError('Missing source metadata or genomic group')
        row=metadata[key[7:]]
        values=tuple(row[column] for column in COLUMNS)
        for column,value in zip(COLUMNS,values):states[column][value]+=1
        if all(value in ('Yes','No') for value in values):
            categories['BW_'+values[0]+'__BL_'+values[1]].append(key)
        else:
            exclusions.append({'key':key,'BW25113_literal':values[0],'BL21_literal':values[1],'group':groups[key]})
    named={category:{'count':len(categories[category]),'keys':sorted(categories[category]),
                      'group_count':len({groups[k] for k in categories[category]}),
                      'groups':sorted({groups[k] for k in categories[category]})}
           for category in ['BW_Yes__BL_Yes','BW_Yes__BL_No','BW_No__BL_Yes','BW_No__BL_No']}
    directions=[category for category in ['BW_Yes__BL_No','BW_No__BL_Yes'] if named[category]['count']]
    return {'records':len(keys),'literal_state_counts':{c:dict(sorted(v.items())) for c,v in states.items()},
            'complete_literal_pairs':sum(v['count'] for v in named.values()),'categories':named,
            'incomplete_or_nonbinary_pairs':exclusions,'mixed_state_directions_present':directions,
            'opposing_host_preferences_observed':len(directions)==2,
            'fixed_host_order_can_satisfy_every_mixed_state_comparison':len(directions)<=1,
            'has_any_mixed_state_comparisons':bool(directions)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/receptor-source-host-contrasts.json')
    args=parser.parse_args()
    policy=json.loads(POLICY.read_text())
    for item in policy['inputs']:
        if sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Bound source changed')
    with (ROOT/policy['inputs'][0]['path']).open(newline='') as stream:
        rows=list(csv.DictReader(stream,delimiter='\t'))
    metadata={r['Phage']:r for r in rows}
    if len(rows)!=260 or len(metadata)!=260:raise ValueError('Cohort metadata cardinality changed')
    roles=json.loads((ROOT/policy['inputs'][1]['path']).read_text())
    groups=json.loads((ROOT/policy['inputs'][2]['path']).read_text())['grouping']['container_groups']
    all_keys=[k for bucket in roles['buckets'].values() for k in bucket]
    if len(set(all_keys))!=260 or len(all_keys)!=260:raise ValueError('Partition accounting differs')
    partitions={bucket:summarize(keys,metadata,groups) for bucket,keys in roles['buckets'].items()}
    train_groups={groups[k] for k in roles['buckets']['candidate_training']}
    test_groups={groups[k] for k in roles['buckets']['candidate_testing']}
    if train_groups & test_groups:raise ValueError('Source split group overlap')
    report={'schema_version':1,'status':'literal_susceptibility_contrast_audit_not_training_admission',
            'script_sha256':sha(__file__),'policy_sha256':sha(POLICY),'inputs':policy['inputs'],
            'partitions':partitions,'cohort_records_accounted_once':260,'candidate_group_overlap':[],
            'labels_used_to_change_partitions':False,'susceptibility_semantics_validated_for_model':False,
            'original_ando_case_file_loaded':False,
            'target_related_source_annotation_bucket_audited_but_not_admitted':True,
            'model_fitted_or_scores_generated':False,
            'new_original_benchmark_coverage':0,'both_unseen_contract_changed':False}
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({k:{c:v['count'] for c,v in r['categories'].items()} for k,r in partitions.items()},indent=2))


if __name__=='__main__':main()
