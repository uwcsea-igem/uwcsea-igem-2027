"""Recompute the complete ECOR group/exclusion report from saved skani rows."""
import argparse
import csv
import json
from pathlib import Path

import importlib.util

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,default=ROOT/'data/processed/ecor-host-reference-groups.json')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/ecor-host-reference-groups-replay.json')
    args=parser.parse_args()
    path=ROOT/'scripts/group-ecor-host-references.py'
    spec=importlib.util.spec_from_file_location('producer',path)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    report=json.loads(args.report.read_text())
    if producer.sha(path)!=report['script_sha256'] or producer.sha(producer.POLICY)!=report['policy_sha256']:
        raise ValueError('Executed producer/policy changed')
    for item in report['inputs']:
        if producer.sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Bound dependency changed')
    for item in report['files']:
        if producer.sha(item['path'])!=item['file_sha256']:raise ValueError('Reference changed')
    matrix=report['matrix']
    if producer.sha(matrix['path'])!=matrix['sha256']:raise ValueError('Raw matrix changed')
    summarizer=producer.load('saved_matrix_summary',ROOT/'scripts/audit-receptor-host-relatedness.py')
    near=producer.load('fixed_near_rule',ROOT/'scripts/audit-genophi-host-overlap.py')
    policy=json.loads(producer.POLICY.read_text())
    with Path(matrix['path']).open() as stream:
        summary=summarizer.summarize(list(csv.DictReader(stream,delimiter='\t')),report['files'],policy['near_reference_rule'],near.classify_pair)
    if summary!=report['comparison_summary']:raise ValueError('Complete comparison summary replay differs')
    prior=json.loads((ROOT/'data/processed/receptor-host-relatedness-audit.json').read_text())
    current={(r['reference'],r['query']):r for r in summary['comparisons']}
    exact=len(prior['comparisons'])==25 and all(current.get((r['reference'],r['query']))==r for r in prior['comparisons'])
    if exact!=report['prior25_control_exact']:raise ValueError('Prior control replay differs')
    partition=producer.partition_components(summary['reference_components'] if exact else None,[f['id'] for f in report['files']])
    if partition!=report['partition']:raise ValueError('Group/exclusion replay differs')
    result={'schema_version':1,'status':'exact_complete_matrix_and_partition_replay',
            'report':{'path':str(args.report),'sha256':producer.sha(args.report)},
            'script_sha256':producer.sha(__file__),'matrix':matrix,'directed_pairs':summary['directed_pair_count'],
            'prior25_control_exact':exact,'candidate_source_hosts':partition['candidate_source_hosts'],
            'target_related_source_hosts':partition['target_related_source_hosts'],
            'candidate_source_group_count':partition.get('candidate_source_group_count'),
            'skani_network_or_model_execution':False}
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({'status':result['status'],'directions':result['directed_pairs'],'candidate_groups':result['candidate_source_group_count']}))


if __name__=='__main__':main()
