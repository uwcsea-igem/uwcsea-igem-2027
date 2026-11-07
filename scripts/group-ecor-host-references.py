"""Group authenticated ECOR and target references under the existing ANI rule."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/ecor-host-relatedness-policy.json'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def partition_components(components,all_ids):
    if components is None:
        return {'status':'HOLD_incomplete_reference_grouping','source_groups':None,
                'target_related_source_hosts':None,'candidate_source_hosts':None}
    flattened=[name for group in components for name in group]
    if len(flattened)!=len(set(flattened)) or set(flattened)!=set(all_ids):
        raise ValueError('Components must cover each reference exactly once')
    source_groups={};excluded=[];candidates=[]
    for group in components:
        members=sorted(group)
        key='ecor-host-group-'+hashlib.sha256(json.dumps(members,separators=(',',':')).encode()).hexdigest()[:16]
        target=any(name.startswith('original:') for name in members)
        source=[name for name in members if name.startswith('ecor:')]
        source_groups[key]={'members':members,'target_containing':target,'source_members':source}
        (excluded if target else candidates).extend(source)
    return {'status':'complete_host_reference_grouping','source_groups':source_groups,
            'target_related_source_hosts':sorted(excluded),'candidate_source_hosts':sorted(candidates),
            'candidate_source_group_count':sum(bool(g['source_members']) and not g['target_containing'] for g in source_groups.values())}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,default=ROOT/'data/raw/genophi-external-validation/ecor-host-groups-v1')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/ecor-host-reference-groups.json')
    args=parser.parse_args()
    policy_bytes=POLICY.read_bytes();policy=json.loads(policy_bytes)
    script_hash=sha(__file__)
    bindings=policy['inputs']+[policy['executable']]
    for item in bindings:
        if sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Bound dependency changed')
    for item in policy['files']:
        if sha(item['path'])!=item['file_sha256']:raise ValueError('Reference FASTA changed')
    if len(policy['files'])!=30 or policy['expected_directed_pairs_including_self']!=900:
        raise ValueError('Expected all30 references and900 directions')
    if args.workdir.exists() or args.output.exists():raise ValueError('Use fresh destinations')
    args.workdir.mkdir(parents=True)
    summarizer=load('host_relatedness',ROOT/'scripts/audit-receptor-host-relatedness.py')
    near=load('near_rule',ROOT/'scripts/audit-genophi-host-overlap.py')
    binary=str(ROOT/policy['executable']['path'])
    if subprocess.check_output([binary,'--version'],text=True).strip()!=policy['executable']['version']:
        raise ValueError('skani version differs')
    files=args.workdir/'files.txt';files.write_text(''.join(str(Path(x['path']).resolve())+'\n' for x in policy['files']))
    matrix=args.workdir/'skani.tsv';p=policy['parameters']
    command=[binary,'dist','--rl',str(files),'--ql',str(files),'-c',str(p['compression']),'-m',str(p['marker_compression']),
             '-s',str(p['screen_percent']),'--min-af',str(p['minimum_reported_af_percent']),'-t',str(p['threads']),'-o',str(matrix)]
    plan={'policy_sha256':hashlib.sha256(policy_bytes).hexdigest(),'script_sha256':script_hash,
          'files':policy['files'],'command':command,'inputs':bindings}
    (args.workdir/'frozen-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    result=subprocess.run(command,capture_output=True,text=True,timeout=policy['timeout_seconds'])
    (args.workdir/'run.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    with matrix.open() as stream:
        summary=summarizer.summarize(list(csv.DictReader(stream,delimiter='\t')),policy['files'],policy['near_reference_rule'],near.classify_pair)
    prior=json.loads((ROOT/'data/processed/receptor-host-relatedness-audit.json').read_text())
    current={(r['reference'],r['query']):r for r in summary['comparisons']}
    differences=[{'prior':row,'current':current.get((row['reference'],row['query']))} for row in prior['comparisons']
                 if current.get((row['reference'],row['query']))!=row]
    controls_pass=not differences and len(prior['comparisons'])==25
    admitted=partition_components(summary['reference_components'] if controls_pass else None,[x['id'] for x in policy['files']])
    report={'schema_version':1,**plan,'comparison_summary':summary,'prior25_control_exact':controls_pass,
            'prior_control_differences':differences,'partition':admitted,
            'matrix':{'path':str(matrix),'sha256':sha(matrix)},'outcomes_loaded':False,'model_fitted_or_scored':False,
            'experimental_stock_identity_verified':False,'final_training_admission':False,'new_original_benchmark_coverage':0,
            'interpretation':'Reference relatedness under the existing buffered rule; all known references retained as bridges. No stock identity, taxonomy, assembly prediction or phenotype is inferred.'}
    for item in bindings:
        if sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Dependency changed during comparison')
    for item in policy['files']:
        if sha(item['path'])!=item['file_sha256']:raise ValueError('Reference changed during comparison')
    if POLICY.read_bytes()!=policy_bytes or sha(__file__)!=script_hash:raise ValueError('Execution policy or producer changed')
    with args.output.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
    print(json.dumps({'comparison_status':summary['status'],'prior25_control_exact':controls_pass,'partition':admitted}),flush=True)


if __name__=='__main__':main()
