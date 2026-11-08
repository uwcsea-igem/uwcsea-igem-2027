"""Observe/resume an existing orphaned search; validate artifacts without inventing an exit code."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/official-variant-recovery-policy.json'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


def observe(pid):
    state=subprocess.run(['ps','-ww','-p',str(pid),'-o','state='],capture_output=True,text=True).stdout.strip()
    if not state:return None
    command=subprocess.check_output(['ps','-ww','-p',str(pid),'-o','command='],text=True).strip()
    started=subprocess.check_output(['ps','-p',str(pid),'-o','lstart='],text=True).strip()
    return {'pid':pid,'state':state,'command':command,'start_description':started}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/official-variant-alignments-recovered.json')
    args=parser.parse_args()
    policy_bytes=POLICY.read_bytes();policy=json.loads(policy_bytes);code_hash=sha(__file__)
    for item in policy['bindings']:
        if sha(item['path'])!=item['sha256']:raise ValueError('Recovery-bound input changed')
    plan=json.loads(Path(policy['frozen_plan']).read_text())
    for item in plan['inputs']+[plan['subjects_fasta'],plan['record_manifest']]:
        if sha(item['path'])!=item['sha256']:raise ValueError('Original scientific input changed')
    work=Path(policy['workdir']);receipt_path=work/'recovery-observation.json'
    if args.output.exists() or receipt_path.exists():raise ValueError('Use a new recovery receipt; never restart/resume blindly')
    before=observe(policy['pid'])
    if before is None or any(before[k]!=policy['process_identity'][k] for k in ['command','start_description']):
        raise ValueError('Original process identity changed or missing before attachment')
    if before['command']!=' '.join(plan['commands'][1]):raise ValueError('Observed executable/arguments differ from frozen command')
    receipt={'policy_sha256':hashlib.sha256(policy_bytes).hexdigest(),'recovery_helper_sha256':code_hash,
             'attached_epoch':time.time(),'before':before,'resume_signal_sent':False,'search_restarted':False,
             'exit_code':None,'exit_status_observable_for_nonchild':False}
    if before['state'].startswith('T'):
        os.kill(policy['pid'],signal.SIGCONT);receipt['resume_signal_sent']=True
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'attached_existing_pid':policy['pid'],'resumed':receipt['resume_signal_sent'],'exit_code':'unobservable_nonchild'}),flush=True)
    while True:
        current=observe(policy['pid'])
        if current is None or current['state'].startswith('Z'):break
        if any(current[k]!=before[k] for k in ['command','start_description']):raise ValueError('PID identity changed')
        print('OBSERVED same existing search',policy['pid'],current['state'],flush=True)
        time.sleep(30)
    receipt.update({'terminal_observed_epoch':time.time(),'terminal_observation':current})
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    for item in policy['bindings']+plan['inputs']+[plan['subjects_fasta'],plan['record_manifest']]:
        if sha(item['path'])!=item['sha256']:raise ValueError('Bound input changed while observing')
    validator=load('original_hsp_validator',ROOT/'scripts/align-receptor-cohort.py')
    batches=load('batch_protocol',ROOT/'scripts/replay-receptor-alignments.py')
    mapping=json.loads(Path(plan['record_manifest']['path']).read_text())
    lengths={r['id']:r['retained_bases'] for r in mapping['records']}
    if len(lengths)!=273:raise ValueError('Recovery requires full273-record inventory')
    raw=work/'blastn.tsv';numeric=work/'recovered-blastn-data.tsv'
    with raw.open() as source,numeric.open('x') as dest:
        completion=batches.validate_batches(source,lengths,dest,validator.validate_hsp)
    first=json.loads((work/'step-1.json').read_text())
    if first['exit_code']!=0 or first['argv']!=plan['commands'][0]:raise ValueError('Original database construction not verified')
    if POLICY.read_bytes()!=policy_bytes or sha(__file__)!=code_hash:raise ValueError('Recovery protocol changed')
    result={**plan,'status':'validated_recovered_all_record_nucleotide_hsps',
            'original_alignment_script_sha256':plan['script_sha256'],
            'recovery':{'policy':{'path':str(POLICY),'sha256':hashlib.sha256(policy_bytes).hexdigest()},
                        'helper':{'path':__file__,'sha256':code_hash},'observation':{'path':str(receipt_path),'sha256':sha(receipt_path)}},
            'execution_status':'orphaned_controller_exit_code_unobservable','search_rerun':False,
            'commands':[first,{'argv':plan['commands'][1],'exit_code':None,'terminal_observed':True}],
            'database_records_verified':273,'query_completion':completion,
            'raw_hsps':{'path':str(raw),'sha256':sha(raw),'bytes':raw.stat().st_size},
            'numeric_hsps':{'path':str(numeric),'sha256':sha(numeric),'bytes':numeric.stat().st_size},
            'artifact_completion_verified':True,'execution_exit_success_claimed':False,
            'phenotype_tables_loaded':False,'training_or_prediction':False,'new_original_benchmark_coverage':0}
    with args.output.open('x') as dest:json.dump(result,dest,indent=2,sort_keys=True);dest.write('\n')
    print(json.dumps({'artifact_status':result['status'],'execution_exit_code':None,'HSPs':completion['reported_hsp_rows']}),flush=True)


if __name__=='__main__':main()
