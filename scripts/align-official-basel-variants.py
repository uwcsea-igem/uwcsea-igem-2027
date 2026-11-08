"""Full common-context BLAST of old references plus three official DNA variants."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/official-variant-alignment-policy.json'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


def prepare(policy):
    import Bio
    if Bio.__version__!='1.86':raise ValueError('Frozen preparation requires Biopython1.86')
    old=json.loads((ROOT/policy['old_alignment_report']).read_text())
    for field in ['record_manifest','subjects_fasta']:
        if sha(old[field]['path'])!=old[field]['sha256']:raise ValueError('Old corpus/map changed')
    original=Path(old['subjects_fasta']['path']).read_bytes()
    if not original.endswith(b'\n'):raise ValueError('Require explicit record boundary at old FASTA end')
    mapping=json.loads(Path(old['record_manifest']['path']).read_text())
    records=list(mapping['records'])
    if len(records)!=270:raise ValueError('Expected full original270-record corpus')
    official=json.loads((ROOT/policy['official_replay']).read_text())
    inputs=json.loads((ROOT/policy['official_authentication']).read_text())
    by_name={e['phage']:e for e in inputs['entries']}
    observations={o['phage']:o for o in official['observations']}
    if official['author_DNA_mismatches']!=policy['official_variants']:raise ValueError('Variant roster changed')
    parser=module('official_source_parser',ROOT/'scripts/authenticate-basel-phage-accessions.py')
    fasta=original
    for name in policy['official_variants']:
        item=next(x for x in by_name[name]['official_payloads'] if x['kind']=='fasta')
        if sha(item['path'])!=item['sha256']:raise ValueError('Official payload changed')
        record=parser.sequence_record(Path(item['path']).read_bytes(),'fasta',by_name[name]['requested_accession'])
        expected=observations[name]
        if record['sequence_sha256']!=expected['official_sequence_sha256']:raise ValueError('Official DNA changed')
        identifier='official_'+name
        records.append({'id':identifier,'source_id':name,'role':'official_reference_version',
                        'source_record_id':record['versioned_accession'],'retained_bases':record['length'],
                        'sequence_sha256':record['sequence_sha256'],
                        'container_key':'official:'+name+':'+record['versioned_accession']})
        fasta+=f'>{identifier}\n{record["sequence"]}\n'.encode()
    if len(records)!=273 or len({r['id'] for r in records})!=273 or len({r['container_key'] for r in records})!=266:
        raise ValueError('Full augmented record/container cardinality differs')
    return fasta,{'schema_version':1,'records':records,'original_records_retained':270,
                  'official_variant_records_added':3,'nucleotide_records':273,'container_nodes':266,
                  'DNA_concatenated_trimmed_or_designed':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,default=ROOT/'data/raw/genophi-external-validation/official-variant-alignments-v1')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/official-variant-alignments.json')
    args=parser.parse_args()
    policy_bytes=POLICY.read_bytes();policy=json.loads(policy_bytes);script_hash=sha(__file__)
    bindings=policy['inputs']+list(policy['executables'].values())
    def verify():
        for item in bindings:
            if sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Bound dependency changed')
        if POLICY.read_bytes()!=policy_bytes or sha(__file__)!=script_hash:raise ValueError('Producer/policy changed')
    verify()
    fasta,mapping=prepare(policy)
    if hashlib.sha256(fasta).hexdigest()!=policy['expected_subjects_sha256']:raise ValueError('Augmented corpus differs from frozen plan')
    if args.workdir.exists() or args.output.exists():raise ValueError('Use fresh output destinations')
    args.workdir.mkdir(parents=True)
    subjects=args.workdir/'subjects.fna';subjects.write_bytes(fasta)
    manifest=args.workdir/'record-manifest.json';manifest.write_text(json.dumps(mapping,indent=2,sort_keys=True)+'\n')
    validators=module('original_hsp_validator',ROOT/'scripts/align-receptor-cohort.py')
    batches=module('batch_protocol_validator',ROOT/'scripts/replay-receptor-alignments.py')
    if policy['blastn_parameters']!=validators.PARAMETERS or policy['blastn_fields']!=validators.FIELDS:
        raise ValueError('Search/HSP schema changed from prior method')
    for name,item in policy['executables'].items():
        if subprocess.check_output([str(ROOT/item['path']),'-version'],text=True).splitlines()[0]!=item['version']:
            raise ValueError('Executable version differs')
    database=args.workdir/'db';raw=args.workdir/'blastn.tsv';numeric=args.workdir/'blastn-data.tsv'
    commands=[[str(ROOT/policy['executables']['makeblastdb']['path']),'-in',str(subjects),'-dbtype','nucl','-parse_seqids','-out',str(database)],
              [str(ROOT/policy['executables']['blastn']['path']),'-query',str(subjects),'-db',str(database),
               *policy['blastn_parameters'],'-outfmt','7 '+' '.join(policy['blastn_fields']),'-out',str(raw)]]
    plan={'schema_version':1,'policy_sha256':hashlib.sha256(policy_bytes).hexdigest(),'script_sha256':script_hash,
          'subjects_fasta':{'path':str(subjects),'sha256':sha(subjects)},
          'record_manifest':{'path':str(manifest),'sha256':sha(manifest)},'commands':commands,
          'inputs':bindings,'status':'frozen_before_database_construction','nucleotide_queries':273,'container_nodes':266}
    (args.workdir/'frozen-plan.json').write_text(json.dumps(plan,indent=2,sort_keys=True)+'\n')
    ledger=[]
    for index,command in enumerate(commands,1):
        verify();log=args.workdir/f'step-{index}.log';start=time.monotonic()
        with log.open('xb') as stream:
            child=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            state={'step':index,'pid':child.pid,'argv':command,'started_epoch':time.time(),'state':'running'}
            (args.workdir/'process.json').write_text(json.dumps(state,indent=2)+'\n')
            print(json.dumps(state),flush=True)
            while True:
                try:code=child.wait(timeout=30);break
                except subprocess.TimeoutExpired:print(f'RUNNING same pid={child.pid}, step={index}, elapsed={time.monotonic()-start:.0f}s',flush=True)
        state.update({'state':'terminal','exit_code':code});(args.workdir/'process.json').write_text(json.dumps(state,indent=2)+'\n')
        ledger.append({'argv':command,'exit_code':code,'elapsed_seconds':time.monotonic()-start,'log':str(log),'log_sha256':sha(log)})
        (args.workdir/f'step-{index}.json').write_text(json.dumps(ledger[-1],indent=2)+'\n')
        if code:raise ValueError('Search command failed; no automatic restart')
        if index==1:
            text=log.read_text()
            if '273 sequences' not in text:raise ValueError('Database count not verified as273')
    verify()
    lengths={r['id']:r['retained_bases'] for r in mapping['records']}
    with raw.open() as stream,numeric.open('x') as target:
        completion=batches.validate_batches(stream,lengths,target,validators.validate_hsp)
    verify()
    result={**plan,'status':'validated_all_record_nucleotide_hsps','commands':ledger,'query_completion':completion,
            'raw_hsps':{'path':str(raw),'sha256':sha(raw),'bytes':raw.stat().st_size},
            'numeric_hsps':{'path':str(numeric),'sha256':sha(numeric),'bytes':numeric.stat().st_size},
            'database_records_verified':273,'old_HSP_rows_reused':False,'similarities_or_groups_calculated':False,
            'phenotype_tables_loaded':False,'training_or_prediction':False,
            'new_original_benchmark_coverage':0}
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({'status':result['status'],'queries':273,'HSPs':completion['reported_hsp_rows']}),flush=True)


if __name__=='__main__':main()
