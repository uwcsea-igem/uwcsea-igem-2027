"""Original-R-checked common-context groups with explicit official-version nodes."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/official-variant-group-policy.json'


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);sys.modules[name]=result;spec.loader.exec_module(result)
    return result


def add_official_records(records,containers,observed,official):
    records=dict(records);containers=dict(containers)
    names=set()
    for row in observed:
        name=row['source_id'];identifier=row['id']
        if name not in official or name in names or identifier in records:
            raise ValueError('Unknown/duplicate official node or record ID')
        source=official[name];key='official:'+name+':'+source['official_versioned_accession']
        if (row['role']!='official_reference_version' or row['container_key']!=key or key in containers
            or row['source_record_id']!=source['official_versioned_accession']
            or row['retained_bases']!=source['official_length']
            or row['sequence_sha256']!=source['official_sequence_sha256']):
            raise ValueError('Official node differs from authenticated genotype')
        names.add(name)
        records[identifier]={**row,'offset':0}
        containers[key]={'container_key':key,'source_id':name,'role':'official_reference_version',
                         'cohort_member':False,'published_model_role':None,'records':[identifier],
                         'length':row['retained_bases'],'ambiguous_bases':source['official_ambiguous_bases']}
    if names!=set(official):raise ValueError('Official genotype inventory incomplete')
    return records,dict(sorted(containers.items()))


def require_execution_evidence(alignment,frozen,policy):
    if (any(alignment[k]!=frozen[k] for k in ['subjects_fasta','record_manifest','policy_sha256','script_sha256'])
        or len(alignment['commands'])!=2 or [c['argv'] for c in alignment['commands']]!=frozen['commands']):
        raise ValueError('Search input/commands differ from frozen plan')
    if alignment['status']=='validated_all_record_nucleotide_hsps':
        if any(c['exit_code']!=0 for c in alignment['commands']):
            raise ValueError('Unknown or failed exit code cannot be treated as zero')
        return 'both_original_command_exit_codes_observed_zero'
    if alignment['status']!='validated_recovered_all_record_nucleotide_hsps' or not policy.get('allow_explicit_orphan_artifact_recovery'):
        raise ValueError('Recovery requires an explicit selected consumer policy')
    if (alignment['commands'][0]['exit_code']!=0 or alignment['commands'][1]['exit_code'] is not None
        or not alignment['commands'][1].get('terminal_observed')
        or alignment.get('execution_exit_success_claimed') is not False
        or alignment.get('artifact_completion_verified') is not True
        or alignment.get('search_rerun') is not False
        or alignment.get('execution_status')!='orphaned_controller_exit_code_unobservable'):
        raise ValueError('Recovered execution/data evidence is inconsistent')
    for field,name in [('helper','recovery_helper'),('policy','recovery_policy')]:
        if alignment['recovery'][field]['sha256']!=policy['inputs'][name]['sha256']:
            raise ValueError('Recovery provenance differs from selected policy')
    return 'artifact_validated_but_search_exit_code_unknown'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--policy',type=Path,default=POLICY)
    parser.add_argument('--alignment-report',type=Path,default=ROOT/'data/processed/official-variant-alignments.json')
    parser.add_argument('--workdir',type=Path,default=ROOT/'data/raw/genophi-external-validation/official-variant-groups-v1')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/official-variant-genome-groups.json')
    parser.add_argument('--matrix',type=Path,default=ROOT/'data/processed/official-variant-genome-similarities.tsv')
    args=parser.parse_args()
    if not args.alignment_report.exists():raise ValueError('Complete accepted alignment report required; never infer completion from a running-state file')
    old=load('existing_genome_aggregation',ROOT/'scripts/aggregate-receptor-genomes.py')
    policy_bytes=args.policy.read_bytes();policy=json.loads(policy_bytes);wrapper=old.bound(__file__)
    for spec in policy['inputs'].values():old.checked_file(spec)
    if old.digest(old.BRIDGE.encode())!=policy['GPL_external_boundary']['bridge_sha256']:
        raise ValueError('Pinned external bridge differs')
    if policy['grouping']!={'similarity_percent':70.0,'boundary_buffer_percentage_points':0.01}:
        raise ValueError('Grouping threshold changed')
    for path in [args.workdir,args.output,args.matrix]:
        if path.exists():raise ValueError('Use fresh destinations')
    if args.output.resolve().is_relative_to(args.workdir.resolve()) or args.matrix.resolve().is_relative_to(args.workdir.resolve()):
        raise ValueError('Public outputs must be outside raw workdir')
    alignment_bytes=args.alignment_report.read_bytes();alignment=json.loads(alignment_bytes)
    if (alignment['status'] not in ['validated_all_record_nucleotide_hsps','validated_recovered_all_record_nucleotide_hsps'] or alignment['nucleotide_queries']!=273
        or alignment['database_records_verified']!=273 or alignment['policy_sha256']!=policy['inputs']['alignment_policy']['sha256']
        or alignment['script_sha256']!=policy['inputs']['alignment_helper']['sha256']):
        raise ValueError('Require accepted273-query common-context search')
    frozen=json.loads(old.checked(policy['inputs']['frozen_alignment_plan']))
    execution_evidence=require_execution_evidence(alignment,frozen,policy)
    if alignment['status']=='validated_recovered_all_record_nucleotide_hsps':
        for spec in alignment['recovery'].values():old.checked_file(spec)
    allocation=json.loads(old.checked(policy['inputs']['streaming_control_report']))
    if (allocation['status']!='passed' or not allocation['exact_pair_array_equality']
        or allocation['pair_count']!=4851 or allocation['self_controls']!=99 or allocation['self_failures']
        or allocation['inputs']['bridge']['sha256']!=old.digest(old.BRIDGE.encode())):
        raise ValueError('External allocation-control evidence differs')
    for spec in alignment['inputs']+[alignment[k] for k in ['subjects_fasta','record_manifest','raw_hsps','numeric_hsps']]:
        old.checked_file(spec)
    record_map=json.loads(Path(alignment['record_manifest']['path']).read_text())
    sources=json.loads(old.checked(policy['inputs']['source_manifest']))['entries']
    targets=[g for g in json.loads(old.checked(policy['inputs']['target_manifest']))['genomes'] if g['role']=='target_reference']
    base=[r for r in record_map['records'] if r['role']!='official_reference_version']
    additions=[r for r in record_map['records'] if r['role']=='official_reference_version']
    records,containers=old.build_mapping(sources,targets,base)
    correspondence=json.loads(old.checked(policy['inputs']['official_replay']))
    official={r['phage']:r for r in correspondence['observations'] if not r['exact_author_DNA_match']}
    if sorted(official)!=policy['official_variants']:raise ValueError('Official variant set changed')
    records,containers=add_official_records(records,containers,additions,official)
    if len(records)!=273 or len(containers)!=266:raise ValueError('Augmented full-corpus counts differ')
    if sorted(c['source_id'] for c in containers.values() if c['ambiguous_bases'])!=policy['known_ambiguous_source_ids']:
        raise ValueError('Ambiguous self-control inventory differs from frozen policy')
    old.verify_fasta(alignment['subjects_fasta']['path'],records)
    validator=load('original_alignment_validation',policy['inputs']['original_HSP_validator']['path'])
    batches=load('batch_aware_validation',policy['inputs']['batch_validator']['path'])
    lengths={key:r['retained_bases'] for key,r in records.items()}
    with Path(alignment['raw_hsps']['path']).open() as raw,Path(alignment['numeric_hsps']['path']).open() as numeric:
        class ExactSink:
            def write(self,line):
                if numeric.readline()!=line:raise ValueError('Numeric rows/order differ from raw HSPs')
        completion=batches.validate_batches(raw,lengths,ExactSink(),validator.validate_hsp)
        if numeric.read(1) or completion!=alignment['query_completion']:raise ValueError('Alignment completion replay differs')
    gp17=json.loads(old.checked(policy['inputs']['gp17_report']))
    near=gp17['direct_near_hit_container_keys']
    broad=sorted(set(gp17['broader_diagnostic_container_keys'])-set(near))
    if not set(near+broad)<=set(containers):raise ValueError('gp17 seed node missing')
    # The input policy freezes the gp17 producer/status acceptance and exact key set
    # via its hash-bound completed receipt before any numeric calculation.
    if gp17['status']!=policy['accepted_gp17_status']:raise ValueError('gp17 evidence has not passed')
    args.workdir.mkdir(parents=True)
    transformed=args.workdir/'container-hsps.tsv'
    with Path(alignment['numeric_hsps']['path']).open() as source,transformed.open('x') as output:
        rows=0
        for line in old.transform_rows(source,records,containers):output.write(line);rows+=1
    if rows!=completion['reported_hsp_rows']:raise ValueError('Transformed row count differs')
    length_file=args.workdir/'container-lengths.json';length_file.write_bytes(old.encoded({k:c['length'] for k,c in containers.items()}))
    record_file=args.workdir/'container-record-map.json';record_file.write_bytes(old.encoded({'records':records,'containers':containers}))
    bridge=args.workdir/'external-calculator.py';bridge.write_text(old.BRIDGE)
    compare=load('frozen_original_R_comparison',policy['inputs']['R_comparison_helper']['path'])
    if old.digest(compare.R_DRIVER.encode())!=policy['GPL_external_boundary']['R_driver_sha256']:
        raise ValueError('Original R driver changed')
    runtime=json.loads(old.checked(policy['inputs']['R_runtime_manifest']))
    verified_runtime=compare.validate_runtime(Path(policy['R_home']),Path(policy['R_library']),runtime)
    rdir=args.workdir/'original-R';rdir.mkdir()
    driver=rdir/'run-original.R';driver.write_text(compare.R_DRIVER)
    (rdir/'lengths.tsv').write_text('id\tlength\n'+''.join(f"{k}\t{c['length']}\n" for k,c in containers.items()))
    adapter_file=args.workdir/'adapter-similarities.json'
    pyargv=[sys.executable,str(bridge),policy['inputs']['adapter_source']['path'],str(transformed),str(length_file),str(adapter_file)]
    rargv=[str(Path(policy['R_home'])/'bin/exec/R'),'--vanilla','--slave','-f',str(driver),'--args',
           policy['inputs']['R_upstream_source']['path'],str(transformed),str(rdir)]
    plan={'policy_sha256':old.digest(policy_bytes),'wrapper':wrapper,'alignment_report':old.bound(args.alignment_report),
          'record_mapping':old.bound(record_file),'container_lengths':old.bound(length_file),'HSP_rows':rows,
          'transformed_HSPs':old.bound(transformed),'adapter_command':pyargv,'R_command':rargv,
          'expected_unordered_pairs':35245,'expected_nonself_directions':70490,'self_controls':266}
    (args.workdir/'frozen-calculation-plan.json').write_bytes(old.encoded(plan))
    commands=[old.command(pyargv,args.workdir/'adapter.log',policy['timeouts']['adapter'])]
    validation={};grouping=None;adapter=None
    if commands[0]['returncode']==0:
        adapter=json.loads(adapter_file.read_bytes());old.validate_pairs(adapter['pairs'],containers)
        comparison_policy={'inputs':{'genome_lengths':old.bound(length_file),'adapter_output':old.bound(adapter_file)},
                           'scope':{'distinct_unordered_pairs':35245},
                           'agreement_requirements':{'unrounded_similarity_absolute_tolerance_percentage_points':1e-12}}
        (args.workdir/'R-comparison-inputs.json').write_bytes(old.encoded(comparison_policy))
        env=os.environ.copy();env.update(R_HOME=str(Path(policy['R_home']).resolve()),R_LIBS_USER=str(Path(policy['R_library']).resolve()),
                                        DYLD_LIBRARY_PATH=str((Path(policy['R_home'])/'lib').resolve()))
        commands.append(old.command(rargv,rdir/'execution.log',policy['timeouts']['original_R'],env=env))
        if commands[-1]['returncode']==0:
            generic=compare.compare(comparison_policy,rdir)
            selves=old.self_checks(adapter,compare.read_tsv(rdir/'directional.tsv'),containers)
            passed=not generic['discrepancies'] and not generic['original_R_errors'] and all(s['passed'] for s in selves)
            validation={'generic_original_helper_result':generic,'dataset_specific_self_controls':selves,
                        'dataset_specific_directional_equivalence_verified':passed,'universal_upstream_equivalence_verified':False}
            if passed:grouping=old.components(adapter['pairs'],containers,near,broad,69.99)
    if adapter is not None:
        with args.matrix.open('x') as stream:
            writer=csv.writer(stream,delimiter='\t',lineterminator='\n');writer.writerow(old.FIELDS)
            for row in adapter['pairs']:
                writer.writerow([row['query'],row['subject'],containers[row['query']]['source_id'],containers[row['subject']]['source_id'],
                                 repr(row['similarity_percent_unrounded']),row['query_nident'],row['subject_nident'],row['query_length'],row['subject_length']])
    panel={r['phage']:('source:'+r['phage'] if r['exact_author_DNA_match'] else 'official:'+r['phage']+':'+r['official_versioned_accession'])
           for r in correspondence['observations']}
    if not set(panel.values())<=set(containers):raise ValueError('Panel genotype key missing')
    for spec in policy['inputs'].values():old.checked_file(spec)
    for spec in alignment['inputs']+[alignment[k] for k in ['subjects_fasta','record_manifest','raw_hsps','numeric_hsps']]:old.checked_file(spec)
    if args.policy.read_bytes()!=policy_bytes or old.bound(__file__)!=wrapper or args.alignment_report.read_bytes()!=alignment_bytes:
        raise ValueError('Executed source changed')
    result={'schema_version':1,'status':'completed_original_R_checked_official_variant_groups' if grouping else 'HOLD_numeric_calculation_or_validation',
            'policy_sha256':old.digest(policy_bytes),'wrapper':wrapper,'inputs':policy['inputs'],'plan':plan,
            'commands':commands,'runtime_content_verification':verified_runtime,'containers':list(containers.values()),
            'numeric_validation':validation,'grouping':grouping,'official_panel_node_map':panel,
            'source_execution_evidence':execution_evidence,
            'numeric_matrix':old.bound(args.matrix) if args.matrix.exists() else None,
            'numeric_matrix_for_grouping_verified':grouping is not None,
            'final_training_admission':False,'phenotype_tables_or_models_used':False,'new_original_benchmark_coverage':0,
            'raw_artifacts':[old.bound(p) for p in sorted(args.workdir.rglob('*')) if p.is_file()]}
    with args.output.open('xb') as stream:stream.write(old.encoded(result))
    print(json.dumps({'status':result['status'],'groups':len(grouping['groups']) if grouping else None}))


if __name__=='__main__':main()
