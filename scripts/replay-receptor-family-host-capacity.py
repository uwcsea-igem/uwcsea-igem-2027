"""Verify and replay the completed family host pilot without sequence searches."""
import argparse
import importlib.util
import json
from pathlib import Path

from tipcraft import protein_features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, default=Path('data/processed/receptor-family-host-capacity.json'))
    parser.add_argument('--output', type=Path, default=Path('data/processed/receptor-family-host-capacity-replay.json'))
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('capacity', 'scripts/check-receptor-family-host-capacity.py')
    pilot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pilot)
    report = json.loads(args.report.read_text())
    if pilot.sha(spec.origin) != report['script_sha256'] or pilot.sha(pilot.POLICY) != report['policy_sha256']:
        raise ValueError('Executed producer or policy changed')
    policy = json.loads(pilot.POLICY.read_text())
    freeze = json.loads(Path(report['vocabulary_freeze']['path']).read_text())
    bindings = report['inputs'] + freeze['artifacts'] + [report['vocabulary_freeze']]
    for item in bindings:
        if pilot.sha(item['path']) != item['sha256']:
            raise ValueError(f"Binding mismatch: {item['path']}")
    work = Path(report['vocabulary_freeze']['path']).parent
    stored = json.loads((work/'vocabulary.json').read_text())
    members = {(family, kmer): feature for family,kmer,feature in stored['members']}
    expected = stored['expected_patterns']
    inventories = {}
    for role,prefix in [('vocabulary_fit_only','s'),('projection_only','q')]:
        names, seqs, genomes = [], {}, {}
        for gi,item in enumerate(f for f in policy['files'] if f['role']==role):
            names.append(item['sample_id'])
            records = []
            name, chunks = None, []
            for line in Path(item['path']).read_text().splitlines():
                if line.startswith('>'):
                    if name is not None:
                        records.append(''.join(chunks).removesuffix('*'))
                    name, chunks = line[1:].split()[0], []
                elif line.strip():
                    chunks.append(line)
            if name is not None:
                records.append(''.join(chunks).removesuffix('*'))
            if len(records)!=item['expected_proteins']:
                raise ValueError('Protein count changed')
            for pi,seq in enumerate(records):
                key=f'{prefix}{gi:03d}p{pi:06d}'
                seqs[key],genomes[key]=seq,item['sample_id']
        filename='training.faa' if prefix=='s' else 'queries.faa'
        regenerated=''.join(f'>{p}\n{s}\n' for p,s in seqs.items()).encode()
        if regenerated!=(work/filename).read_bytes():
            raise ValueError('Normalized FASTA differs from bound inputs')
        inventories[role]=(names,seqs,genomes)
    names,train,genomes=inventories['vocabulary_fit_only']
    lengths={p:len(s) for p,s in train.items()}
    families=protein_features.read_families((work/'clusters.tsv').read_text().splitlines(),train)
    self_hits,_=protein_features.best_hits((work/'self.tsv').read_text().splitlines(),lengths,lengths,families)
    fitted,patterns=pilot.fit(names,train,genomes,self_hits)
    if fitted!=members or {k:list(v) for k,v in patterns.items()}!=expected:
        raise ValueError('Training-only vocabulary replay differs')
    for role,filename,field in [('vocabulary_fit_only','self.tsv','training_projection'),
                                ('projection_only','query.tsv','target_projection')]:
        names,seqs,genomes=inventories[role]
        hits,_=protein_features.best_hits((work/filename).read_text().splitlines(),{p:len(s) for p,s in seqs.items()},lengths,families)
        projected=pilot.project(names,seqs,genomes,hits,members,expected,policy['threshold'])
        if projected!=report[field]:
            raise ValueError(f'Projection replay differs: {field}')
    files=['training.faa','queries.faa','clusters.tsv','self.tsv','query.tsv','vocabulary.json','vocabulary-freeze.json']
    artifacts=[{'path':str(work/name),'sha256':pilot.sha(work/name)} for name in files]
    receipt={'schema_version':1,'status':'exact_training_vocabulary_and_both_projection_replay',
             'report':{'path':str(args.report),'sha256':pilot.sha(args.report)},
             'script_sha256':pilot.sha(__file__),'artifacts':artifacts,
             'sequence_searches_or_models_executed':False,'bound_input_count':len(bindings),
             'all_source_and_target_rows_and_match_fractions_equal':True}
    pilot.dump(args.output,receipt)
    print(json.dumps({'status':receipt['status'],'artifacts':len(artifacts)}))


if __name__=='__main__':
    main()
