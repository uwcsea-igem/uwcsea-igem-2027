"""Label-free family-conditioned host 5-mer capacity pilot; no interaction model."""
import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import time

from tipcraft import protein_features

POLICY = Path('software/tipcraft/benchmarks/receptor-resource/family-host-capacity-policy.json')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    with Path(path).open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def kmers(seq):
    return {seq[i:i+5] for i in range(len(seq)-4)}


def fit(names, sequences, genomes, hits):
    counts, occurrence = Counter(), defaultdict(set)
    for protein, seq in sequences.items():
        family = hits[protein]['family']
        for kmer in kmers(seq):
            key = (family, kmer)
            counts[key] += 1
            occurrence[key].add(genomes[protein])
    patterns, members, expected = {}, {}, {}
    for key in sorted(counts):
        if counts[key] < 2:
            continue
        pattern = tuple(int(n in occurrence[key]) for n in names)
        if pattern not in patterns:
            feature = f'f{len(patterns)}'
            patterns[pattern] = feature
            expected[feature] = pattern
        members[key] = patterns[pattern]
    return members, expected


def project(names, sequences, genomes, hits, members, expected, threshold):
    by_family = defaultdict(list)
    for (family, kmer), feature in members.items():
        by_family[family].append((kmer, feature))
    reached = {g: set() for g in names}
    matched = {g: set() for g in names}
    for protein, hit in hits.items():
        genome = genomes[protein]
        family = hit['family']
        reached[genome].add(family)
        windows = kmers(sequences[protein])
        for kmer, feature in by_family[family]:
            if kmer in windows:
                matched[genome].add((family, kmer, feature))
    rows, diagnostics = {}, {}
    features = list(expected)
    for genome in names:
        denominator = Counter(feature for family in sorted(reached[genome]) for _, feature in by_family[family])
        numerator = Counter(feature for _, _, feature in matched[genome])
        rows[genome] = [int(denominator[f] > 0 and numerator[f] / denominator[f] >= threshold) for f in features]
        diagnostics[genome] = {f: {'matched_members': numerator[f], 'reached_members': denominator[f],
                                  'fraction': numerator[f]/denominator[f] if denominator[f] else None} for f in features}
    return {'features': features, 'matrix': rows, 'member_match_counts': diagnostics,
            'zero_vector_genomes': [g for g in names if not any(rows[g])]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, default=Path('data/raw/receptor-resource/family-host-capacity-v1'))
    parser.add_argument('--output', type=Path, default=Path('data/processed/receptor-family-host-capacity.json'))
    args = parser.parse_args()
    policy = json.loads(POLICY.read_text())
    inputs = policy['inputs'] + policy['files']
    for item in inputs:
        if sha(item['path']) != item['sha256']:
            raise ValueError(f"Input mismatch: {item['path']}")
    if sha(policy['binary']) != policy['binary_sha256']:
        raise ValueError('MMseqs binary mismatch')
    if args.workdir.exists() or args.output.exists():
        raise ValueError('Use fresh destinations')
    args.workdir.mkdir(parents=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    work = args.workdir
    spec = importlib.util.spec_from_file_location('pilot', 'scripts/check-genophi-phage-features.py')
    pilot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pilot)
    bindings = inputs + [{'path': str(POLICY), 'sha256': sha(POLICY)},
                         {'path': __file__, 'sha256': sha(__file__)},
                         {'path': policy['binary'], 'sha256': sha(policy['binary'])}]
    ledger = []
    start = time.monotonic()

    def run(command, stdin=None):
        log = work / f'step-{len(ledger)+1:02d}.log'
        print('START', len(ledger)+1, command[1], flush=True)
        before = time.monotonic()
        with log.open('xb') as out:
            result = subprocess.run(command, input=stdin, stdout=out, stderr=subprocess.STDOUT,
                                    timeout=policy['command_timeout_seconds'])
        ledger.append({'argv': command, 'returncode': result.returncode, 'log': str(log),
                       'seconds': time.monotonic()-before, 'log_sha256': sha(log)})
        dump(work / f'step-{len(ledger):02d}.json', ledger[-1])
        result.check_returncode()

    def load_files(items, prefix):
        seqs, genome_map, names = {}, {}, []
        for gi, item in enumerate(items):
            names.append(item['sample_id'])
            records = pilot.proteins(Path(item['path']).read_bytes())
            if len(records) != item['expected_proteins']:
                raise ValueError('Protein count mismatch')
            for pi, (_, seq) in enumerate(records):
                if not set(seq) <= set('ACDEFGHIKLMNPQRSTVWY'):
                    raise ValueError('This frozen uniform-call pilot requires standard20 inputs')
                name = f'{prefix}{gi:03d}p{pi:06d}'
                seqs[name], genome_map[name] = seq, item['sample_id']
        return names, seqs, genome_map

    train_files = [f for f in policy['files'] if f['role'] == 'vocabulary_fit_only']
    names, train, train_genomes = load_files(train_files, 's')
    fasta = ''.join(f'>{p}\n{s}\n' for p, s in train.items()).encode()
    (work/'training.faa').write_bytes(fasta)
    binary = policy['binary']
    path = lambda x: str(work/x)
    flags = ['-c','0.8','--min-seq-id','0.4','-s','7.5','--threads','4','-v','3']
    tail = ['--threads','4','-v','3']
    run([binary,'createdb','stdin',path('train_db'),'-v','3'], fasta)
    run([binary,'cluster',path('train_db'),path('clusters'),path('tmp_cluster'),*flags])
    run([binary,'createtsv',path('train_db'),path('train_db'),path('clusters'),path('clusters.tsv'),*tail])
    run([binary,'search',path('train_db'),path('train_db'),path('self_hits'),path('tmp_self'),*flags])
    run([binary,'createtsv',path('train_db'),path('train_db'),path('self_hits'),path('self.tsv'),*tail])
    families = protein_features.read_families((work/'clusters.tsv').read_text().splitlines(), train)
    lengths = {p: len(s) for p,s in train.items()}
    hits, self_diagnostics = protein_features.best_hits((work/'self.tsv').read_text().splitlines(), lengths, lengths, families)
    if set(hits) != set(train):
        raise ValueError('Incomplete self-search')
    members, expected = fit(names, train, train_genomes, hits)
    control = project(names, train, train_genomes, hits, members, expected, policy['threshold'])
    expected_rows = {g: [expected[f][i] for f in expected] for i,g in enumerate(names)}
    control_pass = control['matrix'] == expected_rows and bool(expected)
    dump(work/'vocabulary.json', {'members': [[f,k,v] for (f,k),v in members.items()],
                                 'expected_patterns': expected, 'training_names': names})
    freeze = {'target_search_started': False, 'self_projection_pass': control_pass,
              'bindings': bindings, 'artifacts': [{'path': str(work/n), 'sha256': sha(work/n)} for n in
                  ['training.faa','clusters.tsv','self.tsv','vocabulary.json']]}
    dump(work/'vocabulary-freeze.json', freeze)
    report = {'schema_version': 1, 'policy_sha256': sha(POLICY), 'script_sha256': sha(__file__),
              'status': 'self_projection_failed' if not control_pass else 'capacity_pilot_completed',
              'self_projection_pass': control_pass, 'training_projection': control,
              'training_expected_vectors': expected_rows, 'self_search_diagnostics': self_diagnostics,
              'family_count': len(set(families.values())), 'retained_family_kmer_members': len(members),
              'merged_feature_count': len(expected), 'phenotypes_loaded': False, 'model_fitted': False,
              'original_benchmark_rows_admitted': 0}
    if control_pass:
        for item in bindings + freeze['artifacts']:
            if sha(item['path']) != item['sha256']:
                raise ValueError('Pre-target freeze changed')
        target_files = [f for f in policy['files'] if f['role'] == 'projection_only']
        target_names, targets, target_genomes = load_files(target_files, 'q')
        query_fasta = ''.join(f'>{p}\n{s}\n' for p,s in targets.items()).encode()
        (work/'queries.faa').write_bytes(query_fasta)
        run([binary,'createdb','stdin',path('query_db'),'-v','3'],query_fasta)
        run([binary,'search',path('query_db'),path('train_db'),path('query_hits'),path('tmp_query'),*flags])
        run([binary,'createtsv',path('query_db'),path('train_db'),path('query_hits'),path('query.tsv'),*tail])
        query_hits, diagnostics = protein_features.best_hits((work/'query.tsv').read_text().splitlines(),
            {p:len(s) for p,s in targets.items()},lengths,families)
        report['query_search_diagnostics'] = diagnostics
        report['target_projection'] = project(target_names,targets,target_genomes,query_hits,members,expected,policy['threshold'])
        groups = defaultdict(list)
        for genome,row in report['target_projection']['matrix'].items():
            groups[tuple(row)].append(genome)
        report['target_collision_groups'] = [gs for gs in groups.values() if len(gs)>1]
        report['all_three_targets_distinguishable'] = len(groups)==3
    report.update({'commands': ledger, 'inputs': bindings, 'vocabulary_freeze': {'path': str(work/'vocabulary-freeze.json'),
                   'sha256': sha(work/'vocabulary-freeze.json')}, 'elapsed_seconds': time.monotonic()-start})
    for item in bindings + freeze['artifacts']:
        if sha(item['path']) != item['sha256']:
            raise ValueError('Bound input/artifact changed during pilot')
    dump(args.output, report)
    print(json.dumps({'status':report['status'],'output':str(args.output),'self_projection_pass':control_pass}),flush=True)


if __name__ == '__main__':
    main()
