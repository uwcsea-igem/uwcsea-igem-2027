"""Replay official-record identity and author correspondence as separate claims."""
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    import Bio
    if Bio.__version__!='1.86':raise ValueError('Require Biopython1.86')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,default=ROOT/'data/processed/basel-phage-accession-authentication.json')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/basel-phage-accession-replay.json')
    args=parser.parse_args()
    path=ROOT/'scripts/authenticate-basel-phage-accessions.py'
    spec=importlib.util.spec_from_file_location('producer',path)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    sha=lambda p:producer.sha(Path(p).read_bytes())
    report=json.loads(args.report.read_text())
    if sha(path)!=report['script_sha256'] or sha(producer.POLICY)!=report['policy_sha256']:
        raise ValueError('Executed producer or policy differs')
    for item in report['inputs']:
        if sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Bound source changed')
    authority=json.loads((ROOT/report['inputs'][0]['path']).read_text())
    roster=authority['phage_id_to_accession']
    ids=[e['phage'] for e in report['entries']]
    if len(ids)!=len(set(ids)) or set(ids)!=set(roster):raise ValueError('Full52 roster not preserved')
    observations=[];artifacts=[]
    for entry in report['entries']:
        phage=entry['phage'];source=entry['author_container']
        if sha(source['path'])!=source['sha256']:raise ValueError('Author container changed')
        author=producer.sequence_record(Path(source['path']).read_bytes(),'genbank')
        if author['sequence_sha256']!=source['sequence_sha256'] or author['length']!=source['length']:
            raise ValueError('Author DNA summary differs')
        artifacts.append({'phage':phage,'kind':'author_container','path':source['path'],'sha256':source['sha256']})
        formats={}
        for item in entry['official_payloads']:
            if item['kind'] in formats or sha(item['path'])!=item['sha256']:raise ValueError('Official artifact changed/duplicated')
            formats[item['kind']]=producer.sequence_record(Path(item['path']).read_bytes(),item['format'],roster[phage])
            artifacts.append({'phage':phage,'kind':item['kind'],'path':item['path'],'sha256':item['sha256']})
        if set(formats)!={'annotation','fasta'}:raise ValueError('Both official formats required for this replay')
        official=formats['annotation'];fasta=formats['fasta']
        if official['versioned_accession']!=fasta['versioned_accession'] or official['sequence']!=fasta['sequence']:
            raise ValueError('Official formats disagree')
        exact=official['sequence']==author['sequence']
        if exact!=entry['status'].startswith('authenticated_'):raise ValueError('Original correspondence status replay differs')
        if not exact and entry.get('error')!='Official DNA differs from authenticated author container':
            raise ValueError('Original unresolved reason differs')
        observations.append({'phage':phage,'requested_accession':roster[phage],
            'official_versioned_accession':official['versioned_accession'],'official_formats_agree':True,
            'official_length':official['length'],'official_sequence_sha256':official['sequence_sha256'],
            'official_ambiguous_bases':official['ambiguous_bases'],'author_length':author['length'],
            'author_sequence_sha256':author['sequence_sha256'],'exact_author_DNA_match':exact,
            'old_author_sequence_groups_reusable_by_exact_DNA_correspondence_only':exact})
    mismatches=[o['phage'] for o in observations if not o['exact_author_DNA_match']]
    if report['authenticated_phages']!=len(observations)-len(mismatches) or report['unresolved_phages']!=mismatches:
        raise ValueError('Public count/status replay differs')
    result={'schema_version':1,'status':'complete_official_format_and_author_correspondence_replay',
            'report':{'path':str(args.report),'sha256':sha(args.report)},'script_sha256':sha(__file__),
            'parser':'Biopython '+Bio.__version__,'official_record_format_pairs_verified':len(observations),
            'exact_author_DNA_matches':len(observations)-len(mismatches),'author_DNA_mismatches':mismatches,
            'verified_artifact_count':len(artifacts),'artifacts':artifacts,'observations':observations,
            'official_sequence_version_counts':dict(Counter(o['official_versioned_accession'].rsplit('.',1)[1] for o in observations)),
            'network_requests_models_or_sequence_edits_executed':False,
            'interpretation':'All official records can be authenticated even when three correspondences to the old author cache fail. Those three official genomes need their own relatedness checks before old groups can be reused; neither representation is declared erroneous or stock-identical.'}
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['official_record_format_pairs_verified','exact_author_DNA_matches','author_DNA_mismatches','verified_artifact_count']}))


if __name__=='__main__':main()
