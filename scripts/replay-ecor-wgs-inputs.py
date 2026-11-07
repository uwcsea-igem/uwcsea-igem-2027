"""Offline replay of every authenticated ECOR contig and provenance binding."""
import argparse
from collections import defaultdict
import gzip
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,default=ROOT/'data/processed/ecor-wgs-input-authentication.json')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/ecor-wgs-input-replay.json')
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('ecor_input_producer',ROOT/'scripts/authenticate-ecor-wgs-inputs.py')
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    sha=lambda path:producer.sha(Path(path).read_bytes())
    report=json.loads(args.report.read_text())
    if sha(spec.origin)!=report['script_sha256'] or sha(producer.POLICY)!=report['policy_sha256']:
        raise ValueError('Executed producer or policy differs')
    for item in report['inputs']+report['pilot_files']:
        if sha(ROOT/item['path'])!=item['sha256']:raise ValueError('Bound source differs')
    source=json.loads((ROOT/report['inputs'][0]['path']).read_text())
    expected=source['host_id_to_wgs_master_accession']
    hosts=[entry['host'] for entry in report['entries']]
    if len(hosts)!=len(set(hosts)) or set(hosts)!=set(expected):raise ValueError('Host inventory incomplete or duplicate')
    policy=json.loads(producer.POLICY.read_text())
    artifacts=[];fingerprints={};duplicates=defaultdict(list)
    total_contigs=total_bases=ambiguous=0
    for entry in report['entries']:
        host=entry['host']
        if entry['status']!='authenticated_official_wgs_contig_set':raise ValueError('Unresolved host; full replay unavailable')
        for field in ['master_file','fasta_gzip','fasta','contig_manifest']:
            item=entry[field]
            if sha(item['path'])!=item['sha256']:raise ValueError('Artifact changed: '+item['path'])
            artifacts.append({'host':host,'kind':field,'path':item['path'],'sha256':item['sha256']})
        master=producer.parse_master(Path(entry['master_file']['path']).read_bytes(),host,expected[host],policy['expected_project'])
        expanded=gzip.decompress(Path(entry['fasta_gzip']['path']).read_bytes())
        if expanded!=Path(entry['fasta']['path']).read_bytes():raise ValueError('Expanded FASTA differs')
        records=producer.parse_fasta(expanded,master)
        stored=json.loads(Path(entry['contig_manifest']['path']).read_text())
        if stored['master']!=master or stored['records']!=records or stored['DNA_concatenated_trimmed_or_designed']:
            raise ValueError('Full record inventory replay differs')
        bases=sum(r['length'] for r in records);amb=sum(r['ambiguous_bases'] for r in records)
        if (entry['contigs'],entry['total_bases'],entry['ambiguous_bases'])!=(len(records),bases,amb):
            raise ValueError('Public summary differs from records')
        parts=sorted((r['length'],r['sequence_sha256']) for r in records)
        fingerprint=producer.sha(json.dumps(parts,separators=(',',':')).encode())
        fingerprints[host]=fingerprint;duplicates[fingerprint].append(host)
        total_contigs+=len(records);total_bases+=bases;ambiguous+=amb
    result={'schema_version':1,'status':'all25_official_contig_sets_replayed_offline',
            'report':{'path':str(args.report),'sha256':sha(args.report)},'script_sha256':sha(__file__),
            'verified_raw_artifact_count':len(artifacts),'artifacts':artifacts,
            'hosts':len(hosts),'contigs':total_contigs,'total_bases':total_bases,'ambiguous_bases':ambiguous,
            'contig_sequence_multiset_fingerprints':fingerprints,
            'identical_contig_sequence_multiset_groups':[sorted(gs) for gs in duplicates.values() if len(gs)>1],
            'fingerprint_definition':'SHA256 of compact JSON of sorted (length, sequence SHA256) tuples, preserving duplicate contig multiplicity; no DNA concatenation',
            'fingerprint_limit':'Not reverse-complement or fragmentation invariant; no-hit is not biological independence or relatedness clearance',
            'network_requests_or_models_executed':False,'stock_identity_verified':False}
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['status','hosts','contigs','total_bases','ambiguous_bases','verified_raw_artifact_count']}))


if __name__=='__main__':main()
