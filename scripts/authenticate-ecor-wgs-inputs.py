"""Authenticate official ECOR WGS contig sets without reading host-range outcomes."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/ecor-wgs-input-policy.json'


def sha(data):return hashlib.sha256(data).hexdigest()


def parse_master(data,host,accession,project):
    text=data.decode('utf-8')
    match=re.search(r'^ID\s+([A-Z]{4})(\d{2})000000;.*?; (\d+) SQ\.$',text,re.M)
    if not match or match[1]!=accession[:4]:raise ValueError('WGS master identity missing/different')
    if not re.search(r'^PR\s+Project:'+re.escape(project)+r';',text,re.M):raise ValueError('WGS project differs')
    strain='MOD1-'+host
    if f'/strain="{strain}"' not in text:raise ValueError('WGS strain differs')
    ranges=re.findall(r'^WGS\s+([A-Z]{4}\d{8})-([A-Z]{4}\d{8})\s*$',text,re.M)
    if not ranges:raise ValueError('WGS contig ranges missing')
    expected=[]
    for first,last in ranges:
        if first[:6]!=match[1]+match[2] or last[:6]!=first[:6] or int(last[6:])<int(first[6:]):
            raise ValueError('WGS contig range differs from master')
        expected.extend(first[:6]+f'{i:06d}' for i in range(int(first[6:]),int(last[6:])+1))
    if len(expected)!=int(match[3]) or len(set(expected))!=len(expected):raise ValueError('WGS master range count differs')
    samples=re.findall(r'^DR\s+BioSample;\s+(\S+)\.',text,re.M)
    return {'wgs_master':match[1]+match[2]+'000000','wgs_build':match[2],
            'set_id':match[1]+match[2],'declared_contigs':len(expected),'expected_contigs':expected,
            'strain':strain,'project':project,'biosample':samples[0] if len(samples)==1 else None}


def parse_fasta(data,master):
    records=[]
    header=None;parts=[]
    def finish():
        if header is None:return
        match=re.match(r'ENA\|([A-Z]{4}\d{8})\|\1\.(\d+)\s+(.+)',header)
        if not match or f"strain {master['strain']} " not in match[3]:raise ValueError('Unversioned or mismatched FASTA header')
        sequence=''.join(parts).upper()
        if not sequence or not set(sequence)<=set('ACGTRYSWKMBDHVN'):raise ValueError('Empty or invalid DNA contig')
        counts=Counter(sequence)
        records.append({'accession':match[1],'sequence_version':int(match[2]),
                        'versioned_accession':match[1]+'.'+match[2],'length':len(sequence),
                        'sequence_sha256':sha(sequence.encode()),'ambiguous_bases':sum(n for c,n in counts.items() if c not in 'ACGT')})
    for line in data.decode('ascii').splitlines():
        if line.startswith('>'):
            finish();header=line[1:];parts=[]
        elif line:
            if header is None or any(c.isspace() for c in line):raise ValueError('Malformed FASTA')
            parts.append(line)
    finish()
    ids=[r['accession'] for r in records]
    if len(ids)!=len(set(ids)) or set(ids)!=set(master['expected_contigs']):raise ValueError('Missing, duplicate or unexpected WGS contigs')
    return records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,default=ROOT/'data/raw/genophi-external-validation/ecor-wgs-inputs-v1')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/ecor-wgs-input-authentication.json')
    args=parser.parse_args()
    policy=json.loads(POLICY.read_text())
    for item in policy['inputs']+policy['pilot_files']:
        if sha((ROOT/item['path']).read_bytes())!=item['sha256']:raise ValueError('Bound input differs')
    authority=json.loads((ROOT/policy['inputs'][0]['path']).read_text())
    roster=authority['host_id_to_wgs_master_accession']
    with (ROOT/policy['inputs'][1]['path']).open(newline='') as stream:metadata={r['Strain']:r for r in csv.DictReader(stream)}
    if len(roster)!=25 or any(metadata[h]['NCBI accession']!=a for h,a in roster.items()):raise ValueError('Frozen host roster differs')
    if args.workdir.exists() or args.output.exists():raise ValueError('Use fresh workdir and output')
    args.workdir.mkdir(parents=True)
    def fetch(url,path):
        with urlopen(Request(url,headers={'User-Agent':'TipCraft-source-audit/0.1'}),timeout=policy['timeout_seconds']) as response:
            data=response.read(policy['max_download_bytes']+1)
            if len(data)>policy['max_download_bytes']:raise ValueError('Download exceeds fixed byte limit')
            receipt={'url':url,'final_url':response.url,'http_status':response.status,'retrieved_utc':datetime.now(timezone.utc).isoformat()}
        path.write_bytes(data)
        return data,receipt
    def process(pair):
        host,accession=pair
        dest=args.workdir/host;dest.mkdir()
        result={'host':host,'source_wgs_accession':accession,'status':'pending','source_scaffolds':int(metadata[host]['Scaffolds'])}
        try:
            master_url=policy['master_url_template'].format(accession=accession)
            if host=='ECOR12':
                raw=(ROOT/policy['pilot_files'][0]['path']).read_bytes();(dest/'master.embl').write_bytes(raw)
                master_receipt={'url':master_url,'reused_pilot':policy['pilot_files'][0]}
            else:raw,master_receipt=fetch(master_url,dest/'master.embl')
            master=parse_master(raw,host,accession,policy['expected_project'])
            result['master']={k:v for k,v in master.items() if k!='expected_contigs'}
            result['master_file']={'path':str(dest/'master.embl'),'sha256':sha(raw),'retrieval':master_receipt}
            result['source_scaffold_count_matches_master']=result['source_scaffolds']==master['declared_contigs']
            url=policy['fasta_url_template'].format(prefix3=accession[:3].lower(),set_id=master['set_id'])
            if host=='ECOR12':
                compressed=(ROOT/policy['pilot_files'][1]['path']).read_bytes();(dest/'contigs.fasta.gz').write_bytes(compressed)
                retrieval={'url':url,'reused_pilot':policy['pilot_files'][1]}
            else:compressed,retrieval=fetch(url,dest/'contigs.fasta.gz')
            expanded=gzip.decompress(compressed)
            records=parse_fasta(expanded,master)
            (dest/'contigs.fasta').write_bytes(expanded)
            manifest={'host':host,'source_wgs_accession':accession,'master':master,'records':records,
                      'DNA_concatenated_trimmed_or_designed':False}
            manifest_path=dest/'contig-manifest.json';manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
            result.update({'status':'authenticated_official_wgs_contig_set','contigs':len(records),
                           'total_bases':sum(r['length'] for r in records),'ambiguous_bases':sum(r['ambiguous_bases'] for r in records),
                           'sequence_versions':sorted({r['sequence_version'] for r in records}),
                           'fasta_gzip':{'path':str(dest/'contigs.fasta.gz'),'bytes':len(compressed),'sha256':sha(compressed),'retrieval':retrieval},
                           'fasta':{'path':str(dest/'contigs.fasta'),'bytes':len(expanded),'sha256':sha(expanded)},
                           'contig_manifest':{'path':str(manifest_path),'sha256':sha(manifest_path.read_bytes())}})
        except Exception as exc:
            result.update({'status':'unresolved_input','error':str(exc)})
        (dest/'receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        print(host,result['status'],flush=True)
        return result
    with ThreadPoolExecutor(max_workers=policy['workers']) as pool:entries=list(pool.map(process,sorted(roster.items())))
    report={'schema_version':1,'status':'completed_source_authentication_with_explicit_per_host_status',
            'script_sha256':sha(Path(__file__).read_bytes()),'policy_sha256':sha(POLICY.read_bytes()),
            'inputs':policy['inputs'],'pilot_files':policy['pilot_files'],'entries':entries,
            'roster_hosts':len(roster),'authenticated_hosts':sum(e['status']=='authenticated_official_wgs_contig_set' for e in entries),
            'unresolved_hosts':[e['host'] for e in entries if e['status']!='authenticated_official_wgs_contig_set'],
            'source_scaffold_metric_mismatches':[e['host'] for e in entries if e.get('source_scaffold_count_matches_master') is False],
            'host_range_outcomes_loaded':False,'experimental_stock_identity_verified':False,
            'training_admission_or_prediction':False,'new_original_benchmark_coverage':0,
            'GCA_accessions_inferred_for_all_hosts':False,'Picard_collection_downloaded':False}
    for item in policy['inputs']+policy['pilot_files']:
        if sha((ROOT/item['path']).read_bytes())!=item['sha256']:raise ValueError('Bound input changed during execution')
    with args.output.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({'authenticated_hosts':report['authenticated_hosts'],'unresolved_hosts':report['unresolved_hosts']}),flush=True)


if __name__=='__main__':main()
