"""Authenticate all52 published phage accessions against official and author DNA."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import threading
import time
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
POLICY=ROOT/'software/tipcraft/benchmarks/receptor-resource/basel-phage-accession-policy.json'


def sha(data):return hashlib.sha256(data).hexdigest()


def sequence_record(data,fmt,expected_accession=None):
    from Bio import SeqIO
    records=list(SeqIO.parse(io.StringIO(data.decode('utf-8')),fmt))
    if len(records)!=1:raise ValueError('Require exactly one complete sequence record')
    record=records[0]
    sequence=str(record.seq).upper()
    if not sequence or not set(sequence)<=set('ACGTRYSWKMBDHVN'):
        raise ValueError('Empty, undefined or invalid nucleotide sequence')
    identifier=record.id.split('|')[-1]
    if expected_accession is not None and not re.fullmatch(re.escape(expected_accession)+r'\.\d+',identifier):
        raise ValueError('Official accession/version does not match request')
    return {'record_id':record.id,'versioned_accession':identifier if expected_accession is not None else None,
            'sequence':sequence,'sequence_sha256':sha(sequence.encode()),'length':len(sequence),
            'ambiguous_bases':sum(c not in 'ACGT' for c in sequence),'topology':record.annotations.get('topology')}


def compare_records(annotation,fasta,author):
    if annotation['versioned_accession']!=fasta['versioned_accession'] or annotation['sequence']!=fasta['sequence']:
        raise ValueError('Official formats disagree on accession version or DNA')
    if annotation['sequence']!=author['sequence']:
        raise ValueError('Official DNA differs from authenticated author container')
    return True


def main():
    import Bio
    if Bio.__version__!='1.86':raise ValueError('Frozen parser requires Biopython1.86')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,default=ROOT/'data/raw/genophi-external-validation/basel-phage-accessions-v1')
    parser.add_argument('--output',type=Path,default=ROOT/'data/processed/basel-phage-accession-authentication.json')
    args=parser.parse_args()
    policy=json.loads(POLICY.read_text());script_hash=sha(Path(__file__).read_bytes())
    for item in policy['inputs']:
        if sha((ROOT/item['path']).read_bytes())!=item['sha256']:raise ValueError('Bound source changed')
    authority=json.loads((ROOT/policy['inputs'][0]['path']).read_text())
    cohort=json.loads((ROOT/policy['inputs'][1]['path']).read_text())
    pilot=json.loads((ROOT/policy['inputs'][2]['path']).read_text())
    roster=authority['phage_id_to_accession']
    authors={e['source']['primary_id']:e for e in cohort['entries']}
    probes={e['phage']:e for e in pilot['observations']}
    if len(roster)!=52 or set(roster)-authors.keys() or set(probes)!={'Bas01','Bas12','Bas14'}:
        raise ValueError('Frozen roster or probe set differs')
    if args.workdir.exists() or args.output.exists():raise ValueError('Use fresh destinations')
    args.workdir.mkdir(parents=True)
    lock=threading.Lock();last_ncbi=[0.0]
    def fetch(url,path):
        if 'ncbi.nlm.nih.gov/' in url:
            with lock:
                delay=policy['ncbi_minimum_request_spacing_seconds']-(time.monotonic()-last_ncbi[0])
                if delay>0:time.sleep(delay)
                last_ncbi[0]=time.monotonic()
        with urlopen(Request(url,headers={'User-Agent':'TipCraft-source-audit/0.1'}),timeout=policy['timeout_seconds']) as response:
            data=response.read(policy['max_download_bytes']+1)
            if len(data)>policy['max_download_bytes']:raise ValueError('Response exceeds fixed byte limit')
            receipt={'requested_url':url,'final_url':response.url,'http_status':response.status,
                     'retrieved_utc':datetime.now(timezone.utc).isoformat()}
        path.write_bytes(data)
        return data,receipt
    def process(pair):
        phage,accession=pair
        dest=args.workdir/phage;dest.mkdir()
        result={'phage':phage,'requested_accession':accession,'status':'unresolved_input','official_payloads':[]}
        try:
            source=authors[phage];path=Path(source['container_path']);raw=path.read_bytes()
            if sha(raw)!=source['sha256']:raise ValueError('Author container hash differs')
            author=sequence_record(raw,'genbank')
            result['author_container']={'path':str(path),'sha256':sha(raw),'record_id':author['record_id'],
                                        'sequence_sha256':author['sequence_sha256'],'length':author['length']}
            route='RefSeq_NC' if accession.startswith('NC_') else 'PZ_MZ'
            parsed={}
            for kind in ['annotation','fasta']:
                url=policy['routes'][route][kind].format(accession=accession)
                fmt='fasta' if kind=='fasta' else 'genbank' if route=='RefSeq_NC' else 'embl'
                target=dest/(kind+'.txt')
                if phage in probes:
                    prior=next(p for p in probes[phage]['official_payloads'] if ('_fasta' in Path(p['path']).name)==(kind=='fasta'))
                    data=(ROOT/prior['path']).read_bytes()
                    if sha(data)!=prior['file_sha256'] or url!=prior['url']:raise ValueError('Probe bytes or route differ')
                    target.write_bytes(data);retrieval={'requested_url':url,'reused_pilot_payload':prior}
                else:data,retrieval=fetch(url,target)
                result['official_payloads'].append({'kind':kind,'format':fmt,'path':str(target),'sha256':sha(data),
                                                    'bytes':len(data),'retrieval':retrieval})
                parsed[kind]=sequence_record(data,fmt,accession)
            compare_records(parsed['annotation'],parsed['fasta'],author)
            result.update({'status':'authenticated_official_formats_and_author_DNA_exact_match',
                           'versioned_accession':parsed['annotation']['versioned_accession'],
                           'sequence_sha256':author['sequence_sha256'],'length':author['length'],
                           'ambiguous_bases':author['ambiguous_bases'],'official_topology':parsed['annotation']['topology']})
        except Exception as exc:result['error']=str(exc)
        (dest/'receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        print(phage,result['status'],flush=True)
        return result
    with ThreadPoolExecutor(max_workers=policy['workers']) as pool:entries=list(pool.map(process,sorted(roster.items())))
    result={'schema_version':1,'status':'completed_full_roster_accession_authentication_with_per_record_status',
            'script_sha256':script_hash,'policy_sha256':sha(POLICY.read_bytes()),'inputs':policy['inputs'],
            'parser':'Biopython '+Bio.__version__,'entries':entries,'roster_phages':len(roster),
            'authenticated_phages':sum(e['status'].startswith('authenticated_') for e in entries),
            'unresolved_phages':[e['phage'] for e in entries if not e['status'].startswith('authenticated_')],
            'interaction_outcomes_loaded':False,'experimental_stock_identity_verified':False,
            'protein_annotation_equivalence_verified':False,'training_or_prediction':False,
            'new_original_benchmark_coverage':0,'DNA_edited_rotated_or_reverse_complemented':False}
    for item in policy['inputs']:
        if sha((ROOT/item['path']).read_bytes())!=item['sha256']:raise ValueError('Source changed during execution')
    if sha(Path(__file__).read_bytes())!=script_hash:raise ValueError('Executed producer changed')
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({'authenticated_phages':result['authenticated_phages'],'unresolved_phages':result['unresolved_phages']}),flush=True)


if __name__=='__main__':main()
