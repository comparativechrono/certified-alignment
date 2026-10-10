#!/usr/bin/env python3
"""Resolve source objects from the public HPRC bucket, without decoding reads."""
import argparse, concurrent.futures, csv, hashlib, json
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode
import xml.etree.ElementTree as ET
BASE = 'https://s3-us-west-2.amazonaws.com/human-pangenomics'
NS = {'s': 'http://s3.amazonaws.com/doc/2006-03-01/'}

def sha(b): return hashlib.sha256(b).hexdigest()

def resolve(sample):
    hits=[]
    for project in ['HPRC', 'HPRC_PLUS']:
        prefix=f'working/{project}/{sample}/raw_data/Illumina/child/'
        u=BASE+'?'+urlencode({'list-type':'2','prefix':prefix,'max-keys':1000})
        raw=urlopen(u,timeout=90).read()
        root=ET.fromstring(raw)
        if root.findtext('s:IsTruncated',namespaces=NS)=='true': raise ValueError('Unpaginated listing')
        for obj in root.findall('s:Contents',NS):
            key=obj.findtext('s:Key',namespaces=NS)
            if key.endswith('/'+sample+'.final.cram'):
                hits.append({'sample':sample,'cram_url':BASE+'/'+key,'cram_size':int(obj.findtext('s:Size',namespaces=NS)), 'cram_etag':obj.findtext('s:ETag',namespaces=NS), 'cram_last_modified':obj.findtext('s:LastModified',namespaces=NS),'listing_url':u,'listing_sha256':sha(raw),'listing_object_key':key})
    if len(hits)!=1: raise ValueError(f'{sample}: expected exactly one CRAM, found {hits}')
    out=hits[0]
    for suffix in ['.crai','.md5']:
        b=urlopen(out['cram_url']+suffix,timeout=90).read()
        out[suffix[1:]+'_url']=out['cram_url']+suffix
        out[suffix[1:]+'_bytes']=len(b)
        out[suffix[1:]+'_sha256']=sha(b)
        if suffix=='.md5':out['cram_published_md5_text']=b.decode().strip()
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--cases',type=Path,required=True);p.add_argument('--controls',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    samples=set()
    for path in [a.cases,a.controls]:
        with path.open() as f:samples.update(row['sample'] for row in csv.DictReader(f,delimiter='\t'))
    with concurrent.futures.ThreadPoolExecutor(5) as ex: records=list(ex.map(resolve,sorted(samples)))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps({'source':'HPRC public bucket Illumina/child objects linked by Locityper publication data availability','records':records,'caveat':'Public donor-matched Illumina read objects; exact input library equivalence to archived Locityper calls is not asserted.'},indent=2)+'\n')
    print(json.dumps({'samples':len(records),'output':str(a.output)}))
if __name__=='__main__':main()
