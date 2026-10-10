#!/usr/bin/env python3
"""Acquire donor-matched public CRAM regional records, retaining sequence/qualities.

No mapping-quality, proper-pair, duplicate or QC filtering is applied here.
Records outside the declared windows, including globally unmapped reads, are not
searched. This is a regional-ascertainment audit, not complete WGS recruitment.
"""
import argparse, concurrent.futures, csv, datetime, hashlib, json, os, sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import urlopen, Request
import pysam

REFERENCE = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'

def sha_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def download(url,path,expected_sha=None):
    if not path.exists():
        tmp=path.with_suffix(path.suffix+'.part')
        with urlopen(url,timeout=120) as src,tmp.open('wb') as out:
            while b:=src.read(1048576):out.write(b)
        tmp.replace(path)
    digest=sha_file(path)
    if expected_sha and digest!=expected_sha:raise ValueError(f'Checksum mismatch: {path}')
    return digest

def parse_bed(path):
    intervals=[]
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split();intervals.append((f[0],int(f[1]),int(f[2])))
    if not intervals:raise ValueError(f'Empty BED: {path}')
    return intervals

def merge_intervals(intervals,flank,lengths):
    by=defaultdict(list)
    for c,s,e in intervals:
        if c not in lengths or not 0<=s<e<=lengths[c]:raise ValueError(f'Invalid GRCh38 interval {c}:{s}-{e}')
        by[c].append((max(0,s-flank),min(lengths[c],e+flank)))
    out=[]
    for c in by:
        merged=[]
        for s,e in sorted(by[c]):
            if merged and s<=merged[-1][1]:merged[-1]=(merged[-1][0],max(e,merged[-1][1]))
            else:merged.append((s,e))
        out.extend((c,s,e) for s,e in merged)
    return out

def acquire_sample(job):
    sample,rows,source,args=job
    cache=Path(args['cache']);output=Path(args['output']);db=Path(args['db'])
    index=cache/(sample+'.cram.crai')
    indexsha=download(source['crai_url'],index,source['crai_sha256'])
    with urlopen(Request(source['cram_url'],method='HEAD'),timeout=90) as head:
        if head.headers.get('ETag')!=source['cram_etag'] or int(head.headers['Content-Length'])!=source['cram_size']:raise ValueError(f'Source object metadata changed: {sample}')
    outrows=[];errors=[]
    for row in rows:
        locus=row['locus'];name=sample+'__'+locus
        try:
            # A fresh handle per locus prevents a failed network read poisoning later windows.
            with pysam.AlignmentFile(source['cram_url'],'rc',index_filename=str(index),reference_filename=args['reference_fasta'],threads=2) as remote:
                h=remote.header.to_dict();sq={r['SN']:r for r in h['SQ']}; lengths={c:r['LN'] for c,r in sq.items()}
                rg=h.get('RG',[])
                if not rg or any(r.get('SM')!=sample for r in rg):raise ValueError(f'Unexpected donor in read group: {sample}')
                bed=db/'loci'/locus/'ref.bed';unexpanded=parse_bed(bed)
                regions=merge_intervals(unexpanded,args['flank'],lengths)
                for chrom in dict.fromkeys(r[0] for r in regions):
                    expected=args['reference_dictionary'][chrom]
                    if sq[chrom].get('M5')!=expected['M5'] or lengths[chrom]!=int(expected['LN']):raise ValueError(f'CRAM/reference dictionary mismatch: {sample} {chrom}')
                target=output/(name+'.bam');meta=output/(name+'.metadata.json')
                if target.exists() and meta.exists():
                    old=json.loads(meta.read_text())
                    checks=[old.get('acquisition_verified') is True,old.get('counts',{}).get('all_records',0)>0,old['bam_sha256']==sha_file(target),old['ref_bed_sha256']==sha_file(bed),old['flank_bp']==args['flank'],old['cram_url']==source['cram_url'],old['cram_etag']==source['cram_etag'],old['crai_sha256']==indexsha,old['reference_url']==REFERENCE,old['reference_fai_sha256']==args['reference_fai_sha256'],old.get('reference_local_sha256')==args['reference_local_sha256']]
                    if not all(checks):raise ValueError(f'Stale or unverified acquisition: {name}')
                    outrows.append(old);continue
                tmp=target.with_suffix('.bam.part');counts=Counter();qflags=defaultdict(set);seen=set()
                with pysam.AlignmentFile(str(tmp),'wb',header=remote.header) as dst:
                    for chrom,start,end in regions:
                        for read in remote.fetch(chrom,start,end):
                            key=(read.query_name,read.flag,read.reference_id,read.reference_start,read.cigarstring)
                            if key in seen:continue
                            seen.add(key);dst.write(read);counts['all_records']+=1
                            if read.is_unmapped:counts['unmapped_records_in_windows']+=1
                            if read.is_duplicate:counts['duplicate_records']+=1
                            if read.is_qcfail:counts['qcfail_records']+=1
                            if read.is_secondary:counts['secondary_records']+=1
                            if read.is_supplementary:counts['supplementary_records']+=1
                            if not read.is_secondary and not read.is_supplementary:
                                counts['primary_records']+=1
                                if read.is_paired:qflags[read.query_name].add(1 if read.is_read1 else 2 if read.is_read2 else 0)
                if counts['all_records']==0:raise ValueError('unresolved_inadequate_retrieval: zero records; not biological absence')
                tmp.replace(target)
                sortedbam=target.with_suffix('.sorted.bam')
                pysam.sort('-o',str(sortedbam),str(target));sortedbam.replace(target)
                pysam.index(str(target))
                counts['primary_paired_qnames']=len(qflags)
                counts['primary_paired_qnames_with_both_ends']=sum(v=={1,2} for v in qflags.values())
                counts['primary_paired_qnames_without_both_ends']=sum(v!={1,2} for v in qflags.values())
                # Independent sequential re-read verifies the complete local BAM and record count.
                with pysam.AlignmentFile(str(target),'rb') as check:
                    observed=sum(1 for _ in check.fetch(until_eof=True))
                if observed!=counts['all_records']:raise ValueError(f'Local BAM record count mismatch: {name}')
                record={
                    'sample':sample,'locus':locus,'role':row['role'],'status':'verified','acquisition_verified':True,
                    'assembly':'GRCh38_full_analysis_set_plus_decoy_hla','coordinate_convention':'0-based half-open',
                    'historical_ref_bed_intervals':unexpanded,'fetched_intervals':regions,'flank_bp':args['flank'],'ref_bed_sha256':sha_file(bed),
                    'cram_url':source['cram_url'],'cram_published_md5_text':source['cram_published_md5_text'],
                    'whole_cram_checksum_verification':'Publisher sidecar only; not locally recomputed because only indexed slices were downloaded',
                    'cram_etag':source['cram_etag'],'source_etag_and_size_checked_at_acquisition':True,'cram_size':source['cram_size'],'crai_url':source['crai_url'],'crai_sha256':indexsha,
                    'reference_url':REFERENCE,'reference_fai_sha256':args['reference_fai_sha256'],'reference_dict_sha256':args['reference_dict_sha256'],
                    'reference_local_path':args['reference_fasta'],'reference_local_sha256':args['reference_local_sha256'],
                    'reference_dictionary_m5_match':True,'cram_reference_sequences':[sq[c] for c in dict.fromkeys(r[0] for r in regions)],
                    'read_groups':[{k:r[k] for k in ['ID','SM','PL','PM','LB','DS'] if k in r} for r in rg],
                    'bam_path':str(target.resolve()),'bam_sha256':sha_file(target),'bam_bytes':target.stat().st_size,'bai_sha256':sha_file(str(target)+'.bai'),
                    'counts':dict(counts),'mate_retrieval':'No fetch outside windows; paired qnames without both ends explicitly counted. Globally unmapped and elsewhere-mapped reads absent unless their stored coordinates fall inside windows.',
                    'filters':'None during acquisition','pysam_version':pysam.__version__,'htslib_version':pysam.__samtools_version__,
                    'retrieved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()
                }
                meta.write_text(json.dumps(record,indent=2)+'\n');outrows.append(record)
                print(json.dumps({'sample':sample,'locus':locus,'bam_bytes':record['bam_bytes'],'records':counts['all_records'],'status':'verified'}),flush=True)
        except Exception as exc:
            err={'sample':sample,'locus':locus,'error':repr(exc),'status':'unresolved_acquisition_failure'}
            errors.append(err);print(json.dumps(err),flush=True)
    return {'records':outrows,'errors':errors}

def main():
    p=argparse.ArgumentParser();p.add_argument('--cases',type=Path,required=True);p.add_argument('--controls',type=Path,required=True);p.add_argument('--sources',type=Path,required=True);p.add_argument('--db',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--cache',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--flank',type=int,default=1000);p.add_argument('--workers',type=int,default=4);p.add_argument('--protocol-commit',required=True);p.add_argument('--reference-fasta',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);a.cache.mkdir(parents=True,exist_ok=True)
    fhash=download(REFERENCE+'.fai',a.cache/'reference.fa.fai')
    dhash=download(REFERENCE[:-3]+'.dict',a.cache/'reference.dict')
    refdict={}
    for line in (a.cache/'reference.dict').read_text().splitlines():
        if line.startswith('@SQ\t'):
            fields=dict(x.split(':',1) for x in line.split('\t')[1:]);refdict[fields['SN']]=fields
    rows=[]
    for path in [a.cases,a.controls]:
        with path.open() as f:rows.extend(csv.DictReader(f,delimiter='\t'))
    keys=[(r['sample'],r['locus']) for r in rows]
    if len(set(keys))!=len(keys):raise ValueError('Duplicate acquisition keys')
    sources={r['sample']:r for r in json.loads(a.sources.read_text())['records']};by=defaultdict(list)
    for row in rows:by[row['sample']].append(row)
    args={'output':str(a.output),'cache':str(a.cache),'db':str(a.db),'flank':a.flank,'reference_fai_sha256':fhash,'reference_dict_sha256':dhash,'reference_dictionary':refdict,'reference_fasta':str(a.reference_fasta.resolve()),'reference_local_sha256':sha_file(a.reference_fasta)}
    jobs=[(s,sorted(by[s],key=lambda r:r['locus']),sources[s],args) for s in sorted(by)]
    records=[];errors=[]
    with concurrent.futures.ProcessPoolExecutor(a.workers) as ex:
        futures={ex.submit(acquire_sample,j):j[0] for j in jobs}
        for f in concurrent.futures.as_completed(futures):
            try:
                result=f.result();records.extend(result['records']);errors.extend(result['errors'])
            except Exception as e:
                for row in by[futures[f]]:
                    errors.append({'sample':futures[f],'locus':row['locus'],'error':repr(e)});print(json.dumps(errors[-1]),flush=True)
    a.manifest.parent.mkdir(parents=True,exist_ok=True)
    a.manifest.write_text(json.dumps({'status':'complete' if not errors else 'partial_acquisition','acquisition_mode':'fresh_per_locus_handles_verified_local_reference','protocol_commit':a.protocol_commit,'acquisition_script_sha256':sha_file(__file__),'command_argv':[sys.executable]+sys.argv,'records':sorted(records,key=lambda r:(r['sample'],r['locus'])),'errors':errors,'expected_sample_locus_pairs':len(rows),'reference_fai_sha256':fhash,'reference_dict_sha256':dhash,'interpretation':'Regional read audit. No assertion of identical raw libraries to historical caller input. No global unmapped-read search.'},indent=2)+'\n')
    print(json.dumps({'successful_pairs':len(records),'expected_pairs':len(rows),'errors':errors}),flush=True)
    if errors:raise SystemExit(1)
if __name__=='__main__':main()
