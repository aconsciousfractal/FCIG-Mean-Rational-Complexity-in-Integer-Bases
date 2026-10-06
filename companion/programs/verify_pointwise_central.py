"""Distinct all-threshold central check of the reconstructed pointwise premises."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from archive import data
from verify_all_central import enumerate_central, strict_load, fresh_output


def need(ok, message):
    if not ok:raise ValueError(message)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind',required=True,choices=('semiprime','base6','sixth'))
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--input-sha256',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--start',type=int,default=0)
    parser.add_argument('--stop',type=int)
    args=parser.parse_args()
    raw=args.input.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==args.input_sha256,'fresh pointwise premise byte pin')
    payload=strict_load(raw);rows=payload['rows']
    need(payload.get('status') == {'semiprime':'PORTABLE_E032_RECONSTRUCTION_PASS','base6':'E029_PASS','sixth':'E030_PASS'}[args.kind], 'completed pointwise source status')
    need(type(rows) is list and all(type(r) is dict and type(r.get('M')) is int for r in rows), 'exact whole pointwise domain')
    expected=([v['M'] for v in data(28,'SEMIPRIME_DOMAIN.json')['rows']] if args.kind=='semiprime'
              else [6**n for n in range(1,9)] if args.kind=='base6' else [b**6 for b in (10,12,14)])
    need([r['M'] for r in rows]==expected and len(set(expected))==len(expected),'complete original pointwise domain')
    stop=len(rows) if args.stop is None else args.stop
    need(0<=args.start<stop<=len(rows),'nonempty bounded domain chunk')
    need(not args.output.exists() and not args.output.is_symlink() and args.output.parent.is_dir(),'fresh output file')
    results=[]
    for index in range(args.start,stop):
        row=rows[index];M=row['M'];actual=enumerate_central(M);k=actual['k']
        hist=row['upper_histogram'] if args.kind=='semiprime' else row['histogram']
        need(all(type(x) is list and len(x)==2 and all(type(v) is int for v in x) and x[0]>0 and x[1]>0 for x in hist),'exact pointwise histogram')
        need([h for h,c in hist]==sorted({h for h,c in hist}) and sum(c for h,c in hist)==M,'histogram order and mass')
        bins={h:c for h,c in hist if h<=k};prefix=0
        for T in range(1,k+1):
            prefix+=bins.get(T,0)
            need(prefix==actual['cdf'][T],'every pointwise central prefix M='+str(M)+' T='+str(T))
        results.append({'domain_index':index,'M':M,'k':k,'central_states':actual['cdf'][k],
                        'central_vector_sha256':actual['central_vector_sha256'],
                        'central_cdf_sha256':hashlib.sha256(json.dumps(actual['cdf'],separators=(',',':')).encode()).hexdigest(),
                        'counters':actual['counters'],'status':'ALL_POINTWISE_CENTRAL_THRESHOLDS_PASS'})
        if (index+1)%100==0 or index+1==stop:print('POINTWISE '+args.kind+' '+str(index+1)+'/'+str(stop),file=sys.stderr,flush=True)
    result={'status':'COMPLETE_POINTWISE_CENTRAL_CHUNK_PASS','kind':args.kind,'input_sha256':args.input_sha256,
            'start':args.start,'stop':stop,'full_domain_count':len(rows),'rows':results,
            'complete_domain_verified':args.start==0 and stop==len(rows),'all_thresholds_checked':sum(r['k'] for r in results),
            'producer_imported':False,'true_far_minima_computed':False}
    with fresh_output(args.output) as stream:json.dump(result,stream,sort_keys=True,separators=(',',':'));stream.write('\n')
    print(json.dumps({k:result[k] for k in ('status','kind','start','stop','complete_domain_verified','all_thresholds_checked')},sort_keys=True,separators=(',',':')))


if __name__=='__main__':main()
