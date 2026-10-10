"""Witness histogram certification and finite-to-infinite closure."""
import csv,json,hashlib
from pathlib import Path
from math import gcd,prod
from central import *
from secondary_log import contained
from inventory import write_json,write_csv

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def candidate():
    M=210;hs=[];P=1;checked=0
    for s in range(M):
        H=max(1,min(s,M-s))
        for d in range(1,M//2+1):
            if gcd(d,M)!=1:continue
            a=s*d%M;H=min(H,max(min(a,M-a),d))
        for d in range(1,H):
            if gcd(d,M)==1:
                a=s*d%M;need(min(a,M-a)>=H,'candidate minimality');checked+=1
        hs.append(H);P*=H
    need(P**170>210**18690,'candidate separation integer comparison')
    return dict(status='PASS',M=M,heights=hs,product=str(P),minimality_denominators=checked,c0_lower='2/85',integer_gate='P^170 > 210^18690')

def certify(work,out):
    work=Path(work);out=Path(out);rows=[];logs=set();mass=reps=0
    domain=list(csv.DictReader((out/'finite_domain.csv').open()))
    for idx,d in enumerate(domain):
        M=int(d['M']);B=int(d['primitive_base']);g=int(d['max_length'])
        hist=[tuple(map(int,line.split())) for line in (work/f'M{M}.hist').read_text().splitlines()]
        need(hist and sum(n for h,n in hist)==M,'full histogram mass')
        need(all(1<=h<=M//2 and n>0 for h,n in hist),'histogram value range')
        need(all(hist[i][0]<hist[i+1][0] for i in range(len(hist)-1)),'strict histogram order')
        total=sum(n*ticks(h)[1] for h,n in hist);lm=ticks(M)[0]
        upper=F(2*total-M*lm,2*M*LS);target=DELTA*log_lo(B)
        need(upper<target,'FINITE GATE FAIL M='+str(M))
        rows.append({**{k:int(v) for k,v in d.items()},'representatives':M//2+1,'bins':len(hist),'witness_sha256':sha(work/f'M{M}.bin'),'histogram_sha256':sha(work/f'M{M}.hist'),'sum_log_upper_ticks':total,'log_M_lower_ticks':lm,'mu_upper':str(upper),'mu_upper_decimal':out_decimal(upper),'target_lower':str(target)})
        mass+=M;reps+=M//2+1;logs.update(h for h,n in hist);logs.add(M)
        if idx%250==0:print('CERTIFIED',idx,'/',len(domain),flush=True)
    worst=max(rows,key=lambda x:F(x['mu_upper']))
    write_csv(out/'finite_certificate.csv',rows)
    cert=dict(status='PASS',moduli=len(rows),full_state_mass=mass,representatives=reps,all_witness_means_negative=all(F(x['mu_upper'])<0 for x in rows),worst_upper=worst,minimality_claimed=False,distinct_log_inputs=len(logs))
    write_json(out/'finite_summary.json',cert)
    for i,x in enumerate(sorted(logs)):
        lo,hi=ticks(x)
        need(contained(x,lo,hi,LS),'independent log interval at '+str(x))
        if i%5000==0:print('INDEPENDENT LOG',i,'/',len(logs),flush=True)
    write_json(out/'logarithm_checks.json',dict(status='PASS',all_census_log_inputs=len(logs),primary_terms=28,secondary_terms=96,secondary_arithmetic='outward integer fixed-point, scale 10^45',scale=LS,all_inputs_checked=True))
    write_json(out/'candidate.json',candidate())
    routes=[];support_map={}
    for x in csv.DictReader((out/'selected_353.csv').open()):
        Q=int(x['radical']);v=max(F(0),F(x['tail_mu_upper']))/log_lo(Q)
        routes.append((f'tail radical {Q}, M >= {x["tail_cutoff"]}',Q,v))
        support_map[Q]=[v]
    for x in csv.DictReader((out/'analytic_moduli.csv').open()):
        Q=int(x['radical']);v=max(F(0),F(x['mu_upper']))/log_lo(int(x['primitive_base']))
        routes.append((f'analytic M={x["M"]}',Q,v));support_map[Q].append(v)
    for x in rows:
        v=max(F(0),F(x['mu_upper']))/log_lo(x['primitive_base']);Q=x['radical']
        routes.append((f'witness M={x["M"]}',Q,v));support_map[Q].append(v)
    need(all(v<DELTA for label,Q,v in routes),'all full family routes strictly separated')
    label,worst_Q,gamma=max(routes,key=lambda x:x[2]);need(len(support_map)==353,'closed support count')
    closure_rows=[]
    for Q in sorted(support_map):
        v=max(support_map[Q]);closure_rows.append(dict(radical=Q,common_delta_upper=str(v),decimal_upper=out_decimal(v),all_positive_valuations=True,all_lengths=True,strict_gap_to_2_over_85=str(DELTA-v)))
    write_csv(out/'closure_by_support.csv',closure_rows)
    cl=dict(status='PASS',selected_supports_closed=353,all_selected_complements_checked=True,all_positive_valuations=True,all_lengths=True,common_delta_upper=str(gamma),common_delta_upper_decimal=out_decimal(gamma),strict_margin_to_2_over_85=str(DELTA-gamma),worst_majorant_route=label,worst_route_radical=worst_Q,remaining_three_prime_supports_after_G53_and_F6p_import=89,G53_large_prime_result_replayed=False,F6p_replayed=False,entire_three_prime_sector_closed=False,global_Cstar_identified=False)
    write_json(out/'closure.json',cl);print(json.dumps(cl,indent=2),flush=True)
    return cl

if __name__=='__main__':
    import sys
    certify(sys.argv[1],sys.argv[2])
