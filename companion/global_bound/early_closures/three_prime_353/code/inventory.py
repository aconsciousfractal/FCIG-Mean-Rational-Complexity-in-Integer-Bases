"""Exact complete complement: three positive valuations and primitive-base normalization."""
import csv,json,heapq
from pathlib import Path
from itertools import combinations
from math import prod,gcd
from central import *

def write_json(path,obj):Path(path).write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
def write_csv(path,rows):
    need(bool(rows),'nonempty table')
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def power_tuples(primes,cap):
    vals=[(1,())]
    for pr in primes:
        nxt=[]
        for M,es in vals:
            M*=pr;e=1
            while M<cap:nxt.append((M,es+(e,)));M*=pr;e+=1
        vals=nxt
    return sorted(vals)

def heap_moduli(primes,cap):
    q=prod(primes)
    if q>=cap:return []
    pending=[q];seen={q};out=[]
    while pending:
        M=heapq.heappop(pending);out.append(M)
        for p in primes:
            y=M*p
            if y<cap and y not in seen:seen.add(y);heapq.heappush(pending,y)
    return out

def primitive_base(M,primes):
    es=prime_exponents(M,primes);g=gcd(*es)
    B=prod(p**(e//g) for p,e in zip(primes,es))
    need(B**g==M and gcd(*(e//g for e in es))==1,'primitive base identity')
    return B,g,es

def build(out,work,root):
    out=Path(out);work=Path(work);root=Path(root)
    out.mkdir(parents=True,exist_ok=True);work.mkdir(parents=True,exist_ok=True)
    primes=sieve(47);need(len(primes)==15 and all(prime_trial(p) for p in primes),'prime domain')
    supports=[t for t in combinations(primes,3) if t[:2]!=(2,3)]
    alternate=[]
    for Q in range(2,prod(primes[-3:])+1):
        z=Q;ff=[]
        for p in primes:
            if z%p==0:
                ff.append(p);z//=p
        if z==1 and len(ff)==3 and ff[:2]!=[2,3]:alternate.append(tuple(ff))
    need(sorted(alternate)==sorted(supports),'independent squarefree radical scan')
    need(len(supports)==442,'complete input support set')
    old=list(csv.DictReader((root/'imports/frontier_442.csv').open()))
    need({tuple(int(x[k]) for k in ('r','s','p')) for x in old}==set(supports),'imported frontier identity')
    old_tail={int(x['radical']) for x in old if x['tail_from_M']!='not certified'}
    profiles=[];selected=[];unresolved=[];pending=[];analytic=[];cost=[]
    for prs in supports:
        r,s,p=prs;Q=prod(prs);target=DELTA*log_lo(Q);C=coeffs(*prs);lc=central_constants(*prs)
        f=dict(r=r,s=s,p=p,radical=Q,full_unit_density=str(lc['u']),full_primitive_density_factor=str(lc['rho']),full_Dplus=str(lc['dp']),full_Dminus=str(lc['dminus']),full_interval_span=str(lc['ki']),full_floor_K=str(lc['K']),full_harmonic_h=str(lc['h']),leading_upper=str(C[0]),C0=str(C[1]),C1=str(C[2]),C2=str(C[3]),C3=str(C[4]),target_lower=str(target),tail_cutoff='none',tail_mu_upper='none',old_tail_available=Q in old_tail)
        if C[0]>=target:
            unresolved.append(dict(r=r,s=s,p=p,radical=Q,leading_upper=str(C[0]),target_lower=str(target),status='CLO; this full-central majorant has no certified tail'))
            profiles.append(f);continue
        K=Q
        # Scalar profile at K; K need not have the target radical.
        while improved(K,*prs)>=target:K*=2
        tail_value=improved(K,*prs)
        need(tail_value<target,'strict tail certificate')
        f.update(tail_cutoff=K,tail_mu_upper=str(tail_value));profiles.append(f)
        vals=power_tuples(prs,K);heap=heap_moduli(prs,K)
        need([M for M,es in vals]==heap,'independent complete valuation enumeration')
        pc=ac=pm=0;radonly_mass=0
        for M,es in vals:
            B,g,es2=primitive_base(M,prs);need(es==es2,'independent exponent extraction')
            upper=true_upper(M,*prs);rec=dict(M=M,r=r,s=s,p=p,radical=Q,primitive_base=B,max_length=g,alpha=es[0],beta=es[1],gamma=es[2])
            if upper>=target:radonly_mass+=M
            if upper<DELTA*log_lo(B):
                analytic.append({**rec,'mu_upper':str(upper),'target_lower':str(DELTA*log_lo(B))});ac+=1
            else:pending.append(rec);pc+=1;pm+=M
        selected.append(dict(r=r,s=s,p=p,radical=Q,tail_cutoff=K,tail_mu_upper=str(tail_value),target_lower=str(target),complement_moduli=len(vals),analytic_moduli=ac,witness_moduli=pc,witness_state_mass=pm,new_tail_not_in_parent=Q not in old_tail))
        cost.append(dict(radical=Q,full_central_radical_only_mass=radonly_mass,primitive_base_witness_mass=pm,complement_moduli=len(vals),witness_moduli=pc))
    pending.sort(key=lambda x:x['M']);analytic.sort(key=lambda x:x['M'])
    need(len({x['M'] for x in pending+analytic})==len(pending)+len(analytic),'global unique moduli')
    need(len(selected)==353 and len(unresolved)==89,'selected and remaining counts')
    need(old_tail<={x['radical'] for x in selected},'all 340 earlier tail supports included')
    need(max(x['M'] for x in pending)<200000000,'safe census integer domain')
    write_csv(out/'profiles_442.csv',profiles);write_csv(out/'selected_353.csv',selected)
    write_csv(out/'remaining_89.csv',unresolved);write_csv(out/'finite_domain.csv',pending)
    write_csv(out/'analytic_moduli.csv',analytic);write_csv(out/'cost_comparison.csv',cost)
    (work/'domain.txt').write_text(''.join(f"{x['M']} {x['radical']}\n" for x in pending))
    cert=dict(status='PASS',input_supports=len(supports),selected_supports=len(selected),old_tail_supports=len(old_tail),new_tail_supports=len(selected)-len(old_tail),remaining_supports=len(unresolved),tail_cutoff_max=max(x['tail_cutoff'] for x in selected),complement_moduli=len(pending)+len(analytic),analytic_moduli=len(analytic),witness_moduli=len(pending),witness_state_mass=sum(x['M'] for x in pending),representatives=sum(x['M']//2+1 for x in pending),max_witness_modulus=max(x['M'] for x in pending),independent_radical_scan=True,valuation_enumeration_two_routes=True,primitive_base_normalization=True,all_selected_complements_checked=False,global_Cstar_identified=False)
    write_json(out/'inventory.json',cert);print(json.dumps(cert,indent=2),flush=True)
    return cert

if __name__=='__main__':
    import sys
    build(sys.argv[1],sys.argv[2],Path(__file__).resolve().parents[1])
