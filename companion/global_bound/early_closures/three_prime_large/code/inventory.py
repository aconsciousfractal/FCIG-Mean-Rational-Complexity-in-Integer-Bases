"""Complete support reduction. Every decision is a rational comparison."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,prod
from bisect import bisect_left
import csv,json,time
from exact import *
from family import *

def write_json(path,obj):Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def write_csv(path,rows):
    need(bool(rows),'nonempty table')
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def factors(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            a=0
            while n%p==0:n//=p;a+=1
            out.append((p,a))
        p+=1
    if n>1:out.append((n,1))
    return out

def build(out,work):
    out=Path(out);work=Path(work);out.mkdir(parents=True,exist_ok=True);work.mkdir(parents=True,exist_ok=True)
    write_json(out/'constants.json',constants())
    primes=sieve(RAD_CUTOFF//10)
    need(all(prime_trial(p) for p in primes),'all generated primes trial-checked')
    pairs=[]
    for i,r in enumerate(primes):
        if r**3>=RAD_CUTOFF:break
        for s in primes[i+1:]:
            if r*s*max(P_MIN,s+1)>=RAD_CUTOFF:break
            if (r,s)!=(2,3):pairs.append((r,s))
    pairs2=[]
    for q in range(10,(RAD_CUTOFF-1)//P_MIN+1):
        fac=factors(q)
        if len(fac)==2 and all(e==1 for p,e in fac):
            r,s=[p for p,e in fac]
            if (r,s)!=(2,3) and q*max(P_MIN,s+1)<RAD_CUTOFF:pairs2.append((r,s))
    need(sorted(pairs)==sorted(pairs2),'complete ordered pair inventory')
    pairrows=[];rootrows=[];analytic=[];pending=[];root_keys=set();t=time.time();fullsupportcount=0;empty_pairs=0
    for index,(r,s) in enumerate(pairs):
        q=r*s;start=max(P_MIN,s+1)
        def gate(p):return DELTA*log_lo(q*p)-bound(q*p,r,s,p,False)
        lo=start;hi=start
        while gate(hi)<=0:hi*=2
        while lo<hi:
            mid=(lo+hi)//2
            if gate(mid)>0:hi=mid
            else:lo=mid+1
        cut=hi;mar=gate(cut);need(mar>0,'certified pair prime cutoff')
        row=dict(parameters(r,s));row.update({'p_start':start,'p_cut':cut,'cutoff_margin':str(mar)});pairrows.append(row)
        ap=bisect_left(primes,start);bp=bisect_left(primes,(RAD_CUTOFF-1)//q+1)
        fullsupportcount+=max(0,bp-ap)
        if bp<=ap:empty_pairs+=1
        for p in primes[ap:bisect_left(primes,min(cut,(RAD_CUTOFF-1)//q+1))]:
            Q=q*p;target=DELTA*log_lo(Q)
            mc=Q
            while bound(mc,r,s,p,False)>=target:mc*=2
            need((r,s,p) not in root_keys,'unique support');root_keys.add((r,s,p))
            rootrows.append({'r':r,'s':s,'p':p,'radical':Q,'M_cut':mc,'cutoff_margin':str(target-bound(mc,r,s,p,False))})
            a=r
            while a*s*p<mc:
                b=s
                while a*b*p<mc:
                    d=p
                    while a*b*d<mc:
                        M=a*b*d;v=bound(M,r,s,p)
                        rec={'M':M,'r':r,'s':s,'p':p,'radical':Q}
                        if v<target:analytic.append({**rec,'mu_upper':str(v),'target_lower':str(target)})
                        else:pending.append(rec)
                        d*=p
                    b*=s
                a*=r
        if index%100==0:print('PAIR',index,'/',len(pairs),'root exceptions',len(rootrows),'pending',len(pending),'seconds',round(time.time()-t,2),flush=True)
    pending.sort(key=lambda x:x['M']);analytic.sort(key=lambda x:x['M']);rootrows.sort(key=lambda x:x['radical'])
    need(len(set(x['M'] for x in pending+analytic))==len(pending+analytic),'unique moduli')
    maxmc=max(x['M_cut'] for x in rootrows)
    rootmap={(x['r'],x['s'],x['p']):x for x in rootrows}
    # Independent full integer scan, using trial division rather than power tuples.
    pending2=[];an2=[]
    for M in range(2,maxmc):
        f=factors(M)
        if len(f)!=3:continue
        rr,ss,pp=[p for p,e in f]
        key=(rr,ss,pp)
        if key not in rootmap:continue
        row=rootmap[key]
        if M>=row['M_cut']:continue
        if bound(M,rr,ss,pp)<DELTA*log_lo(row['radical']):an2.append(M)
        else:pending2.append(M)
    need(pending2==[x['M'] for x in pending],'full-scan finite pending domain')
    need(an2==[x['M'] for x in analytic],'full-scan analytic finite domain')
    write_csv(out/'pair_profiles.csv',pairrows);write_csv(out/'support_exceptions.csv',rootrows)
    write_csv(out/'analytic_moduli.csv',analytic);write_csv(out/'finite_domain.csv',pending)
    (work/'domain.txt').write_text(''.join(f"{x['M']} {x['radical']}\n" for x in pending))
    cert={'status':'PASS','largest_prime_threshold':P_MIN,'radical_cutoff':RAD_CUTOFF,'pair_families':len(pairrows),'empty_pair_domains':empty_pairs,'ordered_supports_below_radical_cutoff':fullsupportcount,'support_exceptions':len(rootrows),'analytic_valuation_moduli':len(analytic),'finite_moduli':len(pending),'all_state_mass':sum(x['M'] for x in pending),'representatives':sum(x['M']//2+1 for x in pending),'max_finite_M':max(x['M'] for x in pending),'max_M_cut':maxmc,'max_finite_p':max(x['p'] for x in pending),'finite_pair_count':len(set((x['r'],x['s']) for x in pending)),'primes_generated':len(primes),'primes_trial_verified':True,'pair_inventory_two_routes':True,'finite_domain_complete_integer_scan':True,'floating_point_decisions':False,'global_Cstar_identified':False}
    write_json(out/'inventory.json',cert);return cert
if __name__=='__main__':
    import sys
    print(json.dumps(build(sys.argv[1],sys.argv[2]),indent=2),flush=True)
