"""Finite upper witnesses, logarithm consumer, exact family assembly."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,prod,comb
import csv,json,hashlib
from exact import *
from family import *
from inventory import write_json,write_csv

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while b:=f.read(1<<20):h.update(b)
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
                a=s*d%M;need(min(a,M-a)>=H,'candidate independent minimality');checked+=1
        hs.append(H);P*=H
    need(P**170>210**18690,'candidate positive separation')
    return {'status':'PASS','M':210,'heights':hs,'product':str(P),'minimality_denominators':checked,'integer_gate':'P^170 > 210^18690','c0_lower':'2/85'}

def finite(work,out):
    work=Path(work);out=Path(out)
    domain=list(csv.DictReader((out/'finite_domain.csv').open()));rows=[];logs=set();mass=reps=bins=nonnegative=0
    for idx,d in enumerate(domain):
        M=int(d['M']);Q=int(d['radical']);hist=[]
        for line in (work/f'M{M}.hist').read_text().splitlines():
            h,n=map(int,line.split());hist.append((h,n));logs.add(h)
        need(hist and sum(n for h,n in hist)==M,'complete consumed histogram')
        need(all(1<=h<=M//2 and n>0 for h,n in hist),'histogram range')
        need(all(hist[j][0]<hist[j+1][0] for j in range(len(hist)-1)),'strict histogram ordering')
        total=sum(n*ticks(h)[1] for h,n in hist);lm=ticks(M)[0];logs.add(M)
        mu=F(total,M*LS)-F(lm,2*LS);target=DELTA*log_lo(Q)
        need(mu<target,'finite support mean gate M='+str(M))
        if mu>=0:nonnegative+=1
        rows.append({**{k:int(v) for k,v in d.items()},'representatives':M//2+1,'mass':M,'bins':len(hist),'witness_sha256':sha(work/f'M{M}.bin'),'histogram_sha256':sha(work/f'M{M}.hist'),'sum_log_upper_ticks':total,'log_M_lower_ticks':lm,'mu_upper':str(mu),'mu_upper_decimal':out_decimal(mu),'target_lower':str(target)})
        mass+=M;reps+=M//2+1;bins+=len(hist)
        if idx%1000==0:print('CERTIFIED MODULI',idx,'/',len(domain),flush=True)
    write_csv(out/'finite_certificate.csv',rows)
    worst=max(rows,key=lambda x:F(x['mu_upper']))
    result={'status':'PASS','finite_moduli':len(rows),'all_state_mass':mass,'representatives':reps,'histogram_bins':bins,'all_means_negative':nonnegative==0,'nonnegative_witness_means':nonnegative,'worst_upper':worst,'distinct_log_inputs':len(logs),'minimality_claimed':False}
    write_json(out/'finite_summary.json',result)
    for index,x in enumerate(sorted(logs)):
        a,b=alternate_log(x);lo,hi=ticks(x)
        need(F(lo,LS)<=a<=b<=F(hi,LS),'independent logarithm bounds '+str(x))
        if index%5000==0:print('INDEPENDENT LOGARITHMS',index,'/',len(logs),flush=True)
    write_json(out/'logarithm_checks.json',{'status':'PASS','all_census_log_inputs':len(logs),'primary':'exact atanh series 28 terms','independent':'exact negative logarithm series 96 terms','scale':LS})
    return result

def closure(out):
    out=Path(out);records=[]
    constants=json.loads((out/'constants.json').read_text());g=F(constants['generic_upper_at_cutoff'])
    records.append(('radical>=1000000',max(F(0),g/log_lo(RAD_CUTOFF))))
    for d in csv.DictReader((out/'pair_profiles.csv').open()):
        r,s,p=[int(d[k]) for k in ('r','s','p_cut')];Q=r*s*p
        records.append((f'pair {r},{s}, p>={p}',max(F(0),bound(Q,r,s,p,False)/log_lo(Q))))
    for d in csv.DictReader((out/'support_exceptions.csv').open()):
        r,s,p,M,Q=[int(d[k]) for k in ('r','s','p','M_cut','radical')]
        records.append((f'support {r},{s},{p}, M>={M}',max(F(0),bound(M,r,s,p,False)/log_lo(Q))))
    for d in csv.DictReader((out/'analytic_moduli.csv').open()):
        records.append((f"analytic modulus {d['M']}",max(F(0),F(d['mu_upper'])/log_lo(int(d['radical'])))))
    for d in csv.DictReader((out/'finite_certificate.csv').open()):
        records.append((f"finite modulus {d['M']}",max(F(0),F(d['mu_upper'])/log_lo(int(d['radical'])))))
    need(all(v<DELTA for tag,v in records),'strict common family ceiling')
    tag,gamma=max(records,key=lambda x:x[1]);need(gamma<DELTA,'finite-max common strict gap')
    primes=sieve(P_MIN-1);need(all(prime_trial(p) for p in primes),'remaining small primes')
    residual=[];tail=[]
    for i,r in enumerate(primes):
        for j in range(i+1,len(primes)):
            s=primes[j]
            if (r,s)==(2,3):continue
            for p in primes[j+1:]:
                Q=r*s*p;row={'r':r,'s':s,'p':p,'radical':Q,'status':'CLO; no complete valuation closure claimed'};residual.append(row)
                d=parameters(r,s);u=F(d['u'])
                leading=F(d['T0_tick'],CS)+phi_hi(5*u/p)+u*(1+ELL_HI)/p
                target=DELTA*log_lo(Q)
                if leading<target:
                    M=Q
                    while bound(M,r,s,p,False)>=target:M*=2
                    tail.append({'r':r,'s':s,'p':p,'radical':Q,'leading_upper':str(leading),'target_lower':str(target),'tail_from_M':M,'tail_margin':str(target-bound(M,r,s,p,False)),'finite_complement_checked':False})
                else:
                    tail.append({'r':r,'s':s,'p':p,'radical':Q,'leading_upper':str(leading),'target_lower':str(target),'tail_from_M':'not certified','tail_margin':'not certified','finite_complement_checked':False})
    need(len(primes)==15 and len(residual)==comb(15,3)-13==442,'complete residual-support count')
    need(all(not prime_trial(n) for n in range(48,53)),'no prime between 47 and 53')
    write_csv(out/'remaining_442_supports.csv',residual);write_csv(out/'remaining_support_tail_frontier.csv',tail)
    write_json(out/'candidate.json',candidate())
    obj={'status':'PASS','theorem':'all rad(M)=r*s*p, r<s<p, (r,s)!=(2,3), p>=53: mu(M)<(2/85)*log(r*s*p)','largest_prime_threshold':P_MIN,'all_lengths_and_positive_valuations':True,'all_cases_complete':True,'common_delta_upper':str(gamma),'common_delta_upper_decimal':out_decimal(gamma),'strict_margin_to_2_over_85':str(DELTA-gamma),'worst_majorant_route':tag,'candidate_lower':'2/85','remaining_three_prime_supports_with_internal_F6p_import':len(residual),'remaining_supports_analytic_tails':sum(isinstance(x['tail_from_M'],int) for x in tail),'remaining_supports_without_tail_certificate':sum(not isinstance(x['tail_from_M'],int) for x in tail),'remaining_support_finite_complements_checked':False,'entire_three_prime_sector_closed':False,'global_Cstar_identified':False,'F6p_import_replayed':False}
    write_json(out/'family_closure.json',obj);return obj

def run(work,out):
    result=finite(work,out);cl=closure(out);return {'finite':result,'closure':cl}
if __name__=='__main__':
    import sys
    print(json.dumps(run(sys.argv[1],sys.argv[2]),indent=2),flush=True)
