"""Exact certificate assembly from consumed witness histograms."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,isqrt,prod
import csv,json,hashlib
from exact import *

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while b:=f.read(1<<20):h.update(b)
    return h.hexdigest()

def candidate():
    M=210;heights=[];Hprod=1;den_checks=0
    for s in range(M):
        H=max(1,min(s,M-s));best=(min(s,M-s),1)
        for d in range(1,M//2+1):
            if gcd(d,M)!=1:continue
            a=s*d%M;a=min(a,M-a)
            h=max(a,d)
            if h<H:H=h;best=(a,d)
        for d in range(1,H):
            if gcd(d,M)==1:
                a=s*d%M
                need(min(a,M-a)>=H,'candidate lower height')
                den_checks+=1
        heights.append(H);Hprod*=H
    need(Hprod**170>210**18690,'candidate integer product gate')
    need(log_lo(210)>F(16,3),'candidate base logarithm')
    need(log_lo(762)>F(51,8),'log762 rational threshold')
    need(F(3,20)/log_lo(762)<F(2,85),'new-family normalized gate')
    return {'status':'PASS','M':210,'heights':heights,'product':str(Hprod),'minimality_denominator_checks':den_checks,'product_comparison':'P^170 > 210^18690','c0_lower':'2/85','new_family_delta_upper':str(F(3,20)/log_lo(762))}

def constants():
    a=F(2,3);l=F(7,10);beta=F(1,5);r=F(5,7);r2=F(99,70);pi=F(22,7)
    need(F(7,5)**2<2<F(99,70)**2,'sqrt 2 enclosure')
    need(sum(F(1,n*n) for n in range(1,8))>F(3,2),'zeta lower')
    zeta_hi=sum(F(1,n*n) for n in range(1,6))+F(2,11)
    need(zeta_hi==F(65159,39600) and zeta_hi<F(5,3),'zeta upper')
    need(log_hi(2)<l,'log 2 upper')
    alo=1/zeta_hi
    delta=alo*(F(11,16)+l/2)-F(1,3)-5*l/12
    need(delta==F(2885,521272) and delta>F(1,200),'leading margin')
    need(alo/2-F(5,12)<0,'log 2 coefficient sign')
    c0=a+2+8*beta/9+l*(a+F(16,9))+(5+5*a/8)*r+(4+2*l)/3+13*pi*r2/6
    c1=F(2,3)+l/3+8*r/9
    d0=16*l/9+2
    need(c0<20 and c1<F(31,20) and d0<F(13,4),'polynomial coefficients')
    phi127=phi_hi(127)
    need(phi127<F(71,1000),'prime cost scalar')
    need(log_hi(387)<6,'cutoff logarithm')
    cutoff_upper=F(71,1000)-F(1,200)+F(17,3810)+(20+F(31,20)*6)/387+(F(13,4)+F(4,9)*6)/387**2
    need(387**2<=150000,'cutoff square')
    need(cutoff_upper<F(3,20),'infinite-region gate')
    return {'status':'PASS','alpha_interval':[str(alo),'2/3'],'sqrt2_interval':['7/5','99/70'],'log2_upper':'7/10','delta_lower':str(delta),'linear_constant_upper':str(c0),'logarithmic_coefficient_upper':str(c1),'quadratic_constant_upper':str(d0),'coarsened_coefficients':['20','31/20','13/4','4/9'],'phi127_upper':str(phi127),'phi127_upper_decimal':out_decimal(phi127),'cutoff_modulus':150000,'cutoff_upper':str(cutoff_upper),'cutoff_upper_decimal':out_decimal(cutoff_upper),'cutoff_margin':str(F(3,20)-cutoff_upper),'old_Psi127_interval':[str(log_lo(F(133,127))+F(6,127)*log_lo(F(133,6))),str(log_hi(F(133,127))+F(6,127)*log_hi(F(133,6)))],'new_prime_loss_above_cutoff_upper':str(phi127+F(13,6)*pi*r2/387)}

def finite(work,result):
    work=Path(work);result=Path(result)
    domain=[tuple(map(int,l.split())) for l in (work/'domain.txt').read_text().splitlines()]
    rows=[];total_mass=0;total_reps=0;bins=0;log_values=set()
    for M,q in domain:
        hist=[]
        for line in (work/f'M{M}.hist').read_text().splitlines():
            h,n=map(int,line.split());hist.append((h,n));log_values.add(h)
        need(hist and all(1<=h<=M//2 and n>0 for h,n in hist),'histogram range')
        need(all(hist[i][0]<hist[i+1][0] for i in range(len(hist)-1)),'histogram ordering')
        need(sum(n for h,n in hist)==M,'histogram total mass')
        upper=sum(n*ticks(h)[1] for h,n in hist)
        lower_M=ticks(M)[0];log_values.add(M)
        mu=F(upper,M*LS)-F(lower_M,2*LS)
        need(mu<F(3,20),f'finite uniform mean: M={M}')
        need(mu<F(0),f'finite strengthened negativity: M={M}')
        rows.append({'M':M,'p':q//6,'representatives':M//2+1,'mass':M,'bins':len(hist),'witness_sha256':sha(work/f'M{M}.bin'),'histogram_sha256':sha(work/f'M{M}.hist'),'sum_upper_ticks':str(upper),'log_M_lower_ticks':str(lower_M),'mu_upper':str(mu),'mu_upper_decimal':out_decimal(mu)})
        total_mass+=M;total_reps+=M//2+1;bins+=len(hist)
    with (result/'finite_certificate.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    best=max(rows,key=lambda r:F(r['mu_upper']))
    cert={'status':'PASS','moduli':len(rows),'all_state_mass':total_mass,'representatives':total_reps,'total_histogram_bins':bins,'all_finite_means_negative':True,'worst_upper':best,'distinct_integer_log_inputs':len(log_values),'witness_minimality_claimed':False}
    (result/'finite_summary.json').write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n')
    # Independent logarithm series encloses every input actually used by the census.
    count=0
    for x in sorted(log_values):
        lo,hi=alternate_log(x);a,b=ticks(x)
        need(F(a,LS)<=lo<=hi<=F(b,LS),f'independent log ticks x={x}')
        count+=1
    audit={'status':'PASS','all_census_integer_log_inputs':count,'primary':'range-reduced atanh series, 28 terms','secondary':'negative logarithm Taylor series, 96 terms','tick_scale':LS}
    (result/'logarithm_checks.json').write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n')
    return cert

def run(work,result):
    result=Path(result)
    for name,obj in [('candidate',candidate()),('constants',constants())]:
        (result/(name+'.json')).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    return finite(work,result)

if __name__=='__main__':
    import sys
    print(json.dumps(run(sys.argv[1],sys.argv[2]),indent=2),flush=True)
