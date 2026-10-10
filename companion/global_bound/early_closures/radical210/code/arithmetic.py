"""Exact arithmetic, no assertions and no floating-point acceptance gates."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict
from math import gcd,prod
from struct import unpack_from
import csv,json
S=10**9; LS=10**15; TERMS=20; PDEN=4503599627370496
Q=210; T=24; J=96; AG=16; MASS=1<<28

def need(ok:bool,message:str)->None:
    if not ok:raise ArithmeticError(message)
def ceildiv(a:int,b:int)->int:
    need(b>0,'positive denominator');return -((-a)//b)
def atom_log(a:int,b:int):
    need(a>=0 and b>0 and 3*a<=b,'atanh range')
    if not a:return 0,0
    ap=a;bp=b;lo=0
    for j in range(TERMS):
        lo+=(2*LS*ap)//((2*j+1)*bp);ap*=a*a;bp*=b*b
    rn=2*LS*ap*b*b;rd=bp*(2*TERMS+1)*(b*b-a*a)
    return lo,lo+TERMS+ceildiv(rn,rd)
L2,U2=atom_log(1,3)
@lru_cache(None)
def logs(n:int):
    need(n>=1,'positive logarithm input');k=n.bit_length()-1;p=1<<k
    lo,hi=atom_log(n-p,n+p);return lo+k*L2,hi+k*U2

def shape_rows(root):
    rows=[tuple(map(int,l.split())) for l in (root/'data/shape_centres.tsv').read_text().splitlines()]
    need(len(rows)==2688,'128 by 21 centres')
    for i,(k,t,a,z,b,d) in enumerate(rows):
        need((k,t)==(-125+i//21,-10+i%21),'centre ordering')
        need(a>0 and z==0 and d>0,'triangular positive centre')
    return rows

def make_cells(root):
    vals=list(map(int,(root/'data/partition.txt').read_text().split()))
    den=vals[0];rs=vals[1:165]
    need(den==PDEN and len(vals)==186,'partition format')
    need(all(rs[j]<rs[j+1] for j in range(163)),'ordered radial bounds')
    # Gauss gives rho <= 2/sqrt(3); the final upper endpoint covers it.
    need(3*rs[162]**2>=4*den**2,'top Gauss domain covered')
    out=[];corners=0
    for k,t,a,z,b,d in shape_rows(root):
        lo=rs[k+159];hi=rs[k+160];tl=max(-20,2*t-1);th=min(20,2*t+1);worst=S
        for x in (lo,hi):
            for tt in (tl,th):
                A=40*d*x;B=(a*tt-40*b)*x;D=40*a*den;det=A*D;tr=A*A+B*B+D*D
                def ok(K):return K>=S and (K*K+S*S)*det*det>=K*S*(tr*tr-2*det*det)
                left=S;right=4*S;need(ok(right),'condition bounded by 4')
                while left<right:
                    mid=(left+right)//2
                    if ok(mid):right=mid
                    else:left=mid+1
                need(ok(left),'condition polynomial');worst=max(worst,left);corners+=1
        out.append((k,t,lo,hi,tl,th,worst))
    return out,{'cells':len(out),'condition_corners':corners,'largest_condition_upper':str(F(max(r[-1] for r in out),S))}

def encode_cells(rows):return ''.join(' '.join(map(str,r))+'\n' for r in rows)

def uniform_certificate(prefix:Path):
    rho0=None;results=[]
    for k,cap in [(1,F(7,25)),(2,F(27,200)),(3,F(61,500))]:
        cellrows=[list(map(int,l.split())) for l in Path(str(prefix)+f'.cells{k}').read_text().splitlines()]
        need(len(cellrows)==2688,'all cells retained')
        for c,row in enumerate(cellrows):
            need(row[0]==c and 0<=row[2]<576 and 0<=row[3]<33,'cell label and angle witness')
        rr=[list(map(int,l.split())) for l in Path(str(prefix)+f'.radial{k}').read_text().splitlines()]
        need(len(rr)==128,'all radial bins')
        for j,r in enumerate(rr):
            need(r[0]==j-125,'radial ordering')
            group=cellrows[j*21:(j+1)*21]
            need(r[3]==max(c[1] for c in group),'radial majorant exact maximum')
            need(r[4] in [c[0] for c in group if c[1]==r[3]],'radial witness belongs to maximum')
            if j:need(rr[j-1][2]==r[1],'contiguous radial coverage')
        b=[F(r[3],S) for r in rr]
        for j in range(126,-1,-1):b[j]=max(b[j],b[j+1])
        A=F(29,14);rho=F(rr[0][1],PDEN)
        need(rho>F(1,2**20),'cusp cutoff lower bound')
        need(U2<LS*F(7,10),'log2 upper bound')
        # log(209)<6, log(210)<6: e>5/2 and (5/2)^6>210.
        need(F(5,2)**6>210,'elementary exponential comparison')
        t0=13+3*k+4*rho
        bulk=b[-1]+sum((b[j]-b[j+1])*min(F(1),A*F(rr[j][2],PDEN)) for j in range(127))
        cusp=A*rho*(max(t0-b[0],F(0))+F(1,2))
        total=bulk+cusp;need(total<cap,f'uniform cap k={k}')
        results.append({'k':k,'bulk_upper':str(bulk),'cusp_upper':str(cusp),'total_upper':str(total),'rational_cap':str(cap),'strict_comparison':True,'isotonic_radial_numerators':[int(v*S) for v in b]})
        rho0=rho
    # Polynomial division proving pi<22/7 by a positive integral.
    quotient=[4,0,-4,0,5,-4,1];product=[0]*9
    for j,c in enumerate(quotient):product[j]+=c;product[j+2]+=c
    product[0]-=4
    need(product==[0,0,0,0,1,-4,6,-4,1],'pi polynomial identity')
    need(sum(F(c,j+1) for j,c in enumerate(quotient))==F(22,7),'integrated polynomial')
    return {'status':'PASS','claim_level':'CL3 finite constants for the general theorem proved in the dossier','Q':Q,'directions':J,'angle_grid_points':33,'shape_cells':2688,'marks_per_cell':576,'envelope_rows_per_horizon':2688*576,'cusp_start':str(rho0),'short_vector_CDF_coefficient':'29/14','short_vector_source_bound_all_L_ge_1':True,'final_modulus_requires_radical_210':True,'rows':results,'uses_historical_tail35':False,'uses_equdistribution_or_Weil':False}

def moduli():
    ans=[]
    for m in range(Q,150000,Q):
        n=m
        for p in [2,3,5,7]:
            while n%p==0:n//=p
        if n==1:ans.append(m)
    # A second enumeration by exponent vectors, bounded using integer products.
    other=[]
    a=2
    while a*105<150000:
        b=3
        while a*b*35<150000:
            c=5
            while a*b*c*7<150000:
                d=7
                while a*b*c*d<150000:other.append(a*b*c*d);d*=7
                c*=5
            b*=3
        a*=2
    need(ans==sorted(other) and len(set(other))==len(other),'two complete domain enumerations')
    return ans

def finite_certificate(histfile:Path,wdir:Path):
    hist=defaultdict(dict)
    for r in csv.DictReader(histfile.open()):
        m,h,n=int(r['M']),int(r['height']),int(r['count'])
        need(h not in hist[m] and n>0,'distinct histogram height and positive mass');hist[m][h]=n
    need(sorted(hist)==moduli(),'complete finite domain')
    rows=[]
    for m,hh in sorted(hist.items()):
        need(sum(hh.values())==m,'histogram mass')
        lo=sum(n*logs(h)[0] for h,n in hh.items());hi=sum(n*logs(h)[1] for h,n in hh.items());lm,um=logs(m)
        mu=F(2*hi-m*lm,2*m*LS);delta=F(2*hi-m*lm,2*m*lm)
        gate=170*hi<89*m*lm
        if m!=Q:need(gate and delta<F(19,1000),f'finite competitor M={m}')
        rows.append({'M':m,'sum_log_height_lower_num':lo,'sum_log_height_upper_num':hi,'scale':LS,'mu_upper':str(mu),'Delta_upper':str(delta),'below_2_over_85':gate})
    binary=(wdir/'w210.bin').read_bytes();need(len(binary)==106*8,'candidate bytes')
    heights=[]
    for s in range(106):
        a,d=unpack_from('<ii',binary,s*8);need(d>0 and gcd(d,Q)==1 and (a-s*d)%Q==0,'candidate witness')
        heights.append(max(abs(a),d))
    full=[heights[min(s,Q-s)] for s in range(Q)];checks=0
    for s,h in enumerate(full):
        for d in range(1,h):
            if gcd(d,Q)!=1:continue
            checks+=1;a=s*d%Q;need(max(d,min(a,Q-a))>=h,'candidate smaller denominator excluded')
    p=prod(full);need(p**170>Q**18690,'candidate exact lower comparison')
    factors={};v=p
    for prime in range(2,106):
        while v%prime==0:factors[str(prime)]=factors.get(str(prime),0)+1;v//=prime
    need(v==1,'candidate complete factorization')
    n2=next(r for r in rows if r['M']==Q*Q)
    need(F(n2['mu_upper'])<F(32,255),'N2 below candidate mean')
    need(10*logs(150000)[0]>119*LS,'large-base logarithm threshold')
    need(logs(420)[0]>6*LS,'first other base logarithm')
    need(3*logs(210)[0]>16*LS,'candidate base logarithm')
    return {'status':'PASS','claim_level':'CL3','distinct_moduli':len(rows),'full_state_mass':sum(hist),'stored_witnesses':sum(m//2+1 for m in hist),'N1_excluded_count':len(rows)-1,'histograms_are_admissible_witness_heights_not_claimed_minimal':True,'candidate_all_states':Q,'candidate_smaller_unit_denominator_tests':checks,'candidate_P':str(p),'candidate_factorization':factors,'candidate_P170_gt_21018690':True,'N2_modulus_reused_from_N1_census':44100,'N2_mu_upper':n2['mu_upper'],'log_scale':LS,'log_terms':TERMS,'rows':rows}

def conclusion(uniform,finite):
    need([r['rational_cap'] for r in uniform['rows']]==['7/25','27/200','61/500'],'uniform interfaces')
    need(finite['N1_excluded_count']==120 and finite['candidate_P170_gt_21018690'],'finite interfaces')
    need(F(7,25)/F(119,10)==F(2,85),'large-base comparison')
    need(F(27,200)/6==F(9,400)<F(2,85),'other-base N>=2 comparison')
    need(F(61,500)/F(16,3)==F(183,8000)<F(2,85),'base210 N>=3 comparison')
    need(F(finite['N2_mu_upper'])<F(32,255),'N2 interface')
    return {'status':'PASS','claim_level':'CL5 with finite CL3 dependencies','closed_sector':'all integer B>=2 with rad(B)=210, all N>=1','unique_maximizing_pair':[210,1],'sector_maximum_exact':'log(P^2/210^210)/(420*log(210))','all_other_pairs_strictly_below_2_over_85':True,'global_Cstar_value_determined':False,'other_radicals_excluded_here':False,'historical_base210_or_tail35_theorem_used':False,'formalized':False,'independent_specialist_review':False}
