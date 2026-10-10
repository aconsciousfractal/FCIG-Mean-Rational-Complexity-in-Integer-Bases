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
Q=6; T=24; J=96; AG=16; MASS=1<<28

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
    for k,cap in [(1,F(393,1000)),(2,F(187,1000)),(3,F(17,250)),(4,-F(3,500)),(5,-F(19,500))]:
        cellrows=[list(map(int,l.split())) for l in Path(str(prefix)+f'.cells{k}').read_text().splitlines()]
        need(len(cellrows)==2688,'all cells retained')
        for c,row in enumerate(cellrows):
            need(row[0]==c and 0<=row[2]<12 and 0<=row[3]<33,'cell label and angle witness')
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
        # Conservative cusp constants: log(5)<6 and log(6)<6; e>5/2.
        need(F(5,2)**6>6,'elementary exponential comparison')
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
    return {'status':'PASS','claim_level':'CL3 finite constants for the general theorem proved in the dossier','Q':Q,'directions':J,'angle_grid_points':33,'shape_cells':2688,'marks_per_cell':12,'envelope_rows_per_horizon':2688*12,'cusp_start':str(rho0),'short_vector_CDF_coefficient':'29/14','short_vector_source_bound_all_L_ge_1':True,'observable_is_relaxed_mod6':True,'final_modulus_radical_unrestricted':True,'hypothesis': '6^k divides M','rows':results,'uses_historical_tail35':False,'uses_equdistribution_or_Weil':False}
