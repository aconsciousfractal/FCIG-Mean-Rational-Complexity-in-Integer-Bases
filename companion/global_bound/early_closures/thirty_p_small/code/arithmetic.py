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

