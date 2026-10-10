"""Independent consumer checks; exact declared domains, no sampling promoted to exhaustion."""
from pathlib import Path
from struct import iter_unpack
from fractions import Fraction as F
from math import gcd
import arithmetic as ar
import reference_family as ref
import check_routes
S=10**9;T=16;J=64

def witnesses(root,prefix):
    ss=ar.shape_rows(root);marks=[ref.rep(l) for l in ref.LABELS];n=len(marks)
    raw=Path(str(prefix)+'.witnesses').read_bytes();ar.need(len(raw)==len(ss)*n*4,'complete coefficient file')
    largest=0
    for row,(i,j) in enumerate(iter_unpack('<hh',raw)):
        _,_,a,z,b,d=ss[row//n];r=marks[row%n]
        ar.need(gcd(i*r[0]+j*r[1],ref.Q)==1,'full mark at every initial witness')
        ar.need(i*a+j*b!=0 or j*d!=0,'nonzero initial witness')
        largest=max(largest,abs(i),abs(j))
    return {'all_initial_coefficients':len(ss)*n,'largest_absolute_coefficient':largest}

def common_costs(root,prefix,cells):
    shapes=ar.shape_rows(root);ids=[];total=2688*33*J
    for line in Path(str(prefix)+'.cellcosts').read_text().splitlines():
        c,ii,j,cost=map(int,line.split());ids.append((c*33+ii)*J+j)
        _,_,a,z,b,d=shapes[c];_,_,lo,hi,tl,th,K=cells[c]
        u=(T,j-T) if j<2*T else (3*T-j,T);n=ii-16
        best=None
        for x in (lo,hi):
            for tt in (tl,th):
                A=40*d*x;B=(a*tt-40*b)*x;D=40*a*ar.PDEN
                ua=u[0]*A+u[1]*B;ub=u[1]*D
                R=max(abs(16*ua+n*ub),abs(-n*ua+16*ub));num=R*R;den=T*T*A*D*(256+n*n)
                v=F(ar.logs(num)[1]-ar.logs(den)[0],2*ar.LS) if num!=den else F(0)
                best=v if best is None else max(best,v)
        ar.need(best<=F(cost,S),'independent cell logarithmic upper')
    exp=list(range(0,total,997))
    if exp[-1]!=total-1:exp.append(total-1)
    ar.need(ids==exp,'declared common-cost index set');ar.logs.cache_clear()
    return {'all_common_cell_costs_regenerated_by_each_main_engine':total,'independent_common_costs':len(ids),'corners_per_checked_cost':4,'stride':997}

def run_one(root,p,prefix,four=True):
    ref.setup(p,four)
    out={'p':p,'q':ref.Q,'algebra':ref.algebra_checks()}
    out.update(witnesses(root,prefix));out.update(check_routes.check(root,prefix))
    return out
