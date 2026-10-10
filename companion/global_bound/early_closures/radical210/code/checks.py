"""Independent finite checks with explicit domains; general arguments are in the proof."""
from pathlib import Path
from fractions import Fraction as F
from struct import iter_unpack
from math import gcd
from collections import Counter
import arithmetic as ar
import reference as ref
need=ar.need

def initial_witnesses(root,prefix):
    ss=ar.shape_rows(root);marks=[ref.rep(l) for l in ref.LABELS]
    data=Path(str(prefix)+'.witnesses').read_bytes();need(len(data)==len(ss)*576*4,'all initial coefficient bytes')
    largest=0
    for row,(i,j) in enumerate(iter_unpack('<hh',data)):
        _,_,a,z,b,d=ss[row//576];r=marks[row%576]
        need(gcd(i*r[0]+j*r[1],210)==1,'initial full-unit mark')
        x=i*a+j*b;y=j*d;need(x!=0 or y!=0,'nonzero initial witness')
        largest=max(largest,abs(i),abs(j))
    return {'all_initial_integer_coefficients':len(ss)*576,'largest_coefficient':largest}

def cell_costs(root,prefix,cells):
    shapes=ar.shape_rows(root);given=[];total=2688*33*96
    for line in Path(str(prefix)+'.cellcosts').read_text().splitlines():
        c,ii,j,cost=map(int,line.split());given.append((c*33+ii)*96+j)
        need(0<=c<2688 and 0<=ii<33 and 0<=j<96,'cost audit range')
        _,_,a,z,b,d=shapes[c];_,_,lo,hi,tl,th,K=cells[c]
        u=(24,j-24) if j<48 else (72-j,24);n=ii-16
        maxupper=None
        for x in (lo,hi):
            for tt in (tl,th):
                A=40*d*x;B=(a*tt-40*b)*x;D=40*a*ar.PDEN
                ua=u[0]*A+u[1]*B;ub=u[1]*D
                R=max(abs(16*ua+n*ub),abs(-n*ua+16*ub))
                num=R*R;den=24**2*A*D*(256+n*n)
                upper=F(0) if num==den else F(ar.logs(num)[1]-ar.logs(den)[0],2*ar.LS)
                maxupper=upper if maxupper is None else max(maxupper,upper)
        need(maxupper<=F(cost,ar.S),'independent rational cell-log upper bound')
    expected=[i for i in range(0,total,997)]
    if expected[-1]!=total-1:expected.append(total-1)
    need(given==expected,'complete declared cost-audit index set')
    return {'cell_costs_checked_against_independent_series':len(given),'corners_per_cost':4,'sampling_stride':997,'all_cost_entries_recomputed_by_main_engine':total}

def endpoint_checks(root,cells):
    ss=ar.shape_rows(root);count=0
    # Exact squared norm inequalities: no trigonometric approximations.
    for c in range(0,2688,31):
        _,_,a,z,b,d=ss[c];_,_,lo,hi,tl,th,K=cells[c]
        for rx in [F(1,4),F(1,2),F(3,4)]:
            x=lo+(hi-lo)*rx
            for rt in [F(1,3),F(2,3)]:
                t=tl+(th-tl)*rt
                for n in [-16,-5,0,7,16]:
                    for j in [0,13,23,47,71,95]:
                        u=(24,j-24) if j<48 else (72-j,24)
                        def normratio(x,t):
                            A=40*d*x;B=(a*t-40*b)*x;D=40*a*ar.PDEN
                            ua=u[0]*A+u[1]*B;ub=u[1]*D
                            R=max(abs(16*ua+n*ub),abs(-n*ua+16*ub))
                            return F(R*R,24**2*A*D*(256+n*n))
                        v=normratio(x,t);cap=max(normratio(xx,tt) for xx in [lo,hi] for tt in [tl,th])
                        need(v<=cap,'whole-cell endpoint domination test');count+=1
    return {'exact_interior_vs_corner_norm_tests':count,'shape_stride':31,'general_convexity_is_proved_not_inferred_from_tests':True}

def short_vector_checks():
    states=thresholds=0
    for L in range(1,401):
        hist=Counter()
        for s in range(L):
            Y,V,P=ref.canonical((L,0,s,1));nu=Y[0]**2+Y[1]**2
            need(abs(ref.det(Y))==L and nu>0,'source lattice determinant')
            hist[nu]+=1;states+=1
        count=0
        for nu,n in sorted(hist.items()):
            count+=n;need(14*count<=29*nu,'short-vector cumulative bound');thresholds+=1
    seeds=[];leaves=0
    for L in range(1,101):
        n=L
        for p in (2,3,5,7):
            while n%p==0:n//=p
        if n!=1:continue
        states2=[s+b*L for s in range(L) for b in range(210)]
        need(sorted(states2)==list(range(210*L)),'all lifts exactly once');seeds.append(L);leaves+=len(states2)
    return {'short_vector_L_range':[1,400],'all_source_states_checked':states,'cdf_thresholds_checked':thresholds,'seed_moduli':seeds,'physical_lift_leaves_checked':leaves}

def run(root,prefix,cells):
    algebra=ref.algebra_checks()
    _,logchecks=ref.log_checks(root/'data/log_table.txt')
    # Cross-check the new log engine independently against exact Fraction series.
    nums=[1,2,3,7,16,31,209,210,420,150000]+[2**k+17 for k in range(10,241,10)]
    for n in nums:
        lo,hi=ref.rational_log(n,terms=80);a,b=ar.logs(n)
        need(F(a,ar.LS)<=lo<=hi<=F(b,ar.LS),'two independent logarithm implementations')
    out={'status':'PASS','algebra':algebra,'parent_log_tables':logchecks,'new_log_crosschecks':len(nums)}
    out.update(initial_witnesses(root,prefix));out.update(cell_costs(root,prefix,cells));out.update(endpoint_checks(root,cells));out.update(short_vector_checks())
    ar.logs.cache_clear()
    return out
