"""Recompute finite geometric profiles and all base/length covering gates exactly."""
from pathlib import Path
from fractions import Fraction as F
import arithmetic as ar
from finite import primes_to_113,upper_decimal
S=10**9

def configurations():
    return [(p,3 if p==5 else 2 if p in (7,11,13) else 1) for p in primes_to_113()]
def rational_caps(p):
    return {5:[F(285,1000),F(1,8),F(31,500)],7:[F(107,500),F(67,1000)],11:[F(129,1000),F(19,1000)],13:[F(21,200),F(7,1000)]}.get(p,[F(7,100)])

def profile(p,k,prefix,cap):
    q=6*p;labels=12*(p+1)
    cc=[list(map(int,l.split())) for l in Path(str(prefix)+f'.cells{k}').read_text().splitlines()]
    rr=[list(map(int,l.split())) for l in Path(str(prefix)+f'.radial{k}').read_text().splitlines()]
    ar.need(len(cc)==2688 and len(rr)==128,'complete cells and radial slices')
    for c,x in enumerate(cc):ar.need(x[0]==c and 0<=x[2]<labels and 0<=x[3]<33,'cell record and extremizer')
    for j,r in enumerate(rr):
        group=cc[j*21:(j+1)*21]
        ar.need(r[0]==j-125 and r[3]==max(v[1] for v in group),'radial maximum')
        ar.need(any(v[0]==r[4] and v[1]==r[3] for v in group),'radial extremizer')
        if j:ar.need(rr[j-1][2]==r[1],'radial cover has no gaps')
    b=[F(r[3],S) for r in rr]
    for j in range(126,-1,-1):b[j]=max(b[j],b[j+1])
    A=F(29,14);rho=F(rr[0][1],ar.PDEN)
    # Entire cusp, not a finite radial sample. Use full upper logs for q-1,q and -log rho.
    t=F(ar.logs(q-1)[1],ar.LS)+F(k*ar.logs(q)[1],2*ar.LS)+F(ar.logs(ar.PDEN)[1]-ar.logs(rr[0][1])[0],2*ar.LS)+4*rho
    bulk=b[-1]+sum((b[j]-b[j+1])*min(F(1),A*F(rr[j][2],ar.PDEN)) for j in range(127))
    cusp=A*rho*(max(t-b[0],F(0))+F(1,2));total=bulk+cusp
    ar.need(total<cap,f'full uniform cap p={p}, k={k}')
    return {'p':p,'q':q,'k':k,'labels':labels,'bulk_upper':str(bulk),'cusp_upper':str(cusp),'total_upper':str(total),'rational_cap':str(cap),'upper_decimal':upper_decimal(total),'radial_isotonic_numerators':[int(v*S) for v in b]}

def covering_gates(rows,census):
    mu={(r['p'],r['k']):F(r['rational_cap']) for r in rows};eta=F(117,5000)
    # Each row describes a full infinite region of base--length pairs.
    specs=[(5,1,200000,'N=1, B>=200000'),(5,2,210,'N=2, B>=210'),(5,3,30,'N>=3'),(7,1,10000,'N=1, B>=10000'),(7,2,42,'N>=2'),(11,1,250,'N=1, B>=250'),(11,2,66,'N>=2'),(13,1,90,'N=1, B>=90'),(13,2,78,'N>=2')]
    specs.extend((p,1,6*p,'N>=1') for p in primes_to_113() if p>=17)
    gates=[]
    for p,k,Bmin,domain in specs:
        lo=F(ar.logs(Bmin)[0],ar.LS);upper=mu[p,k]/lo
        ar.need(upper<eta,f'infinite family gate p={p}, k={k}')
        gates.append({'p':p,'k':k,'B_min':Bmin,'domain':domain,'mu_cap':str(mu[p,k]),'delta_upper':str(upper),'margin_to_117_over_5000':str(eta-upper)})
    worstfinite=F(census['worst_delta_upper']);ar.need(worstfinite<eta,'all finite cases under common uniform cap')
    worst=max([worstfinite]+[F(r['delta_upper']) for r in gates])
    ar.need(worst<eta<F(2,85),'common strict separation')
    return {'status':'PASS','primes':primes_to_113(),'prime_count':28,'all_lengths_and_positive_valuations':True,'mu_bounds_are_on_marked_means_when_other_primes_occur_in_M':True,'gate_rows':gates,'finite_cases':census['counts'],'uniform_delta_upper':str(worst),'uniform_delta_upper_decimal':upper_decimal(worst),'clean_rational_cap':'117/5000','original_candidate_lower':'2/85','global_Cstar_identified':False}
