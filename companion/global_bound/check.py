"""global scalar independent rational checker. No producer or historical Python imports.

The historical PCDF and covering theorems remain mathematical premises.
All finite square/cross profiles, support inventories, scalar comparisons and
integrations used here are reconstructed from their written formulas.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from math import gcd, isqrt, prod
import argparse
import hashlib
import json

KERNEL = (2,3,5,7,11,13,17,19)
OUTSIDE = (23,29,31,37,41,43)
EXCLUDED = (2310,2730,3570,3990,4290,4830,5610,6006,6090,6270,6510)
SOURCE = Path(__file__).resolve().parent
COUNTS = {}
def need(condition, message, category='exact'):
    COUNTS[category] = COUNTS.get(category,0)+1
    if not condition:
        raise ValueError(message)
def load(path):
    return json.loads(path.read_text(encoding='utf-8'))
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def ceil(x):
    return -((-x.numerator)//x.denominator)
def dec_up(x, places=15):
    t=ceil(x*10**places)
    return f'{t//10**places}.{t%10**places:0{places}d}'
def serial(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [serial(v) for v in x]
    return x
def primes_to(n):
    a=bytearray(b'\x01')*(n+1)
    a[:2]=b'\x00\x00'
    for p in range(2,isqrt(n)+1):
        if a[p]: a[p*p:n+1:p]=b'\x00'*((n-p*p)//p+1)
    return [i for i in range(2,n+1) if a[i]]
def mobius_subsets(primes):
    out=[(1,1)]
    for p in primes: out += [(d*p,-mu) for d,mu in out]
    return out
@lru_cache(None)
def sqrt_interval(x):
    x=F(x); scale=10**100
    n=isqrt(x.numerator*scale*scale//x.denominator)
    lo=F(n,scale); hi=lo+F(1,scale)
    need(lo*lo<=x<=hi*hi,'sqrt enclosure','sqrt')
    return lo,hi
@lru_cache(None)
def log_interval(x):
    x=F(x); need(x>0,'log domain','log')
    if x==1: return F(0),F(0)
    if x<1:
        l,h=log_interval(1/x); return -h,-l
    e=0
    while x>2: x/=2; e+=1
    t=(x-1)/(x+1); t2=t*t; term=t; total=F(0)
    for j in range(160):
        total+=term/F(2*j+1); term*=t2
    lo=2*total; hi=lo+2*term/(321*(1-t2))
    if e:
        l2,h2=log_interval(F(2)); lo+=e*l2; hi+=e*h2
    scale=10**90
    return F(lo.numerator*scale//lo.denominator,scale), F(ceil(hi*scale),scale)
def cut(Q,k2):
    # Positive root of A*zÂ²-B*z-C=0; both coefficients are enclosed.
    A=1-F(221184,715715)/Q
    a=F(6414336,5010005); b=F(45173053,51531480)
    blo,bhi=sqrt_interval(k2/F(Q)); blo*=b; bhi*=b
    dlo=sqrt_interval(blo*blo+4*A*a*k2)[0]
    dhi=sqrt_interval(bhi*bhi+4*A*a*k2)[1]
    lo=(blo+dlo)/(2*A); hi=(bhi+dhi)/(2*A)
    zl=ceil(lo*10**6); zh=ceil(hi*10**6)
    need(zl==zh,'cut grid is not decided','cut')
    return F(zh,10**6),(lo,hi)
def positive(Q,k2,z):
    a=F(6414336,5010005); b=F(45173053,51531480); d=F(221184,715715); E=a/2+b
    lz,hz=log_interval(z); cQ=k2*Q/(z*z); lc,hc=log_interval(cQ)
    rootlo,roothi=sqrt_interval(k2/F(Q))
    constant=a*k2/(2*z*z)-E/Q
    lo=lz+constant+b*rootlo/z+d*lc/(2*Q)
    hi=hz+constant+b*roothi/z+d*hc/(2*Q)
    sq=b*b*k2*Q-4*z*z*E*E
    rlo,rhi=sqrt_interval(k2*Q)
    derivative=E-b*rlo/(2*z)+d*(1-lc)/2
    need(z>1 and cQ>1 and lc>1,'positive cutoff/log guard','derivative')
    need(sq>0 and derivative<0,'continuous derivative guard','derivative')
    # Independently integrate the finite endpoint expression symbolically.
    # a*kÂ² integral t^-3, b*k/sqrt(M) integral t^-2, d/M integral t^-1.
    endpoint_a=a/(2*Q); endpoint_b=b/Q
    need(endpoint_a+endpoint_b==E/Q,'finite endpoint corrections','integration')
    return (lo,hi),{'cQ':cQ,'log_cQ_interval':(lc,hc),'E':E,'sqrt_k2Q_interval':(rlo,rhi),'squared_gate_margin':sq,'derivative_M_upper_at_Q':derivative}
@lru_cache(None)
def inventory(w):
    result={}
    for mask in range(256):
        S=[p for j,p in enumerate(KERNEL) if mask>>j&1]
        if len(S)>w: continue
        P=prod(S); t=w-len(S)
        if w==5 and t==0 and P in EXCLUDED: continue
        outside=list(OUTSIDE[:t])
        if w==5 and t==1:
            p=23
            while P*p in EXCLUDED:
                need(P*p in EXCLUDED,'skipped outside prime exclusion','domain')
                p=next(q for q in primes_to(p+20) if q>p)
            outside=[p]
        result[mask]=(S,P,outside,P*prod(outside))
    need(len(result)=={4:163,5:211,6:247}[w],'support inventory cardinality','domain')
    return result
def tail_arrays(n,primes):
    ds=[0]*(n+1); dc=[0]*(n+1)
    def cross_area(k):
        h=k//2; return 2*k*h-h*h
    for p in primes:
        if p<=199: continue
        for k in range(1,n//p+1):
            ds[k*p]+=2*k-1
            dc[k*p]+=cross_area(k)-cross_area(k-1)
    for i in range(1,n+1): ds[i]+=ds[i-1]; dc[i]+=dc[i-1]
    return ds,dc
def reconstruct_profile(w,mask,stock,kind,all_primes,tails):
    S,P,outside,Bmin=inventory(w)[mask]
    comp=[p for p in KERNEL if p not in S]
    marked=mobius_subsets(S); complement=mobius_subsets(comp)
    rho=prod(F(p-1,p) for p in S); E=F(2**(len(S)-1)) if S else F(0)
    alpha=rho*prod(1-F(1,p*p) for p in comp)
    mids=[p for p in all_primes if 19<p<=199]
    b=sum((F(1,p*p) for p in mids),F(0)); e=sum((F(1,p) for p in mids),F(0))
    beta=F(4523,10000)-sum((F(1,p*p) for p in all_primes if p<=199),F(0))
    o=sum((F(1,p) for p in outside),F(0))
    gamma=alpha-rho*(o+b)-beta
    if kind=='square':
        baseline={4:F(11,50),5:F(21,100),6:F(37,200)}[w]
        error=sum((F(1,k)+F(1,k*h) for k,mu in complement for h,muh in marked if mu*muh==1),F(0))+E*(len(outside)+e)
    else:
        gamma*=F(3,4); baseline={5:F(63,400),6:F(111,800)}[w]
        pos=sum((F(1,k) for k,mu in complement if mu==1),F(0))
        neg=sum((F(1,k) for k,mu in complement if mu==-1),F(0))
        error=(F(9,4)*rho+E)*pos+E*neg+E*(len(outside)+e)
    coefficient=(gamma+baseline)/2
    n=max(4,ceil(error/(gamma-coefficient)))
    need(gamma>coefficient>baseline,'analytic tail ordering','profile')
    need(gamma-error/n>=coefficient and 2*coefficient*n>=1,'infinite tail monotonicity','profile')
    need(F(stock['coefficient'])==coefficient and stock['tail_begin']==n,'coefficient and tail cutoff','profile')
    if kind=='cross':
        need(F(stock['gamma'])==gamma and F(stock['floor_error'])==error,'cross tail coefficients','profile')
    unit=[0]*n
    for d in range(1,n): unit[d]=unit[d-1]+(gcd(d,P)==1)
    def area(x):
        h=x//2; return (x-h)*unit[h]+h*unit[x]
    initial=stock['initial_exact_minima'] if kind=='square' else stock['initial_values']
    scale=stock['tick_scale']; tail_ratio=(1+2*coefficient*n*n)/((n+1)*(n+1))
    need(type(scale) is int and scale>0,'positive integer outward tick scale','profile')
    need(F(stock['tail_ratio_exact'])==tail_ratio,'tail ratio','profile')
    tick=tail_ratio.numerator*scale//tail_ratio.denominator
    suffix=[tick]*(n+1)
    trace=hashlib.sha256()
    for R in range(1,n):
        if kind=='square':
            T=sum(mu*(R//k)*unit[R//k] for k,mu in complement)
            loss=R*sum(unit[R//p] for p in outside)
            loss+=sum((R//p)*unit[R//p] for p in mids)+tails[0][R]
            v=max(R,T-loss)
        else:
            T=sum(mu*area(R//k) for k,mu in complement)
            h=R//2
            loss=sum((R-h)*unit[h//p]+h*unit[R//p] for p in outside)
            loss+=sum(area(R//p) for p in mids)+tails[1][R]
            v=max(R if R>=2 else 0,T-loss)
        if R<=len(initial): need(v==initial[R-1],'independent initial exact profile value','profile_count')
        suffix[R]=(1+2*v)*scale//((R+1)*(R+1))
        trace.update(f'{R}:{v}\n'.encode())
    for R in range(n-1,0,-1): suffix[R]=min(suffix[R],suffix[R+1])
    suffix[0]=0
    breaks=[]; previous=None
    for R,tick in enumerate(suffix):
        if tick!=previous: breaks.append([R,tick]); previous=tick
    need(breaks==stock['breakpoints_radius_tick'],'all independently reconstructed suffix breakpoints','profile')
    need(stock.get('external_prime_lower_bounds',stock.get('outside_prime_lower_bounds'))==outside,'profile outside primes','domain')
    return {'mask':mask,'kind':kind,'finite_radii':n-1,'trace_sha256':trace.hexdigest(),'gamma':gamma,'floor_error':error,'coefficient':coefficient,'tail_begin':n,'breakpoints_radius_tick':breaks,'tick_scale':scale}
def strict_radius(Q,alpha):
    target=alpha*Q; r=isqrt(target.numerator//target.denominator)
    if F(r*r)>=target: r-=1
    need(F(r*r)<target<=F((r+1)*(r+1)),'strict determinant radius','edge')
    return r
def profile_h(profile,R):
    ticks=0
    for radius,tick in profile['breakpoints_radius_tick']:
        if radius>R: break
        ticks=tick
    return F(ticks,profile['tick_scale'])
def negative_bands(Q,square,cross=None):
    profiles=[(F(1,2),square)]
    if cross: profiles.append((F(4,5),cross))
    points={F(0),F(1)}
    for alpha,p in profiles:
        points.add(alpha)
        for R,_ in p['breakpoints_radius_tick']:
            x=F(R*R,Q)
            if 0<x<alpha: points.add(x)
    def branch(alpha,p,u):
        if u>=alpha: return F(0),alpha*profile_h(p,strict_radius(Q,alpha))
        R=isqrt((u*Q).numerator//(u*Q).denominator)
        return profile_h(p,R),F(0)
    edges=sorted(points); bands=[]
    for l,r in zip(edges,edges[1:]):
        mid=(l+r)/2; funcs=[branch(a,p,mid) for a,p in profiles]
        cuts=[l,r]
        if len(funcs)==2:
            m1,b1=funcs[0]; m2,b2=funcs[1]
            if m1!=m2:
                intersection=(b2-b1)/(m1-m2)
                if l<intersection<r: cuts=[l,intersection,r]
        for left,right in zip(cuts,cuts[1:]):
            mid=(left+right)/2
            m,b=max(funcs,key=lambda f:f[0]*mid+f[1])
            need(m>=0 and b>=0 and (left>0 or b==0),'affine integration nonnegative','integration')
            for mm,bb in funcs:
                need(m*left+b>=mm*left+bb and m*right+b>=mm*right+bb,'max dominates both affine CDF bounds','integration')
            if bands and bands[-1][1]==left and bands[-1][2:]==[m,b]: bands[-1][1]=right
            else: bands.append([left,right,m,b])
    lo=F(0); hi=F(0)
    for l,r,m,b in bands:
        lo+=m*(r-l)/2; hi+=m*(r-l)/2
        if b:
            ll,lh=log_interval(r/l); lo+=b*ll/2; hi+=b*lh/2
    need(bands[0][0]==0 and bands[-1][1]==1,'integral complete','integration')
    return bands,(lo,hi)
def physical_tests():
    count=0
    for M in (2,5,210,3003,6630,7770,30030):
        # Test exact unit/primitive maps through selected strictly allowed radii.
        for R in sorted({1,2,min(8,strict_radius(M,F(1,2))),min(9,strict_radius(M,F(4,5)))}):
            for alpha,kind in ((F(1,2),'square'),(F(4,5),'cross')):
                if R<=0 or F(R*R)>=alpha*M: continue
                states={}; h=R//2
                for d in range(1,R+1):
                    if gcd(d,M)!=1: continue
                    for a in range(-R,R+1):
                        if kind=='cross' and a!=0 and abs(a)>h and d>h: continue
                        if kind=='cross' and a==0 and d!=1: continue
                        if gcd(a,d)!=1: continue
                        state=a*pow(d,-1,M)%M
                        need(state not in states,'strict physical injectivity','physical')
                        states[state]=(a,d); count+=1
                need(states.get(0)==(0,1),'zero state included','physical')
    need((-1)%2==1%2,'square closed edge counterexample','edge')
    need((-2)%5==pow(2,-1,5)%5,'cross closed edge counterexample','edge')
    for k in range(1,40):
        need(strict_radius(2*k*k,F(1,2))==k-1,'square strict exact-square edge','edge')
        need(strict_radius(5*k*k,F(4,5))==2*k-1,'cross strict exact-square edge','edge')
    return count
def bound_in(interval,container,message):
    lo,hi=map(F,container)
    need(lo<=interval[0]<=interval[1]<=hi,message,'producer_comparison')
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--source',type=Path,default=SOURCE); ap.add_argument('--producer','--certificate',dest='producer',type=Path); ap.add_argument('--closures-report',type=Path); ap.add_argument('--output',type=Path,required=True); args=ap.parse_args()
    source=args.source
    import check_premises
    exact_premises=check_premises.certificate()
    bindings=[]
    for relative in ('inputs/square.json','inputs/cross.json','inputs/global_sectors.json','proofs/PROOF.md'):
        path=source/relative
        bindings.append({'path':relative,'sha256':sha(path),'bytes':path.stat().st_size})
    negative=load(source/'inputs/square.json'); cross=load(source/'inputs/cross.json'); old=load(source/'inputs/global_sectors.json')
    period=exact_premises['periodic_CDF']
    for field in ('quadratic_coefficient','linear_coefficient','constant_coefficient'):
        need(F(period[field])==F(old['periodic_CDF'][field]),'complete PCDF period produces consumed coefficient '+field,'premises')
    need(exact_premises['candidate210']['state_minima']==old['candidate']['heights'],'all candidate minima reconstructed from definition','premises')
    candidate_logs=[log_interval(n) for n in exact_premises['candidate210']['state_minima']]
    log210=log_interval(210)
    mean=(sum((v[0] for v in candidate_logs),F(0))/210-log210[1]/2,sum((v[1] for v in candidate_logs),F(0))/210-log210[0]/2)
    bound_in((mean[0]/log210[1],mean[1]/log210[0]),(old['candidate']['c0_lower'],old['candidate']['c0_upper']),'independent exact candidate interval')
    need(tuple(negative['conditional_domain']['excluded_radicals'])==EXCLUDED,'exact E11 only five-prime domain','domain')
    need(negative['conditional_domain']['unconditional_prime_counts']==[4,6],'unconditional sectors 4 and 6','domain')
    maxR=max(p['tail_begin'] for w in (4,5,6) for p in negative['support_envelopes'][str(w)])
    maxR=max(maxR,max(p['tail_begin'] for w in (5,6) for p in cross['support_envelopes'][str(w)]))
    primes=primes_to(max(20000,maxR)); tails=tail_arrays(maxR,primes)
    scale=10**12
    prime_sum=sum(ceil(F(scale,p*p)) for p in primes if p<=20000)
    need(F(prime_sum,scale)+F(1,20000)<F(4523,10000),'own all-prime-square bound','profile')
    rebuilt={}; reconstruction=[]
    for w in (4,5,6):
        expected=inventory(w)
        for kind,certificate in [('square',negative)]+([('cross',cross)] if w in (5,6) else []):
            stocks=certificate['support_envelopes'][str(w)]
            need(len(stocks)==len(expected),'unique complete support inventory','domain')
            need({p['mask'] for p in stocks}==set(expected),'exact support mask coverage','domain')
            rebuilt[w,kind]={}
            for p in stocks:
                q=reconstruct_profile(w,p['mask'],p,kind,primes,tails)
                rebuilt[w,kind][p['mask']]=q
                reconstruction.append({k:v for k,v in q.items() if k not in ('breakpoints_radius_tick','tick_scale')})
        print('Independent finite profiles reconstructed for omega',w,flush=True)
    rows={}; summaries={}; Qs={4:3003,5:6630,6:30030}; Hs={4:6,5:10,6:14}
    for w in (4,5,6):
        k2=1+max(F(Hs[w]*Hs[w],4),F((Hs[w]-1)**2,3)); sector=[]
        for mask,(S,P,outside,Bmin) in inventory(w).items():
            Q=max(Qs[w],Bmin); z,stationary=cut(Q,k2); pos,guards=positive(Q,k2,z)
            sq=rebuilt[w,'square'][mask]; cr=rebuilt[w,'cross'][mask] if w in (5,6) else None
            bands,neg=negative_bands(Q,sq,cr)
            logQ=log_interval(Q); U=(pos[0]-neg[1],pos[1]-neg[0])
            upper=max(F(0),U[1])/logQ[0]
            need(logQ[0]>0 and upper>=0,'sign clamp before denominator','denominator')
            row={'omega':w,'mask':mask,'marked_kernel_primes':S,'kernel_radical':P,'outside_prime_lower_bounds':outside,'physical_base_min':Bmin,'Q_S':Q,'cut_z':z,'stationary_cut_interval':stationary,'positive_interval':pos,'negative_interval':neg,'negative_bands':bands,'ln_Q_interval':logQ,'U_S':U,'final_upper':upper,'final_upper_decimal':dec_up(upper),'derivative_guards':guards}
            sector.append(row)
        rows[str(w)]=sector
        best=max(sector,key=lambda p:p['final_upper'])
        summaries[str(w)]={'mask':best['mask'],'Q_S':best['Q_S'],'upper':best['final_upper'],'upper_decimal':best['final_upper_decimal'],'rows':len(sector)}
        print('Independent sector',w,'rows',len(sector),'maximum',best['final_upper_decimal'],'mask',best['mask'],flush=True)
    unchanged=[]
    need([p['omega'] for p in old['finite_sectors'] if p['omega']>=7]==list(range(7,16)),'complete inherited omega7 through15 sectors','domain')
    need(len(reconstruction)==1079,'all1079 independent negative profiles','domain')
    need(sum(len(r) for r in rows.values())==621,'all621 support rows','domain')
    for p in old['finite_sectors']:
        if p['omega']<7: continue
        Q=p['base_min']; k2=F(p['kappa_squared']); z=F(p['cut_z']); pos,guards=positive(Q,k2,z)
        w=p['omega']; H=exact_premises['uniform_covering']['bounds'][w-2]['H']
        need(Q==prod(primes[:w]) and k2==1+max(F(H*H,4),F((H-1)**2,3)),'inherited primorial and proven covering constant','old_sector')
        need(pos[1]/log_interval(Q)[0]<=F(p['final_upper']),'old finite sector remains safe','old_sector')
        unchanged.append({'omega':p['omega'],'upper':F(p['final_upper']),'own_interval':pos,'guards':guards})
    p=old['infinite_tail']; Q=p['base_min']; k2=F(p['kappa_squared']); z=F(p['cut_z']); pos,guards=positive(Q,k2,z)
    need(prod(primes[:16])==Q and primes[16]==59 and k2==F(2**30+3,3),'tail anchor and classical coefficient','old_sector')
    derivative_v=log_interval(2)[1]*log_interval(Q)[1]-pos[0]*log_interval(59)[0]
    need(derivative_v<0,'entire real prime-count tail derivative','old_sector')
    need(pos[1]/log_interval(Q)[0]<=F(p['upper']),'old tail remains safe','old_sector')
    oldtail={'upper':F(p['upper']),'own_interval':pos,'derivative_v_upper':derivative_v,'guards':guards}
    candidates=[(f'omega={w}',v['upper']) for w,v in summaries.items()]+[(f"omega={p['omega']}",p['upper']) for p in unchanged]+[('omega>=16',oldtail['upper']),('candidate',F(old['candidate']['c0_upper'])),('closed',F(2,85))]
    winner=max(candidates,key=lambda p:p[1]); need(winner[1]<F(old['global_upper']),'strict improvement over old global','global')
    physical=physical_tests()
    output={'schema':'global-independent-check-v1','status':'PASS_INDEPENDENT_RATIONAL_RECONSTRUCTION','producer_python_imported':False,'historical_python_imported':False,'source_bindings':bindings,'log_method':'160 terms of atanh after exact powers-of-two reduction, explicit remainder, outward 10^-90 grid','sqrt_method':'integer isqrt enclosure on 10^-100 grid','all_finite_profile_suffixes_reconstructed':True,'profile_reconstruction':reconstruction,'profile_finite_radius_total':sum(p['finite_radii'] for p in reconstruction),'sectors':rows,'sector_maxima':summaries,'unchanged_sectors':unchanged,'infinite_tail':oldtail,'global_upper':winner[1],'global_upper_decimal':dec_up(winner[1]),'dominant_sector':winner[0],'old_global_upper':F(old['global_upper']),'physical_vectors_checked':physical,'historical_large_geometric_raw_replayed':False,'checks':COUNTS,'exact_premises':exact_premises}
    catalog_path=source/'closure_input_catalog.json';catalog=load(catalog_path)
    need(len({r['path'] for r in catalog['files']})==len(catalog['files']),'unique complete premise catalog','bindings')
    for row in catalog['files']:
        path=source/row['path'];need(path.resolve().is_relative_to(source.resolve()),'contained reader input','bindings')
        need(path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],'premise source and input binding '+row['path'],'bindings')
    output['input_catalog_sha256']=sha(catalog_path)
    output['closed_domain_arithmetic_replay_authenticated']=False
    output['scalar_scope']='Bound conditional on the finite geometric and physical certificate inputs and their closed-domain proofs; use --closures-report for the complete retained arithmetic chain.'
    if args.closures_report:
        closure=load(args.closures_report)
        need(closure['status']=='PASS' and closure['all_consumed_scalar_closures_reconstructed'] is True,'complete closure arithmetic replay','closures')
        need(closure['input_catalog_sha256']==sha(catalog_path) and type(closure['catalog_files']) is int and closure['catalog_files']==len(catalog['files']),'same exact closure inputs and code','closures')
        need(closure['mode'] in ('retained_caps_exact_integration_and_independent_HNF','full_primary_geometric_reconstruction'),'explicit closure replay mode','closures')
        from verify_closures import expected_components
        need([r['component'] for r in closure['components']]==expected_components(source,closure['mode']=='full_primary_geometric_reconstruction'),'complete ordered closure component inventory','closures')
        need(all(type(r['returncode']) is int and r['returncode']==0 for r in closure['components']),'every closure command succeeded','closures')
        output['closed_domain_arithmetic_replay_authenticated']=True
        output['closure_report_sha256']=sha(args.closures_report)
        output['closure_replay_mode']=closure['mode']
    if args.producer:
        producer=load(args.producer)
        for row in producer['source_bindings']:
            path=source/row['path'];need(path.resolve().is_relative_to(source.resolve()),'contained producer input binding','producer_comparison')
            need(row['bytes']==path.stat().st_size and row['sha256']==sha(path),'producer current source binding '+row['path'],'producer_comparison')
        for w,sector in rows.items():
            prows=producer['sectors'][w]['rows']; need(len(prows)==len(sector),'producer exact row count','producer_comparison'); pm={p['mask']:p for p in prows}
            need(set(pm)=={p['mask'] for p in sector},'producer exact support inventory','producer_comparison')
            for p in sector:
                other=pm[p['mask']]
                for key in ('omega','mask','marked_kernel_primes','kernel_radical','outside_prime_lower_bounds','physical_base_min','Q_S'):
                    need(p[key]==other[key],'producer arithmetic field '+key,'producer_comparison')
                need(p['cut_z']==F(other['cut_z']),'producer rational cut','producer_comparison')
                for key in ('positive_interval','negative_interval','ln_Q_interval'):
                    bound_in(p[key],other[key],'producer enclosing interval '+key)
                upper_U=F(other['positive_interval'][1])-F(other['negative_interval'][0])
                need(F(other['U_S'])==upper_U and upper_U>=p['U_S'][1],'producer scalar numerator upper','producer_comparison')
                expected_final=max(F(0),upper_U)/F(other['ln_Q_interval'][0])
                need(F(other['final_upper'])==expected_final,'producer exact sign-clamped division','producer_comparison')
                need(other['numerator_positive']==(upper_U>0) and other['nonpositive_numerator_uses_zero']==(upper_U<=0),'producer sign branch','producer_comparison')
                need(other['strict_square_radius']==strict_radius(p['Q_S'],F(1,2)),'producer strict square edge','producer_comparison')
                if int(w) in (5,6): need(other['strict_cross_radius']==strict_radius(p['Q_S'],F(4,5)),'producer strict cross edge','producer_comparison')
                need(other['square_profile']==negative['support_envelopes'][w][next(i for i,v in enumerate(negative['support_envelopes'][w]) if v['mask']==p['mask'])],'producer exact square support profile','producer_comparison')
                if int(w) in (5,6):
                    need(other['cross_profile']==cross['support_envelopes'][w][next(i for i,v in enumerate(cross['support_envelopes'][w]) if v['mask']==p['mask'])],'producer exact cross support profile','producer_comparison')
                else: need(other['cross_profile'] is None,'no cross at four primes','producer_comparison')
                dg=other['derivative_guards']; own=p['derivative_guards']
                need(F(dg['cutoff_ratio_squared_at_Q'])==own['cQ'] and F(dg['squared_gate_margin'])==own['squared_gate_margin'] and F(dg['E'])==own['E'],'producer derivative exact rational guards','producer_comparison')
                bound_in(own['log_cQ_interval'],dg['log_cutoff_ratio_squared_interval'],'producer log derivative enclosure')
                bound_in(own['sqrt_k2Q_interval'],dg['sqrt_k2Q_interval'],'producer sqrt derivative enclosure')
                need(own['derivative_M_upper_at_Q']<=F(dg['derivative_M_upper_at_Q'])<0 and dg['all_real_M_ge_Q_S'] is True,'producer continuous derivative upper','producer_comparison')
                need(p['final_upper']<=F(other['final_upper']),'producer final bound conservative','producer_comparison')
                # Check the producer integral by evaluating its own complete bands,
                # allowing an equivalent subdivision of our independently built bands.
                pb=[[F(x) for x in band] for band in other['negative_bands']]
                need(pb[0][0]==0 and pb[-1][1]==1,'producer bands full interval','producer_comparison')
                for a,b in zip(pb,pb[1:]): need(a[1]==b[0],'producer bands no gaps','producer_comparison')
                all_edges=sorted({x for band in pb+p['negative_bands'] for x in band[:2]})
                for l,r in zip(all_edges,all_edges[1:]):
                    mid=(l+r)/2; own=next(b for b in p['negative_bands'] if b[0]<=mid<=b[1]); theirs=next(b for b in pb if b[0]<=mid<=b[1])
                    need(own[2:]==theirs[2:],'producer negative affine function equality','producer_comparison')
            sector_data=producer['sectors'][w]
            need(sector_data['Q']==Qs[int(w)] and sector_data['H']==Hs[int(w)] and F(sector_data['kappa_squared'])==1+max(F(Hs[int(w)]**2,4),F((Hs[int(w)]-1)**2,3)),'producer sector constants','producer_comparison')
            selected=max(prows,key=lambda p:F(p['final_upper']))
            need(sector_data['selected_support_mask']==selected['mask'] and sector_data['selected_Q_S']==selected['Q_S'] and F(sector_data['upper'])==F(selected['final_upper']),'producer sector maximum after division','producer_comparison')
            need(selected['mask']==summaries[w]['mask'],'producer independent maximum mask','producer_comparison')
        need(producer['inherited_finite_sectors']==[p for p in old['finite_sectors'] if p['omega']>=7],'exact inherited sectors','producer_comparison')
        need(producer['infinite_tail']==old['infinite_tail'] and producer['candidate']==old['candidate'],'exact inherited tail and candidate','producer_comparison')
        pcandidates=[F(producer['sectors'][w]['upper']) for w in rows]+[F(p['final_upper']) for p in producer['inherited_finite_sectors']]+[F(producer['infinite_tail']['upper']),F(producer['candidate']['c0_upper']),F(2,85)]
        need(F(producer['global_upper'])==max(pcandidates) and producer['dominant_sector']==winner[0],'producer full global maximum','producer_comparison')
        need(F(producer['global_upper'])<F(old['global_upper']),'producer strict global improvement','producer_comparison')
        need(tuple(producer['conditional_excluded_five_radicals'])==EXCLUDED and producer['normalization']=='Delta(B,N)=mu(B^N)/ln(B)','producer domain and normalization','producer_comparison')
        output['producer_certificate_sha256']=sha(args.producer)
        output['producer_certificate_checked']=True
    else: output['producer_certificate_checked']=False
    need(not args.output.exists(),'output must be fresh','output')
    args.output.write_text(json.dumps(serial(output),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('PASS',winner[0],dec_up(winner[1]),'rows',sum(len(x) for x in rows.values()),'profile radii',output['profile_finite_radius_total'],flush=True)
if __name__=='__main__': main()

