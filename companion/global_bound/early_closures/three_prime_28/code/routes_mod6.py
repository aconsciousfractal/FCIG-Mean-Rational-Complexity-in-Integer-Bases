"""Independent all-edge basis/mark/condition check plus four directions per edge.
The discrete directions left unchecked here are all regenerated in the C++ replay.
"""
from pathlib import Path
from fractions import Fraction as F
import arithmetic as ar
import reference as ref
need=ar.need
S=10**9;FQ=1<<16;LQ=1<<15;T=24;J=96


def gauss(B):
    Y=tuple(B);V=(1,0,0,1)
    for _ in range(201):
        nu=Y[0]**2+Y[1]**2;nv=Y[2]**2+Y[3]**2
        if nu>nv:
            A=(0,1,1,0);Y=ref.mm(A,Y);V=ref.mm(A,V);nu=nv
        t=(2*(Y[0]*Y[2]+Y[1]*Y[3])+nu)//(2*nu)
        if not t:break
        A=(1,0,-t,1);Y=ref.mm(A,Y);V=ref.mm(A,V)
    else:raise ArithmeticError('independent reduction termination')
    if ref.det(Y)<0:
        A=(1,0,0,-1);Y=ref.mm(A,Y);V=ref.mm(A,V)
    need(ref.mm(V,B)==Y and ref.det(V)==1,'oriented unimodular basis transport')
    return Y,V


def rounddir(a,b):
    R=max(abs(a),abs(b));need(R>0,'nonzero direction')
    if abs(a)>=abs(b):
        if a<0:a,b=-a,-b
        n=T*(a+b);idx,n=divmod(n,a)
    else:
        if b<0:a,b=-a,-b
        n=T*(b-a);i,n=divmod(n,b);idx=2*T+i
    return idx%J,(n*FQ)//R,R


def check(root,prefix):
    shapes=[r[2:] for r in ar.shape_rows(root)]
    den,rb,_=ref.load_partition(root/'data/partition.txt')
    seen=0;ordinary=fallback=dc=0;maxcond=0;independent_log_checks=0
    for line in Path(str(prefix)+'.audit').read_text().splitlines():
        v=list(map(int,line.split()));e,d,l=v[:3];Y=tuple(v[3:7]);dest,K,base,error,di,fr=v[7:13]
        need(e==seen,'all routes in canonical order');seen+=1
        H,r=ref.BR[e%12];raw=ref.mm(H,shapes[e//12]);yy,V=gauss(raw)
        need(yy==Y and l==ref.idx(ref.mark(ref.mv(V,r))),'basis and mark independently recovered')
        D=ref.det(Y);nu=Y[0]**2+Y[1]**2
        # Integer comparison to all exact partition boundaries.
        idx=next(i for i,b in enumerate(rb) if nu*den<=b*D);kr=idx-160
        shear=(40*(Y[0]*Y[2]+Y[1]*Y[3])+nu)//(2*nu)
        target=-1 if kr< -125 else (kr+125)*21+shear+10
        need(target==d and -10<=shear<=10,'declared partition and coverage')
        if d<0:
            fallback+=1;need(dest==-1 and len(v)==13,'fallback representation')
            ii,ff,R=rounddir(-Y[1],Y[0]);need((ii,ff)==(di,fr),'fallback direction')
            a,b=ar.logs(5);c,z=ar.logs(D);a2,b2=ar.logs(max(abs(Y[0]),abs(Y[1])));a3,b3=ar.logs(nu)
            need(F(b+ar.ceildiv(z,2)+b2-a3,ar.LS)<=F(base,S),'fallback logarithmic base')
            x=fr if fr>=FQ//2 else fr+1
            bend=F(x*(FQ-x),FQ*FQ*T*T)
            need(F(error,S)>=F(4*nu,D)+bend+F(1,FQ)+F(1,1<<28),'fallback curvature and both quantizations')
            continue
        ordinary+=1;need(dest==12*d+l and len(v)==29,'ordinary route representation')
        C=shapes[d]
        N=(C[3]*Y[0]-C[1]*Y[2],C[3]*Y[1]-C[1]*Y[3],-C[2]*Y[0]+C[0]*Y[2],-C[2]*Y[1]+C[0]*Y[3])
        nd=ref.det(N);tr=sum(x*x for x in N)
        need(nd>0 and K>=S and (K*K+S*S)*nd*nd>=K*S*(tr*tr-2*nd*nd),'condition certificate polynomial')
        maxcond=max(maxcond,K)
        for z in range(13,29,4):
            j,ii,ff,cost=v[z:z+4];need(j==[0,23,47,95][(z-13)//4],'selected direction audit index')
            u=(T,j-T) if j<2*T else (3*T-j,T)
            a=u[0]*N[0]+u[1]*N[2];b=u[0]*N[1]+u[1]*N[3]
            i2,f2,R=rounddir(a,b);need((i2,f2)==(ii,ff),'exact direction bin/fraction')
            la,ua=ar.logs(R*R);lb,ub=ar.logs(T*T*D*ref.det(C))
            need(F(ua-lb,2*ar.LS)<=F(cost,LQ),'independent log correction bound')
            independent_log_checks+=1;dc+=1
    need(seen==2688*12,'all geometric edges checked')
    return {'all_geometric_edges':seen,'ordinary_edges':ordinary,'fallback_edges':fallback,
            'condition_polynomials_checked':ordinary,'max_condition_num':maxcond,'condition_scale':S,
            'direction_routes_checked':dc,'directions_per_ordinary_edge':4,
            'log_corrections_independent_series':independent_log_checks,
            'all_direction_routes_are_not_independently_rechecked':True}
