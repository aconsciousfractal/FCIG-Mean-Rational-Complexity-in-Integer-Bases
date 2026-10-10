"""Full-unit central density with a two-prime marked Farey tail.
The scalar profile is evaluable off the actual support; theorem use is separate.
"""
from pair_parameters import *
from math import prod

@lru_cache(None)
def central_constants(r,s,p):
    need(r<s<p and all(prime_trial(a) for a in (r,s,p)), 'three primes')
    R=r*s*p;ph=(r-1)*(s-1)*(p-1);u=F(ph,R);rho=F(R,(r+1)*(s+1)*(p+1))
    U=0;dm=dp=0
    for n in range(R):
        if n and gcd(n,R)==1:U+=1
        e=R*U-ph*n;dp=max(dp,e);dm=min(dm,e)
    dp=F(dp,R);dmin=F(dm,R);dn=u-dmin;ki=dp+dn;K=max(dp,2*u-dmin,u)
    H=F(0);un=0;hh=F(0)
    for n in range(1,17):
        if gcd(n,R)==1:H+=F(1,n);un+=1
        hh=max(hh,H-u*log_lo(n))
    hh=max(hh,H-u*log_lo(16)+(16*u-un+dp)/16)
    hh=F(upper_tick(hh,1000),1000)
    return dict(R=R,u=u,rho=rho,dp=dp,dminus=dn,ki=ki,K=K,h=hh)

@lru_cache(None)
def coeffs(r,s,p):
    d=parameters(r,s);f=central_constants(r,s,p)
    u=f['u'];rho=f['rho'];rho0=F(d['rho']);ec0=2*(f['K']*f['h']+u);ec1=2*f['K']*u
    c=2*AH*rho;theta=AH*(1-rho0);ki=f['ki'];b=F(d['bad_incidence_b']);u0=F(d['u'])
    e0=F(d['deleted_discrepancy'])*F(d['harmonic_h'])+2*F(d['deleted_c']);e1=F(d['deleted_discrepancy'])*u0
    lc=F(1,2)-AL*rho-F(3,4)*u*u
    le=ELL_HI if lc>=0 else ELL_LO
    leading=(1-rho0)/2-3*AL*(1-rho0)/8-AL*rho+F(3,4)*u*u+lc*le+phi_hi(5*u0/p)
    C0=c*(1+ELL_HI)+ec0*(1+ELL_HI/2)+2*u*ki*(3-2*SQ_LO)+ELL_HI*u*ki+(1+5*theta/4+4*e0)*SQ_HI/2+b*PI_HI*SQ_HI
    C1=ec1*(1+ELL_HI/2)+2*e1*SQ_HI;C2=ki*ki*ELL_HI+2*e0;C3=2*e1
    need(C0>=C1 and 2*C2>=C3,'decreasing improved errors')
    out=[F(upper_tick(x),CS) for x in (leading,C0,C1,C2,C3)]
    need(out[1]>=out[2] and 2*out[3]>=out[4],"rounded coefficients retain monotonicity")
    return out

def improved(M,r,s,p):
    need(isinstance(M,int) and M>=2,"scalar profile argument >=2")
    a=isqrt(M);l=log_hi(a);T,c0,c1,c2,c3=coeffs(r,s,p)
    return T+(c0+c1*l)/a+(c2+c3*l)/(a*a)


def prime_exponents(M, primes):
    need(isinstance(M,int) and M>=2,'integer modulus >=2')
    z=M; es=[]
    for p in primes:
        e=0
        while z%p==0:z//=p;e+=1
        need(e>0,'all support primes present')
        es.append(e)
    need(z==1,'radical must be exactly the three displayed primes')
    return tuple(es)

def true_upper(M,r,s,p):
    prime_exponents(M,(r,s,p))
    return improved(M,r,s,p)
