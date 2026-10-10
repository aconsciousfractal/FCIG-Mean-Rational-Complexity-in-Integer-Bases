"""Two-prime marked Farey bound; exact parameter generation and family gates."""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd,isqrt
from exact import *
AL,AH=alpha_interval()
ELL_LO=log_lo(2);ELL_HI=log_hi(2)
SQ_LO=F(isqrt(2*CS*CS),CS);SQ_HI=SQ_LO+F(1,CS)
PI_HI=F(22,7)
DELTA=F(2,85)
RAD_CUTOFF=10**6
P_MIN=53

@lru_cache(None)
def parameters(r,s):
    need(prime_trial(r) and prime_trial(s) and r<s,'two distinct ordered primes')
    need((r,s)!=(2,3),'family 6p excluded from the new domain')
    q=r*s;ph=(r-1)*(s-1);u=F(ph,q);rho=F(q,(r+1)*(s+1))
    cn=(r-1)*s*s+(s-1)*r*r-(q-1);c=F(cn,q*q)
    U=S=0;dmin=dmax=0;zmin=zmax=0;b=F(0)
    for n in range(q):
        if n and gcd(n,q)==1:U+=1;S+=n
        dn=q*U-ph*n;dmin=min(dmin,dn);dmax=max(dmax,dn)
        if n:
            b=max(b,F(3*dn,q),F(2*q*(2*n*U+S)-5*ph*n*n,2*q*n))
        zn=q*s*(r-1)*(n//r)+q*r*(s-1)*(n//s)-q*(q-1)*(n//q)-cn*n
        zmin=min(zmin,zn-cn);zmax=max(zmax,zn)
    dp=F(dmax,q);dm=u-F(dmin,q);span=dp+dm
    K=max(dp,2*u-F(dmin,q),u);dd=F(max(-zmin,zmax),q*q)
    H=F(0);U16=0;hv=F(0)
    for n in range(1,17):
        if gcd(n,q)==1:H+=F(1,n);U16+=1
        hv=max(hv,H-u*log_lo(n))
    hv=max(hv,H-u*log_lo(16)+(16*u-U16+dp)/16)
    h=F(upper_tick(hv,1000),1000)
    need(dp<=2 and dm<=2 and span<=4 and K<=3 and b<=8,'universal local bounds')
    need(0<c<F(1,2) and dd<=2,'universal deleted-density bounds')
    ec0=2*(K*h+u+dp);ec1=2*K*u
    de0=dd*h+2*c;de1=dd*u
    cdf=2*AH*rho;theta=AH*(1-rho)
    beta=3-2*SQ_LO
    # T0 decreases with alpha; choose the endpoint of log(2) according to its sign.
    lc=F(1,2)-AL*rho-F(3,4)*u*u
    le=ELL_HI if lc>=0 else ELL_LO
    t0=F(1,2)-3*AL/8-rho*(5*AL/8+F(1,2))+F(3,4)*u*u+lc*le
    C0=cdf*(1+ELL_HI)+ec0*(1+ELL_HI/2)+2*u*span*beta+ELL_HI*u*span+(1+5*theta/4+4*de0)*SQ_HI/2+b*PI_HI*SQ_HI
    C1=ec1*(1+ELL_HI/2)+2*de1*SQ_HI
    C2=span*span*ELL_HI+2*de0;C3=2*de1
    need(C0>=C1 and 2*C2>=C3,'decreasing errors for A>=1')
    d={'r':r,'s':s,'q':q,'u':str(u),'rho':str(rho),'Dplus':str(dp),'Dminus':str(dm),'interval_span':str(span),'primitive_floor_K':str(K),'bad_incidence_b':str(b),'deleted_c':str(c),'deleted_discrepancy':str(dd),'harmonic_h':str(h),'T0_tick':upper_tick(t0),'C0_tick':upper_tick(C0),'C1_tick':upper_tick(C1),'C2_tick':upper_tick(C2),'C3_tick':upper_tick(C3)}
    return d

def bound(M,r,s,p,check_domain=True):
    need(M>=2 and p>=7 and p>s and r<s,'three-prime analytic domain')
    if check_domain:
        need(prime_trial(p),'third parameter prime')
        z=M
        for t in (r,s,p):
            need(z%t==0,'every displayed prime must divide the modulus')
            while z%t==0:z//=t
        need(z==1,'radical must be exactly the three displayed primes')
    d=parameters(r,s);u=F(d['u']);a=isqrt(M);l=log_hi(a)
    return (F(d['T0_tick'],CS)+phi_hi(5*u/p)+u*(1+ELL_HI)/p
            +(F(d['C0_tick'],CS)+F(d['C1_tick'],CS)*l)/a
            +(F(d['C2_tick'],CS)+F(d['C3_tick'],CS)*l)/(a*a))

def generic(M,p):
    need(M>=2 and p>=P_MIN,'generic domain')
    a=isqrt(M);l=log_hi(a)
    return F(parameters(2,5)['T0_tick'],CS)+phi_hi(F(2,p))+F(2,5)*(1+ELL_HI)/p+(70+14*l)/a+(18+4*l)/(a*a)

def constants():
    need(SQ_LO*SQ_LO<2<SQ_HI*SQ_HI,'sqrt2 rational enclosure')
    need(F(3,5)<AL<AH<F(2,3),'zeta bounds')
    need(F(2,3)<ELL_LO<=ELL_HI<F(7,10),'log2 bounds')
    # Monotonicity of T0 in inverse primes on [0,1/2]^2.
    Klo=F(3,5)*(F(5,8)+F(2,3))+F(1,2)
    need(Klo>F(16,27),'T0 monotonicity coefficient')
    t10=F(parameters(2,5)['T0_tick'],CS)
    need(t10<-F(3,50),'T0 upper outside pair (2,3)')
    C0=F(4,3)*F(17,10)+12*F(27,20)+F(8,5)+F(14,5)+(1+F(5,12)+12)*F(5,7)+8*F(22,7)*F(99,70)
    C1=6*F(27,20)+8*F(5,7)
    need(C0<70 and C1<14 and 16*F(7,10)+6<18,'generic coefficient ceilings')
    Ksharp=AL*(F(5,8)+ELL_LO)+F(1,2)
    Chi=F(3,4)*(1-ELL_LO)
    derivative_lower=4*Ksharp/9-Chi
    excision_derivative_upper=F(5,P_MIN)*log_hi(1+F(P_MIN,2))+(1+ELL_HI)/P_MIN
    need(8*Ksharp/27>Chi,'decreasing auxiliary derivative on inverse-prime box')
    need(derivative_lower>excision_derivative_upper,'leading mean monotonicity including prime repair')
    g=generic(RAD_CUTOFF,P_MIN);tar=DELTA*log_lo(RAD_CUTOFF)
    need(g<tar,'all radical>=1000000 families excluded')
    return {'status':'PASS','alpha_lo':str(AL),'alpha_hi':str(AH),'sqrt2_lo':str(SQ_LO),'sqrt2_hi':str(SQ_HI),'pi_upper':'22/7','T0_q10_upper':str(t10),'monotonicity_coefficient_lower':str(Klo),'generic_C0_upper':str(C0),'generic_C1_upper':str(C1),'largest_prime_threshold':P_MIN,'leading_derivative_lower':str(derivative_lower),'repair_derivative_upper':str(excision_derivative_upper),'derivative_margin':str(derivative_lower-excision_derivative_upper),'radical_cutoff':RAD_CUTOFF,'generic_upper_at_cutoff':str(g),'generic_upper_decimal':out_decimal(g),'target_lower_at_cutoff':str(tar),'strict_margin':str(tar-g)}
