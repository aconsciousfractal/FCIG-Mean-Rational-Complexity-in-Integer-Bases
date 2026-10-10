# Supportwise global bound

Throughout, states s modulo M are uniform and include zero. Define

    R_M(s)=min{max(|a|,d): d>0, gcd(d,M)=1, a=sd (mod M)},
    Z_M(s)=R_M(s)/sqrt(M), mu(M)=E log Z_M,
    Delta(B,N)=mu(B^N)/log B, C*=sup_(B>=2,N>=1) Delta(B,N).

The norm defining R is the infinity norm. All logarithms are natural. Numerator zero is admitted. Every theorem below concerns all positive prime valuations and all integer N>=1. Write P(M)=E(log Z)_+ and L(M)=E(-log Z)_+, so mu=P-L exactly.

## Primitive Euclidean distribution

For Lambda_s={(a,d):a=sd (mod M)} let lambda_s be its shortest nonzero Euclidean length. Its determinant is M. Set Y=sqrt(Mx). For 1<=Y<M a shortest vector has positive d and is primitive in Lambda_s. In coordinates (a,d)=(Mm+sd,d), this means gcd(m,d)=1; it does not mean gcd(a,d)=1. Put h=gcd(d,M), r=d/h, A=sqrt(Y^2-d^2). The exact number of primitive incidences with denominator d is

    c_d = phi(d)/phi(r) * #{b integer: |b|<=A/h, gcd(b,r)=1}.

Indeed a=hb and b=rs+(M/h)m, with r invertible modulo M/h. The h lifts of the state traverse one class modulo r in m and its lifts modulo d. The reduction of units modulo d onto units modulo r is surjective, each fibre having phi(d)/phi(r) elements. For r=1 use phi(1)=1 and include b=0. For every X>=0 the symmetric unit count is at most phi(r)(2X/r+1): complete periods contribute exactly 2phi(r), and the remaining units pair as u,r-u. Consequently, with f_Y(t)=2sqrt(Y^2-t^2)+t,

    M Pr(lambda_s<=Y) <= sum_(1<=d<=floorY) [phi(d)/d] f_Y(d).

Use P=30030 and w_d=prod_(p|P)(1-1_(p|d)/p), which bounds phi(d)/d for every d. Its mean alpha and cumulative discrepancy extrema ell,h are computed over the entire period:

    alpha=442368/715715,
    ell=-134597/1431430, h=744773/1431430.

Periodicity extends these bounds to every integer index. Concavity of 2sqrt(Y^2-t^2) and summing the lower trapezoids gives

    sum_d f_Y(d) <= (pi+1)Y^2/2-Y/2+1/2.

For n=floorY, the trapezoid bound gives the error sqrt(Y^2-n^2); use sqrt(u)<=(u+1)/2 and n+1<=Y+1. The function f_Y is unimodal, starts at 2Y and has maximum sqrt(5)Y. Summation by parts with ell<=D<=h therefore bounds its weighted discrepancy by [(h-ell)sqrt5+2ell]Y. Thus the incidence sum is at most

    aY^2+bY+d,
    a=6414336/5010005, b=45173053/51531480, d=221184/715715.

Here pi<22/7 follows from the positive integral of t^4(1-t)^4/(1+t^2) over [0,1], equal to 22/7-pi, and sqrt5<161/72 follows by squaring. After division by M this proves

    Pr(lambda_s^2/M<=x) <= min(1,a*x+b*sqrt(x/M)+d/M).

If Y<1 the event is empty; if Y>=M then a*x>=a>1, so the same bound is valid in those cases too, including M=1. The additional 3x/2 bound used by finite closure profiles has a separate finite-incidence and tail proof in the premise documentation.

## Uniform covering and coupled shear

For at most k distinct primes, arbitrary forbidden residue classes are simultaneously avoided in every H_k consecutive integers, where

    H_1,...,H_14 = 2,4,6,10,14,22,26,44,58,66,98,148,204,360.

Here is the exhaustive reduction checked by the covering component. Separate the presence of 2, and let S be the actual subset of {3,5,7,11}. If 2 is present pass to the H/2 allowed parity positions, whose indices have one forbidden class for each remaining odd prime. Otherwise retain length H. CRT translates forbidden small-prime classes to zero classes. For all interval starts modulo prod(S), inclusion-exclusion independently counts the small-sieve survivors. Each of the t outside primes removes at most ceil(L/p_j) positions, where p_j is the jth prime greater than 11. Thus it suffices that

    min_start survivors(S,L)-sum_(j<=t) ceil(L/p_j) >=1.

All 399 admissible (k,parity,S) cases satisfy this. Missing primes can be added to extend a smaller support; the conclusion for the larger support implies it for the smaller. This does not assume primorial extremality.

Let u be a shortest nonzero vector in Lambda_s, complete it to a lattice basis u,v and reduce v=theta*u+v_perp to |theta|<=1/2. Then lambda^4(1-theta^2)<=M^2 by minimality of u. The vertical projection of Lambda_s is Z, hence gcd(u_d,v_d)=1. If u_d is a unit, the direct witness is covered by the bound below. Otherwise at least one prime dividing u_d is automatically repaired by every v+j*u. Each remaining prime forbids one class of j. There are at most omega(M)-1 such primes.

Take H=2h from the covering table and t=|theta|. Allowed j have successive gaps at most H. The nearest permitted integer to -theta is at distance at most h-t: if zero is allowed use t<=h-t; otherwise the neighbouring permitted integers are -l,r, l+r<=2h. If l<=h use l-t; if l>=h+1 then r<=h-1 and r+t<=h-t. Hence

    R_M(s) <= kappa(H) M/lambda_s,
    kappa(H)^2=1+max(h^2,(2h-1)^2/3).

To cover the whole real shear interval, for h=1,2,3 use

    h^2(1-t^2)-(h-t)^2=t[2h-(h^2+1)t]>=0;

for h>=4 use

    (2h-1)^2(1-t^2)-3(h-t)^2
      =(1-2t)[2(h^2-h+1)t+(h^2-4h+1)]>=0.

Both endpoints supply the stated maximum. The direct unit branch has cost at most (2/sqrt3)M/lambda, whose square 4/3 is below kappa^2. No argument assumes squarefreeness of M.

## Positive part, finite endpoint and continuous monotonicity

For a rational cut z>1 let k2=kappa^2 and E=a/2+b. Since lambda>=1, the repaired witness has no tail beyond kappa*sqrt(M). Integrating the primitive CDF only to that finite endpoint gives

    P(M) <= C(M,k2,z)
    =log z+a*k2/(2z^2)+b*sqrt(k2/M)/z
      +d*log(k2*M/z^2)/(2M)-E/M.

In particular the constant d/M is never integrated to infinity. The subtraction E/M includes both endpoint corrections. If, at M=Q,

    log(k2*Q/z^2)>1, b^2*k2*Q>4z^2*E^2,

then M^2*dC/dM = E-b*sqrt(k2*M)/(2z)+d[1-log(k2*M/z^2)]/2 is negative for every real M>=Q. Both nonconstant terms become smaller as M increases. Each support row verifies these conditions at its actual rational cut. The positive stationary root is rounded outward to a grid of denominator 10^6; exact optimality of the cut is unnecessary.

## Negative square and cross profiles

Fix K={2,3,5,7,11,13,17,19} and S=K intersect supp(M). This same S is kept for the entire integral. Each stored profile is a lower bound derived by exact Mobius counting and ordered outside-prime loss bounds; the independent checker recomputes every finite radius and every suffix minimum, together with the analytical infinite tail. The formulas and their tail inequalities are given in NEGATIVE_PROFILES.md.

For a square radius R, primitive vectors with denominator unit and |a|<=R, 0<d<=R have pairwise determinants of absolute value at most 2R^2. If 2R^2<M, distinct primitive vectors represent distinct states. The zero vector state is represented by (0,1); all other vectors occur with both numerator signs. The resulting exact count is 1+2V(R). In the cross, with r=floor(R/2), use the union of |a|<=R,0<d<=r and |a|<=r,0<d<=R. Its determinant bound is R^2+r^2<=5R^2/4; strict injectivity holds for 5R^2<4M. Inclusion-exclusion removes the overlap; the physical norm remains max(|a|,d).

Let h_S(R) be the nondecreasing rational suffix envelope below the primitive count divided by (R+1)^2. For M>=Q, the square lower CDF in u is u*h_S(floor sqrt(Q*u)) for u<1/2, continued by h_S(r_strict)/2 thereafter. The cross is the analogous expression with threshold 4/5. The strict radius is the largest integer R with R^2<alpha*Q, including the exact-square boundaries. Since these are lower bounds for the same CDF, take their pointwise maximum, without adding overlapping vector counts.

The negative expectation is half the integral of this lower CDF divided by u on [0,1]. Every breakpoint and square/cross intersection is rational. On each band [l,r], the maximum is m*u+b and its integral is [m(r-l)+b*log(r/l)]/2. If l=0 then b=0, and no logarithm is taken there. The checker verifies coverage of [0,1], equality of the independent affine functions and conservative interval arithmetic.

## Domain partition and supportwise division

Five-prime negative profiles exclude exactly E={2310,2730,3570,3990,4290,4830,5610,6006,6090,6270,6510}; four- and six-prime profiles are unconditional. Every excluded domain has a separate all-length, all-valuation closure certificate. Four-prime families have (core,new-prime minimum)

    (30,7),(42,11),(70,11),(66,13),(105,11),(78,17),
    (110,13),(102,19),(154,13),(165,13),(130,17),(114,23).

The 30p family includes the candidate and has a bound by c0; the other families and the eleven radicals are below 2/85. Two independent exhaustive radical inventories show that 3003 and 6630 are the first unclosed four- and five-prime supports. The <=3 domains are handled by the early-closure component, with the existing prime-power and two-prime theorems. Six primes require at least 30030.

For each residual cardinality w=4,5,6 and mask S, the outside primes are bounded coordinatewise by the first primes greater than 19. When w=5 and there is only one outside prime, skip precisely those whose product with prod(S) belongs to E. When none is outside, omit precisely the masks whose radical is in E. This gives 163,211,247 rows respectively. Define

    Q_S=max(Q_w,prod(S)*prod(outside_lower_bounds)).

For every base in that row, B>=Q_S and M=B^N>=Q_S, regardless of valuations or N. Keep S fixed in C(Q_S), the negative integral L_S(Q_S), and the divisor log Q_S. Set U_S=C(Q_S)-L_S(Q_S). If U_S>0 then Delta<=U_S/log Q_S; if U_S<=0 then Delta<=0. The safe interval upper is therefore max(0,U_hi)/(log Q_S)_lo. The sign test precedes denominator substitution.

For w=7,...,15 integrate the positive bound at the corresponding primorial and use L>=0. For w>=16 use the classical covering H<=2^(w-1) and the weaker shear constant (H^2+3)/3. Write v=w-16, P16=prod(first16primes), and next prime=59. The coefficient is dominated by kappa16*2^v, the cut scales by 2^v, and B>=P16*59^v. The upper has the form

    [C(P16,kappa16,z16)+v*log2]/[log P16+v*log59].

Its derivative has sign log2*log P16-C(P16)*log59, checked strictly negative by rational logarithm intervals. Thus its anchor bounds the entire real v>=0 tail. Finite sectors end at 15 and the tail starts at 16. The coupled shear coefficient is not assumed to scale this way.

Finally maximize all support rows, the nine remaining finite sectors, the infinite tail and all closed/candidate domains. This yields C*<0.241933621228105. Exact C*=c0, attainment, and a classification of the sign of every offset are not asserted.

The scalar checker is independent of the producer's implementation. Exact rational logarithm and square-root enclosures are numerical proof tools, not formalization of the written mathematical lemmas. Full reconstruction of finite geometric premises is separately indexed and cannot be replaced by sampled audits or a receipt saying PASS.
