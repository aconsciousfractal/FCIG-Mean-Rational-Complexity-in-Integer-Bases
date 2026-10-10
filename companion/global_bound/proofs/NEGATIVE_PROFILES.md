# Exact support profiles and their infinite tails

Use the kernel K={2,3,5,7,11,13,17,19}, its product P, a fixed physical sub-support S, P_S=prod(S), and t=w-|S|. Let U_S(n)=#{1<=d<=n:gcd(d,P_S)=1}, rho=prod_(p inS)(1-1/p), and E=0 for empty S and 2^(|S|-1) otherwise. For integer n>=0 the Mobius expansion gives |U_S(n)-rho*n|<=E: the positive and negative floor errors each use half of the 2^|S| terms. Write C=K\S, and let ell_j be the ordered outside-prime lower bounds of this row.

For the square define

    T_S(R)=sum_(k|prodC) mu(k) floor(R/k) U_S(floor(R/k)),
    D_S(R)=T_S(R)-R*sum_j U_S(floor(R/ell_j))
      -sum_(19<p<=199) floor(R/p)*U_S(floor(R/p))
      -sum_(199<p<=R) floor(R/p)^2,
    v_S(R)=max(R,D_S(R)).

The Mobius sum removes common kernel prime factors, while retaining the marked denominator-unit condition. An outside physical prime q removes at most R*U_S(floor(R/q)) pairs; replacing q by ell_j increases this loss. Common nonkernel primes use the other two loss sums. Union bounds allow overlaps safely. The baseline R is supplied by (a,1), 1<=a<=R. Hence there are at least 1+2v_S(R) primitive, unit-denominator signed square vectors, including (0,1).

The exact finite T_S is reconstructed by the checker for every radius below its own tail begin. Expanding U_S and bounding floors gives

    T_S(R)>=alpha*R^2-C_S*R,
    alpha=rho*prod_(p inC)(1-1/p^2),
    C_S=sum_(k|prodC,h|P_S,mu(k)mu(h)=1) (1/k+1/(kh)).

For positive terms use floor(x)floor(y)>=xy-x-y; for negative terms use floor(x)floor(y)<=xy. Let b=sum_(19<p<=199)1/p^2, e=sum_(19<p<=199)1/p, beta=4523/10000-sum_(p<=199)1/p^2 and o=sum_j1/ell_j. Then

    D_S(R)>=gamma*R^2-F_S*R,
    gamma=alpha-rho*(o+b)-beta,
    F_S=C_S+E*(t+e).

The prime reciprocal-square bound 4523/10000 is checked from every prime <=20000 with upper rounding to denominator10^12 and an integer-tail integral1/20000. Thus it bounds all primes without imported numerical asymptotics. Each square row uses baseline c_w=11/50,21/100,37/200 for w=4,5,6, respectively; c=(gamma+c_w)/2 and n=max(4,ceil(F_S/(gamma-c))). Explicit rational checks give gamma>c>c_w and gamma-F_S/n>=c, hence v_S(R)>=cR^2 for every R>=n.

For the cross put h=floor(R/2) and f_S(A)=(A-floor(A/2))*U_S(floor(A/2))+floor(A/2)*U_S(A). Then

    T_cross(R)=sum_(k|prodC) mu(k) f_S(floor(R/k)),
    D_cross(R)=T_cross(R)
      -sum_j[(R-h)*U_S(floor(h/ell_j))+h*U_S(floor(R/ell_j))]
      -sum_(19<p<=199) f_S(floor(R/p))
      -sum_(199<p<=R) f_empty(floor(R/p)),
    v_cross(R)=max(R*1_(R>=2),D_cross(R)).

Dividing both coordinates by a common kernel factor preserves the cross because floor(floor(R/2)/k)=floor(floor(R/k)/2). The two arms are counted by splitting at d=h, so there is no double counting. With E as above, for every real x>=0,

    (3rho/4)*x^2-(9rho/4+E)*x
      <=f_S(floor x)<=(3rho/4)*x^2+E*x.

To see this first apply |U_S(n)-rho n|<=E to integer A, noting 3A^2/4-3A/4<=2A floor(A/2)-floor(A/2)^2<=3A^2/4. Then use floor(x)^2>=x^2-2x and floor(x)<=x. Consequently

    gamma_cross=(3/4)*(alpha-rho*(o+b)-beta),
    F_cross=(9rho/4+E)*sum_(k|prodC,mu(k)=1)1/k
      +E*sum_(k|prodC,mu(k)=-1)1/k+E*(t+e).

Use c_w=63/400,111/800 for w=5,6, and the same c,n construction and explicit inequalities as for the square. This proves the infinite cross tail, rather than extrapolating a finite scan.

For either shape the infinite ratio (1+2cR^2)/(R+1)^2 increases when 2cR>=1, checked at n and therefore throughout R>=n. Its value at n is a lower bound for every later radius. The suffix envelope at r is the minimum of this tail value and all finite (1+2v_S(R))/(R+1)^2 for r<=R<n, rounded downward to its recorded tick scale. Set the value at r=0 to zero. This construction is nondecreasing and bounds the entire infinite suffix, not merely the stored radii.

For M>=Q and u below a shape's strict determinant threshold, let R=floor sqrt(Mu). Then R+1>sqrt(Mu), so the primitive count divided by M exceeds u*h_S(R), which is at least u*h_S(floor sqrt(Qu)) by monotonicity. At the strict threshold alpha use the limiting value from below and retain it as a plateau for larger u, by monotonicity of the physical CDF. This proves exactly the lower CDF used in PROOF.md, including exact-square threshold cases and zero.
