# Mean rational complexity in integer bases: a uniform bound and a two-prime theorem

**Oleksiy Babanskyy**

This repository accompanies the [paper](paper/Mean-Rational-Complexity-in-Integer-Bases.pdf)
and its [standalone LaTeX source](paper/Mean-Rational-Complexity-in-Integer-Bases.tex).

For $M=B^N$, $B\ge2$, $N\ge1$, define
$$
R_{B,N}(s)=\min\{\max(|a|,d):a\equiv sd\pmod M,\ d>0,\ \gcd(d,M)=1\}.
$$
All $M$ residues, including zero, have uniform weight; logarithms are real and unrounded.

The paper proves
$$
\mathbb E\log_B R_{B,N}<\frac N2+0.241933621228105
$$
for every integer base and every length. It strengthens this to
$\mathbb E\log_B R_{B,N}\le N/2$ for prime-power bases and to
$\mathbb E\log_B R_{B,N}<N/2$ for bases with exactly two distinct prime factors,
including arbitrary positive valuations. The binary theorem is imported from
the [earlier binary paper](https://github.com/aconsciousfractal/FCIG-Exact-Mean-Bound-for-the-Rational-Complexity-of-Binary-Strings).

The zero correction is false in general: $B=30,N=2$ is an exact integer-product
counterexample. If $\operatorname{rad}(B)=510510$ and $B^N\ge510510^2$,
the centered mean exceeds $1/(1000\log B)$.
For the optimal additive correction,
$$
0.023672606561106<C^*=\sup_{B,N}\left(\mathbb E\log_B R_{B,N}-N/2\right)<0.241933621228105.
$$
Two simpler independent bounds remain available. A corollary using Kanold's classical unit-gap bound as reported by
[Costello–Watts, arXiv:1208.5342v2, pp. 1–2](https://arxiv.org/pdf/1208.5342v2), gives
$$
C^*\le\frac{\log67+1}{2\log30}<\frac{5037}{6400}<\frac45.
$$
Thus the mean is also strictly below $N/2+4/5$ in every base and length.
The original $49/50$ proof remains self-contained. The sharper bound combines supportwise positive and negative log contributions. It is not claimed optimal.

The positive per-base correction tends to zero for large bases; the global
supremum reduces to finitely many bases. A maximizing base has at least three
primes and is not a perfect power. Its value and a maximizing finite length
are not determined.

The general bound addresses [Lin–Xiao–Chen, Question 3.2](https://doi.org/10.3934/math.2026823)
in this explicitly defined full-unit rational-height model. Additional FCSR
connection conditions and rounded register lengths require separate comparison.

## Verify

Python >=3.10 for the original two-prime checks and the scalar computation;
Python >=3.11 for the marked-70 component, including the new global retained chain.
These arithmetic checks use the standard library.
The new full geometric reconstruction requires GNU g++ with C++17, OpenMP and
`__int128`, on little-endian Linux. From the repository root:

~~~text
python -B scripts/check_manifest.py
python -B scripts/reproduce.py smoke
python -B scripts/reproduce.py smoke --optimized
python -B scripts/reproduce.py retained
python -B scripts/reproduce.py retained --optimized
python -B scripts/reproduce.py semiprime
python -B scripts/reproduce.py semiprime --optimized
python -B -m unittest discover -s tests -v
python -B -O -m unittest discover -s tests -v
~~~

The original two-prime proof-premise reconstruction passed normally and with
optimization: all 37 dependency stages, all 4,596 semiprime histograms and log
ledgers, and distinct checks of every central threshold and gcd stratum in
the nine aggregate batches. Those numerical premises are preserved in this
revision; the sharper global bound has its own reconstruction chain below.

~~~text
python -B scripts/check_complete_evidence.py
python -B -O scripts/check_complete_evidence.py
python -B scripts/reconstruct_all.py --work full_normal
python -B scripts/reconstruct_all.py --work full_optimized --optimized
~~~

The first two commands validate the distributed evidence and recompute all
semiprime log bounds. The last two regenerate the complete computational
chain in fresh work directories. The shorter defaults above retain their
original, narrower scopes.
[Reproduction](REPRODUCE.md) states the exact scopes.

[Reviewer path](README_REVIEWER.md) · [Claim map](docs/CLAIM_MAP.md) ·
[Sources](docs/SOURCES.md) · [Limits](docs/KNOWN_LIMITS.md) ·
[Provenance](docs/PROVENANCE.md)

## Sharper global bound

The new self-contained companion is in
[companion/global_bound](companion/global_bound/README.md). Outputs must be new;
the following commands deliberately write outside the distributed tree.

~~~text
python -B companion/global_bound/producer.py --build --output ../new-global-certificate.json
python -B companion/global_bound/verify_closures.py --output ../new-closure-check
python -B companion/global_bound/check.py --certificate ../new-global-certificate.json --closures-report ../new-closure-check/report.json --output ../new-global-check.json
~~~

The producer and distinct scalar checker evaluate all 621 support sectors and
the remaining cardinality tail. Retained closure checks verify their finite
premises and bind their complete-domain report to the scalar check. Without
that report, the scalar check explicitly remains conditional on those
closure premises. Full C++ reconstruction commands and exact scopes are in
[REPRODUCE.md](REPRODUCE.md).
The printed proof supplies the infinite transfer and tail inequalities.
The earlier 37-stage two-prime reconstruction remains a separate chain.

## Scope and attribution

The manuscript supplies the infinite estimates and exact finite proof
dependencies. Classical lattice, CRT, Möbius, Farey and moment tools are credited.
No optimal constant, complete mixed law, effective convergence rate,
three-prime sign classification, universal FCSR convention or Collatz theorem
is claimed. Chen–Winterhof source comparisons use the accessible preprint;
the final journal body was not read. Priority is not certified.

[AI_USE.md](AI_USE.md) discloses substantive assistance and shared review exposure.

Version **0.2.0**. [Changes](CHANGELOG.md).
Code, author-generated certificates and documentation: [MIT](LICENSE).
Original manuscript source/PDF: [CC BY 4.0](LICENSE_MANUSCRIPT.md).
[License scopes](LICENSE_SCOPE.md) · [Third-party notices](THIRD_PARTY_NOTICES.md).
