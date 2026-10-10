# Complete marked profiles

`profiles_plan.json` lists all 18 profiles: two marked-70 anchors and the
sixteen full marks `70p` for primes from 11 through 71. The grids, angular
nodes, cap integrals, branches and tensor summaries are compact inputs.

From the parent component directory, run `python -B verify.py --quick` for
exact stored reintegration, complete independent HNF algebra and the 96
stratified geometric witnesses per profile. Run `python3 -B verify.py
--rebuild` with GNU C++17/OpenMP for the complete producer/consumer route.
The latter checks all raw geometric records. Sampled audits do not replace
that exhaustive predicate or the analytical continuous-coverage proof.

`PROOF.md` and `imports/DIRECT_METHOD_SOURCE.tex` give the transfer and
continuous-coverage contracts. Code, numerical data and companion proof
documents use the MIT license of the enclosing repository.
