# Marked-70 family certificates

This component proves the uniform four-prime family bound for radical `70p`,
every prime `p>=11`, every positive valuation and every length `N>=1`.
It also provides the marked-70 anchor profiles and the `770ell` insertion
tail used by the article. The physical correction is always
`Delta(B,N)=mu(B^N)/log(B)`, with signed numerator, positive full-unit
denominator, infinity norm and all states including zero. It does not
identify the global optimizer.

Copy this entire directory anywhere. Python 3.11 or later and its standard
library suffice for the compact check:

```
python -B verify.py --quick --output-dir runs/quick
python -B -O verify.py --quick --output-dir runs/optimized
```

Run from this directory and choose new output directories. Each command
checks the complete public manifest, all logarithm inputs, all 18 exact
stored profile integrations, complete HNF branch enumeration, repair and
insertion predicates, infinite gates, the finite domain and the current
closure assembly. The additional geometric audit checks 96 stratified
witness samples per profile, including every angular node and cell vertex
for those samples. It does not check every geometric witness or recreate
the 663 physical finite streams. These narrower scopes appear in its JSON
reports. `RECONSTRUCTION_HISTORY.json` records the separately observed full
reconstruction of all primary records on 10 October 2026.

For the complete route on a little-endian GNU/Linux host with GNU C++17
and OpenMP:

```
python3 -B verify.py --rebuild --output-dir runs/full
```

The public generator and independently compiled consumer produce and
check every geometric witness in all 18 profiles, compare every complete
cell output and full witness digest, rebuild the tensor algebra and
perform the separate Python audits. They then generate and consume every
finite physical stream, compare all finite bounds and stream digests,
and execute six rejection checks. Jobs run sequentially with four threads;
each inner stage has a 600-second timeout. The enclosing geometric route
has a 3,600-second timeout. Failure logs and generated streams remain
available under the new output directory, including failed attempts. The raw geometric and finite streams occupy
2,195,472,384 and 540,868,424 bytes respectively, before auxiliary files.
No network or external solver is required. UBSan is enabled for the finite
physical runner; the geometric runner uses optimized C++17/OpenMP.

The manifest covers compact scientific inputs, code and documents.
Generated `runs/`, `bin/`, `work/` and `logs/` are local execution material,
not part of that payload. The supplied streams are hashes and numerical
summaries; the complete generator recreates raw streams rather than
publishing gigabytes of transient output.

The transfer and continuous-coverage arguments are in `profiles/PROOF.md`
and `profiles/imports/DIRECT_METHOD_SOURCE.tex`, with the article's
computational appendix. Repair and insertion proofs are in `checks/`;
the current exact closure is `results/closure_certificate.json`.
Functional `review/check_final_closure.py` checks that closure by exact
regeneration; its directory name does not assert human review.
Numerical checks depend on the stated analytical lemmas and are not a
formal proof-assistant verification.

Code, numerical data and companion documentation use the repository MIT
license. The article source and PDF use CC BY 4.0, as specified by the
repository `LICENSE_SCOPE.md` and `LICENSE_MANUSCRIPT.md`. Authorship is
given by the repository `CITATION.cff`.
