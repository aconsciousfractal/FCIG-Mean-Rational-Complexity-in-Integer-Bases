# Repair and insertion predicates

The four `verify_*.py` scripts regenerate exact repair, gate and sharpened
insertion certificates. `imports/q70_k1.json`, `q70_k2.json` and `q770.json`
are byte-identical copies of the corresponding profile results and are
joined by the closure checker. The mathematical proofs are adjacent.

The scripts' finite model regressions are separate from their analytical
uniform conclusions. The top-level `verify.py` runs all four scripts and
the full family closure assembly. Their regenerated certificates are
checked against the public manifest after the run. These repair checks
do not alone close small-prime profiles or the finite physical complement.

Code, numerical certificates and companion proof documents use the MIT
license of the enclosing repository.
