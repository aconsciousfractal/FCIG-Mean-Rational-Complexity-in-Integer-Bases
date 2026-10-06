# Reviewer path

## First reading

Read the [definition and theorems](paper/Mean-Rational-Complexity-in-Integer-Bases.pdf),
[claim map](docs/CLAIM_MAP.md), [sources](docs/SOURCES.md) and [limits](docs/KNOWN_LIMITS.md).
Keep all-residue weight, positive full-unit denominator, max-height and unrounded
logarithm fixed. A connection denominator constrained to -1 modulo the base is
a different convention.

## Proof checks

Check the reduced lattice basis, unit-wheel completion and shortest-vector
incidence bound. The numerical 49/50 argument takes common base-independent
caps before the supremum. The prime-power and two-prime arguments do not use
that strengthening, so the dependency order is acyclic.
Corollary 7.2 separately imports Kanold's bound as reported by
Costello–Watts v2, Introduction pp. 1–2. Check the equality of the cyclic-gap
and consecutive-interval conventions and the common cap before taking the
suprema. The resulting $4/5$ bound does not replace the self-contained proof
or determine the optimum.

For exact central counts, inspect primitive cancellation, state-gcd layers,
same-sign injection, opposite-sign closed collision phases and the separate
zero. Signed inverse witnesses have numerator ±g; a missing-prime repair
reverses the sign. The child modulus2 is one atom.

The early base-six/sixth-power/squarefree histograms use a pointwise envelope,
explained in the paper separately from the later aggregate CDF calendar.
The negative outward-log sign of each finite histogram is a premise to check,
not a consequence of domination alone. Wider early denominator caps and later
refined caps are both justified and distinguished.

Check the two-prime Farey difference witness, positive repaired denominators,
open hole endpoints, one integer per hole and strict far threshold.
Then check all central/middle/far compensation errors and the radical-six gain.
The five analytical bounds pair each error with its own simultaneous gains.

The coverage order is: large-modulus estimates and sixth powers; four complete
radical families; fourth powers/cubes/squares; then the length-one inventory.
Its 19,239 profiles split 7,857 ratio + 6,972 analytic + 4,410 finite; the final branches
contain 1,087 + 1,422 + 1,901 profiles. Reindex at the same physical modulus only.

For the lower/positive results, verify the primitive union upper, logarithmic
integral, exact B30,N2 product and radical510510 constants.
The optimal-correction argument proves finite bases, not finite-length
attainment. Non-perfect-power does not mean squarefree.

## Reproduction and identity

Follow [REPRODUCE.md](REPRODUCE.md). The complete source reconstruction
has passed in both modes. The new evidence includes all 4,596 semiprime
histograms and consumed ledgers, 4,919 aggregate profiles, and every central
threshold and proper gcd stratum. The 135 validation comparisons are counted
separately. A distinct signed legal-pair enumeration checks the entire central
range without importing the producer arithmetic. It also checks the complete
semiprime, base-six and sixth-power central prefixes.

The complete-evidence command checks identities, domains, coverage and distinct
log bounds; it does not repeat the central enumeration. Use the full fresh
reconstruction command for that. Analytic arguments and far-envelope domination
remain written proof dependencies. No true far-minimum vector is claimed for
aggregate upper distributions.

[MANIFEST_SHA256.txt](MANIFEST_SHA256.txt) binds the reader tree.
The manifest checker works without Git and ignores generated build/work trees.
The replay driver refuses reused work directories before scientific launch,
and verifies public stored/uncompressed identities during materialization.

Public sources and retained metadata are derived by documented relative-path,
comment and digest-binding transformations; provenance identifies original
and exported identities. Such transformations do not imply historical jobs
were run with the exported sources. [AI-use disclosure](AI_USE.md).
