# Reproduction

Python >=3.10 and the standard library suffice for the original two-prime
chain and the scalar checks. The marked-70 component, including the new global
retained chain, requires Python >=3.11. Full geometric reconstruction of the new
global bound additionally requires the GNU C++ environment described below.
No Git, network, internal checkout or external Python package is required.
Run from the repository root.

## Complete distributed evidence

~~~text
python -B scripts/check_manifest.py
python -B scripts/check_complete_evidence.py
python -B -O scripts/check_complete_evidence.py
~~~

The evidence checker verifies stored and uncompressed identities, exact ordered
domains, all reported central coverage and source bindings. It recomputes the
distinct log enclosures for all 4,596 reconstructed semiprime histograms and
their 705,257 consumed ledger entries, using 6,929 distinct integers. It also
rechecks the retained arithmetic and complete final 4,410-profile closure.

It validates records of central enumeration; it does not rerun that enumeration
or the original producers. Its output states these boundaries explicitly.

## Complete fresh reconstruction

~~~text
python -B scripts/reconstruct_all.py --work full_normal
python -B scripts/reconstruct_all.py --work full_optimized --optimized
~~~

Each command requires a new simple work name. Existing work directories,
symlinks and escaped paths are refused. Sources, adapters, inputs and the
reconstruction plan are checked before scientific use. Generated outputs feed
the subsequent consumers in the same workspace; the inventory certificate is
derived from the two newly generated inventories. Duplicate JSON keys,
nonfinite values, changed types and checkpoint byte changes are rejected.

The full route runs 37 ordered producer/checker stages, exports every semiprime
histogram and ledger, recomputes the distinct log checks, and distinctly
enumerates all aggregate and pointwise central prefixes. Only the full route
reports complete reconstruction. Optional --stage reconstruction or --stage
central reports partial work and still regenerates the prerequisites.

All constituent scientific stages were actually executed, normally and with
optimization, with identical scientific output bytes. The reader driver
packages the qualified orchestration and is separately tested with control
fixtures and a fresh real-source prefix; the entire new outer command was
not additionally rerun after assembling it. The observed complete dependency
stages took about 18 minutes per mode on the local machine, plus about six
minutes for the nine distinct aggregate batches and about one minute for the
pointwise central checks. These are execution observations, not portable
runtime or RAM guarantees. Gzip byte identity across Python versions is not
promised.

## Exact observed coverage

- 37 dependency stages per mode, including prime-power bridges, mixed finite
  cases, all semiprime envelopes, analytical bounds, two distinct inventory
  methods, branch calendars/aggregates, checkers and physical controls.
- 4,596 semiprime histograms and ledgers, representing 41,585,055 upper-envelope
  entries. Every reconstructed outward logarithmic sum is strictly negative.
- 4,919 distinct aggregate profiles, all 883,401 central thresholds and
  17,931,560 proper-stratum threshold cells, including forced zeros.
- 135 separately counted validation comparisons on 103 distinct moduli,
  with 2,191 central thresholds and 19,349 proper-stratum threshold cells.
  Validation domains can overlap primary domains; do not add them as new moduli.
- Complete semiprime/base-six/sixth-power pointwise central prefixes:
  407,316 / 2,187 / 5,472 thresholds.
- 355 analytic + 99 finite square profiles; 6,972 analytic + 4,410 finite
  length-one profiles. The last union is 1,087 + 1,422 + 1,901 with no remainder.
- 25,073,766 signed physical-control types reexecuted across the three branches.

The distinct central algorithm visits all positive a,d <= floor(sqrt(M)),
requires gcd(d,M)=1, inserts both signs and retains the smallest max(a,d).
It includes nonprimitive pairs, treats zero separately, and checks every
threshold and proper gcd stratum. It imports neither the producer arithmetic
nor the interval oracle. Far upper distributions remain justified by the
written domination argument; exact true far-minimum vectors are not required.

## Short baseline and tests

~~~text
python -B scripts/reproduce.py smoke
python -B scripts/reproduce.py smoke --optimized
python -B scripts/reproduce.py retained
python -B scripts/reproduce.py retained --optimized
python -B scripts/reproduce.py semiprime
python -B scripts/reproduce.py semiprime --optimized
python -B -m unittest discover -s tests -v
python -B -O -m unittest discover -s tests -v
~~~

Smoke covers 1,326 defining-minimum states and the counterexample at 900.
Retained checks the original certificates and partitions. Semiprime reconstructs
four cases by default; --all expands its mathematical domain but printing to
stdout alone is not the persisted full evidence export used by the full driver.
The original scope flags in these outputs remain accurate for those commands.
[expected_replay.json](companion/expected_replay.json) records their unchanged
canonical identities. An actual fixed source/oracle qualifier is available:

~~~text
python -B scripts/reproduce.py replay --batch 39 --script qualify_oracle.py --work oracle_normal
python -B scripts/reproduce.py replay --batch 39 --script qualify_oracle.py --work oracle_optimized --optimized
~~~

Generated build and replay_work directories are excluded from the manifest.
The distributed tree preserves exact bytes through the declared Git filters.
Large full reconstructions are optional and are not part of default hosted CI.
Directory/file symlink creation tests can be skipped on Windows without the
required privilege. The guards remain enabled.

## Build the paper

With an existing pdflatex installation:

~~~text
python -B scripts/build_pdf.py
~~~

Set PDFLATEX if needed. Shell escape and MiKTeX package installation are disabled.
The distributed PDF was compiled and rendered for visual inspection.
Recompilation can change PDF metadata/fonts; the manifest binds the distributed
artifact. No TeX packages are bundled.

These checks do not establish originality, human specialist review, formal
proof verification, universal FCSR equivalence or a Collatz theorem.

## Version 0.2.0 global correction

All defaults under `companion/global_bound` are relative and may be relocated.
Write outputs outside the distributed tree, or in a fresh ignored work directory.
No command silently accepts an old output as a fresh reconstruction.

~~~text
python -B companion/global_bound/producer.py --build --output ../new-global-certificate.json
python -B companion/global_bound/verify_closures.py --output ../new-retained-closures
python -B -O companion/global_bound/verify_closures.py --output ../new-retained-closures-optimized
python -B companion/global_bound/check.py --certificate ../new-global-certificate.json --closures-report ../new-retained-closures/report.json --output ../new-global-check.json
python -B -O companion/global_bound/check.py --certificate ../new-global-certificate.json --closures-report ../new-retained-closures-optimized/report.json --output ../new-global-check-optimized.json
~~~

The scalar checker recomputes all negative profiles and infinite profile
tails, all 621 sector integrals and the scalar/cardinality tail.
The premise checker also recomputes the finite incidence period, covering
predicates, shear comparisons and candidate minimum histogram.
Retained closure checking validates the other numerical finite premises,
exact finite-domain partitions and support joins.
The scalar check authenticates the complete retained report against its
component domain and the distributed closure inventory. Omitting
`--closures-report` produces an explicitly conditional scalar result.

For full geometric reconstruction, use a little-endian Linux environment
with GNU g++ C++17, `__int128` and OpenMP. No third-party Python library is
required. This longer path produces and consumes the complete witness
streams and finite residuals, rather than inferring their cost predicates
from retained profile maxima or samples.

~~~text
python3 -B companion/global_bound/reconstruct_closures.py --output ../new-full-closures --threads 4
~~~

Consult the global companion README for individual early closure routes,
thread limits, output refusal rules and the exact current execution record.
The new chain does not rerun the unchanged original two-prime companion.
Its older complete route remains `scripts/reconstruct_all.py`.
