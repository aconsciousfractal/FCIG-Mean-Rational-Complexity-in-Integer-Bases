# Global bound companion

This directory computes the supportwise bound C* < 0.241933621228105 for every integer base B >= 2 and every N >= 1. Python 3.10 or later and the standard library suffice for the producer and scalar checker. The complete retained closure chain includes the marked70 carrier and requires Python 3.11 or later. Copy the entire directory anywhere; all default inputs are relative to these scripts.

Run the commands below from the repository root. All outputs are outside the entire distributed repository. For a standalone copy, replace `companion/global_bound/` with the path to that copy and keep the outputs outside it.

```
python -B companion/global_bound/producer.py --build --output ../new-certificate.json
python -B companion/global_bound/check.py --certificate ../new-certificate.json --output ../new-check.json
python -B -O companion/global_bound/check.py --certificate ../new-certificate.json --output ../new-optimized-check.json
python -B companion/global_bound/verify_closures.py --output ../new-closure-checks
python -B companion/global_bound/check.py --certificate ../new-certificate.json --closures-report ../new-closure-checks/report.json --output ../new-complete-check.json
python -B companion/global_bound/reconstruct_closures.py --output ../new-full-reconstruction --threads 4
```

Every output file or directory must be new. Failures retain their logs. The checker does not import the producer. It reconstructs all 1,079 finite square/cross profiles by independent Mobius counts, includes their infinite tails, then independently evaluates all 621 supportwise rows and the remaining scalar sectors. It reconstructs the primitive period, covering cases, integer gap predicates and the 210 candidate minima.

`verify_closures.py` recomputes exact radial integrals, full HNF algebra, input logarithm tables, insertion predicates, all-valuation family inventories, infinite gates and eleven radical closures from retained numerical inputs. Its geometric cost audits are sampled; retained cell caps are imported finite scientific certificates. This command does not recreate the primary geometric witnesses or the primary shortest-vector streams.

`reconstruct_closures.py` recreates the primary geometric witnesses for the later family and radical components in `premises`, and the 70p component in `marked70`, with a C++ producer and separately compiled consumer. It compares their cells, branches, tensors, integrals and witness digests, and regenerates their finite physical residuals. The early closures are checked from retained outputs by this command; their full reconstruction routes are listed below. GNU C++17, OpenMP and a little-endian GNU/Linux host are required. UBSan is enabled by the drivers that support it. Jobs run sequentially with at most four threads and a 3,600 second limit per stage; timed-out process groups are terminated. Allow several gigabytes of temporary disk space and up to a few hours. Failed work and logs are preserved.

The eight early closures have separate full comparison drivers. Run these sequentially from the repository root on the same GNU/Linux host; omit every development `--generate` option. Their complete reconstructions were observed before this revision with the same numerical kernels, inputs and expected outputs. This revision does not claim a fresh repeat of those expensive computations. The early drivers require the caller to supply fresh report paths; they do not all reject an existing report file.

```
python -B companion/global_bound/early_closures/three_prime_family/verify_all.py --sanitize --keep-work ../new-three-prime-work
python -B companion/global_bound/early_closures/three_prime_28/verify_all.py --sanitize
python -B companion/global_bound/early_closures/three_prime_89/verify_all.py --sanitize --threads 4 --keep-work --report ../new-three-prime-89.json
python -B companion/global_bound/early_closures/three_prime_353/verify_all.py --sanitize --keep-work --report ../new-three-prime-353.json
python -B companion/global_bound/early_closures/three_prime_large/verify_all.py --sanitize --keep-work --report ../new-three-prime-large.json
python -B companion/global_bound/early_closures/radical210/verify_all.py --sanitize
python -B companion/global_bound/early_closures/thirty_p_small/verify_all.py --sanitize --threads 4 --keep-work --report ../new-thirty-p-small.json
python -B companion/global_bound/early_closures/thirty_p_large/verify_all.py --sanitize --threads 4 --keep-work --work-dir ../new-thirty-p-large-work --report ../new-thirty-p-large.json
```

The early closures for at most three primes and the 30p family are in `early_closures`; the reconstructed 70p family is in `marked70`. The numerical premise catalog is `closure_input_catalog.json`; `finite_certificate_index.json` maps full marks and transfer indices to grids, caps and reconstruction components. The mathematical implications of these finite predicates, including all positive valuations and all N, are proved in `proofs/PROOF.md` and the article's computational appendix. Square/cross finite-count and analytic-tail formulas are proved in `proofs/NEGATIVE_PROFILES.md`.

Code and companion documentation are covered by the repository MIT license. The article source and PDF are covered by CC BY 4.0, as stated in LICENSE_MANUSCRIPT.md. Authorship is as stated in CITATION.cff.
