# Early closure certificates

These eight cases supply the finite and profile premises for the three-prime
sector and the family rad(B)=30p. Mathematical all-length and arbitrary-valuation
arguments are in the main manuscript. The data retain upper bounds, rather than
claiming exact minima unless the case explicitly provides them.

`verify_early.py` checks the retained inputs, exact joins and domain partitions.
Each case also includes its full comparison reconstruction driver `verify_all.py`.
Those longer routes require GNU g++ (C++17, __int128, OpenMP where used) and run on
a little-endian Linux system. Invoke them with Python -B, optionally -O, and
`--sanitize`; use the individual driver's --help for retained-work arguments.
Do not use --generate when verifying distributed evidence.

The numerical kernels and inputs were fully reconstructed on 8 October 2026.
The present revision preserves those numerical identities. It omits historical
working notes and redirects their theorem statements to the main manuscript.
Runtime metadata and an imported-parent flag are not mathematical certificates:
each parent's separate case and the joint checker remain necessary.

MIT licensing follows the repository's LICENSE_SCOPE.md. The main manuscript
contains the mathematical narrative under its separate manuscript license.
