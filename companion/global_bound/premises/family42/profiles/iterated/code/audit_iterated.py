from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json, struct, hashlib
import arithmetic as ar, reference_family as ref, independent as ind, check_routes as cr
from exact import need, log_lo, log_hi, alternate_log, out_decimal
ROOT = Path(__file__).resolve().parents[1]
PREF = ROOT / 'results/q462'

def main():
    global PREF
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--prefix', type=Path, default=PREF)
    ap.add_argument('--reference-dir', type=Path, default=ROOT / 'results')
    ap.add_argument('--output-dir', type=Path, default=ROOT / 'results')
    args = ap.parse_args()
    PREF = args.prefix
    cells, cellmeta = ar.make_cells(ROOT)
    need(ar.encode_cells(cells) == (ROOT / 'data/cells.txt').read_text(), 'all closed-cell condition bounds')
    _, logmeta = ref.log_checks(ROOT / 'data/log_table.txt')
    ind.T = cr.T = 8
    ind.J = cr.J = 32
    algebra_geometry = ind.run_one(ROOT, 11, PREF, True)
    common = ind.common_costs(ROOT, PREF, cells)
    arrays = {}
    logs = set()
    states = 0
    for k in (1, 2):
        M = 462 ** k
        path = args.reference_dir / f'M{M}_q462.min'
        raw = path.read_bytes()
        need(len(raw) == 4 * M, 'complete definition bank')
        R = struct.unpack('<' + 'I' * M, raw)
        need(all((1 <= x <= M // 2 for x in R)), 'height domain')
        need(all((R[s] == R[-s % M] for s in range(M))), 'inversion symmetry including zero')
        arrays[k] = R
        logs.update(R)
        logs.add(M)
        states += M
    hi = {h: log_hi(h) for h in logs}
    lo = {h: log_lo(h) for h in logs}
    for h in sorted(logs):
        a, b = alternate_log(h)
        need(lo[h] <= a <= b <= hi[h], 'every physical integer logarithm independently enclosed')
    physical = []
    points = conditional = 0
    for line in Path(str(PREF) + '.points').read_text().splitlines():
        k, b = map(int, line.split())
        if k == 0:
            need(b >= 0, 'identity observable')
            continue
        M = 462 ** k
        mean = sum((hi[h] for h in arrays[k])) / M - lo[M] / 2
        need(mean <= F(b, ar.S), 'complete root minimum mean versus iterated envelope')
        physical.append({'k': k, 'M': M, 'mean_upper': str(mean), 'profile_root_upper': str(F(b, ar.S))})
    for line in Path(str(PREF) + '.bank').read_text().splitlines():
        k, s, b = map(int, line.split())
        need(0 <= s < 462, 'point bankstate')
        if k == 0:
            mean = hi[arrays[1][s]] - lo[462] / 2
            points += 1
        elif k == 1:
            mean = sum((hi[arrays[2][s + 462 * t]] for t in range(462))) / 462 - lo[462 ** 2] / 2
            conditional += 1
        else:
            raise ArithmeticError('bank horizon')
        need(mean <= F(b, ar.S), 'complete conditional physical mean versus iterated point bound')
    need(points == 462 and conditional == 462, 'every declared physical parent checked')
    fmean = sum((hi[h] for h in arrays[1])) / 462 - lo[462] / 2
    ratio = max(F(0), fmean) / lo[462]
    need(ratio < F(2, 85), 'finite base462 excluded exactly')
    inventory = []

    def rec(j, n, ee):
        if j == 4:
            inventory.append({'B': n, 'exponents': ee})
            return
        p = (2, 3, 7, 11)[j]
        n *= p
        e = 1
        while n < 924:
            rec(j + 1, n, ee + [e])
            n *= p
            e += 1
    rec(0, 1, [])
    scan = []
    for n in range(462, 924, 462):
        z = n
        for p in (2, 3, 7, 11):
            while z % p == 0:
                z //= p
        if z == 1:
            scan.append(n)
    need(sorted((r['B'] for r in inventory)) == scan == [462], 'exhaustive finite base inventory below924 by two methods')
    finite = {'status': 'PASS', 'M': 462, 'B': 462, 'N': 1, 'base_cutoff': 924, 'inventory': inventory, 'exact_minimum_histogram': sorted(Counter(arrays[1]).items()), 'mean_upper': str(fmean), 'mean_upper_decimal': out_decimal(fmean), 'normalized_upper': str(ratio), 'normalized_upper_decimal': out_decimal(ratio), 'states': 462, 'representatives': 232, 'primitive_denominator_enumeration': True, 'independent_log_arguments': len(logs)}
    (args.output_dir / 'finite_base462.json').write_text(json.dumps(finite, indent=2) + '\n', newline='\n')
    result = {'status': 'PASS', 'cell_inputs': cellmeta, 'log_table': logmeta, 'initial_coefficients_and_routes': algebra_geometry, 'common_costs': common, 'physical': {'complete_states': states, 'pointwise_checks': points, 'conditional_checks': conditional, 'rows': physical}, 'independent_integer_logarithms': len(logs), 'second_geometric_audit_exhaustive': False}
    (args.output_dir / 'independent_checks.json').write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps({'status': 'PASS', 'finite_base462_mu_upper': out_decimal(fmean), 'physical_states': states, 'second_geometric_audit_exhaustive': False}), flush=True)
if __name__ == '__main__':
    main()
