"""Concrete rejection and necessity gates for the explicit mark (3,5,11,p)."""
from pathlib import Path
from math import gcd, prod
import argparse, tempfile, subprocess, os, shutil, json, struct, hashlib
from fractions import Fraction as F
from exact import need
from assemble import grid
from core import mark, MIN_NEW_PRIME, CORE_PRIMES
ROOT = Path(__file__).resolve().parents[1]
EXPECTED_CASES = 16

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--consumer', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--log', type=Path, required=True)
    a = ap.parse_args()
    env = {**os.environ, 'CSTAR_THREADS': '1'}
    tmp = Path(tempfile.mkdtemp(prefix='p72-165-negative-'))
    rows = []
    logs = []
    q, nl = mark(MIN_NEW_PRIME)
    expected = 2688 * nl * 8
    try:
        cases = [('truncated_stream', -8, MIN_NEW_PRIME, 1, 'truncated witness file'), ('extra_record', 8, MIN_NEW_PRIME, 1, 'extra witness data'), ('nonunit_zero_coefficients', 0, MIN_NEW_PRIME, 1, 'full mark unit'), ('composite_parameter', 0, 25, 1, 'prime input')] + [(f'forbidden_core_prime_{pr}', 0, pr, 1, 'prime parameter range') for pr in CORE_PRIMES] + [('unsupported_core_index', 0, 1, 2, 'marked-165 divisor index multiplier')]
        for name, delta, prime, mult, token in cases:
            pref = tmp / name
            with Path(str(pref) + '.witness').open('wb') as f:
                f.truncate(expected + delta)
            z = subprocess.run([str(a.consumer), 'consume', str(prime), str(mult), str(ROOT / 'data'), str(pref)], env=env, text=True, capture_output=True)
            need(z.returncode != 0 and token in z.stderr, 'reject ' + name + ' with specific reason')
            rows.append({'name': name, 'status': 'REJECTED', 'expected_diagnostic': token, 'nonzero_exit': True})
            logs.append(name + '\n' + z.stdout + z.stderr)
            Path(str(pref) + '.witness').unlink()
        ps = CORE_PRIMES + (MIN_NEW_PRIME,)
        branches = [list(map(int, line.split())) for line in (ROOT / 'results' / f'q{q}.branches').read_text().splitlines()]
        need(len(branches) == nl, 'full-mark fixture branch enumeration')

        def crt(xs):
            return sum((x * (q // p) * pow(q // p, -1, p) for x, p in zip(xs, ps))) % q
        for bad in ps:
            name = f'only_prime_{bad}_nonunit'
            pref = tmp / name
            block = bytearray()
            checked = 0
            for _, aa, bb, dd, r0, r1 in branches:
                aa_local = []
                bb_local = []
                for pr in ps:
                    target = 0 if pr == bad else 1
                    if r0 % pr:
                        aa_local.append(target * pow(r0, -1, pr) % pr)
                        bb_local.append(0)
                    else:
                        need(r1 % pr != 0, 'primitive local fixture image')
                        aa_local.append(0)
                        bb_local.append(target * pow(r1, -1, pr) % pr)
                ci, cj = (crt(aa_local), crt(bb_local))
                need(gcd(ci * r0 + cj * r1, q) == bad, 'only designated prime rejects every fixture record')
                checked += 1
                block.extend(struct.pack('<ii', ci, cj))
            with Path(str(pref) + '.witness').open('wb') as f:
                for _ in range(2688):
                    f.write(block)
            need(Path(str(pref) + '.witness').stat().st_size == expected, 'exact synthetic marked stream length')
            z = subprocess.run([str(a.consumer), 'consume', str(MIN_NEW_PRIME), '1', str(ROOT / 'data'), str(pref)], env=env, text=True, capture_output=True)
            need(z.returncode != 0 and 'full mark unit' in z.stderr, 'specific missing-prime denominator rejection')
            rows.append({'name': name, 'status': 'REJECTED', 'expected_diagnostic': 'full mark unit', 'nonzero_exit': True, 'marked_primes': list(ps), 'unique_nonunit_prime': bad, 'coefficient_pairs_checked': checked, 'block_sha256': hashlib.sha256(block).hexdigest(), 'stream_records': 2688 * nl})
            logs.append(name + '\n' + z.stdout + z.stderr)
            Path(str(pref) + '.witness').unlink()
        gr = tmp / 'gapped'
        (gr / 'data').mkdir(parents=True)
        raw = (ROOT / 'data/grid.tsv').read_text().splitlines()
        rs = [list(map(int, s.split())) for s in raw[1:]]
        for r in rs[:21]:
            r[2] -= 1
        (gr / 'data/grid.tsv').write_text(raw[0] + '\n' + ''.join((' '.join(map(str, r)) + '\n' for r in rs)), newline='\n')
        try:
            grid(gr)
        except ArithmeticError as exc:
            need('radial continuity' in str(exc), 'specific radial gap failure')
            rows.append({'name': 'radial_gap', 'status': 'REJECTED', 'expected_diagnostic': 'radial continuity'})
            logs.append('radial_gap\n' + str(exc) + '\n')
        else:
            raise ArithmeticError('gapped radial grid accepted')
        squared = [max((4 + 3 * n) ** 2, (-n + 12) ** 2) * F(1, 16 + n * n) for n in range(-4, 5)]
        need(max(squared) < 10, 'finite angle nodes miss the continuous maximum')
        rows.append({'name': 'angular_penalty_necessity', 'status': 'PASS', 'vector': [1, 3], 'sampled_squared_max': str(max(squared)), 'continuous_squared_max': '10', 'node_tangents': 'n/4 for -4<=n<=4'})
        need(F(1) / F(2, 3) == F(3, 2), 'linear primitive-CDF constant is sharp')
        rows.append({'name': 'cdf_below_three_halves_fails', 'status': 'REJECTED', 'example_M': 3, 'x': '2/3', 'probability': '1', 'constant_rejected': '1499/1000'})
        positive = []
        samples = [list(map(int, line.split())) for line in (ROOT / 'results' / f'q{q}.audit').read_text().splitlines()]
        for sample in samples:
            cell, h, ci, cj, u, v = sample[:6]
            _, aa, bb, dd, r0, r1 = branches[h]
            denominator_mark = ci * r0 + cj * r1
            if denominator_mark % 2 == 0:
                need(gcd(denominator_mark, q) == 1 and (u, v) == (ci * aa + cj * bb, cj * dd), 'actual even denominator mark admitted by the odd-core consumer')
                positive.append(dict(cell=cell, branch=h, ci=ci, cj=cj, denominator_mark=denominator_mark))
        need(len(positive) >= 8, 'positive allowed2 fixture has concrete records from the fully consumed first small profile')
        rows.append({'name': 'allowed_prime2_actual_records', 'status': 'PASS', 'prime_allowed': 2, 'q': q, 'records': positive, 'scope': 'actual sampled coefficient records from the complete producer/consumer stream; 2 is not forbidden by the odd mark'})
        need(F(1499, 1000) * F(2, 3) < 1, 'explicit CDF counterexample')
        need(len(rows) == EXPECTED_CASES, 'all declared fixture gates')
        result = {'status': 'PASS', 'cases': rows, 'isolated_synthetic_files': True, 'positive_geometry_gate': 'all selected profiles are separately produced and consumed in full'}
        a.output.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
        a.log.write_text('\n'.join(logs), newline='\n')
        print('PASS negative profile gates', len(rows), flush=True)
    finally:
        shutil.rmtree(tmp)
if __name__ == '__main__':
    main()
