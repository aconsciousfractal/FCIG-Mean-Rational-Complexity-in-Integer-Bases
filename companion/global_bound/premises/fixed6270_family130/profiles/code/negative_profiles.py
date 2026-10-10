"""Execute rejection gates and verify accepted records divisible by unmarked7/13."""
from pathlib import Path
from math import gcd
import argparse, tempfile, subprocess, os, shutil, json, struct, hashlib
from fractions import Fraction as F
from exact import need
from assemble import grid
from core import mark, CORE_PRIMES
ROOT = Path(__file__).resolve().parents[1]
EXPECTED_CASES = 15

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--consumer', type=Path, required=True)
    ap.add_argument('--positive-prefix', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--log', type=Path, required=True)
    a = ap.parse_args()
    env = {**os.environ, 'CSTAR_THREADS': '1'}
    tmp = Path(tempfile.mkdtemp(prefix='p72-6270-negative-'))
    rows = []
    logs = []
    q, nl = mark(1)
    expected = 2688 * nl * 8
    try:
        cases = [('truncated_stream', -8, 1, 1, 'truncated witness file'), ('extra_record', 8, 1, 1, 'extra witness data'), ('nonunit_zero_coefficients', 0, 1, 1, 'full mark unit'), ('unsupported_mark_selector', 0, 2, 1, 'fixed6270 mark selector'), ('unsupported_mixed_index', 0, 1, 3, 'marked6270 index multiplier')]
        for name, delta, selector, mult, token in cases:
            pref = tmp / name
            with Path(str(pref) + '.witness').open('wb') as f:
                f.truncate(expected + delta)
            a.consumer.chmod(448)
            z = subprocess.run([str(a.consumer), 'consume', str(selector), str(mult), str(ROOT / 'data'), str(pref)], env=env, text=True, capture_output=True)
            need(z.returncode != 0 and token in z.stderr, 'reject ' + name + ' with specific reason')
            rows.append({'name': name, 'status': 'REJECTED', 'expected_diagnostic': token, 'nonzero_exit': True})
            logs.append(name + '\n' + z.stdout + z.stderr)
            Path(str(pref) + '.witness').unlink()
        ps = CORE_PRIMES
        branches = [list(map(int, line.split())) for line in Path(str(a.positive_prefix) + '.branches').read_text().splitlines()]
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
                need(gcd(ci * r0 + cj * r1, q) == bad, 'exactly designated prime rejects every fixture record')
                checked += 1
                block.extend(struct.pack('<ii', ci, cj))
            with Path(str(pref) + '.witness').open('wb') as f:
                for _ in range(2688):
                    f.write(block)
            need(Path(str(pref) + '.witness').stat().st_size == expected, 'exact synthetic stream length')
            a.consumer.chmod(448)
            z = subprocess.run([str(a.consumer), 'consume', '1', '1', str(ROOT / 'data'), str(pref)], env=env, text=True, capture_output=True)
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
        squared = [max((8 + 3 * n) ** 2, (-n + 24) ** 2) * F(1, 64 + n * n) for n in range(-8, 9)]
        need(max(squared) < 10, 'finite angle nodes miss the continuous maximum')
        rows.append({'name': 'angular_penalty_necessity', 'status': 'PASS', 'vector': [1, 3], 'sampled_squared_max': str(max(squared)), 'continuous_squared_max': '10', 'node_tangents': 'n/8 for -8<=n<=8'})
        need(F(1499, 1000) * F(2, 3) < 1, 'linear primitive-CDF counterexample')
        rows.append({'name': 'cdf_below_three_halves_fails', 'status': 'REJECTED', 'example_M': 3, 'x': '2/3', 'probability': '1', 'constant_rejected': '1499/1000'})
        witness = Path(str(a.positive_prefix) + '.witness')
        need(witness.stat().st_size == expected, 'real positive stream length')
        digest = hashlib.sha256()
        with witness.open('rb') as f:
            for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
                digest.update(block)
        for unmarked in (7, 13):
            need(gcd(unmarked, q) == 1, 'accepted fixture prime is outside full mark')
            examples = []
            with witness.open('rb') as f:
                for i in range(2688 * nl):
                    rec = f.read(8)
                    need(len(rec) == 8, 'positive record available')
                    ci, cj = struct.unpack('<ii', rec)
                    h = i % nl
                    r0, r1 = branches[h][4:]
                    val = ci * r0 + cj * r1
                    if val and val % unmarked == 0 and (val % q % unmarked == 0) and (gcd(val, q) == 1):
                        examples.append({'record_index': i, 'cell': i // nl, 'branch': h, 'ci': ci, 'cj': cj, 'mark_linear_form_integer': val, 'reduced_mark_linear_form': val % q})
                        if len(examples) == 32:
                            break
            need(len(examples) == 32, 'actual original records with unmarked prime accepted by completed full consumer')
            rows.append({'name': f'unmarked_prime{unmarked}_accepted_in_primary_stream', 'status': 'PASS', 'unmarked_prime': unmarked, 'actual_records': examples, 'witness_sha256': digest.hexdigest(), 'acceptance_evidence': 'same whole stream is consumed successfully by the primary consumer; these32records are individually checked here'})
        need(len(rows) == EXPECTED_CASES, 'all declared fixture gates')
        a.output.write_text(json.dumps({'status': 'PASS', 'cases': rows, 'isolated_synthetic_files': True, 'positive_geometry_gate': 'all selected original streams are separately produced and consumed in full'}, indent=2) + '\n', newline='\n')
        a.log.write_text('\n'.join(logs), newline='\n')
        print('PASS negative profile gates', len(rows), flush=True)
    finally:
        shutil.rmtree(tmp)
if __name__ == '__main__':
    main()
