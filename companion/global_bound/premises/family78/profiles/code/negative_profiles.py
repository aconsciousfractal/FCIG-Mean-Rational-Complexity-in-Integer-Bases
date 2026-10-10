"""Exercise concrete rejection gates on isolated synthetic geometric inputs."""
from pathlib import Path
from fractions import Fraction as F
import argparse, tempfile, subprocess, os, shutil, json
from exact import need
from assemble import grid
from core import mark, MIN_NEW_PRIME
ROOT = Path(__file__).resolve().parents[1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--consumer', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--log', type=Path, required=True)
    a = ap.parse_args()
    env = {**os.environ, 'CSTAR_THREADS': '1'}
    tmp = Path(tempfile.mkdtemp(prefix='p72-78-negative-'))
    rows = []
    logs = []
    expected = 2688 * mark(MIN_NEW_PRIME)[1] * 8
    try:
        cases = [('truncated_stream', -8, MIN_NEW_PRIME, 1, 'truncated witness file'), ('extra_record', 8, MIN_NEW_PRIME, 1, 'extra witness data'), ('nonunit_coefficients', 0, MIN_NEW_PRIME, 1, 'full mark unit'), ('composite_parameter', 0, 25, 1, 'prime input'), ('forbidden_new_prime', 0, 13, 1, 'prime parameter range'), ('unsupported_core_index', 0, 1, 3, 'marked-78 index multiplier')]
        for name, delta, prime, mult, token in cases:
            pref = tmp / name
            with Path(str(pref) + '.witness').open('wb') as f:
                f.truncate(expected + delta)
            z = subprocess.run([str(a.consumer), 'consume', str(prime), str(mult), str(ROOT / 'data'), str(pref)], env=env, text=True, capture_output=True)
            need(z.returncode != 0 and token in z.stderr, 'reject ' + name + ' with specific reason')
            rows.append({'name': name, 'status': 'REJECTED', 'expected_diagnostic': token, 'nonzero_exit': True})
            logs.append(name + '\n' + z.stdout + z.stderr)
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
        need(F(1499, 1000) * F(2, 3) < 1, 'explicit CDF counterexample')
        result = {'status': 'PASS', 'cases': rows, 'isolated_synthetic_files': True, 'positive_geometry_gate': 'all nine selected profiles, separately consumed in full'}
        a.output.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
        a.log.write_text('\n'.join(logs), newline='\n')
        print('PASS negative profile gates', len(rows), flush=True)
    finally:
        shutil.rmtree(tmp)
if __name__ == '__main__':
    main()
