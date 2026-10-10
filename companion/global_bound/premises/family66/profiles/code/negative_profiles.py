"""Reject malformed geometric streams and a rational hole in grid coverage."""
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
    ap.add_argument('--output', type=Path, default=ROOT / 'results/negative_checks.json')
    ap.add_argument('--log', type=Path, default=ROOT / 'logs/negative_checks.log')
    a = ap.parse_args()
    env = {**os.environ, 'CSTAR_THREADS': '2'}
    tmp = Path(tempfile.mkdtemp(prefix='p72-66p-negative-'))
    rows = []
    logs = []
    expected = 2688 * mark(MIN_NEW_PRIME)[1] * 8
    try:
        for name, delta, prime, token in [('truncated_stream', -8, MIN_NEW_PRIME, 'truncated witness file'), ('extra_record', 8, MIN_NEW_PRIME, 'extra witness data'), ('nonunit_coefficients', 0, MIN_NEW_PRIME, 'full mark unit'), ('composite_parameter', 0, 15, 'prime input')]:
            pref = tmp / name
            with Path(str(pref) + '.witness').open('wb') as f:
                f.truncate(expected + delta)
            z = subprocess.run([str(a.consumer), 'consume', str(prime), str(ROOT / 'data'), str(pref)], env=env, text=True, capture_output=True)
            need(z.returncode != 0 and token in z.stderr, 'reject ' + name + ' with correct reason')
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
            need('radial continuity' in str(exc), 'specific radial coverage failure')
            rows.append({'name': 'radial_gap', 'status': 'REJECTED', 'expected_diagnostic': 'radial continuity'})
            logs.append('radial_gap\n' + str(exc) + '\n')
        else:
            raise ArithmeticError('gapped radial grid accepted')
        squared = [max((4 + 3 * n) ** 2, (-n + 12) ** 2) * F(1, 16 + n * n) for n in range(-4, 5)]
        need(max(squared) < 10, 'finite angles miss exact continuous maximum for vector(1,3)')
        rows.append({'name': 'angular_penalty_necessity', 'status': 'PASS', 'vector': [1, 3], 'nodes_tangent_n_over_4': list(range(-4, 5)), 'sampled_squared_max': str(max(squared)), 'continuous_squared_max': '10'})
        result = {'status': 'PASS', 'cases': rows, 'synthetic_files_isolated': True, 'expected_positive_consumer_is_exhaustive_profile_q858': True}
        a.output.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
        a.log.write_text('\n'.join(logs), newline='\n')
        print('PASS negative profile checks', len(rows), flush=True)
    finally:
        shutil.rmtree(tmp)
if __name__ == '__main__':
    main()
