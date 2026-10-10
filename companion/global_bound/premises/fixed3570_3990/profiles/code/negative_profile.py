"""Isolated adversarial streams and malformed continuous-domain records."""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, tempfile, subprocess, os, shutil
from integrate import integrate
from core import mark
from exact import need
R = Path(__file__).resolve().parents[1]

def run(consumer, p, m, style, ag, prefix):
    expected = 2667 * mark(p)[1] * m * 8
    cases = []
    logs = []
    env = {**os.environ, 'CSTAR_THREADS': '2'}
    with tempfile.TemporaryDirectory(prefix='p72-3570-3990-negative-') as tmp:
        tmp = Path(tmp)
        for name, delta, token in [('truncated', -8, 'truncated witness file'), ('extra', 8, 'extra witness data'), ('nonunit', 0, 'full mark unit')]:
            pre = tmp / name
            with Path(str(pre) + '.witness').open('wb') as f:
                f.truncate(expected + delta)
            z = subprocess.run([str(consumer), 'consume', str(p), style, str(R / 'data'), str(pre)], env=env, text=True, capture_output=True)
            need(z.returncode and token in z.stderr, 'specific rejection ' + name)
            cases.append(dict(name=name, status='REJECTED', diagnostic=token))
            logs.append(z.stdout + z.stderr)
        z = subprocess.run([str(consumer), 'consume', '15', style, str(R / 'data'), str(tmp / 'nonunit')], env=env, text=True, capture_output=True)
        need(z.returncode and 'declared new-prime domain' in z.stderr, 'unsupported parameter rejected')
        cases.append(dict(name='unsupported_prime_parameter', status='REJECTED', diagnostic='declared new-prime domain'))
        logs.append(z.stdout + z.stderr)
        z = subprocess.run([str(consumer), 'consume', str(p), 'hybrid' if style in ('coarse', 'angular_tail') else 'coarse', str(R / 'data'), str(tmp / 'nonunit')], env=env, text=True, capture_output=True)
        need(z.returncode and 'style and index contract' in z.stderr, 'wrong grid-index pair rejected')
        cases.append(dict(name='grid_index_mismatch', status='REJECTED', diagnostic='style and index contract'))
        logs.append(z.stdout + z.stderr)
        d = tmp / 'badcrop'
        d.mkdir()
        shutil.copy2(R / 'data/log_table.txt', d / 'log_table.txt')
        lines = (R / 'data/grid.tsv').read_text().splitlines()
        lines[0] = str(int(lines[0]) + 1)
        (d / 'grid.tsv').write_text('\n'.join(lines) + '\n', newline='\n')
        z = subprocess.run([str(consumer), 'consume', str(p), style, str(d), str(tmp / 'nonunit')], env=env, text=True, capture_output=True)
        need(z.returncode and 'exact Gauss ceiling' in z.stderr, 'crop domain rejected')
        cases.append(dict(name='wrong_gauss_scale', status='REJECTED', diagnostic='exact Gauss ceiling'))
        logs.append(z.stdout + z.stderr)
        lines = Path(str(prefix) + '.cells').read_text().splitlines()
        row = lines[0].split()
        row[3] = str(int(row[3]) + 1)
        lines[0] = ' '.join(row)
        bad = tmp / 'gap.cells'
        bad.write_text('\n'.join(lines) + '\n', newline='\n')
        try:
            integrate(p, m, style, bad, R, ag)
        except ArithmeticError as e:
            need('exact profile geometry' in str(e), 'specific radial hole rejection')
        else:
            raise ArithmeticError('radial hole accepted')
        cases.append(dict(name='radial_gap', status='REJECTED', diagnostic='exact profile geometry'))
        for invalid_m, invalid_style in ((2, 'hybrid'), (3, 'coarse'), (3, 'hybrid'), (6, 'coarse'), (1, 'angular_tail'), (3, 'angular_tail'), (6, 'angular_tail')):
            try:
                integrate(p, invalid_m, invalid_style, str(prefix) + '.cells', R, ag)
            except ArithmeticError as e:
                need('declared index and grid' in str(e), 'specific unsupported pair rejection')
            else:
                raise ArithmeticError('unsupported index-grid pair accepted')
            cases.append(dict(name=f'unsupported_index_grid_{invalid_m}_{invalid_style}', status='REJECTED', diagnostic='declared index and grid'))
        squared = [max((4 + 3 * n) ** 2, (-n + 12) ** 2) * F(1, 16 + n * n) for n in range(-4, 5)]
        need(max(squared) < 10, 'finite nodes omit a continuous angular maximum')
        cases.append(dict(name='angular_interpolation_necessity', status='PASS', vector=[1, 3], sampled_squared_max=str(max(squared)), continuous_squared_max='10'))
    return (dict(status='PASS', q=210 * p, transfer_index=210 * p * m, cases=cases, synthetic_files_isolated=True), '\n'.join(logs))
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--consumer', type=Path, required=True)
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--multiplier', type=int, required=True)
    ap.add_argument('--style', required=True)
    ap.add_argument('--angular', type=int, default=8)
    ap.add_argument('--prefix', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--log', type=Path, required=True)
    a = ap.parse_args()
    r, log = run(a.consumer, a.p, a.multiplier, a.style, a.angular, a.prefix)
    a.output.write_text(json.dumps(r, indent=2) + '\n', newline='\n')
    a.log.write_text(log, newline='\n')
    print(json.dumps(r, indent=2))
