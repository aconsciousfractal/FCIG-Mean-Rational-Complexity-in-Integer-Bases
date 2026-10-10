from pathlib import Path
from hashlib import sha256
import argparse, json, tempfile, subprocess, os, shutil
from assemble_refined import R, assemble
from exact import need

def run(consumer, prefix):
    expected = 257789952
    cases = []
    log = []
    env = {**os.environ, 'CSTAR_THREADS': '2'}
    with tempfile.TemporaryDirectory(prefix='p72-2310-refined-negative-') as tmp:
        tmp = Path(tmp)
        for name, delta, token in [('truncated', -8, 'truncated witness file'), ('extra', 8, 'extra witness data'), ('nonunit', 0, 'full mark unit')]:
            pre = tmp / name
            with Path(str(pre) + '.witness').open('wb') as f:
                f.truncate(expected + delta)
            z = subprocess.run([consumer, 'consume', str(R / 'data'), str(pre)], env=env, text=True, capture_output=True)
            need(z.returncode and token in z.stderr, 'specific rejection ' + name)
            cases.append(dict(name=name, status='REJECTED', diagnostic=token))
            log.append(z.stdout + z.stderr)
        d = tmp / 'badcrop'
        d.mkdir()
        shutil.copy2(R / 'data/log_table.txt', d / 'log_table.txt')
        lines = (R / 'data/grid.tsv').read_text().splitlines()
        lines[0] = str(int(lines[0]) + 1)
        (d / 'grid.tsv').write_text('\n'.join(lines) + '\n', newline='\n')
        z = subprocess.run([consumer, 'consume', str(d), str(tmp / 'nonunit')], env=env, text=True, capture_output=True)
        need(z.returncode and 'exact Gauss ceiling' in z.stderr, 'crop requires exact domain coverage')
        cases.append(dict(name='wrong_gauss_scale', status='REJECTED', diagnostic='exact Gauss ceiling'))
        log.append(z.stdout + z.stderr)
        lines = Path(str(prefix) + '.cells').read_text().splitlines()
        r = lines[0].split()
        r[3] = str(int(r[3]) + 1)
        lines[0] = ' '.join(r)
        bad = tmp / 'gap'
        Path(str(bad) + '.cells').write_text('\n'.join(lines) + '\n', newline='\n')
        try:
            assemble(bad)
        except ArithmeticError as e:
            need('subcell output exact geometry' in str(e), 'precise refined gap rejection')
        else:
            raise ArithmeticError('refined radial hole accepted')
        cases.append(dict(name='subrectangle_gap', status='REJECTED', diagnostic='subcell output exact geometry'))
    return (dict(status='PASS', cases=cases, synthetic_files_isolated=True), '\n'.join(log))
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--consumer', type=str, required=True)
    ap.add_argument('--prefix', type=Path, default=R / 'results/profile')
    ap.add_argument('--output', type=Path, default=R / 'results/negative.json')
    ap.add_argument('--log', type=Path, default=R / 'logs/negative.log')
    a = ap.parse_args()
    z, log = run(a.consumer, a.prefix)
    a.output.write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    a.log.write_text(log, newline='\n')
    print(json.dumps(z, indent=2))
