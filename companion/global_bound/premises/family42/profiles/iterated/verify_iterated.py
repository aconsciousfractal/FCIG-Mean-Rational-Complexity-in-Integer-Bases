"""Rebuild the q462 two-step envelope and the exact residual base462."""
from pathlib import Path
from fractions import Fraction as F
import sys, os, tempfile, subprocess, json, argparse, hashlib, shutil, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need, log_lo

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--threads', type=int, default=4, choices=(1, 2, 3, 4))
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--keep-work', action='store_true')
    ap.add_argument('--report', type=Path)
    a = ap.parse_args()
    need(sys.byteorder == 'little', 'little-endian int16 record domain')
    temp = Path(tempfile.mkdtemp(prefix='p72-462-iterator-'))
    start = time.monotonic()
    env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    compared = 0

    def run(cmd):
        print('RUN', ' '.join(map(str, cmd)), flush=True)
        subprocess.run(list(map(str, cmd)), env=env, check=True)

    def same(x, y, label):
        nonlocal compared
        need(Path(x).read_bytes() == Path(y).read_bytes(), label)
        compared += 1
    flags = ['-O3', '-std=c++17', '-fopenmp']
    if a.sanitize:
        flags += ['-fsanitize=undefined', '-fno-sanitize-recover=all']
    try:
        for mode in ('produce', 'consume'):
            binary = temp / mode
            ext = ['-DFOUR_PRIMES=1', '-DPRIME=11', '-DDIR_T=8', '-DAUDIT_SPARSE=1'] + (['-DBUILD_CONSUMER=1'] if mode == 'consume' else [])
            run(['g++', *flags, *ext, ROOT / 'code/uniform.cpp', '-o', binary])
            pref = temp / ('q462' if mode == 'produce' else 'q462_checked')
            cmd = [binary, ROOT / 'data/log_table.txt', ROOT / 'data/partition.txt', ROOT / 'data/shape_centres.tsv', pref, 2, ROOT / 'data/cells.txt']
            if mode == 'consume':
                cmd += [temp / 'q462.witnesses']
            run(cmd)
        for ext in ('', '.audit', '.cellcosts', '.points', '.bank', '.cells1', '.radial1', '.cells2', '.radial2'):
            same(str(temp / 'q462') + ext, str(temp / 'q462_checked') + ext, 'producer-consumer output' + ext)
            same(str(temp / 'q462') + ext, str(ROOT / 'results/q462') + ext, 'immutable expected output' + ext)
        same(temp / 'q462.witnesses', temp / 'q462_checked.witnesses', 'all initial coefficients')
        manifest = json.loads((ROOT / 'results/initial_coefficient_manifest.json').read_text())
        wf = temp / 'q462.witnesses'
        need(wf.stat().st_size == manifest['bytes'] and sha(wf) == manifest['sha256'], 'initial witness hash and length')
        run(['g++', *flags, ROOT / 'code/reference_min.cpp', '-o', temp / 'reference'])
        for k in (1, 2):
            M = 462 ** k
            name = f'M{M}_q462.min'
            run([temp / 'reference', M, 462, temp / name])
            same(temp / name, ROOT / 'results' / name, 'complete definition minima ' + name)
        run([sys.executable, '-B', ROOT / 'code/integrate.py', '--prefix', temp / 'q462', '--output', temp / 'integrated.json'])
        same(temp / 'integrated.json', ROOT / 'results/integrated.json', 'integrated exact bounds')
        run([sys.executable, '-B', ROOT / 'code/audit_iterated.py', '--prefix', temp / 'q462', '--reference-dir', temp, '--output-dir', temp])
        for name in ('independent_checks.json', 'finite_base462.json'):
            same(temp / name, ROOT / 'results' / name, 'independent output ' + name)
        rows = json.loads((temp / 'integrated.json').read_text())
        bound = next((z for z in rows if z['k'] == 2))
        need(F(bound['total_upper']) < F(2, 85) * log_lo(462), 'uniform N>=2 gate')
        finite = json.loads((temp / 'finite_base462.json').read_text())
        need(F(finite['normalized_upper']) < F(2, 85), 'exact base462 gate')
        report = {'status': 'PASS', 'compared_files': compared, 'initial_witnesses': 2688 * 1152, 'uniform_horizons_regenerated': [1, 2], 'selected_uniform_horizon': 2, 'exact_minimum_states': 462 + 462 ** 2, 'finite_base462': True, 'threads': a.threads, 'sanitize': a.sanitize, 'python_optimized': not __debug__, 'runtime_seconds': round(time.monotonic() - start, 3), 'second_geometric_audit_exhaustive': False}
        if a.report:
            a.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
        print(json.dumps(report, indent=2), flush=True)
    finally:
        if a.keep_work:
            print('Retained', temp, flush=True)
        else:
            shutil.rmtree(temp)
if __name__ == '__main__':
    main()
