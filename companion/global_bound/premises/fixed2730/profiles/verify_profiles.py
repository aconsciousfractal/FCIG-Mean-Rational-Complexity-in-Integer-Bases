"""Read-only verification and complete regeneration of the declared profiles."""
from pathlib import Path
from hashlib import sha256
import argparse, json, os, shutil, subprocess, sys, tempfile, time
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'code'))
from exact import need
from integrate import integrate
from checks import input_checks
from integrate_alpha import assemble_alpha
from fractions import Fraction as F

def digest(p):
    h = sha256()
    with Path(p).open('rb') as f:
        while (b := f.read(8 * 1024 * 1024)):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--quick', action='store_true')
    group.add_argument('--rebuild', action='store_true')
    ap.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=2)
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--report', type=Path)
    ap.add_argument('--keep-work', action='store_true')
    a = ap.parse_args()
    start = time.monotonic()
    need(sys.version_info >= (3, 10), 'Python3.10+ with assertions enabled')
    if a.report:
        need(not a.report.resolve().is_relative_to(R), 'write reports outside the frozen component')
    manifest = json.loads((R / 'MANIFEST_SHA256.json').read_text())
    actual = {p.relative_to(R).as_posix() for p in R.rglob('*') if p.is_file() and p.name != 'MANIFEST_SHA256.json'}
    need(actual == set(manifest), 'complete component manifest')
    for name, h in manifest.items():
        need((R / name).resolve().is_relative_to(R) and digest(R / name) == h, 'input hash ' + name)
    plan = json.loads((R / 'profiles_plan.json').read_text())
    anchors = []
    count = 0
    need(input_checks(R) == json.loads((R / 'results/common_input_checks.json').read_text()), 'independent logarithm inputs')
    for row in plan:
        stem = row['stem']
        pre = R / 'results' / stem
        z = integrate(row['p'], row['multiplier'], row['style'], str(pre) + '.cells', R, row['angular'])
        need(z == json.loads(Path(str(pre) + '.json').read_text()), 'every exact integrated field')
        alt = assemble_alpha(row['p'], row['multiplier'], row['style'], row['angular'], pre, F(3, 2))
        need(alt == json.loads(Path(str(pre) + '.alpha3over2.json').read_text()), 'conditional alternative CDF assembly')
        raw = json.loads(Path(str(pre) + '.witness_metadata.json').read_text())
        need(raw['bytes'] == z['witness_bytes'] and len(raw['sha256']) == 64, 'raw contract identity')
        anchors.append(z)
    if a.rebuild:
        need(shutil.which('g++') and sys.byteorder == 'little', 'GNU C++17 and little-endian replay environment')
        tmp = Path(tempfile.mkdtemp(prefix='p72-2730-profile-replay-'))
        env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
        flags = ['-O3', '-std=c++17', '-fopenmp'] + (['-fsanitize=undefined', '-fno-sanitize-recover=all'] if a.sanitize else [])

        def run(cmd):
            print('RUN', ' '.join(map(str, cmd)), flush=True)
            subprocess.run(list(map(str, cmd)), check=True, env=env)

        def same(x, y):
            nonlocal count
            need(digest(x) == digest(y), 'full output identity ' + str(y))
            count += 1
        try:
            for row in plan:
                stem = row['stem']
                pre = tmp / stem
                old = R / 'results' / stem
                defs = [f"-DINDEX_MULTIPLIER={row['multiplier']}", f"-DFINE_ANGULAR={row['angular']}"]
                algebra = tmp / (stem + '_algebra')
                run(['g++', *flags, *defs, R / 'code/algebra.cpp', '-o', algebra])
                run([algebra, row['p'], 1, R / 'data', str(pre) + '.tensor'])
                for mode in ('produce', 'consume'):
                    binary = tmp / (stem + '_' + mode)
                    run(['g++', *flags, *defs, *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), R / 'code/profile.cpp', '-o', binary])
                    run([binary, mode, row['p'], row['style'], R / 'data', pre])
                metadata = json.loads(Path(str(old) + '.witness_metadata.json').read_text())
                need(Path(str(pre) + '.witness').stat().st_size == metadata['bytes'] and digest(str(pre) + '.witness') == metadata['sha256'], 'every regenerated witness byte')
                count += 1
                for ext in ('.cells', '.audit', '.branches', '.tensor'):
                    same(str(pre) + ext, str(old) + ext)
                run([sys.executable, '-B', R / 'code/integrate.py', row['p'], row['multiplier'], row['angular'], row['style'], pre, str(pre) + '.json'])
                same(str(pre) + '.json', str(old) + '.json')
                run([sys.executable, '-B', R / 'code/integrate_alpha.py', row['p'], row['multiplier'], row['style'], '--angular', row['angular'], '--prefix', pre, '--output', str(pre) + '.alpha3over2.json'])
                same(str(pre) + '.alpha3over2.json', str(old) + '.alpha3over2.json')
                run([sys.executable, '-B', R / 'code/audit_profile.py', row['p'], row['multiplier'], row['style'], '--angular', row['angular'], '--prefix', pre, '--output', str(pre) + '.independent.json'])
                same(str(pre) + '.independent.json', str(old) + '.independent.json')
                run([sys.executable, '-B', R / 'code/negative_profile.py', '--consumer', tmp / (stem + '_consume'), '--p', row['p'], '--multiplier', row['multiplier'], '--style', row['style'], '--angular', row['angular'], '--prefix', pre, '--output', str(pre) + '.negative.json', '--log', tmp / 'negative.log'])
                same(str(pre) + '.negative.json', str(old) + '.negative.json')
                Path(str(pre) + '.witness').unlink()
        finally:
            if a.keep_work:
                print('Retained work', tmp)
            else:
                shutil.rmtree(tmp)
    for name, h in manifest.items():
        need(digest(R / name) == h, 'frozen inputs unchanged')
    report = dict(status='PASS', mode='full_witness_rebuild' if a.rebuild else 'stored_cells_and_exact_assembly', profiles=len(plan), manifest_files=len(manifest), compared_files=count, geometric_witness_records=sum((z['geometric_witness_records'] for z in anchors)), geometric_corner_evaluations=sum((z['geometric_corner_evaluations'] for z in anchors)), geometric_witnesses_replayed=a.rebuild, threads=a.threads if a.rebuild else None, sanitize=a.sanitize if a.rebuild else False, runtime_seconds=round(time.monotonic() - start, 3), alternative_geometric_audit_exhaustive=False, historical_numerical_profiles_used=False)
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report, indent=2), flush=True)
if __name__ == '__main__':
    main()
