"""Portable refined-component replay; the declared coarse component is separate."""
from pathlib import Path
from hashlib import sha256
import argparse, json, os, shutil, subprocess, sys, tempfile, time
from assemble_refined import R, assemble
from exact import need
from checks import input_checks

def sha(path):
    h = sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--rebuild', action='store_true')
    ap.add_argument('--threads', type=int, default=4, choices=(1, 2, 3, 4))
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--baseline-root', type=Path)
    ap.add_argument('--report', type=Path)
    ap.add_argument('--keep-work', action='store_true')
    a = ap.parse_args()
    need(a.quick != a.rebuild, 'exactly one execution mode')
    start = time.monotonic()
    manifest = json.loads((R / 'MANIFEST_SHA256.json').read_text())
    for name, d in manifest.items():
        need((R / name).resolve().is_relative_to(R.resolve()) and sha(R / name) == d, 'manifest ' + name)
    need(input_checks(R) == json.loads((R / 'results/common_input_checks.json').read_text()), 'independent input checks')
    expected = json.loads((R / 'results/profile.json').read_text())
    need(assemble() == expected, 'exact hybrid integration')
    if a.baseline_root:
        for suffix in ['cells', 'json', 'branches', 'tensor', 'independent.json']:
            name = 'q2310_K13860_AG4.' + suffix
            need((a.baseline_root / 'results' / name).read_bytes() == (R / 'imports' / name).read_bytes(), 'external baseline bridge ' + suffix)
    compared = 0
    if a.rebuild:
        need(sys.version_info >= (3, 10) and sys.byteorder == 'little' and shutil.which('g++'), 'replay environment')
        work = Path(tempfile.mkdtemp(prefix='p72-2310-refined-'))
        pre = work / 'profile'
        env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
        flags = ['-O3', '-std=c++17', '-fopenmp'] + (['-fsanitize=undefined', '-fno-sanitize-recover=all'] if a.sanitize else [])

        def run(cmd):
            print('RUN', ' '.join(map(str, cmd)), flush=True)
            subprocess.run(list(map(str, cmd)), env=env, check=True)

        def same(x, y):
            nonlocal compared
            need(Path(x).read_bytes() == Path(y).read_bytes(), 'output byte identity ' + str(y))
            compared += 1
        try:
            for mode in ['produce', 'consume']:
                binary = work / mode
                run(['g++', *flags, *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), R / 'code/refine.cpp', '-o', binary])
                run([binary, mode, R / 'data', pre])
                if mode == 'produce':
                    w = json.loads((R / 'results/witness_metadata.json').read_text())
                    raw = Path(str(pre) + '.witness')
                    need(raw.stat().st_size == w['bytes'] and sha(raw) == w['sha256'], 'entire refined raw SHA256 and size')
                    compared += 1
            for ext in ['cells', 'audit', 'branches']:
                same(str(pre) + '.' + ext, R / 'results' / ('profile.' + ext))
            run(['g++', *flags, '-DINDEX_MULTIPLIER=6', R / 'code/algebra.cpp', '-o', work / 'algebra'])
            run([work / 'algebra', 11, 1, R / 'data', str(pre) + '.tensor'])
            same(str(pre) + '.tensor', R / 'results/profile.tensor')
            run([sys.executable, '-B', R / 'assemble_refined.py', pre, str(pre) + '.json'])
            same(str(pre) + '.json', R / 'results/profile.json')
            run([sys.executable, '-B', R / 'audit_refined.py', '--prefix', pre, '--output', str(pre) + '.independent.json'])
            same(str(pre) + '.independent.json', R / 'results/profile.independent.json')
            run([sys.executable, '-B', R / 'negative_refined.py', '--consumer', work / 'consume', '--prefix', pre, '--output', work / 'negative.json', '--log', work / 'negative.log'])
            same(work / 'negative.json', R / 'results/negative.json')
        finally:
            if a.keep_work:
                print('Retained', work)
            else:
                shutil.rmtree(work)
    for name, d in manifest.items():
        need(sha(R / name) == d, 'frozen input unchanged ' + name)
    report = dict(status='PASS', mode='quick_exact_hybrid_assembly' if a.quick else 'full_refined_witness_rebuild', manifest_files=len(manifest), compared_outputs_or_raw_hashes=compared, refined_witnesses_replayed=a.rebuild, coarse_witnesses_replayed_by_this_component=False, external_baseline_bridge_verified=bool(a.baseline_root), coarse_component_required=True, refined_witness_records=32223744, refined_corner_evaluations=17014136832, threads=a.threads if a.rebuild else None, sanitize=a.sanitize if a.rebuild else False, runtime_seconds=round(time.monotonic() - start, 3), secondary_geometry_exhaustive=False)
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report, indent=2), flush=True)
if __name__ == '__main__':
    main()
