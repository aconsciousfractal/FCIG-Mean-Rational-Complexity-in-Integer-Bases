"""Verify the three baseline marked2310 anchors; no whole-family closure claim."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import argparse, json, os, shutil, subprocess, sys, tempfile, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
from assemble import profile
from checks import input_checks

def assembled():
    rows = []
    for z in json.loads((ROOT / 'profiles_plan.json').read_text()):
        got = profile(z['p'], ROOT / 'results' / (z['stem'] + '.cells'), ROOT, z['k'], z['ag'], z['multiplier'])
        stored = json.loads((ROOT / 'results' / (z['stem'] + '.json')).read_text())
        need(set(stored) == set(got) | {'witness_bytes', 'witness_sha256'}, 'profile schema')
        need(all((stored[key] == value for key, value in got.items())), 'every integrated field')
        need(stored['witness_bytes'] == 8 * got['witnesses'] and len(stored['witness_sha256']) == 64, 'raw metadata')
        need(got['q'] == 2310 and got['transfer_index'] == z['index'] and (F(got['total_upper']) < F(z['simple_cap'])), 'anchor domain and simple cap')
        rows.append(dict(stem=z['stem'], mark=2310, divisor_of_M=z['index'], exact_upper=got['total_upper'], upper_decimal=got['upper_decimal'], simple_cap=z['simple_cap'], valid_for_every_positive_multiple_including_extra_primes=True, conditional_normalized_upper_for_B_ge_2310=got['normalized_upper'], delta_gate_at_base2310_under_divisibility=got['direct_radical_gate']))
    return dict(status='PASS', mark=2310, anchors=rows, whole_radical2310_closed_by_this_component=False, historical_numeric_caps_assumed=False, normalization='log B', all_residues_including_zero=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--rebuild', action='store_true')
    ap.add_argument('--threads', type=int, default=2, choices=(1, 2, 3, 4))
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--keep-work', action='store_true')
    ap.add_argument('--report', type=Path)
    a = ap.parse_args()
    need(a.quick != a.rebuild, 'choose exactly one of --quick or --rebuild')
    start = time.monotonic()
    manifest = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text())
    for name, digest in manifest.items():
        p = ROOT / name
        need(p.resolve().is_relative_to(ROOT.resolve()) and sha256(p.read_bytes()).hexdigest() == digest, 'manifest ' + name)
    inp = input_checks(ROOT)
    need(inp == json.loads((ROOT / 'results/common_input_checks.json').read_text()), 'independent logs')
    anchors = assembled()
    need(anchors == json.loads((ROOT / 'results/anchors.json').read_text()), 'exact anchor assembly')
    plan = json.loads((ROOT / 'profiles_plan.json').read_text())
    compared = 0
    if a.rebuild:
        need(sys.version_info >= (3, 10) and sys.byteorder == 'little' and shutil.which('g++'), 'replay environment')
        tmp = Path(tempfile.mkdtemp(prefix='p72-2310-baseline-'))
        env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
        flags = ['-O3', '-std=c++17', '-fopenmp'] + (['-fsanitize=undefined', '-fno-sanitize-recover=all'] if a.sanitize else [])

        def run(cmd):
            print('RUN', ' '.join(map(str, cmd)), flush=True)
            subprocess.run(list(map(str, cmd)), check=True, env=env)

        def same(x, y):
            nonlocal compared
            need(Path(x).read_bytes() == Path(y).read_bytes(), 'byte identity ' + str(y))
            compared += 1
        try:
            for row in plan:
                stem = row['stem']
                pre = tmp / stem
                expected = ROOT / 'results' / stem
                m = row['multiplier']
                index_flags = [] if m is None else [f'-DINDEX_MULTIPLIER={m}']
                binaries = {}
                for mode in ('produce', 'consume'):
                    binary = tmp / (stem + '_' + mode)
                    binaries[mode] = binary
                    run(['g++', *flags, *index_flags, f"-DDEGREE_EXP={row['k']}", f"-DANGULAR={row['ag']}", *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), ROOT / 'code/angle.cpp', '-o', binary])
                    run([binary, mode, row['p'], ROOT / 'data', pre])
                for ext in ('.cells', '.branches', '.audit'):
                    same(str(pre) + ext, str(pre) + '.checked' + ext)
                    same(str(pre) + ext, str(expected) + ext)
                binary = tmp / (stem + '_algebra')
                run(['g++', *flags, *index_flags, ROOT / 'code/algebra.cpp', '-o', binary])
                run([binary, row['p'], row['k'], ROOT / 'data', str(pre) + '.tensor'])
                same(str(pre) + '.tensor', str(expected) + '.tensor')
                extra = [] if m is None else ['--multiplier', m]
                run([sys.executable, '-B', ROOT / 'code/emit_profile.py', row['p'], row['ag'], row['k'], pre, str(pre) + '.json', *extra])
                same(str(pre) + '.json', str(expected) + '.json')
                run([sys.executable, '-B', ROOT / 'code/audit_one.py', row['p'], row['ag'], '--k', row['k'], *extra, '--prefix', pre, '--output', str(pre) + '.independent.json'])
                same(str(pre) + '.independent.json', str(expected) + '.independent.json')
                run([sys.executable, '-B', ROOT / 'code/negative_profiles.py', '--consumer', binaries['consume'], '--multiplier', m or 1, '--output', tmp / 'negative.json', '--log', tmp / 'negative.log'])
                same(tmp / 'negative.json', ROOT / 'results/negative_direct2310.json')
                Path(str(pre) + '.witness').unlink()
        finally:
            if a.keep_work:
                print('Retained', tmp)
            else:
                shutil.rmtree(tmp)
    for name, digest in manifest.items():
        need(sha256((ROOT / name).read_bytes()).hexdigest() == digest, 'unchanged frozen input')
    report = dict(status='PASS', mode='stored_cells_and_exact_assembly' if a.quick else 'full_witness_rebuild', profiles=3, manifest_files=len(manifest), geometric_witnesses_replayed=a.rebuild, compared_files=compared, geometric_witnesses=167215104, geometric_corner_evaluations=6019743744, threads=a.threads if a.rebuild else None, sanitize=a.sanitize if a.rebuild else False, whole_radical2310_closed_by_this_component=False, runtime_seconds=round(time.monotonic() - start, 3), alternative_geometric_audit_exhaustive=False)
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report, indent=2), flush=True)
if __name__ == '__main__':
    main()
