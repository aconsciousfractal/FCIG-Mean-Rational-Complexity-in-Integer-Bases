"""Portable fixed4290 profiles: exact quick replay or complete primary reconstruction."""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
from checks import input_checks
from pcdf_check import check as pcdf_check
from assemble_profiles import assemble
from core import CORE_LABELS
from negative_profiles import EXPECTED_CASES

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def manifest_check():
    z = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text())
    present = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and p.name != 'MANIFEST_SHA256.json'}
    need(set(z) == present, 'exact complete manifest inventory')
    for name, d in z.items():
        need((ROOT / name).resolve().is_relative_to(ROOT) and sha(ROOT / name) == d, 'manifest SHA256 ' + name)
    return len(z)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--quick', action='store_true')
    group.add_argument('--rebuild', action='store_true')
    ap.add_argument('--threads', type=int, default=4, choices=(1, 2, 3, 4))
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--keep-work', action='store_true')
    ap.add_argument('--report', type=Path)
    ap.add_argument('--no-manifest', action='store_true', help=argparse.SUPPRESS)
    a = ap.parse_args()
    need(sys.version_info >= (3, 10), 'Python3.10+ without -O')
    if a.report:
        a.report = a.report.resolve()
        need(not a.report.is_relative_to(ROOT), 'report outside immutable component')
        a.report.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    count = 0 if a.no_manifest else manifest_check()
    plan = json.loads((ROOT / 'profiles_plan.json').read_text())
    need(input_checks(ROOT) == json.loads((ROOT / 'results/common_input_checks.json').read_text()), 'independent logarithm input checks')
    pcdf = pcdf_check()
    need(pcdf == json.loads((ROOT / 'results/pcdf_certificate.json').read_text()), 'exact primitive-CDF envelope and tail')
    custody = {'kind': 'new_reader_inputs'}
    assembly = assemble()
    need(assembly == json.loads((ROOT / 'results/profile_certificate.json').read_text()), 'exact profile integration and marked scope')
    neg = json.loads((ROOT / 'results/negative_checks.json').read_text())
    need(neg['status'] == 'PASS' and len(neg['cases']) == EXPECTED_CASES, 'stored fixture schema')
    temp = Path(tempfile.mkdtemp(prefix='p72-4290-profile-replay-'))
    env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    compared = states = samples = samplelogs = 0
    commands = []
    success = False

    def run(cmd):
        t = time.monotonic()
        print('RUN', ' '.join(map(str, cmd)), flush=True)
        r = subprocess.run(list(map(str, cmd)), env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(r.stdout, end='', flush=True)
        commands.append({'argv': list(map(str, cmd)), 'returncode': r.returncode, 'seconds': round(time.monotonic() - t, 3)})
        need(r.returncode == 0, 'subprocess succeeded')

    def same(x, y, label):
        nonlocal compared
        need(Path(x).read_bytes() == Path(y).read_bytes(), label)
        compared += 1
    try:
        bins = {}
        if a.rebuild:
            need(sys.byteorder == 'little' and shutil.which('g++') is not None, 'little endian and GNU C++17 compiler')
            largest = max((8 * 2688 * CORE_LABELS * z['multiplier'] for z in plan))
            need(shutil.disk_usage(temp).free >= 2 * largest + 64 * 1024 * 1024, 'sufficient temporary disk for stream and sequential fixtures')
            flags = ['-O3', '-std=c++17', '-fopenmp'] + (['-fsanitize=undefined', '-fno-sanitize-recover=all'] if a.sanitize else [])
            algebra = temp / 'algebra'
            run(['g++', *flags, ROOT / 'code/algebra.cpp', '-o', algebra])
            for ag in sorted({z['ag'] for z in plan}):
                bins[ag] = {}
                for mode in ('produce', 'consume'):
                    target = temp / f'{mode}{ag}'
                    run(['g++', *flags, f'-DANGULAR={ag}', *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), ROOT / 'code/angle.cpp', '-o', target])
                    bins[ag][mode] = target
        for z in plan:
            expected = ROOT / 'results' / z['stem']
            pre = temp / z['stem'] if a.rebuild else expected
            if a.rebuild:
                for mode in ('produce', 'consume'):
                    run([bins[z['ag']][mode], mode, z['p'], z['multiplier'], ROOT / 'data', pre])
                for ext in ('.cells', '.branches', '.audit'):
                    same(str(pre) + ext, str(pre) + '.checked' + ext, 'complete producer-consumer agreement ' + ext)
                    same(str(pre) + ext, str(expected) + ext, 'frozen geometry bytes ' + ext)
                run([algebra, z['p'], z['multiplier'], ROOT / 'data', str(pre) + '.tensor'])
                same(str(pre) + '.tensor', str(expected) + '.tensor', 'complete tensor algebra')
                run([sys.executable, '-B', ROOT / 'code/emit_profile.py', z['p'], z['ag'], z['multiplier'], pre, str(pre) + '.json'])
                same(str(pre) + '.json', str(expected) + '.json', 'exact integration and full raw hash')
            audit = temp / (z['stem'] + '.independent.json')
            run([sys.executable, '-B', ROOT / 'code/audit_one.py', z['p'], z['ag'], '--k', z['multiplier'], '--prefix', pre, '--output', audit])
            same(audit, str(expected) + '.independent.json', 'independent HNF/tensor, stratified geometry and complete physical banks')
            rec = json.loads(audit.read_text())
            states += rec['physical']['state_modulus_incidences']
            samples += rec['stratified_geometry']['witnesses']
            samplelogs += rec['stratified_geometry']['independent_corner_logarithms']
            if a.rebuild:
                if z['multiplier'] == 1:
                    run([sys.executable, '-B', ROOT / 'code/negative_profiles.py', '--consumer', bins[8]['consume'], '--positive-prefix', pre, '--output', temp / 'negative.json', '--log', temp / 'negative.log'])
                    same(temp / 'negative.json', ROOT / 'results/negative_checks.json', 'all fixed4290 fixture gates including accepted7')
                Path(str(pre) + '.witness').unlink()
        report = {'status': 'PASS', 'mode': 'full_witness_rebuild' if a.rebuild else 'stored_cells_with_exact_integration_and_secondary_audits', 'manifest_files': count, 'profiles': len(plan), 'completed_generation_custody': custody, 'marked_radical': 4290, 'arbitrary_positive_valuations': True, 'whole_radical_directly_closed_here': assembly['whole_radical_directly_closed_here'], 'cdf_alpha': '3/2', 'cdf_exact_checks': pcdf['checks'], 'geometric_witnesses': assembly['geometric_witnesses'], 'geometric_corner_evaluations': assembly['geometric_corner_evaluations'], 'primary_geometric_witnesses_replayed': a.rebuild, 'fixture_tests_reexecuted': a.rebuild, 'secondary_geometric_audit_exhaustive': False, 'secondary_geometric_sample_witnesses': samples, 'secondary_corner_logarithms': samplelogs, 'secondary_complete_physical_states': states, 'compared_files': compared, 'threads': a.threads, 'sanitize': a.sanitize if a.rebuild else None, 'runtime_seconds': round(time.monotonic() - start, 3), 'commands': commands}
        if not a.no_manifest:
            need(manifest_check() == count, 'manifest unchanged after verification')
        if a.report:
            a.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
        print(json.dumps({k: v for k, v in report.items() if k != 'commands'}, indent=2), flush=True)
        success = True
    finally:
        if a.keep_work or not success:
            print('Retained', temp, flush=True)
        else:
            shutil.rmtree(temp)
if __name__ == '__main__':
    main()
