"""Portable verification of the new 165 profiles and their declared scope."""
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

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def manifest_check():
    z = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text())
    present = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and p.name != 'MANIFEST_SHA256.json'}
    need(set(z) == present, 'exact complete manifest inventory')
    for name, digest in z.items():
        p = ROOT / name
        need(p.resolve().is_relative_to(ROOT), 'contained manifest path')
        need(sha(p) == digest, 'manifest SHA256 ' + name)
    return len(z)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--quick', action='store_true')
    group.add_argument('--rebuild', action='store_true')
    ap.add_argument('--threads', type=int, default=3, choices=(1, 2, 3, 4))
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--keep-work', action='store_true')
    ap.add_argument('--report', type=Path)
    ap.add_argument('--no-manifest', action='store_true', help=argparse.SUPPRESS)
    args = ap.parse_args()
    need(sys.version_info >= (3, 10), 'Python3.10 required')
    if args.report:
        args.report = args.report.resolve()
        need(not args.report.is_relative_to(ROOT), 'report outside immutable component')
        args.report.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    count = 0 if args.no_manifest else manifest_check()
    plan = json.loads((ROOT / 'profiles_plan.json').read_text())
    common = input_checks(ROOT)
    need(common == json.loads((ROOT / 'results/common_input_checks.json').read_text()), 'all independent logarithm table checks')
    pcdf = pcdf_check()
    need(pcdf == json.loads((ROOT / 'results/pcdf_certificate.json').read_text()), 'complete primitive-CDF finite envelope and tail')
    custody = {'kind': 'new_reader_inputs'}
    assembly = assemble()
    need(assembly == json.loads((ROOT / 'results/profile_certificate.json').read_text()), 'all profile integrals and scope reassembled')
    neg = json.loads((ROOT / 'results/negative_checks.json').read_text())
    need(neg['status'] == 'PASS' and len(neg['cases']) == EXPECTED_CASES, 'stored negative report schema')
    temp = Path(tempfile.mkdtemp(prefix='p72-165-profile-replay-'))
    env = {**os.environ, 'CSTAR_THREADS': str(args.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    compared = 0
    states = 0
    sample_witnesses = 0
    sample_logs = 0
    commands = []

    def run(cmd):
        t = time.monotonic()
        print('RUN', ' '.join(map(str, cmd)), flush=True)
        r = subprocess.run(list(map(str, cmd)), env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(r.stdout, end='', flush=True)
        commands.append({'argv': list(map(str, cmd)), 'returncode': r.returncode, 'seconds': round(time.monotonic() - t, 3)})
        need(r.returncode == 0, 'subprocess succeeded')

    def same(a, b, label):
        nonlocal compared
        need(Path(a).read_bytes() == Path(b).read_bytes(), label)
        compared += 1
    success = False
    try:
        mixed = temp / 'mixed_index_certificate.json'
        run([sys.executable, '-B', ROOT / 'code/mixed_indices.py', '--results', ROOT / 'results', '--output', mixed])
        same(mixed, ROOT / 'results/mixed_index_certificate.json', 'all local quotient classes and every signed tensor row for four unequal indices')
        mixed_report = json.loads(mixed.read_text())
        bins = {}
        if args.rebuild:
            need(sys.byteorder == 'little', 'little endian record format')
            need(shutil.which('g++') is not None, 'GNU C++17 compiler')
            largest = max((8 * 2688 * CORE_LABELS * (z['p'] + 1 if z['p'] > 1 else z['multiplier']) for z in plan))
            need(shutil.disk_usage(temp).free >= 2 * largest + 64 * 1024 * 1024, 'sufficient temporary disk for largest stream')
            flags = ['-O3', '-std=c++17', '-fopenmp'] + (['-fsanitize=undefined', '-fno-sanitize-recover=all'] if args.sanitize else [])
            algebra = temp / 'algebra'
            run(['g++', *flags, ROOT / 'code/algebra.cpp', '-o', algebra])
            for ag in sorted({z['ag'] for z in plan}):
                b = {}
                for mode in ('produce', 'consume'):
                    target = temp / f'{mode}{ag}'
                    extra = [f'-DANGULAR={ag}'] + (['-DBUILD_CONSUMER=1'] if mode == 'consume' else [])
                    run(['g++', *flags, *extra, ROOT / 'code/angle.cpp', '-o', target])
                    b[mode] = target
                bins[ag] = b
        for z in plan:
            expected = ROOT / 'results' / z['stem']
            pre = temp / z['stem'] if args.rebuild else expected
            if args.rebuild:
                for mode in ('produce', 'consume'):
                    run([bins[z['ag']][mode], mode, z['p'], z['multiplier'], ROOT / 'data', pre])
                for ext in ('.cells', '.branches', '.audit'):
                    same(str(pre) + ext, str(pre) + '.checked' + ext, 'complete producer-consumer agreement ' + z['stem'] + ext)
                    same(str(pre) + ext, str(expected) + ext, 'frozen geometry bytes ' + z['stem'] + ext)
                run([algebra, z['p'], z['multiplier'], ROOT / 'data', str(pre) + '.tensor'])
                same(str(pre) + '.tensor', str(expected) + '.tensor', 'complete tensor algebra')
                run([sys.executable, '-B', ROOT / 'code/emit_profile.py', z['p'], z['ag'], z['multiplier'], pre, str(pre) + '.json'])
                same(str(pre) + '.json', str(expected) + '.json', 'exact integration and full stream hash')
            audit = temp / (z['stem'] + '.independent.json')
            run([sys.executable, '-B', ROOT / 'code/audit_one.py', z['p'], z['ag'], '--k', z['multiplier'], '--prefix', pre, '--output', audit])
            same(audit, str(expected) + '.independent.json', 'independent algebra, stratified geometry and complete physical banks')
            rec = json.loads(audit.read_text())
            states += rec['physical']['state_modulus_incidences']
            sample_witnesses += rec['stratified_geometry']['witnesses']
            sample_logs += rec['stratified_geometry']['independent_corner_logarithms']
            if args.rebuild:
                Path(str(pre) + '.witness').unlink()
        if args.rebuild:
            run([sys.executable, '-B', ROOT / 'code/negative_profiles.py', '--consumer', bins[4]['consume'], '--output', temp / 'negative.json', '--log', temp / 'negative.log'])
            same(temp / 'negative.json', ROOT / 'results/negative_checks.json', 'all explicit165 negative and necessity gates')
        expected_comparisons = len(plan) * (9 if args.rebuild else 1) + 1 + int(args.rebuild)
        need(compared == expected_comparisons, 'dynamic exact comparison contract for every profile, mixed-index certificate and fixtures')
        report = {'status': 'PASS', 'mode': 'full_witness_rebuild' if args.rebuild else 'stored_cells_with_exact_integration_and_secondary_audits', 'manifest_files': count, 'profiles': len(plan), 'mixed_index_exact_checks': mixed_report['checks'], 'completed_generation_custody': custody, 'small_primes': assembly['closed_primes'], 'all_positive_prime_exponents': True, 'all_N_positive': True, 'cdf_alpha': '3/2', 'cdf_exact_checks': pcdf['checks'], 'geometric_witnesses': assembly['geometric_witnesses'], 'geometric_corner_evaluations': assembly['geometric_corner_evaluations'], 'primary_geometric_witnesses_replayed': args.rebuild, 'negative_tests_reexecuted': args.rebuild, 'secondary_geometric_audit_exhaustive': False, 'secondary_geometric_sample_witnesses': sample_witnesses, 'secondary_corner_logarithms': sample_logs, 'secondary_complete_physical_states': states, 'compared_files': compared, 'expected_comparisons': expected_comparisons, 'threads': args.threads, 'sanitize': args.sanitize if args.rebuild else None, 'runtime_seconds': round(time.monotonic() - start, 3), 'commands': commands}
        if args.report:
            args.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
        print(json.dumps({k: v for k, v in report.items() if k != 'commands'}, indent=2), flush=True)
        success = True
    finally:
        if args.keep_work or not success:
            print('Retained', temp, flush=True)
        else:
            shutil.rmtree(temp)
if __name__ == '__main__':
    main()
