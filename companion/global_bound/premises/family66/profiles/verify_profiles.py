"""Regenerate every selected full-mark66p and marked66 profile, then consume it."""
from pathlib import Path
from fractions import Fraction as F
import argparse, subprocess, sys, os, tempfile, json, time, shutil
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need, sieve, log_lo
from core import MIN_NEW_PRIME
from verify_stored import verify_manifest

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--rebuild', action='store_true')
    ap.add_argument('--threads', type=int, default=4, choices=(1, 2, 3, 4))
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--keep-work', action='store_true')
    ap.add_argument('--report', type=Path)
    args = ap.parse_args()
    if args.quick:
        need(not args.rebuild, 'choose quick or rebuild')
        subprocess.run([sys.executable, '-B', str(ROOT / 'verify_stored.py'), *(['--report', str(args.report)] if args.report else [])], check=True)
        return
    need(sys.version_info >= (3, 10), 'Python3.10 required')
    need(sys.byteorder == 'little', 'little-endian record format')
    need(shutil.which('g++') is not None, 'GNU C++ required')
    manifest_count = verify_manifest()
    cfg = json.loads((ROOT / 'profiles_plan.json').read_text())
    full = [r for r in cfg if r['p'] > 1]
    contract = json.loads((ROOT / 'profile_contract.json').read_text())
    need(sorted((r['p'] for r in full)) == [p for p in sieve(contract['last_small_prime']) if p >= MIN_NEW_PRIME], 'every declared small prime exactly once')
    temp = Path(tempfile.mkdtemp(prefix='p72-66p-profile-'))
    env = {**os.environ, 'CSTAR_THREADS': str(args.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    start = time.monotonic()
    bins = {}
    compared = 0
    rows = []
    second_states = 0

    def run(cmd):
        print('RUN', ' '.join(map(str, cmd)), flush=True)
        subprocess.run(list(map(str, cmd)), check=True, env=env)

    def same(a, b, msg):
        nonlocal compared
        need(Path(a).read_bytes() == Path(b).read_bytes(), msg)
        compared += 1
    try:
        flags = ['-O3', '-std=c++17', '-fopenmp'] + (['-fsanitize=undefined', '-fno-sanitize-recover=all'] if args.sanitize else [])
        algebra = temp / 'algebra'
        run(['g++', *flags, ROOT / 'code/algebra.cpp', '-o', algebra])
        for z in cfg:
            key = (z['ag'], z['k'])
            if key not in bins:
                b = {}
                for mode in ('produce', 'consume'):
                    target = temp / ('bin_' + str(len(bins)) + '_' + mode)
                    ext = [f"-DANGULAR={z['ag']}", f"-DDEGREE_EXP={z['k']}"] + (['-DBUILD_CONSUMER=1'] if mode == 'consume' else [])
                    run(['g++', *flags, *ext, ROOT / 'code/angle.cpp', '-o', target])
                    b[mode] = target
                bins[key] = b
            b = bins[key]
            name = z['stem']
            pre = temp / name
            expected = ROOT / 'results' / name
            for mode in ('produce', 'consume'):
                run([b[mode], mode, z['p'], ROOT / 'data', pre])
            for ext in ('.cells', '.branches', '.audit'):
                same(str(pre) + ext, str(pre) + '.checked' + ext, 'producer-consumer agreement ' + name + ext)
                same(str(pre) + ext, str(expected) + ext, 'expected certificate bytes ' + name + ext)
            run([algebra, z['p'], z['k'], ROOT / 'data', str(pre) + '.tensor'])
            same(str(pre) + '.tensor', str(expected) + '.tensor', 'full algebra bytes ' + name)
            run([sys.executable, '-B', ROOT / 'code/emit_profile.py', z['p'], z['ag'], z['k'], pre, str(pre) + '.json'])
            same(str(pre) + '.json', str(expected) + '.json', 'profile integration and witness SHA ' + name)
            row = json.loads(Path(str(pre) + '.json').read_text())
            if z['p'] > 1:
                need(F(row['total_upper']) < F(2, 85) * log_lo(row['q']), 'strict normalized direct route ' + name)
            else:
                need(F(row['total_upper']) < F(z['simple_cap']), 'marked66 simple cap ' + name)
            run([sys.executable, '-B', ROOT / 'code/audit_one.py', z['p'], z['ag'], '--k', z['k'], '--prefix', pre, '--output', str(pre) + '.independent.json'])
            same(str(pre) + '.independent.json', str(expected) + '.independent.json', 'second-audit output ' + name)
            second_states += json.loads(Path(str(pre) + '.independent.json').read_text())['physical']['state_modulus_incidences']
            if z['p'] == MIN_NEW_PRIME:
                run([sys.executable, '-B', ROOT / 'code/negative_profiles.py', '--consumer', b['consume'], '--output', temp / 'negative.json', '--log', temp / 'negative.log'])
                same(temp / 'negative.json', ROOT / 'results/negative_checks.json', 'negative geometric gates')
            rows.append(row)
            Path(str(pre) + '.witness').unlink()
        script = "import sys,json;from pathlib import Path;sys.path.insert(0,sys.argv[1]);from checks import input_checks;Path(sys.argv[3]).write_text(json.dumps(input_checks(Path(sys.argv[2])),indent=2)+'\\n')"
        run([sys.executable, '-B', '-c', script, ROOT / 'code', ROOT, temp / 'input.json'])
        same(temp / 'input.json', ROOT / 'results/common_input_checks.json', 'all log table and seed checks')
        run([sys.executable, '-B', ROOT / 'code/assemble_small.py', '--check'])
        report = {'status': 'PASS', 'mode': 'full_witness_rebuild', 'manifest_files': manifest_count, 'profiles': len(rows), 'full_mark_primes': len(full), 'all_prime_exponents_and_lengths': True, 'compared_files': compared, 'geometric_witnesses': sum((z['witnesses'] for z in rows)), 'geometric_corner_evaluations': sum((z['witnesses'] * z['angular_nodes'] * 4 for z in rows)), 'second_definition_based_states': second_states, 'threads': args.threads, 'sanitize': args.sanitize, 'python_optimized': not __debug__, 'runtime_seconds': round(time.monotonic() - start, 3), 'second_geometric_audit_exhaustive': False}
        if args.report:
            args.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
        print(json.dumps(report, indent=2), flush=True)
    finally:
        if args.keep_work:
            print('Retained', temp, flush=True)
        else:
            shutil.rmtree(temp)
if __name__ == '__main__':
    main()
