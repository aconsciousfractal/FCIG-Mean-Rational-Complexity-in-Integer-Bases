"""Regenerate the selected full-mark 42p and marked-42 profiles.
Expected files are immutable; all witnesses are generated and consumed in a temporary
directory. A second Python audit is stratified on geometric costs, complete on HNF
branches and initial arithmetic tables. See PROOF.md for the finite-to-general step.
"""
from pathlib import Path
from fractions import Fraction
import argparse, subprocess, sys, os, tempfile, json, time, shutil
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need, sieve

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true', help='check stored caps and route assembly; no geometric witness replay')
    ap.add_argument('--rebuild', action='store_true', help='explicit full replay, also the default')
    ap.add_argument('--threads', type=int, default=3, choices=(1, 2, 3, 4))
    ap.add_argument('--sanitize', action='store_true')
    ap.add_argument('--keep-work', action='store_true')
    ap.add_argument('--report', type=Path)
    args = ap.parse_args()
    if args.quick:
        need(not args.rebuild, 'choose quick or rebuild')
        cmd = [sys.executable, '-B', str(ROOT / 'verify_stored.py')] + (['--report', str(args.report)] if args.report else [])
        subprocess.run(cmd, check=True)
        return
    need(sys.version_info >= (3, 10), 'Python3.10 required')
    need(sys.byteorder == 'little', 'little-endian record format')
    need(shutil.which('g++') is not None, 'GNU C++ required')
    cfg = json.loads((ROOT / 'profiles_plan.json').read_text())
    full = [r for r in cfg if r['p'] > 1]
    need(sorted((r['p'] for r in full)) == [p for p in sieve(103) if p >= 11], 'every declared prime11..103 exactly once')
    temp = Path(tempfile.mkdtemp(prefix='p72-42p-profile-'))
    env = {**os.environ, 'CSTAR_THREADS': str(args.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    start = time.monotonic()
    bins = {}
    compared = 0
    rows = []

    def run(cmd):
        print('RUN', ' '.join(map(str, cmd)), flush=True)
        subprocess.run(list(map(str, cmd)), check=True, env=env)

    def same(a, b, msg):
        nonlocal compared
        need(Path(a).read_bytes() == Path(b).read_bytes(), msg)
        compared += 1
    try:
        for z in cfg:
            folder = ROOT / z['folder']
            key = (z['folder'], z['ag'], z['k'])
            if key not in bins:
                flags = ['-O3', '-std=c++17', '-fopenmp']
                if args.sanitize:
                    flags += ['-fsanitize=undefined', '-fno-sanitize-recover=all']
                b = {}
                for mode in ('produce', 'consume', 'algebra'):
                    target = temp / ('bin_' + str(len(bins)) + '_' + mode)
                    ext = [] if mode == 'algebra' else [f"-DANGULAR={z['ag']}", f"-DDEGREE_EXP={z['k']}"]
                    if mode == 'consume':
                        ext += ['-DBUILD_CONSUMER=1']
                    run(['g++', *flags, *ext, folder / 'code' / ('algebra.cpp' if mode == 'algebra' else 'angle.cpp'), '-o', target])
                    b[mode] = target
                bins[key] = b
            b = bins[key]
            name = z['stem']
            pre = temp / name
            expected = folder / 'results' / name
            for mode in ('produce', 'consume'):
                run([b[mode], mode, z['p'], folder / 'data', pre])
            for ext in ('.cells', '.branches', '.audit'):
                same(str(pre) + ext, str(pre) + '.checked' + ext, 'producer-consumer agreement ' + name + ext)
                same(str(pre) + ext, str(expected) + ext, 'expected certificate bytes ' + name + ext)
            run([b['algebra'], z['p'], z['k'], folder / 'data', str(pre) + '.tensor'])
            same(str(pre) + '.tensor', str(expected) + '.tensor', 'full algebra bytes ' + name)
            run([sys.executable, '-B', folder / 'code/emit_profile.py', z['p'], z['ag'], z['k'], pre, str(pre) + '.json'])
            same(str(pre) + '.json', str(expected) + '.json', 'profile integration and witness SHA ' + name)
            row = json.loads(Path(str(pre) + '.json').read_text())
            if z['p'] > 1:
                from exact import log_lo
                need(Fraction(row['total_upper']) < Fraction(2, 85) * log_lo(z.get('base_min', row['q'])), 'strict normalized direct route ' + name)
            else:
                need(Fraction(row['total_upper']) < Fraction(z['simple_cap']), 'marked-42 simple cap ' + name)
            run([sys.executable, '-B', folder / 'code/audit_one.py', z['p'], z['ag'], '--k', z['k'], '--prefix', pre, '--output', str(pre) + '.independent.json'])
            same(str(pre) + '.independent.json', str(expected) + '.independent.json', 'second-audit output ' + name)
            if z['p'] == 13:
                run([sys.executable, '-B', folder / 'code/negative_profiles.py', '--consumer', b['consume'], '--output', temp / 'negative.json', '--log', temp / 'negative.log'])
                same(temp / 'negative.json', folder / 'results/negative_checks.json', 'negative geometric gates')
            rows.append(row)
            Path(str(pre) + '.witness').unlink()
        script = "import sys,json;from pathlib import Path;sys.path.insert(0,sys.argv[1]);from checks import input_checks;Path(sys.argv[3]).write_text(json.dumps(input_checks(Path(sys.argv[2])),indent=2)+'\\n')"
        run([sys.executable, '-B', '-c', script, ROOT / 'code', ROOT, temp / 'input.json'])
        same(temp / 'input.json', ROOT / 'results/common_input_checks.json', 'all log table and seed checks')
        run([sys.executable, '-B', ROOT / 'iterated/verify_iterated.py', '--threads', args.threads, *(['--sanitize'] if args.sanitize else []), *(['--keep-work'] if args.keep_work else [])])
        run([sys.executable, '-B', ROOT / 'code/assemble_small.py', '--check'])
        report = dict(status='PASS', mode='full_witness_rebuild', profiles=len(rows), full_mark_primes=len(full), all_prime_exponents_and_lengths=True, compared_files=compared, geometric_witnesses=sum((z['witnesses'] for z in rows)), geometric_corner_evaluations=sum((z['witnesses'] * z['angular_nodes'] * 4 for z in rows)), threads=args.threads, sanitize=args.sanitize, python_optimized=not __debug__, runtime_seconds=round(time.monotonic() - start, 3), second_geometric_audit_exhaustive=False)
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
