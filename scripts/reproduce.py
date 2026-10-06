"""Portable command entry point; Python standard library only."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('smoke', 'retained', 'semiprime', 'replay'))
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--all', action='store_true', help='reconstruct all 4596 semiprime envelopes explicitly')
    parser.add_argument('--batch', type=int)
    parser.add_argument('--script')
    parser.add_argument('--work')
    args = parser.parse_args()
    jobs = {'smoke': 'smoke.py', 'retained': 'verify_retained.py', 'semiprime': 'reconstruct_semiprime.py', 'replay': 'replay.py'}
    argv = [sys.executable, '-B'] + (['-O'] if args.optimized else [])
    argv += [str(ROOT / 'companion/programs' / jobs[args.mode])]
    if args.all:
        if args.mode != 'semiprime':
            parser.error('--all belongs to semiprime mode')
        argv.append('--all')
    if args.mode == 'replay':
        if args.batch is None or not args.script or not args.work:
            parser.error('replay requires --batch, --script and --work')
        argv += [str(args.batch), args.script, '--work', args.work]
        if args.optimized:
            argv.append('--optimized')
    elif any(x is not None for x in (args.batch, args.script, args.work)):
        parser.error('--batch, --script and --work belong to replay mode')
    raise SystemExit(subprocess.run(argv, cwd=ROOT / 'companion', check=False).returncode)


if __name__ == '__main__':
    main()
