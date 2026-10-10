"""Fresh direct marked-6 anchor, producer and consumer compiled with UBSan."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, tempfile, time
ROOT = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def need(b, m):
    if not b:
        raise ArithmeticError(m)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--compare', action='store_true')
    ap.add_argument('--threads', type=int, default=1, choices=(1, 2, 3, 4))
    args = ap.parse_args()
    work = Path(tempfile.mkdtemp(prefix='p72-q6-anchor-'))
    start = time.monotonic()
    flags = ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all']
    env = {**os.environ, 'CSTAR_THREADS': str(args.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    inputs = {str(p.relative_to(ROOT)): sha(p) for folder in ['code', 'data'] for p in sorted((ROOT / folder).rglob('*')) if p.is_file()}
    commands = []

    def run(cmd):
        cmd = list(map(str, cmd))
        commands.append(cmd)
        with (work / 'run.log').open('a') as log:
            log.write('RUN ' + ' '.join(cmd) + '\n')
            log.flush()
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    producer, consumer, algebra = [work / n for n in ['produce', 'consume', 'algebra']]
    for mode, target in [('produce', producer), ('consume', consumer)]:
        run(['g++', *flags, '-DANGULAR=16', '-DDEGREE_EXP=1', *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), ROOT / 'code/angle.cpp', '-o', target])
    run(['g++', *flags, ROOT / 'code/algebra.cpp', '-o', algebra])
    pre = work / 'q6'
    run([producer, 'produce', 1, ROOT / 'data', pre])
    run([consumer, 'consume', 1, ROOT / 'data', pre])
    for ext in ['.cells', '.branches', '.audit']:
        need(Path(str(pre) + ext).read_bytes() == Path(str(pre) + '.checked' + ext).read_bytes(), 'producer/consumer byte identity ' + ext)
    run([algebra, 1, 1, ROOT / 'data', str(pre) + '.tensor'])
    run([sys.executable, '-B', ROOT / 'code/emit_profile.py', 1, 16, 1, pre, str(pre) + '.json'])
    profile = json.loads(Path(str(pre) + '.json').read_text())
    need(profile['q'] == 6 and profile['transfer_index'] == 6 and (profile['labels'] == 12) and (profile['witnesses'] == 32256) and (profile['witness_bytes'] == 258048), 'q6 dimensions')
    destinations = []
    for ext in ['.cells', '.branches', '.audit', '.tensor', '.json']:
        out = ROOT / 'results' / ('q6' + ext)
        data = Path(str(pre) + ext).read_bytes()
        if args.compare:
            need(out.read_bytes() == data, 'frozen output ' + ext)
        else:
            out.write_bytes(data)
        destinations.append({'path': str(out.relative_to(ROOT)), 'sha256': sha(out), 'bytes': out.stat().st_size})
    need(all((sha(ROOT / name) == h for name, h in inputs.items())), 'unchanged source and input')
    result = {'status': 'PASS', 'q': 6, 'k': 1, 'angular_nodes': 33, 'cells': 2688, 'labels': 12, 'geometric_witnesses': 32256, 'corner_evaluations_per_pass': 4257792, 'sanitize': True, 'threads': args.threads, 'compare': args.compare, 'all_positive_multiples_of_6': True, 'source_inputs': inputs, 'outputs': destinations, 'profile': profile, 'runtime_seconds': round(time.monotonic() - start, 3), 'raw_work_directory': str(work)}
    (ROOT / 'logs' / ('replay_receipt.json' if args.compare else 'generation_receipt.json')).write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    (ROOT / 'logs' / ('replay.log' if args.compare else 'generation.log')).write_bytes((work / 'run.log').read_bytes())
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_inputs', 'outputs', 'profile']}, indent=2))
    print('q6 cap', profile['total_upper'], profile['upper_decimal'], flush=True)
if __name__ == '__main__':
    main()
