"""Fresh UBSan producer/consumer generation; accepts no previous numeric caps."""
from pathlib import Path
from hashlib import sha256
import argparse, json, os, subprocess, sys, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
ap = argparse.ArgumentParser()
ap.add_argument('--plan', type=Path, default=ROOT / 'initial_plan.json')
ap.add_argument('--work', type=Path, default=ROOT / 'work_initial')
ap.add_argument('--bin-dir', type=Path, default=ROOT / 'bin')
ap.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=2)
ap.add_argument('--report', type=Path, default=ROOT / 'logs/initial_generation_report.json')
a = ap.parse_args()
plan = json.loads(a.plan.read_text())
a.work.mkdir(parents=True, exist_ok=True)
a.bin_dir.mkdir(parents=True, exist_ok=True)
env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
start = time.monotonic()
built = set()
rows = []
flags = ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all']

def run(cmd):
    print('RUN', ' '.join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), check=True, env=env)
run(['g++', *flags, ROOT / 'code/algebra.cpp', '-o', a.bin_dir / 'algebra'])
for z in plan:
    p, k, ag, stem = (z['p'], z['k'], z['ag'], z['stem'])
    pre = a.work / stem
    need(not (ROOT / 'results' / (stem + '.json')).exists(), 'do not overwrite accepted ' + stem)
    one_start = time.monotonic()
    for mode in ('produce', 'consume'):
        target = a.bin_dir / f'{mode}_{k}_{ag}'
        if (mode, k, ag) not in built:
            run(['g++', *flags, f'-DDEGREE_EXP={k}', f'-DANGULAR={ag}', *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), ROOT / 'code/angle.cpp', '-o', target])
            built.add((mode, k, ag))
        run([target, mode, p, ROOT / 'data', pre])
    for ext in ('.cells', '.branches', '.audit'):
        need(Path(str(pre) + ext).read_bytes() == Path(str(pre) + '.checked' + ext).read_bytes(), 'producer-consumer identical ' + stem + ext)
    run([a.bin_dir / 'algebra', p, k, ROOT / 'data', str(pre) + '.tensor'])
    run([sys.executable, '-B', ROOT / 'code/emit_profile.py', p, ag, k, pre, str(pre) + '.json'])
    row = json.loads(Path(str(pre) + '.json').read_text())
    print('OBSERVED', json.dumps({key: row[key] for key in ('q', 'k', 'upper_decimal', 'direct_radical_gate', 'total_upper', 'witnesses', 'witness_sha256')}), flush=True)
    run([sys.executable, '-B', ROOT / 'code/audit_one.py', p, ag, '--k', k, '--prefix', pre, '--output', str(pre) + '.independent.json'])
    for ext in ('.cells', '.branches', '.audit', '.tensor', '.json', '.independent.json'):
        (ROOT / 'results' / (stem + ext)).write_bytes(Path(str(pre) + ext).read_bytes())
    rec = {**z, 'upper': row['total_upper'], 'upper_decimal': row['upper_decimal'], 'direct_radical_gate': row['direct_radical_gate'], 'seconds': round(time.monotonic() - one_start, 3), 'witness_bytes': row['witness_bytes'], 'witness_sha256': row['witness_sha256']}
    rows.append(rec)
    a.report.write_text(json.dumps(dict(status='IN_PROGRESS', completed=rows, threads=a.threads, sanitize=True, flags=flags), indent=2) + '\n', newline='\n')
    print('ACCEPTED', json.dumps(rec), flush=True)
report = dict(status='PASS', completed=rows, threads=a.threads, sanitize=True, flags=flags, runtime_seconds=round(time.monotonic() - start, 3), plan_sha256=sha256(a.plan.read_bytes()).hexdigest())
a.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
print(json.dumps(report, indent=2), flush=True)
