"""One primitive HNF transfer of a predeclared index with full mark2310."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
import argparse, json, os, resource, subprocess, sys, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need, log_lo, out_decimal
from core import mark, transfer_index
ap = argparse.ArgumentParser()
ap.add_argument('--contract', type=Path, required=True)
a = ap.parse_args()
cfg = json.loads(a.contract.read_text())
mult = cfg['index_multiplier']
p = cfg['p_parameter']
ag = cfg['angular_parameter']
q, nl = mark(p)
K = transfer_index(p, 1, mult)
need(q == 2310 and K == cfg['transfer_index'] and (nl == cfg['labels'] == 6912), 'full marked index contract')
need(cfg['primitive_branches'] == nl * mult and cfg['admissible_branches_per_parent'] == K, 'multiplicities')
need(cfg['witnesses'] == 2688 * nl * mult and cfg['witness_bytes'] == 8 * cfg['witnesses'], 'stream size')
need(cfg['threads'] == 2 and cfg['sanitize'] and (K != q * q), 'bounded initial compute contract')
stem = f'q{q}_K{K}_AG{ag}'
pre = ROOT / 'work_mixed' / stem
(ROOT / 'work_mixed').mkdir(exist_ok=True)
bindir = ROOT / f'bin_m{mult}_a{ag}'
bindir.mkdir(exist_ok=True)
need(not (ROOT / 'results' / (stem + '.json')).exists(), 'never overwrite an accepted profile')
env = {**os.environ, 'CSTAR_THREADS': '2', 'PYTHONDONTWRITEBYTECODE': '1'}
flags = ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all', f'-DINDEX_MULTIPLIER={mult}']
start = time.monotonic()
deadline = start + 900
inputs = {str(x.relative_to(ROOT)): sha256(x.read_bytes()).hexdigest() for dr in ('code', 'data') for x in sorted((ROOT / dr).glob('*')) if x.is_file()}
contract_hash = sha256(a.contract.read_bytes()).hexdigest()

def run(cmd):
    print('RUN', ' '.join(map(str, cmd)), flush=True)
    remaining = deadline - time.monotonic()
    need(remaining > 0, 'initial15-minute limit')
    subprocess.run(list(map(str, cmd)), check=True, env=env, timeout=remaining)
for mode in ('produce', 'consume'):
    target = bindir / mode
    run(['g++', *flags, '-DDEGREE_EXP=1', f'-DANGULAR={ag}', *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), ROOT / 'code/angle.cpp', '-o', target])
    run([target, mode, p, ROOT / 'data', pre])
for ext in ('.cells', '.branches', '.audit'):
    need(Path(str(pre) + ext).read_bytes() == Path(str(pre) + '.checked' + ext).read_bytes(), 'producer-consumer identical ' + ext)
run(['g++', *flags, ROOT / 'code/algebra.cpp', '-o', bindir / 'algebra'])
run([bindir / 'algebra', p, 1, ROOT / 'data', str(pre) + '.tensor'])
run([sys.executable, '-B', ROOT / 'code/emit_profile.py', p, ag, 1, pre, str(pre) + '.json', '--multiplier', mult])
z = json.loads(Path(str(pre) + '.json').read_text())
need(z['transfer_index'] == K and z['transfer_multiplier'] == mult and (z['k'] is None), 'explicit non-power metadata')
print('OBSERVED', json.dumps({key: z[key] for key in ('q', 'transfer_index', 'transfer_multiplier', 'upper_decimal', 'total_upper', 'direct_radical_gate', 'normalized_upper')}), flush=True)
run([sys.executable, '-B', ROOT / 'code/audit_one.py', p, ag, '--k', 1, '--multiplier', mult, '--prefix', pre, '--output', str(pre) + '.independent.json'])
need(inputs == {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in inputs}, 'frozen source/data bytes')
need(contract_hash == sha256(a.contract.read_bytes()).hexdigest(), 'frozen index contract')
peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024
need(peak < cfg['memory_limit_bytes'], 'measured child peak below2GiB')
for ext in ('.cells', '.branches', '.audit', '.tensor', '.json', '.independent.json'):
    (ROOT / 'results' / (stem + ext)).write_bytes(Path(str(pre) + ext).read_bytes())
report = dict(status='PASS', stem=stem, mark=q, transfer_index=K, transfer_multiplier=mult, angular_parameter=ag, cap=z['total_upper'], upper_decimal=z['upper_decimal'], strict_delta_gate=z['direct_radical_gate'], required_divisor_of_M=K, all_positive_multiples_including_extra_primes=True, witnesses=z['witnesses'], witness_bytes=z['witness_bytes'], witness_sha256=z['witness_sha256'], threads=2, sanitize=True, compile_flags=flags, peak_child_memory_bytes=peak, runtime_seconds=round(time.monotonic() - start, 3), contract_sha256=contract_hash, sources_and_data_unchanged=True, input_sha256=inputs)
(ROOT / 'logs' / (stem + '_generation_report.json')).write_text(json.dumps(report, indent=2) + '\n', newline='\n')
print('ACCEPTED', json.dumps(report), flush=True)
