"""Reconsume every K6084 witness at AG8; coefficient search is unchanged."""
from pathlib import Path
import hashlib, json, os, subprocess, sys, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
from generate import sha
pre = ROOT / 'raw/q78_K6084'
out = ROOT / 'results/q78_K6084'
old = json.loads((ROOT / 'diagnostics/ag4_K6084/q78_K6084.json').read_text())
need(sha(str(pre) + '.witness') == old['witness_sha256'], 'unchanged full integer stream')
contract = {'q': 78, 'index': 6084, 'angular_grid': 8, 'witness_sha256': old['witness_sha256'], 'witness_bytes': old['witness_bytes'], 'source_sha256': {str(p.relative_to(ROOT)): sha(p) for folder in ('code', 'data') for p in sorted((ROOT / folder).iterdir()) if p.is_file()}, 'old_profile': 'diagnostics/ag4_K6084/q78_K6084.json', 'new_penalty': '1/256', 'old_penalty': '1/64', 'all_records_consumed_again': True}
path = ROOT / 'logs/q78_K6084_refinement_CONTRACT.json'
path.write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
commands = []
start = time.monotonic()
env = {**os.environ, 'CSTAR_THREADS': '3', 'PYTHONDONTWRITEBYTECODE': '1'}

def run(cmd):
    t = time.monotonic()
    r = subprocess.run(list(map(str, cmd)), env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(r.stdout, end='', flush=True)
    commands.append({'cmd': list(map(str, cmd)), 'exit_code': r.returncode, 'seconds': round(time.monotonic() - t, 3), 'transcript': r.stdout})
    need(r.returncode == 0, 'refinement command')
run([ROOT / 'bin/consume8', 'consume', 1, 78, ROOT / 'data', pre])
for ext in ('.cells', '.branches', '.audit'):
    data = Path(str(pre) + '.checked' + ext).read_bytes()
    Path(str(pre) + ext).write_bytes(data)
    Path(str(out) + ext).write_bytes(data)
run([sys.executable, '-B', ROOT / 'code/emit_profile.py', 1, 8, 78, pre, str(out) + '.json'])
run([sys.executable, '-B', ROOT / 'code/audit_one.py', 1, 8, '--k', 78, '--prefix', pre, '--output', str(out) + '.independent.json'])
v = json.loads(Path(str(out) + '.json').read_text())
need(v['witness_sha256'] == old['witness_sha256'], 'exact same coefficient stream')
receipt = {'status': 'PASS', 'mode': 'full_angular_reconsumption_of_existing_integer_stream', 'contract_sha256': sha(path), 'commands': commands, 'seconds': round(time.monotonic() - start, 3), 'witness_bytes': v['witness_bytes'], 'witness_sha256': v['witness_sha256'], 'total_upper': v['total_upper'], 'upper_decimal': v['upper_decimal'], 'producer_search_does_not_use_angular_grid': True, 'fresh_producer_replay_pending': True}
(ROOT / 'logs/q78_K6084_refinement.json').write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
print(json.dumps(receipt, indent=2), flush=True)
