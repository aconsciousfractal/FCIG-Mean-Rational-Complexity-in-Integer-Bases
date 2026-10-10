"""Record a precontract and actual execution for the sixteen geometric gates."""
from pathlib import Path
import subprocess, json, time, hashlib, sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
BINDIR = Path(json.loads((ROOT / 'logs/build_RECEIPT.json').read_text())['binary_directory'])
names = ['code/negative_profiles.py', 'code/assemble.py', 'code/core.py', 'code/exact.py', 'code/kernel.hpp', 'code/angle.cpp', 'data/grid.tsv', 'data/log_table.txt', 'results/q2622.branches', 'results/q2622.audit', 'run_negative.py']
src = {name: sha(ROOT / name) for name in names}
c = {'source_sha256': src, 'consumer_binary_sha256': sha(BINDIR / 'consume4'), 'ubsan': True, 'expected_cases': 18}
cp = ROOT / 'logs/negative_CONTRACT.json'
cp.write_text(json.dumps(c, indent=2) + '\n', newline='\n')
cmd = [sys.executable, '-B', str(ROOT / 'code/negative_profiles.py'), '--consumer', str(BINDIR / 'consume4'), '--output', str(ROOT / 'results/negative_checks.json'), '--log', str(ROOT / 'logs/negative_checks.log')]
t = time.monotonic()
r = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
need(r.returncode == 0, 'all negative gates completed')
print(r.stdout, end='', flush=True)
for name, digest in src.items():
    need(sha(ROOT / name) == digest, 'negative executed input unchanged')
rec = {'status': 'PASS', 'contract_sha256': sha(cp), 'argv': cmd, 'exit_code': r.returncode, 'stdout': r.stdout, 'seconds': round(time.monotonic() - t, 3), 'source_sha256': src, 'result_sha256': sha(ROOT / 'results/negative_checks.json'), 'diagnostics_sha256': sha(ROOT / 'logs/negative_checks.log')}
(ROOT / 'logs/negative_RECEIPT.json').write_text(json.dumps(rec, indent=2) + '\n', newline='\n')
