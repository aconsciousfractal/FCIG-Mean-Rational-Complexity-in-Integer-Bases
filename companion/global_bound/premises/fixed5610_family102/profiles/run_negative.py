"""Freeze all fixture inputs before exercising the actual compiled consumer."""
from pathlib import Path
import json, hashlib, subprocess, sys, time, os
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'code'))
from negative_profiles import EXPECTED_CASES
from exact import need

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def main():
    names = ['code/negative_profiles.py', 'code/assemble.py', 'code/core.py', 'code/exact.py', 'data/grid.tsv', 'data/log_table.txt', 'results/q5610_K5610.branches', 'run_negative.py']
    src = {s: sha(R / s) for s in names}
    binpath = R / 'bin/consume8'
    c = {'source_sha256': src, 'consumer_binary_sha256': sha(binpath), 'positive_stream_generation_receipt_sha256': sha(R / 'logs/q5610_K5610_generation.json'), 'expected_cases': EXPECTED_CASES, 'ubsan': True}
    cp = R / 'logs/negative_CONTRACT.json'
    need(not cp.exists(), 'fresh negative precontract')
    cp.write_text(json.dumps(c, indent=2) + '\n', newline='\n')
    t = time.monotonic()
    z = subprocess.run([sys.executable, '-B', str(R / 'code/negative_profiles.py'), '--consumer', str(binpath), '--positive-prefix', str(R / 'raw/q5610_K5610'), '--output', str(R / 'results/negative_checks.json'), '--log', str(R / 'logs/negative_checks.log')], text=True, capture_output=True, env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
    print(z.stdout + z.stderr, flush=True)
    need(z.returncode == 0, 'all fixture gates passed')
    need(all((sha(R / s) == h for s, h in src.items())), 'fixture inputs unchanged')
    (R / 'logs/negative_RECEIPT.json').write_text(json.dumps({'status': 'PASS', 'exit_code': z.returncode, 'contract_sha256': sha(cp), 'source_sha256': src, 'result_sha256': sha(R / 'results/negative_checks.json'), 'diagnostics_sha256': sha(R / 'logs/negative_checks.log'), 'seconds': round(time.monotonic() - t, 3), 'transcript': z.stdout + z.stderr}, indent=2) + '\n', newline='\n')
if __name__ == '__main__':
    main()
