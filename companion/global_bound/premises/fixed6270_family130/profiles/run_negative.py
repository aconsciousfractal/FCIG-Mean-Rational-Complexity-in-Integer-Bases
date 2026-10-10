"""Freeze all fixture inputs before exercising the actual compiled consumer."""
from pathlib import Path
import argparse, json, hashlib, subprocess, sys, time, os
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
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw-work', type=Path, required=True)
    ap.add_argument('--consumer-binary', type=Path)
    args = ap.parse_args()
    args.raw_work = args.raw_work.resolve()
    need(args.raw_work.is_dir() and (not args.raw_work.is_relative_to(R)), 'existing external temporary directory')
    prefix = args.raw_work / 'q6270_K6270'
    need(Path(str(prefix) + '.branches').read_bytes() == (R / 'results/q6270_K6270.branches').read_bytes(), 'fixture raw branch identity')
    names = ['code/negative_profiles.py', 'code/assemble.py', 'code/core.py', 'code/exact.py', 'data/grid.tsv', 'data/log_table.txt', 'results/q6270_K6270.branches', 'run_negative.py']
    src = {s: sha(R / s) for s in names}
    build = json.loads((R / 'logs/build_RECEIPT.json').read_text())
    binpath = (args.consumer_binary or Path(build['binary_directory']) / 'consume8').resolve()
    need(sha(binpath) == json.loads((R / 'logs/build_RECEIPT.json').read_text())['binaries']['consume8'], 'exact compiled consumer bytes')
    c = {'source_sha256': src, 'consumer_binary_sha256': sha(binpath), 'positive_stream_generation_receipt_sha256': sha(R / 'logs/q6270_K6270_generation.json'), 'expected_cases': EXPECTED_CASES, 'ubsan': True}
    cp = R / 'logs/negative_CONTRACT.json'
    need(not cp.exists(), 'fresh negative precontract')
    cp.write_text(json.dumps(c, indent=2) + '\n', newline='\n')
    t = time.monotonic()
    z = subprocess.run([sys.executable, '-B', str(R / 'code/negative_profiles.py'), '--consumer', str(binpath), '--positive-prefix', str(prefix), '--output', str(R / 'results/negative_checks.json'), '--log', str(R / 'logs/negative_checks.log')], text=True, capture_output=True, env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
    print(z.stdout + z.stderr, flush=True)
    need(z.returncode == 0, 'all fixture gates passed')
    need(all((sha(R / s) == h for s, h in src.items())), 'fixture inputs unchanged')
    (R / 'logs/negative_RECEIPT.json').write_text(json.dumps({'status': 'PASS', 'exit_code': z.returncode, 'contract_sha256': sha(cp), 'source_sha256': src, 'result_sha256': sha(R / 'results/negative_checks.json'), 'diagnostics_sha256': sha(R / 'logs/negative_checks.log'), 'seconds': round(time.monotonic() - t, 3), 'transcript': z.stdout + z.stderr}, indent=2) + '\n', newline='\n')
if __name__ == '__main__':
    main()
