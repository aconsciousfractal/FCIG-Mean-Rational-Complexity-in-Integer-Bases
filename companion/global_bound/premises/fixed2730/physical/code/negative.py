"""Re-execute tiny adversarial fixtures against the actual physical consumer."""
from pathlib import Path
import argparse, hashlib, json, subprocess
ROOT = Path(__file__).resolve().parents[1]

def run(work):
    work.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'code/negative_harness.cpp'
    binary = work / 'negative'
    compile = ['g++', '-O2', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all', str(source), '-o', str(binary)]
    inputs = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in ('code/negative_harness.cpp', 'code/negative.py', 'code/finite.cpp', 'data/log_table.txt', 'data/log_seeds.tsv')}
    contract = dict(status='FROZEN_BEFORE_RUN', inputs=inputs, compile=compile)
    (work / 'CONTRACT.json').write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    c = subprocess.run(compile, capture_output=True, text=True)
    (work / 'compile.log').write_text(c.stdout + c.stderr, newline='\n')
    c.check_returncode()
    r = subprocess.run([str(binary), str(ROOT / 'data'), str(work / 'fixtures')], capture_output=True, text=True)
    output = r.stdout + r.stderr
    (work / 'run.log').write_text(output, newline='\n')
    r.check_returncode()
    assert output.count('REJECT ') == 15
    assert 'PASS all15 bounded negative checks and2 positive consumer fixtures' in output
    assert inputs == {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in inputs}
    result = dict(status='PASS', negative_fixtures=15, positive_consumer_fixtures=2, valid_modulus=2730, valid_records=1366, sanitize=True, inputs_unchanged=True, large_raw_replay=False)
    (work / 'RECEIPT.json').write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--work-dir', type=Path, required=True)
    a = p.parse_args()
    run(a.work_dir.resolve())
