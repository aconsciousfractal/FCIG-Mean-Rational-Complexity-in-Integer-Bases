"""Adversarial fixtures for both full marks, using the actual physical consumer."""
from pathlib import Path
import argparse, hashlib, json, subprocess
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def run(work):
    assert not work.is_relative_to(ROOT)
    work.mkdir(parents=True, exist_ok=False)
    inputs = {p: sha(ROOT / p) for p in ['code/negative_harness.cpp', 'code/negative.py', 'code/finite.cpp', 'data/log_table.txt', 'data/log_seeds.tsv']}
    command = ['g++', '-O2', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all', str(ROOT / 'code/negative_harness.cpp'), '-o', str(work / 'negative')]
    contract = dict(status='FROZEN_BEFORE_RUN', inputs=inputs, compile=command, expected_negative=37, expected_positive=4)
    (work / 'CONTRACT.json').write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    with (work / 'compile.log').open('w') as f:
        subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True)
    r = subprocess.run([str(work / 'negative'), str(ROOT / 'data'), str(work / 'fixtures')], capture_output=True, text=True)
    output = r.stdout + r.stderr
    (work / 'run.log').write_text(output, newline='\n')
    r.check_returncode()
    assert output.count('REJECT ') == 37
    assert 'PASS all37 bounded negative checks and4 positive consumer fixtures' in output
    assert inputs == {p: sha(ROOT / p) for p in inputs}
    result = dict(status='PASS', negative_fixtures=37, positive_consumer_fixtures=4, radicals=[3570, 3990], valid_records=[1786, 1996], sanitize=True, inputs_unchanged=True, contract_sha256=sha(work / 'CONTRACT.json'), run_log_sha256=sha(work / 'run.log'))
    (work / 'RECEIPT.json').write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--work-dir', type=Path, required=True)
    a = p.parse_args()
    run(a.work_dir.resolve())
