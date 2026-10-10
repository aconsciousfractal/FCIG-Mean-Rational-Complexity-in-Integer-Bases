"""Produce and consume a fixed-radical2310 physical domain, with exact bounds.

Rows specify M and the base B. Normalization always uses log(B).
The N=1 exploratory inventory is a proven complement once the q2310 cap
and its large-base gate have been accepted by the final assembler.
"""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import argparse, csv, hashlib, json, os, subprocess, sys, time
sys.dont_write_bytecode = True
from exact import need, log_lo, out_decimal
ROOT = Path(__file__).resolve().parent.parent
PRIMES = (2, 3, 5, 7, 11)
Q = 2310
BMAX = 100000

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def bases(limit=BMAX):
    rows = []

    def rec(i, m, ee):
        if i == 5:
            rows.append((m, ee))
            return
        p = PRIMES[i]
        x = m * p
        e = 1
        while x < limit:
            rec(i + 1, x, ee + [e])
            x *= p
            e += 1
    rec(0, 1, [])
    rows.sort()
    independent = []
    for m in range(Q, limit, Q):
        x = m
        for p in PRIMES:
            while x % p == 0:
                x //= p
        if x == 1:
            independent.append(m)
    need([m for m, e in rows] == independent, 'two complete all-valuation inventories')
    need(len(rows) == len(set((m for m, e in rows))), 'no duplicate base')
    return rows

def inventory():
    rows = bases()
    need(len(rows) == 30 and sum((m for m, e in rows)) == 1284360, 'fixed N1 inventory')
    records = [dict(M=m, B=m, N=1, exponents=e) for m, e in rows]
    return dict(status='PASS', radical=Q, base_upper_exclusive=BMAX, all_positive_valuations=True, rows=records, states=sum((r['M'] for r in records)), representatives=sum((r['M'] // 2 + 1 for r in records)), maximum_modulus=max((r['M'] for r in records)), analytic_large_base_gate_required=True, all_N_theorem_claimed=False)

def run(work, threads, compare):
    work.mkdir(parents=True, exist_ok=False)
    rows = inventory()['rows']
    domain = ''.join((f"{r['M']} {r['B']}\n" for r in rows))
    need(domain == (ROOT / 'results/finite_domain.txt').read_text(), 'frozen physical domain')
    inputs = ['code/finite.cpp', 'code/finite2310.py', 'results/finite_domain.txt', 'data/log_table.txt', 'data/log_seeds.tsv']
    flags = ['-O2', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all']
    command = ['g++', *flags, str(ROOT / 'code/finite.cpp'), '-o', str(work / 'finite2310')]
    execute = [str(work / 'finite2310'), 'run', str(ROOT / 'results/finite_domain.txt'), str(ROOT / 'data'), str(work), str(work / 'finite_bounds.csv'), str(threads)]
    contract = dict(status='FROZEN_BEFORE_RUN', inputs={n: sha(ROOT / n) for n in inputs}, compile=command, run=execute, compare=compare, normalization='logB', expected_bounds_sha256=sha(ROOT / 'results/finite_bounds.csv') if compare else None)
    (work / 'CONTRACT.json').write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    with (work / 'compile.log').open('w') as f:
        subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True)
    start = time.monotonic()
    with (work / 'run.log').open('w') as f:
        subprocess.run(execute, stdout=f, stderr=subprocess.STDOUT, check=True)
    elapsed = time.monotonic() - start
    got = [{k: int(v) for k, v in r.items()} for r in csv.DictReader((work / 'finite_bounds.csv').open())]
    need([(r['M'], r['B']) for r in got] == [(r['M'], r['B']) for r in rows], 'all rows consumed')
    raw = []
    means = []
    normalized = []
    for r in got:
        M, B = (r['M'], r['B'])
        n = M // 2 + 1
        need(M % 2 == 0 and r['states'] == M and (r['representatives'] == n), 'full even-modulus mass')
        need(F(r['log_M_lower_ticks'], 10 ** 9) <= log_lo(M) and F(r['log_B_lower_ticks'], 10 ** 9) <= log_lo(B), 'independent logarithms for modulus and base')
        num = 2 * r['sum_log_upper_ticks'] - M * r['log_M_lower_ticks']
        need(num == r['mean_numerator'], 'exact mean numerator')
        upper = F(max(0, num), 2 * M * r['log_B_lower_ticks'])
        need(upper < F(2, 85), 'physical upper bound below2/85')
        means.append((F(num, 2 * M * 10 ** 9), M))
        normalized.append((upper, M))
        path = work / f'M{M}.bin'
        need(path.stat().st_size == 8 * n, 'complete binary length')
        raw.append(dict(M=M, B=B, records=n, bytes=8 * n, sha256=sha(path)))
    need(contract['inputs'] == {n: sha(ROOT / n) for n in inputs}, 'scientific input unchanged')
    result = dict(status='PASS', radical=Q, domain='all rad2310 bases B<100000, N=1', normalized_upper=str(max(normalized)[0]), normalized_upper_decimal=out_decimal(max(normalized)[0]), largest_bound_base=max(normalized)[1], worst_mean_upper=str(max(means)[0]), states=sum((r['M'] for r in got)), representatives=sum((r['representatives'] for r in got)), moduli=len(got), refined_moduli=sum((r['stencil'] for r in got)), all_N_theorem_claimed=False, binary_stream_manifest=raw, output_sha256=sha(work / 'finite_bounds.csv'))
    if compare:
        need((work / 'finite_bounds.csv').read_bytes() == (ROOT / 'results/finite_bounds.csv').read_bytes(), 'replayed physical table byte identity')
        need(result == json.loads((ROOT / 'results/finite_certificate.json').read_text()), 'replayed all raw hashes and bounds')
    else:
        (ROOT / 'results/finite_bounds.csv').write_bytes((work / 'finite_bounds.csv').read_bytes())
        (ROOT / 'results/finite_certificate.json').write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    receipt = dict(status='PASS', seconds=round(elapsed, 3), inputs_unchanged=True, compared_with_frozen_output=compare, all_raw_hashes_compared=compare, contract_sha256=sha(work / 'CONTRACT.json'), output_sha256=result['output_sha256'], moduli=len(got), states=result['states'], representatives=result['representatives'], binary_sha256=sha(work / 'finite2310'), sanitize=True)
    (work / 'RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
    print(json.dumps({**receipt, 'normalized_upper_decimal': result['normalized_upper_decimal']}, indent=2), flush=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--inventory', action='store_true')
    p.add_argument('--work-dir', type=Path)
    p.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=2)
    p.add_argument('--compare', action='store_true')
    a = p.parse_args()
    if a.inventory:
        z = inventory()
        (ROOT / 'results/finite_inventory.json').write_text(json.dumps(z, indent=2) + '\n', newline='\n')
        (ROOT / 'results/finite_domain.txt').write_text(''.join((f"{r['M']} {r['B']}\n" for r in z['rows'])), newline='\n')
        print(json.dumps({k: v for k, v in z.items() if k != 'rows'}, indent=2))
        return
    need(a.work_dir is not None, 'work-dir required')
    run(a.work_dir.resolve(), a.threads, a.compare)
if __name__ == '__main__':
    main()
