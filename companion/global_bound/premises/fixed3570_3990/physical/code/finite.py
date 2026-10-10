"""Exact physical complements for radical3570 and3990; normalization is log(B).

Generation may populate only this new working component. Comparison is read-only
and writes all regenerated streams and logs into a fresh external directory.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, csv, hashlib, json, subprocess, sys, time
sys.dont_write_bytecode = True
from exact import need, log_lo, out_decimal
ROOT = Path(__file__).resolve().parents[1]
DOMAINS = {3570: ((2, 3, 5, 7, 17), 10000, [3570, 7140]), 3990: ((2, 3, 5, 7, 19), 6000, [3990])}

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def inventory():
    domains = []
    rows = []
    for q, (primes, limit, expected) in DOMAINS.items():
        found = []

        def visit(i, m, ee):
            if i == len(primes):
                found.append((m, ee))
                return
            p = primes[i]
            v = m * p
            e = 1
            while v < limit:
                visit(i + 1, v, ee + [e])
                v *= p
                e += 1
        visit(0, 1, [])
        found.sort()
        scan = []
        for B in range(q, limit, q):
            m = B
            for p in primes:
                while m % p == 0:
                    m //= p
            if m == 1:
                scan.append(B)
        need([B for B, ee in found] == scan == expected, 'two exhaustive valuation inventories')
        rr = [dict(M=B, B=B, N=1, radical=q, exponents=ee) for B, ee in found]
        rows.extend(rr)
        domains.append(dict(radical=q, prime_support=list(primes), base_upper_exclusive=limit, all_positive_valuations=True, rows=rr))
    need(sum((r['M'] for r in rows)) == 14700, 'complete physical states')
    need(sum((r['M'] // 2 + 1 for r in rows)) == 7353, 'complete compressed representatives')
    return dict(status='PASS', domains=domains, rows=rows, states=14700, representatives=7353, all_N_theorem_claimed=False, large_base_profile_gates_required=True)

def run(work, threads, compare):
    need(not work.is_relative_to(ROOT), 'raw work outside component')
    work.mkdir(parents=True, exist_ok=False)
    inv = inventory()
    rows = inv['rows']
    need(inv == json.loads((ROOT / 'results/finite_inventory.json').read_text()), 'frozen inventory')
    domain = ''.join((f"{r['M']} {r['B']} {r['radical']}\n" for r in rows))
    need(domain == (ROOT / 'results/finite_domain.txt').read_text(), 'frozen domain order')
    inputs = ['code/finite.cpp', 'code/finite.py', 'code/exact.py', 'results/finite_domain.txt', 'results/finite_inventory.json', 'data/log_table.txt', 'data/log_seeds.tsv']
    command = ['g++', '-O2', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all', str(ROOT / 'code/finite.cpp'), '-o', str(work / 'finite')]
    execute = [str(work / 'finite'), 'run', str(ROOT / 'results/finite_domain.txt'), str(ROOT / 'data'), str(work), str(work / 'finite_bounds.csv'), str(threads)]
    contract = dict(status='FROZEN_BEFORE_RUN', inputs={p: sha(ROOT / p) for p in inputs}, compile=command, run=execute, compare=compare, normalization='logB', expected_bounds_sha256=sha(ROOT / 'results/finite_bounds.csv') if compare else None, radicals=[3570, 3990], expected_states=14700, expected_representatives=7353)
    (work / 'CONTRACT.json').write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    with (work / 'compile.log').open('w') as f:
        subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True)
    started = time.monotonic()
    with (work / 'run.log').open('w') as f:
        subprocess.run(execute, stdout=f, stderr=subprocess.STDOUT, check=True)
    elapsed = time.monotonic() - started
    got = [{k: int(v) for k, v in r.items()} for r in csv.DictReader((work / 'finite_bounds.csv').open())]
    need([(r['M'], r['B'], r['radical']) for r in got] == [(r['M'], r['B'], r['radical']) for r in rows], 'complete row identities and order')
    raw = []
    bounds = {q: [] for q in DOMAINS}
    physical_rows = []
    for r in got:
        M, B, q = (r['M'], r['B'], r['radical'])
        n = M // 2 + 1
        need(M % 2 == 0 and r['states'] == M and (r['representatives'] == n), 'full state mass')
        need(F(r['log_M_lower_ticks'], 10 ** 9) <= log_lo(M) and F(r['log_B_lower_ticks'], 10 ** 9) <= log_lo(B), 'separate rational logM and logB')
        num = 2 * r['sum_log_upper_ticks'] - M * r['log_M_lower_ticks']
        need(num == r['mean_numerator'], 'exact mean numerator')
        upper = F(max(0, num), 2 * M * r['log_B_lower_ticks'])
        need(upper < F(2, 85), 'finite normalized gate')
        bounds[q].append((upper, M))
        path = work / f'M{M}.bin'
        need(path.stat().st_size == 8 * n, 'complete stream byte count')
        raw.append(dict(M=M, B=B, radical=q, records=n, bytes=8 * n, sha256=sha(path)))
        physical_rows.append(dict(M=M, B=B, radical=q, normalized_upper=str(upper), normalized_upper_decimal=out_decimal(upper), mean_upper=str(F(num, 2 * M * 10 ** 9)), stencil=r['stencil']))
    sectors = [dict(radical=q, base_upper_exclusive=DOMAINS[q][1], moduli=len(bounds[q]), normalized_upper=str(max(bounds[q])[0]), normalized_upper_decimal=out_decimal(max(bounds[q])[0]), largest_bound_base=max(bounds[q])[1]) for q in DOMAINS]
    need(contract['inputs'] == {p: sha(ROOT / p) for p in inputs}, 'scientific inputs unchanged')
    result = dict(status='PASS', radicals=[3570, 3990], domains=sectors, rows=physical_rows, moduli=3, states=14700, representatives=7353, refined_moduli=sum((r['stencil'] for r in got)), all_N_theorem_claimed=False, normalization='mu(B^N)/log(B)', binary_stream_manifest=raw, output_sha256=sha(work / 'finite_bounds.csv'))
    if compare:
        need((work / 'finite_bounds.csv').read_bytes() == (ROOT / 'results/finite_bounds.csv').read_bytes(), 'physical table byte identity')
        need(result == json.loads((ROOT / 'results/finite_certificate.json').read_text()), 'all raw hashes and exact bounds regenerated')
    else:
        (ROOT / 'results/finite_bounds.csv').write_bytes((work / 'finite_bounds.csv').read_bytes())
        (ROOT / 'results/finite_certificate.json').write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    receipt = dict(status='PASS', seconds=round(elapsed, 3), inputs_unchanged=True, compared_with_frozen_output=compare, all_raw_hashes_compared=compare, contract_sha256=sha(work / 'CONTRACT.json'), output_sha256=result['output_sha256'], moduli=3, states=14700, representatives=7353, binary_sha256=sha(work / 'finite'), sanitize=True)
    (work / 'RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
    print(json.dumps({**receipt, 'domains': sectors}, indent=2), flush=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--inventory', action='store_true')
    p.add_argument('--work-dir', type=Path)
    p.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=1)
    p.add_argument('--compare', action='store_true')
    a = p.parse_args()
    if a.inventory:
        z = inventory()
        (ROOT / 'results/finite_inventory.json').write_text(json.dumps(z, indent=2) + '\n', newline='\n')
        (ROOT / 'results/finite_domain.txt').write_text(''.join((f"{r['M']} {r['B']} {r['radical']}\n" for r in z['rows'])), newline='\n')
        print(json.dumps(z, indent=2))
        return
    need(a.work_dir is not None, 'work-dir required')
    run(a.work_dir.resolve(), a.threads, a.compare)
if __name__ == '__main__':
    main()
