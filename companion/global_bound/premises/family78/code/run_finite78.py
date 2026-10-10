"""Rebuild and consume every witness in a frozen inventory, retaining SHA-256 receipts."""
from pathlib import Path
import argparse, csv, hashlib, json, os, subprocess, sys, time

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work-dir', required=True)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--compare', action='store_true')
    a = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    work = Path(a.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=False)
    if not 1 <= a.threads <= 4:
        raise ValueError('threads must be 1..4')
    domain = root / 'results/finite_domain.txt'
    source = root / 'code/finite.cpp'
    expected = root / 'results/finite_bounds.csv'
    binary = work / 'finite78'
    output = work / 'finite_bounds.csv'
    flags = ['-O2', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all']
    compile_command = ['g++', *flags, str(source), '-o', str(binary)]
    run_command = [str(binary), 'run', str(domain), str(root / 'data'), str(work), str(output), str(a.threads)]
    paths = [Path(__file__).resolve(), source, domain, root / 'data/log_table.txt', root / 'data/log_seeds.tsv']
    contract = dict(source_sha256=sha(source), inputs={str(p.relative_to(root)): sha(p) for p in paths}, compile=compile_command, run=run_command, compare_expected=a.compare, expected_sha256=sha(expected) if a.compare else None, source_bounds=dict(modulus_min=1326, modulus_exclusive_max=1000000, prime_min=17, prime_exclusive_max=350), generated_streams_retained=True)
    (work / 'CONTRACT.json').write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    with (work / 'compile.log').open('w') as log:
        subprocess.run(compile_command, stdout=log, stderr=subprocess.STDOUT, check=True)
    start = time.monotonic()
    with (work / 'run.log').open('w') as log:
        p = subprocess.run(run_command, stdout=log, stderr=subprocess.STDOUT)
    elapsed = time.monotonic() - start
    receipt = dict(status='PASS' if p.returncode == 0 else 'FAIL', returncode=p.returncode, elapsed_seconds=elapsed, binary_sha256=sha(binary), contract_sha256=sha(work / 'CONTRACT.json'))
    if p.returncode:
        (work / 'RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
        print((work / 'run.log').read_text(), flush=True)
        raise RuntimeError('finite witness computation failed; logs retained')
    rr = list(csv.DictReader(output.open()))
    specified = [tuple(map(int, l.split())) for l in domain.read_text().splitlines()]
    if [(int(r['M']), int(r['p'])) for r in rr] != specified:
        raise ValueError('output domain mismatch')
    if len(set((m for m, p in specified))) != len(specified):
        raise ValueError('duplicate modulus')
    manifest = []
    for r in rr:
        M = int(r['M'])
        p = int(r['p'])
        n = M // 2 + 1
        path = work / f'M{M}.bin'
        if path.stat().st_size != 8 * n:
            raise ValueError('stream length')
        if int(r['representatives']) != n or int(r['states']) != M:
            raise ValueError('stream mass')
        if int(r['mean_numerator']) != 2 * int(r['sum_log_upper_ticks']) - M * int(r['log_M_lower_ticks']):
            raise ValueError('exact mean assembly')
        if 85 * max(0, int(r['mean_numerator'])) >= 4 * M * int(r['log_M_lower_ticks']):
            raise ValueError('finite gate')
        manifest.append(dict(M=M, p=p, records=n, bytes=8 * n, sha256=sha(path)))
    unchanged = {str(x.relative_to(root)): sha(x) for x in paths}
    if unchanged != contract['inputs']:
        raise ValueError('scientific input changed during computation')
    receipt.update(moduli=len(rr), states=sum((int(r['M']) for r in rr)), representatives=sum((int(r['representatives']) for r in rr)), stream_bytes=sum((r['bytes'] for r in manifest)), all_upper_means_negative=all((int(r['mean_numerator']) < 0 for r in rr)), refined_moduli=sum((int(r['stencil']) for r in rr)), output_sha256=sha(output), inputs_unchanged=True)
    if a.compare:
        if output.read_bytes() != expected.read_bytes():
            raise ValueError('replayed finite bounds differ from expected')
        receipt['expected_byte_match'] = True
        frozen = json.loads((root / 'results/finite_stream_manifest.json').read_text())
        if manifest != frozen:
            raise ValueError('replayed raw stream hashes differ from expected')
        receipt['all_stream_sha256_match'] = True
    else:
        expected.write_bytes(output.read_bytes())
        (root / 'results/finite_stream_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', newline='\n')
        (root / 'results/finite_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
        for n in ['CONTRACT.json', 'compile.log', 'run.log']:
            (root / 'logs' / ('finite_' + n)).write_bytes((work / n).read_bytes())
    (work / 'RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
    print(json.dumps(receipt, indent=2), flush=True)
if __name__ == '__main__':
    main()
