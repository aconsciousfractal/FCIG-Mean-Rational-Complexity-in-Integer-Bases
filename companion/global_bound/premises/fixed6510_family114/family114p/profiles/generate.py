"""Generate and independently consume a new114 marked profile in a fresh prefix."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time, tempfile
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stem', required=True)
    ap.add_argument('--threads', type=int, default=4)
    args = ap.parse_args()
    need(args.threads in (1, 2, 3, 4), 'one to four threads')
    cfg = json.loads((ROOT / 'profiles_plan.json').read_text())
    z = next((v for v in cfg if v['stem'] == args.stem))
    pre = Path(tempfile.mkdtemp(prefix='p72-114-production-')) / z['stem']
    out = ROOT / 'results' / z['stem']
    need(not Path(str(pre) + '.witness').exists(), 'fresh output prefix')
    executed = ['code/angle.cpp', 'code/kernel.hpp', 'code/search.hpp', 'code/algebra.cpp', 'code/emit_profile.py', 'code/audit_one.py', 'code/checks.py', 'code/assemble.py', 'code/core.py', 'code/exact.py', 'data/grid.tsv', 'data/log_table.txt', 'data/log_seeds.tsv', 'generate.py']
    source = {name: sha(ROOT / name) for name in executed}
    build = json.loads((ROOT / 'logs/build_RECEIPT.json').read_text())
    need(build['status'] == 'PASS', 'compiled binaries receipt')
    BINDIR = Path(build['binary_directory'])
    for name, digest in build['source_sha256'].items():
        need(sha(ROOT / name) == digest, 'compiled current C++ source ' + name)
    names = ['algebra', f"produce{z['ag']}", f"consume{z['ag']}"]
    for name in names:
        need(sha(BINDIR / name) == build['binaries'][name], 'compiled binary hash ' + name)
    contract = {'profile': z, 'source_sha256': source, 'build_receipt_sha256': sha(ROOT / 'logs/build_RECEIPT.json'), 'binary_sha256': {name: sha(BINDIR / name) for name in names}, 'commands_use_ubsan': True, 'threads': args.threads, 'raw_prefix': str(pre), 'scope': 'all M divisible by q*multiplier; explicit marked denominator gcd(d,q)=1; CDF3/2 at seed L'}
    cp = ROOT / 'logs' / (z['stem'] + '_CONTRACT.json')
    need(not cp.exists(), 'fresh precontract')
    cp.write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    env = {**os.environ, 'CSTAR_THREADS': str(args.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    commands = []
    start = time.monotonic()

    def run(cmd):
        tick = time.monotonic()
        exe = Path(str(cmd[0]))
        if exe.is_file() and exe.parent == BINDIR:
            exe.chmod(448)
        r = subprocess.run(list(map(str, cmd)), env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(r.stdout, end='', flush=True)
        commands.append({'cmd': list(map(str, cmd)), 'exit_code': r.returncode, 'seconds': round(time.monotonic() - tick, 3), 'transcript': r.stdout})
        need(r.returncode == 0, 'command succeeded')
    try:
        run([BINDIR / f"produce{z['ag']}", 'produce', z['p'], z['multiplier'], ROOT / 'data', pre])
        witness_path = Path(str(pre) + '.witness')
        with witness_path.open('rb') as synced:
            os.fsync(synced.fileno())
        producer_sha = sha(witness_path)
        producer_bytes = witness_path.stat().st_size
        run([BINDIR / f"consume{z['ag']}", 'consume', z['p'], z['multiplier'], ROOT / 'data', pre])
        need(sha(witness_path) == producer_sha and witness_path.stat().st_size == producer_bytes, 'raw unchanged throughout complete consumption')
        for ext in ('.cells', '.branches', '.audit'):
            a = Path(str(pre) + ext)
            b = Path(str(pre) + '.checked' + ext)
            need(a.read_bytes() == b.read_bytes(), 'complete consumer byte agreement ' + ext)
            Path(str(out) + ext).write_bytes(a.read_bytes())
        run([BINDIR / 'algebra', z['p'], z['multiplier'], ROOT / 'data', str(pre) + '.tensor'])
        Path(str(out) + '.tensor').write_bytes(Path(str(pre) + '.tensor').read_bytes())
        run([sys.executable, '-B', ROOT / 'code/emit_profile.py', z['p'], z['ag'], z['multiplier'], pre, str(out) + '.json'])
        run([sys.executable, '-B', ROOT / 'code/audit_one.py', z['p'], z['ag'], '--k', z['multiplier'], '--prefix', pre, '--output', str(out) + '.independent.json'])
        for name, digest in source.items():
            need(sha(ROOT / name) == digest, 'executed input unchanged ' + name)
        v = json.loads(Path(str(out) + '.json').read_text())
        print(json.dumps({k: v[k] for k in ('p', 'q', 'transfer_index', 'upper_decimal', 'total_upper', 'witnesses', 'witness_bytes', 'witness_sha256')}, indent=2), flush=True)
        receipt = {'status': 'PASS', 'profile': z, 'contract_sha256': sha(cp), 'witness_bytes': v['witness_bytes'], 'witness_sha256': v['witness_sha256'], 'all_geometry_recomputed': True, 'producer_consumer_agreement': True, 'executed_inputs_unchanged': True, 'seconds': round(time.monotonic() - start, 3), 'commands': commands}
        (ROOT / 'logs' / (z['stem'] + '_generation.json')).write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
    except BaseException:
        (ROOT / 'logs' / (z['stem'] + '_FAILED.json')).write_text(json.dumps({'commands': commands, 'seconds': round(time.monotonic() - start, 3)}, indent=2) + '\n', newline='\n')
        raise
if __name__ == '__main__':
    main()
