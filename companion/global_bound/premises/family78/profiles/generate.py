"""Generate and independently consume new 78 marked profiles in a fresh prefix."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time
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
    ap.add_argument('--threads', type=int, default=3)
    args = ap.parse_args()
    cfg = json.loads((ROOT / 'profiles_plan.json').read_text())
    z = next((v for v in cfg if v['stem'] == args.stem))
    pre = ROOT / 'raw' / z['stem']
    out = ROOT / 'results' / z['stem']
    need(not Path(str(pre) + '.witness').exists(), 'fresh output prefix')
    source = {str(p.relative_to(ROOT)): sha(p) for folder in ('code', 'data') for p in sorted((ROOT / folder).iterdir()) if p.is_file()}
    contract = {'profile': z, 'source_sha256': source, 'commands_use_ubsan': True, 'threads': args.threads, 'scope': 'all M divisible by q*multiplier; marked denominator gcd(d,q)=1; CDF 3/2 at seed L'}
    (ROOT / 'logs' / (z['stem'] + '_CONTRACT.json')).write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    env = {**os.environ, 'CSTAR_THREADS': str(args.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    commands = []
    start = time.monotonic()

    def run(cmd):
        tick = time.monotonic()
        r = subprocess.run(list(map(str, cmd)), env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(r.stdout, end='', flush=True)
        commands.append({'cmd': list(map(str, cmd)), 'exit_code': r.returncode, 'seconds': round(time.monotonic() - tick, 3), 'transcript': r.stdout})
        need(r.returncode == 0, 'command succeeded')
    try:
        for mode in ('produce', 'consume'):
            run([ROOT / 'bin' / f"{mode}{z['ag']}", mode, z['p'], z['multiplier'], ROOT / 'data', pre])
        for ext in ('.cells', '.branches', '.audit'):
            a = Path(str(pre) + ext)
            b = Path(str(pre) + '.checked' + ext)
            need(a.read_bytes() == b.read_bytes(), 'complete consumer byte agreement ' + ext)
            Path(str(out) + ext).write_bytes(a.read_bytes())
        run([ROOT / 'bin/algebra', z['p'], z['multiplier'], ROOT / 'data', str(pre) + '.tensor'])
        Path(str(out) + '.tensor').write_bytes(Path(str(pre) + '.tensor').read_bytes())
        run([sys.executable, '-B', ROOT / 'code/emit_profile.py', z['p'], z['ag'], z['multiplier'], pre, str(out) + '.json'])
        run([sys.executable, '-B', ROOT / 'code/audit_one.py', z['p'], z['ag'], '--k', z['multiplier'], '--prefix', pre, '--output', str(out) + '.independent.json'])
        v = json.loads(Path(str(out) + '.json').read_text())
        print(json.dumps(v, indent=2), flush=True)
        receipt = {'status': 'PASS', 'profile': z, 'contract_sha256': sha(ROOT / 'logs' / (z['stem'] + '_CONTRACT.json')), 'witness_bytes': v['witness_bytes'], 'witness_sha256': v['witness_sha256'], 'all_geometry_recomputed': True, 'producer_consumer_agreement': True, 'seconds': round(time.monotonic() - start, 3), 'commands': commands}
        (ROOT / 'logs' / (z['stem'] + '_generation.json')).write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
    except BaseException:
        (ROOT / 'logs' / (z['stem'] + '_FAILED.json')).write_text(json.dumps({'commands': commands, 'seconds': round(time.monotonic() - start, 3)}, indent=2) + '\n', newline='\n')
        raise
if __name__ == '__main__':
    main()
