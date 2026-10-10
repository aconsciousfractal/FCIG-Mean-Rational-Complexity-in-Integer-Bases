"""Read-only portable replay of the shaped-CDF supplement."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys
ROOT = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def check():
    p = ROOT / 'MANIFEST.json'
    m = json.loads(p.read_text())
    expected = {r['path']: r for r in m['files']}
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and p.name != 'MANIFEST.json'}
    if actual != set(expected):
        raise AssertionError('complete supplement inventory differs')
    for name, item in expected.items():
        p = ROOT / name
        if p.stat().st_size != item['bytes'] or sha(p) != item['sha256']:
            raise AssertionError('supplement hash differs: ' + name)
    return (sha(ROOT / 'MANIFEST.json'), expected)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--theory', type=Path, default=ROOT.parent / 'theory')
    a = p.parse_args()
    o = a.output.resolve()
    if o.exists() or o.is_relative_to(ROOT):
        raise AssertionError('fresh output outside supplement required')
    mh, items = check()
    o.mkdir(parents=True)
    (o / 'CONTRACT.json').write_text(json.dumps(dict(manifest_sha256=mh, files=len(items), new_geometric_flows=0, operation='regenerate shaped-CDF rational assembly from frozen theory'), indent=2) + '\n', newline='\n')
    result = o / 'shape_certificate.json'
    run = subprocess.run([sys.executable, '-B', str(ROOT / 'verify_shape.py'), '--theory', str(a.theory.resolve()), '--output', str(result)], capture_output=True, text=True)
    (o / 'run.log').write_text(run.stdout + run.stderr, newline='\n')
    if run.returncode or result.read_bytes() != (ROOT / 'shape_certificate.json').read_bytes():
        raise AssertionError('shaped-CDF replay failed or output differs')
    if check() != (mh, items):
        raise AssertionError('supplement changed during replay')
    z = dict(status='PASS', manifest_sha256=mh, component_file_count=len(items), result_sha256=sha(result), output_byte_identical=True, sources_unchanged=True, new_geometric_flows=0, baseline_component_replay_obligations_preserved=True)
    (o / 'REPORT.json').write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    print(json.dumps(z, indent=2))
if __name__ == '__main__':
    main()
