"""Read-only replay of the Stevens deduction's integer certificate."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys, tempfile
HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    for row in manifest['files']:
        raw = (HERE / row['path']).read_bytes()
        assert len(raw) == row['bytes'] and hashlib.sha256(raw).hexdigest() == row['sha256'], row['path']
    with tempfile.TemporaryDirectory(prefix='stevens2000_replay_') as temp:
        target = Path(temp) / 'certificate.json'
        subprocess.run([sys.executable, str(HERE / 'verify_stevens2000.py'), '--output', str(target)], check=True, stdout=subprocess.PIPE, text=True)
        assert target.read_bytes() == (HERE / 'stevens2000_certificate.json').read_bytes()
        certificate = json.loads(target.read_text())
    result = dict(status='PASS', project='P72_STEVENS_2000', checks=certificate['checks'], certificate_byte_identical=True, classical_import_reproved=False, geometry_replayed=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
