"""Portable, read-only replay of the two light arithmetic certificates."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys, tempfile
HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rows = []
    for row in json.loads((HERE / 'MANIFEST.json').read_text())['files']:
        raw = (HERE / row['path']).read_bytes()
        assert len(raw) == row['bytes'], row['path']
        assert hashlib.sha256(raw).hexdigest() == row['sha256'], row['path']
    with tempfile.TemporaryDirectory(prefix='p2310_theory_replay_') as directory:
        for script, certificate in [('verify_unequal_transfer.py', 'unequal_transfer_certificate.json'), ('verify_global_third.py', 'global_third_certificate.json')]:
            output = Path(directory) / certificate
            subprocess.run([sys.executable, str(HERE / script), '--output', str(output)], check=True, stdout=subprocess.PIPE, text=True)
            assert output.read_bytes() == (HERE / certificate).read_bytes(), certificate
            data = json.loads(output.read_text())
            rows.append(dict(script=script, certificate=certificate, status=data['status'], checks=data['checks'], byte_identical=True))
    result = dict(status='PASS', project='P72_2310_THEORY', certificates=rows, geometry_replayed=False, historical_theorem_replays_performed=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
