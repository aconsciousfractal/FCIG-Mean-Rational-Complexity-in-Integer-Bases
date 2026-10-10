"""Verify the frozen component and replay certificates into a fresh directory."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time
ROOT = Path(__file__).resolve().parent

def need(condition, label):
    if not condition:
        raise RuntimeError(label)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--root', type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    need(not output.is_relative_to(ROOT), 'output must be external to the frozen component')
    need(not output.exists() or not any(output.iterdir()), 'output directory must be empty')
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    manifest = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text())
    entries = manifest['files']
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()}
    need(actual == set(entries) | {'MANIFEST_SHA256.json'}, 'exact component inventory')
    for relative, expected in entries.items():
        raw = (ROOT / relative).read_bytes()
        need(len(raw) == expected['bytes'] and sha256(raw).hexdigest() == expected['sha256'], 'component hash: ' + relative)
    bridges = []
    if args.root is not None:
        mainroot = args.root.resolve()
        provenance = json.loads((ROOT / 'imports/closure/PROVENANCE.json').read_text())
        for row in provenance['sources']:
            relative = row['source'].split('/', 1)[1]
            raw = (mainroot / relative).read_bytes()
            need(len(raw) == row['bytes'] and sha256(raw).hexdigest() == row['sha256'], 'main-project snapshot bridge: ' + relative)
            bridges.append(relative)
    reports = []
    commands = [('insertion78_certificate.json', [sys.executable, '-B', str(ROOT / 'code/verify_insertion78.py'), '--output', str(output / 'insertion78_certificate.json')]), ('gates78_certificate.json', [sys.executable, '-B', str(ROOT / 'code/verify_gates78.py'), '--inputs', str(ROOT / 'imports/closure'), '--output', str(output / 'gates78_certificate.json')])]
    for name, command in commands:
        run = subprocess.run(command, capture_output=True, text=True)
        (output / (name + '.stdout.txt')).write_text(run.stdout, newline='\n')
        (output / (name + '.stderr.txt')).write_text(run.stderr, newline='\n')
        need(run.returncode == 0, 'checker failed: ' + name + '; see external logs')
        generated = (output / name).read_bytes()
        need(generated == (ROOT / 'results' / name).read_bytes(), 'exact certificate replay: ' + name)
        report = json.loads(generated)
        need(report['status'] == 'PASS', 'checker certificate status')
        reports.append(dict(file=name, checks=report['checks'], bytes=len(generated), sha256=sha256(generated).hexdigest()))
    result = dict(status='PASS', payload_files=len(entries), manifest_verified=True, source_bridges=bridges, reports=reports, elapsed_seconds=round(time.monotonic() - started, 6), profile_geometry_replayed_here=False, finite_physical_witnesses_replayed_here=False)
    (output / 'REPLAY.json').write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
