"""Bind a complete profile replay to all frozen component bytes."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time
ROOT = Path(__file__).resolve().parent.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def run(component, output, threads):
    component = component.resolve()
    output = output.resolve()
    if output.exists():
        raise RuntimeError('fresh output directory required')
    output.mkdir(parents=True)
    mpath = component / 'MANIFEST_SHA256.json'
    manifest = json.loads(mpath.read_text())

    def check():
        actual = {str(p.relative_to(component)) for p in component.rglob('*') if p.is_file() and p != mpath}
        if actual != set(manifest):
            raise RuntimeError('complete component manifest')
        for name, h in manifest.items():
            p = component / name
            if not p.resolve().is_relative_to(component) or p.is_symlink() or sha(p) != h:
                raise RuntimeError('frozen component input ' + name)
    check()
    mhash = sha(mpath)
    command = [sys.executable, '-B', str(component / 'verify_profiles.py'), '--rebuild', '--sanitize', '--threads', str(threads), '--report', str(output / 'REPORT.json')]
    contract = dict(status='FROZEN_BEFORE_RUN', export_manifest_sha256=mhash, source_hashes=manifest, run=command, threads=threads, sanitize=True, complete_geometric_witness_rebuild_requested=True)
    cp = output / 'CONTRACT.json'
    cp.write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    start = time.monotonic()
    with (output / 'run.log').open('w') as f:
        p = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
    check()
    if sha(mpath) != mhash:
        raise RuntimeError('component manifest changed')
    if p.returncode:
        raise RuntimeError('profile rebuild failed; ' + str(output / 'run.log'))
    report = json.loads((output / 'REPORT.json').read_text())
    if report['status'] != 'PASS' or not report['sanitize']:
        raise RuntimeError('sanitized replay report required')
    binding = dict(status='PASS', exit_code=p.returncode, export_manifest_sha256=mhash, contract_sha256=sha(cp), receipt_sha256=sha(output / 'REPORT.json'), log_sha256=sha(output / 'run.log'), every_export_hash_verified_before_and_after=True, export_mutations_detected=False, sanitize=True, threads=threads, compared_files=report['compared_files'], profiles=report['profiles'], geometric_witnesses=report['geometric_witnesses'], geometric_corner_evaluations=report['geometric_corner_evaluations'], seconds=round(time.monotonic() - start, 3))
    (output / 'BINDING.json').write_text(json.dumps(binding, indent=2) + '\n', newline='\n')
    print(json.dumps(binding, indent=2), flush=True)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--component', type=Path, default=ROOT / 'profiles')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=2)
    a = p.parse_args()
    run(a.component, a.output, a.threads)
