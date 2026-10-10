"""Bind the complete new refined replay to both frozen geometric components."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time
ROOT = Path(__file__).resolve().parent.parent

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def run(component, baseline, output, threads):
    component = component.resolve()
    baseline = baseline.resolve()
    output = output.resolve()
    if output.exists():
        raise RuntimeError('fresh output directory required')
    inputs = {}
    for label, folder in [('refined', component), ('baseline', baseline)]:
        mp = folder / 'MANIFEST_SHA256.json'
        inputs[label] = dict(manifest_sha256=sha(mp), files=json.loads(mp.read_text()))

    def check():
        for label, folder in [('refined', component), ('baseline', baseline)]:
            mp = folder / 'MANIFEST_SHA256.json'
            entry = inputs[label]
            actual = {str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and p != mp}
            if actual != set(entry['files']) or sha(mp) != entry['manifest_sha256']:
                raise RuntimeError('complete frozen component ' + label)
            for name, h in entry['files'].items():
                p = folder / name
                if not p.resolve().is_relative_to(folder) or p.is_symlink() or sha(p) != h:
                    raise RuntimeError('frozen input ' + label + '/' + name)
    check()
    output.mkdir(parents=True)
    command = [sys.executable, '-B', str(component / 'verify_refined.py'), '--rebuild', '--sanitize', '--baseline-root', str(baseline), '--threads', str(threads), '--report', str(output / 'REPORT.json')]
    contract = dict(status='FROZEN_BEFORE_RUN', components=inputs, run=command, threads=threads, sanitize=True, refined_raw_rebuild_requested=True, coarse_raw_rebuild_requested=False, coarse_dependency_paid_by_separate_baseline_replay=True)
    cp = output / 'CONTRACT.json'
    cp.write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    start = time.monotonic()
    with (output / 'run.log').open('w') as f:
        p = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
    check()
    if p.returncode:
        raise RuntimeError('refined replay failed: ' + str(output / 'run.log'))
    r = json.loads((output / 'REPORT.json').read_text())
    if not (r['status'] == 'PASS' and r['sanitize'] and r['refined_witnesses_replayed'] and r['external_baseline_bridge_verified'] and (r['compared_outputs_or_raw_hashes'] == 8) and (r['refined_witness_records'] == 32223744) and (r['refined_corner_evaluations'] == 17014136832)):
        raise RuntimeError('incomplete refined replay scope')
    binding = dict(status='PASS', exit_code=0, refined_manifest_sha256=inputs['refined']['manifest_sha256'], baseline_manifest_sha256=inputs['baseline']['manifest_sha256'], contract_sha256=sha(cp), receipt_sha256=sha(output / 'REPORT.json'), log_sha256=sha(output / 'run.log'), every_component_hash_verified_before_and_after=True, component_mutations_detected=False, sanitize=True, threads=threads, compared_outputs_or_raw_hashes=r['compared_outputs_or_raw_hashes'], refined_witness_records=r['refined_witness_records'], refined_corner_evaluations=r['refined_corner_evaluations'], baseline_raw_replayed_by_this_invocation=False, seconds=round(time.monotonic() - start, 3))
    (output / 'BINDING.json').write_text(json.dumps(binding, indent=2) + '\n', newline='\n')
    print(json.dumps(binding, indent=2), flush=True)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--component', type=Path, default=ROOT / 'profiles_refined')
    p.add_argument('--baseline', type=Path, default=ROOT / 'profiles')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=4)
    a = p.parse_args()
    run(a.component, a.baseline, a.output, a.threads)
