"""Quick consistency check of stored cells, exact integrals and exhaustive routes.
This deliberately does not certify geometric cell costs from raw witnesses. The
full verifier rebuilds and consumes all witnesses; see generation logs for the
full computation actually executed when this package was prepared.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, subprocess, sys, tempfile, time, struct
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
from assemble_small import assemble

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for z in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(z)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--report', type=Path)
    ap.add_argument('--no-manifest', action='store_true', help=argparse.SUPPRESS)
    a = ap.parse_args()
    start = time.monotonic()
    checked = 0
    if not a.no_manifest:
        files = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text())
        for name, digest in files.items():
            need(sha(ROOT / name) == digest, 'manifest ' + name)
            checked += 1
    plan = json.loads((ROOT / 'profiles_plan.json').read_text())
    for z in plan:
        folder = ROOT / z['folder']
        pref = folder / 'results' / z['stem']
        subprocess.run([sys.executable, '-B', str(folder / 'code/stored_profile.py'), str(z['p']), str(z['ag']), str(z['k']), str(pref)], check=True)
    with tempfile.TemporaryDirectory(prefix='p72-42p-stored-') as d:
        dst = Path(d) / 'integrated.json'
        subprocess.run([sys.executable, '-B', str(ROOT / 'iterated/code/integrate.py'), '--output', str(dst)], check=True)
        need(dst.read_bytes() == (ROOT / 'iterated/results/integrated.json').read_bytes(), 'stored iterated profiles reintegrated')
    from checks import direct, input_checks
    R, tests = direct(462, 462)
    raw = (ROOT / 'iterated/results/M462_q462.min').read_bytes()
    need(tuple(R) == struct.unpack('<' + 'I' * 462, raw), 'fresh exact base462 minimum array')
    inp = input_checks(ROOT)
    need(inp == json.loads((ROOT / 'results/common_input_checks.json').read_text()), 'independent arithmetic input checks')
    closure = assemble(ROOT)
    need((json.dumps(closure, indent=2) + '\n').encode() == (ROOT / 'results/small_prime_closure.json').read_bytes(), 'complete exact route assembly')
    report = {'status': 'PASS', 'mode': 'quick_stored_consistency', 'manifest_files_checked': checked, 'stored_direct_profiles_reintegrated': len(plan), 'stored_iterated_horizons_reintegrated': 2, 'small_base_states_reenumerated': 462, 'small_base_denominator_tests': tests, 'closed_small_primes': closure['closed_primes'], 'common_upper': closure['common_upper'], 'common_upper_decimal': closure['common_upper_decimal'], 'marked42_caps': closure['marked42_caps'], 'geometric_witnesses_replayed': False, 'geometric_cell_costs_recomputed': False, 'uniform_claim_relies_on_recorded_full_generation_and_consumer': True, 'runtime_seconds': round(time.monotonic() - start, 3)}
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report, indent=2), flush=True)
if __name__ == '__main__':
    main()
