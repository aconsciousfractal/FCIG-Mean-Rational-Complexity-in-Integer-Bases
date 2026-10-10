"""Check stored exact integrals and route assembly, without geometric witness replay."""
from pathlib import Path
import argparse, hashlib, json, sys, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
from checks import input_checks
from assemble_small import assemble

def verify_manifest(root=ROOT):
    manifest = json.loads((root / 'MANIFEST_SHA256.json').read_text())
    for name, expected in manifest.items():
        p = root / name
        need(p.is_file() and p.resolve().is_relative_to(root.resolve()), 'contained manifest file ' + name)
        need(hashlib.sha256(p.read_bytes()).hexdigest() == expected, 'manifest SHA256 ' + name)
    return len(manifest)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--report', type=Path)
    ap.add_argument('--no-manifest', action='store_true', help=argparse.SUPPRESS)
    args = ap.parse_args()
    start = time.monotonic()
    checked = 0 if args.no_manifest else verify_manifest()
    inp = input_checks(ROOT)
    need(inp == json.loads((ROOT / 'results/common_input_checks.json').read_text()), 'independent logarithm-table audit')
    z = assemble(ROOT)
    need(z == json.loads((ROOT / 'results/small_prime_closure.json').read_text()), 'stored closure assembly')
    neg = json.loads((ROOT / 'results/negative_checks.json').read_text())
    need(neg['status'] == 'PASS' and len(neg['cases']) == 6, 'stored negative-check report only')
    report = {'status': 'PASS', 'mode': 'stored_cells_and_exact_assembly', 'manifest_files': checked, 'direct_profiles_reintegrated': z['direct_profiles'], 'full_mark_primes': len(z['closed_primes']), 'all_positive_prime_exponents': True, 'all_lengths_N_ge_1': True, 'geometric_witnesses_replayed': False, 'geometric_costs_recomputed': False, 'negative_tests_reexecuted': False, 'independent_logarithm_checks_reexecuted': True, 'common_upper': z['common_upper'], 'runtime_seconds': round(time.monotonic() - start, 3)}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report, indent=2), flush=True)
if __name__ == '__main__':
    main()
