"""Final radical2310 closure, layered over the preserved coarse checkpoint.

The new all-N>=2 gate is strict below c0. It is deliberately tested
separately from the stronger 2/85 gate, which this profile does not pass.
"""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import argparse, hashlib, json
from exact import need, log_lo, out_decimal
from assemble2310 import certificate as baseline_certificate
ROOT = Path(__file__).resolve().parent.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def certificate(root=ROOT):
    old = baseline_certificate(root)
    need(old == json.loads((root / 'results/2310_certificate.json').read_text()), 'preserved coarse checkpoint and all finite N1 results')
    p = root / 'profiles_refined/results/profile.json'
    z = json.loads(p.read_text())
    need(z['status'] == 'PASS' and z['q'] == 2310 and (z['transfer_index'] == 13860) and (z['transfer_multiplier'] == 6) and (z['required_divisor_of_M'] == 13860) and z['valid_for_every_positive_multiple_including_extra_primes'], 'refined full-mark transfer and exact domain')
    valuations = (2, 2, 1, 1, 1)
    primes = (2, 3, 5, 7, 11)
    need(prod((p ** e for p, e in zip(primes, valuations))) == 13860 and max(valuations) == 2 and (2310 ** 2 % 13860 == 0), 'every B^N with radical2310 and N>=2 is an admitted multiple')
    need(z['labels'] == 6912 and z['branches'] == 41472 and (z['new_source_cells'] == 777) and (z['new_subcells'] == 3108) and (z['coarse_cells_retained'] == 1890) and (z['hybrid_cells'] == 4998), 'hybrid certificate geometry counts')
    A = F(z['total_upper'])
    gamma2 = max(A, F(0)) / log_lo(2310)
    need(gamma2 == F(z['normalized_upper']), 'normalization by logB, with B>=2310')
    c0lo = F(old['c0_lower'])
    delta = F(2, 85)
    gamma1 = F(old['N1_uniform_upper'])
    gamma = max(gamma1, gamma2)
    need(gamma1 < delta < c0lo and gamma2 < c0lo and (gamma == gamma2), 'all regions strictly below the candidate lower bound')
    need(gamma2 > delta, 'this accepted profile does not prove the stronger 2/85 statement')
    need(old['full_fine_exclusion_below_2_over_85'] is False and old['full_exclusion_below_c0'] is False, 'coarse checkpoint remains accurately labeled')
    return dict(status='PASS', radical=2310, prime_support=list(primes), theorem='For every B with rad(B)=2310 and every N>=1, Delta(B,N)<=Gamma2310<c0', normalization='Delta(B,N)=mu(B^N)/ln(B)', norm='infinity', logarithms='natural', all_residues_including_zero=True, all_positive_valuations=True, all_N_ge_1=True, uniform_upper=str(gamma), uniform_upper_decimal=out_decimal(gamma), full_exclusion_below_c0=True, full_fine_exclusion_below_2_over_85=False, c0_lower=str(c0lo), strict_margin_to_c0_lower=str(c0lo - gamma), strict_margin_to_2_over_85=str(delta - gamma), regions=[dict(region='N=1, B<100000', upper=old['regions'][0]['upper'], finite_moduli=30, all_finite_states=1284360, all_finite_representatives=642210), dict(region='B>=100000, all N>=1', upper=old['regions'][1]['upper']), dict(region='N>=2, every B with rad(B)=2310', upper=str(gamma2), mean_upper=str(A), required_divisor=13860, minimum_base=2310)], strongest_accepted_route='refined N>=2', observed_maximizer_or_exact_value_claimed=False, geometric_component=dict(path=str(p.relative_to(root)), sha256=sha(p), manifest_sha256=sha(root / 'profiles_refined/MANIFEST_SHA256.json'), coarse_component_manifest_sha256=sha(root / 'profiles/MANIFEST_SHA256.json'), transfer_index=13860, prime_valuations=list(valuations)), coarse_checkpoint=dict(path='results/2310_certificate.json', sha256=sha(root / 'results/2310_certificate.json'), superseded_only_for_final_N_ge_2_bound=True), counting=dict(baseline_geometric_records=167215104, refined_records_reconsumed=32223744, refined_records_are_reused_baseline_subset=True, refined_subrectangle_applications=128894976, baseline_corner_evaluations_per_pass=6019743744, refined_corner_evaluations_per_pass=17014136832, total_corner_evaluations_per_complete_pass=23033880576, physical_and_geometric_counts_not_conflated=True), exact_global_Cstar_proved=False, all_five_prime_supports_excluded=False, all_offsets_sign_classified=False, physical_minimality_claimed=False, continuous_proof_formalized=False, alternative_geometric_audit_exhaustive=False)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    z = certificate()
    s = json.dumps(z, indent=2) + '\n'
    p = ROOT / 'results/final_2310_certificate.json'
    if a.check:
        need(p.read_text() == s, 'frozen final2310 closure')
    else:
        p.write_text(s, newline='\n')
    print(json.dumps({k: z[k] for k in ['status', 'uniform_upper_decimal', 'full_exclusion_below_c0', 'full_fine_exclusion_below_2_over_85']}, indent=2))
