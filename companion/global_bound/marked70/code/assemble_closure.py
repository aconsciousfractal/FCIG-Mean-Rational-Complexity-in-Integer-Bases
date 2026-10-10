#!/usr/bin/env python3
"""Exact assembly of the frozen 70p partition and the 770ell tail.

The continuous arguments are in the proofs, not inferred from sampled points.
The quick assembly checks stored aggregates; raw streams have separate replays.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import argparse, csv, hashlib, io, json
from exact import need, log_lo, phi_hi, out_decimal, sieve
from candidate210 import certificate as candidate_certificate

ROOT = Path(__file__).resolve().parent.parent
DELTA = F(2, 85)
ERR70 = F(971503, 44247)
ERR770 = F(1755677, 43120)

def load(path):
    return json.loads(path.read_text())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def frac(row):
    return F(row['numerator'], row['denominator'])

@lru_cache(None)
def sharp_penalty(p):
    return phi_hi(F(12, 7*p)) + F(9, 14*p)

def sharp_bound(p, B, M, cap):
    need(B > 1 and M >= B and cap >= 0, 'positive normalization domain')
    return (cap + sharp_penalty(p) + ERR70/isqrt(M))/log_lo(B)

def route_csv(rows):
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue()

def certificate(root=ROOT):
    small = load(root/'profiles/results/small_prime_closure.json')
    meta = load(root/'results/inventory_certificate.json')
    theory = load(root/'checks/sharpened_pe_certificate.json')
    for data in (small, meta, theory):
        need(data['status'] == 'PASS', 'component status')
    need(frac(theory['sharp70']['sqrtM_coefficient']) == ERR70,
         'exact sharp70 coefficient bridge')
    need(frac(theory['sharp770']['sqrtM_coefficient']) == ERR770,
         'exact sharp770 coefficient bridge')
    need(frac(theory['linear_remainder_constant']['rational_cap']) == F(3, 16),
         'exact integrated linear remainder cap')
    need(not theory['inventory_or_cutoff_changed'], 'frozen partition retained')

    # These are new profiles from this delivery; raw replay is separate.
    # Historical profile caps are not imported into this closure.
    anchors = []
    for k, cap in ((1, F(7, 50)), (2, F(0))):
        path = root/f'profiles/results/q70_k{k}.json'
        row = load(path)
        need((root/f'checks/imports/q70_k{k}.json').read_bytes() == path.read_bytes(),
             'theory anchor and fresh profile byte identity')
        need(row['q'] == 70 and row['k'] == k and row['transfer_index'] == 70**k,
             'marked anchor and transfer divisor')
        value = F(row['total_upper'])
        need(value < cap, 'fresh anchor below relaxed cap')
        advertised = next(z for z in small['marked70_caps'] if z['k'] == k)
        need(advertised['divisor'] == 70**k and F(advertised['exact_upper']) == value
             and F(advertised['simple_cap']) == cap, 'small closure anchor bridge')
        anchors.append(dict(k=k, divisor_of_M=70**k, exact_upper=str(value),
                            relaxed_cap=str(cap), source_sha256=sha(path),
                            every_multiple_of_divisor=True))

    primes = [p for p in sieve(71) if p >= 11]
    need(small['closed_primes'] == primes and len(primes) == 16,
         'all sixteen small primes exactly once')
    need(small['all_positive_prime_exponents'] and small['all_lengths_N_ge_1']
         and not small['historical_profile_caps_assumed'], 'small family quantifiers')
    small_routes = []
    for row in small['routes']:
        p, q = row['p'], row['q']
        need(q == 70*p and row['base_min'] == q and row['required_divisor_of_M'] == q,
             'full-mark small route domain')
        profile = load(root/f'profiles/results/q{q}.json')
        cap = F(profile['total_upper'])
        upper = max(F(0), cap)/log_lo(q)
        need(profile['q'] == q and profile['k'] == 1 and profile['transfer_index'] == q,
             'full-mark small transfer divisor')
        need(F(row['mean_upper']) == cap and F(row['normalized_upper']) == upper
             and upper < DELTA, 'small normalized bound with correct sign handling')
        small_routes.append(dict(p=p, radical=q, mean_upper=str(cap),
                                 normalized_upper=str(upper)))
    need([r['p'] for r in small_routes] == primes, 'small route list complete')
    small_upper = max(F(r['normalized_upper']) for r in small_routes)
    need(small_upper == F(small['common_upper']), 'common small-prime bound')

    inventory = list(csv.DictReader((root/'results/finite_inventory.csv').open()))
    analytic = [r for r in inventory if r['route'] == 'analytic']
    witnesses = [r for r in inventory if r['route'] == 'witness']
    need(len(inventory) == meta['all_moduli'] == 3518
         and len(analytic) == meta['analytic_moduli'] == 2855
         and len(witnesses) == meta['witness_moduli'] == 663,
         'frozen finite partition counts')
    need(meta['base_upper_exclusive'] == 3_000_000 and meta['prime_range'] == [73, 1000]
         and meta['second_profile_used_at_N1'] == 163, 'finite domain and anchor choices')
    bounds_path = root/'results/finite_bounds.csv'
    bounds = [{k: int(v) for k, v in r.items()} for r in csv.DictReader(bounds_path.open())]
    pairs = [(int(r['M']), int(r['p'])) for r in witnesses]
    need([(r['M'], r['p']) for r in bounds] == pairs, 'every finite witness row present')
    need(len(set(m for m, p in pairs)) == 663, 'finite modulus uniqueness')
    mean_uppers = []
    for row in bounds:
        M = row['M']
        num = 2*row['sum_log_upper_ticks'] - M*row['log_M_lower_ticks']
        need(row['states'] == M and row['representatives'] == M//2+1,
             'zero and midpoint mass included')
        need(row['mean_numerator'] == row['first_numerator'] == num < 0
             and row['stencil'] == 0, 'negative exact finite aggregate, no refinement')
        need(F(row['log_M_lower_ticks'], 10**9) <= log_lo(M), 'finite log direction')
        mean_uppers.append((F(num, 2*M*10**9), M, row['p']))
    receipt = load(root/'results/finite_receipt.json')
    streams = load(root/'results/finite_stream_manifest.json')
    need(receipt['status'] == 'PASS' and receipt['inputs_unchanged']
         and receipt['output_sha256'] == sha(bounds_path), 'primary finite receipt hash')
    need([(r['M'], r['p']) for r in streams] == pairs
         and all(r['records'] == r['M']//2+1 and r['bytes'] == 8*r['records']
                 and len(r['sha256']) == 64 for r in streams), 'stream metadata domain')
    need(sum(r['states'] for r in bounds) == receipt['states'] == 135215780
         and sum(r['representatives'] for r in bounds) == receipt['representatives'] == 67608553
         and sum(r['bytes'] for r in streams) == receipt['stream_bytes'] == 540868424,
         'finite certified mass')

    # Keep the original route choices. Apply the stronger analytic bound only
    # to the same analytic rows; the other 663 rows retain their negative bound.
    sharpened = []
    basic_rows = []
    for row in analytic:
        M, p = int(row['M']), int(row['p'])
        cap = F(0) if M % 4900 == 0 else F(7, 50)
        need(int(row['profile_k']) == (2 if cap == 0 else 1), 'finite anchor divisibility')
        upper = sharp_bound(p, M, M, cap)
        old = F(row['normalized_upper'])
        need(0 < upper < old < DELTA, 'every analytic route improves strictly')
        sharpened.append(dict(B=M, M=M, p=p, profile_k=int(row['profile_k']),
                              source_route='analytic', basic_upper=str(old),
                              sharp_upper=str(upper), sharp_upper_decimal=out_decimal(upper)))
        basic_rows.append((old, 'finite_analytic', M, p))
    regions = []
    specs = [('large_p', 1000, 70000, 70000, F(7, 50)),
             ('N_ge_2', 73, 5110, 26112100, F(0)),
             ('large_base', 73, 3000000, 3000000, F(7, 50))]
    need(len(meta['infinite_regions']) == 3, 'three frozen infinite regions')
    for stored, (name, p, B, M, cap) in zip(meta['infinite_regions'], specs):
        need((stored['name'], stored['p_min'], stored['base_min'], stored['M_min'],
              F(stored['relaxed_cap'])) == (name, p, B, M, cap), 'region corner identity')
        upper = sharp_bound(p, B, M, cap)
        old = F(stored['normalized_upper'])
        need(0 < upper < old < DELTA, 'same infinite region strictly improves')
        regions.append(dict(name=name, p_scalar_min=p, base_min=B, M_min=M,
                            relaxed_cap=str(cap), basic_upper=str(old),
                            sharp_upper=str(upper), sharp_upper_decimal=out_decimal(upper)))
        basic_rows.append((old, name, B, p))
    basic_rows.extend([(small_upper, 'small_profiles', 770, 11), (F(0), 'finite_negative', 0, 0)])
    basic_max = max(basic_rows)
    sharp_choices = [(F(r['sharp_upper']), 'finite_analytic', r['B'], r['p']) for r in sharpened]
    sharp_choices += [(F(r['sharp_upper']), r['name'], r['base_min'], r['p_scalar_min']) for r in regions]
    sharp_choices += [(small_upper, 'small_profiles', 770, 11), (F(0), 'finite_negative', 0, 0)]
    final = max(sharp_choices)
    gamma = final[0]
    need(basic_max == (F(25221644582229033312307, 1072177116035751725286301),
                       'finite_analytic', 189140, 193), 'basic bound exact maximum')
    need(final == (F(503673291285648397571314, 23233333961830267469605761),
                   'finite_analytic', 65030, 929), 'sharp bound exact maximum')
    need(gamma < basic_max[0] < DELTA, 'strict uniform separation')

    # A second insertion preserves all four old primes, using radius six.
    anchor_path = root/'profiles/results/q770.json'
    q770 = load(anchor_path)
    need((root/'checks/imports/q770.json').read_bytes() == anchor_path.read_bytes(),
         'new five-prime anchor byte identity')
    need(q770['q'] == q770['transfer_index'] == 770 and q770['k'] == 1
         and F(q770['total_upper']) < F(9, 80), 'five-prime anchor for all multiples')
    ell, B = 300, 231000
    sharp5 = (F(9, 80)+phi_hi(F(180, 77*ell))+F(135, 154*ell)+ERR770/isqrt(B))/log_lo(B)
    basic5 = (F(9, 80)+phi_hi(F(360, 77*ell))+F(453, 10*isqrt(B)))/log_lo(B)
    need(0 < sharp5 < basic5 < DELTA, 'unchanged five-prime tail gate')
    need(next(p for p in sieve(400) if p >= ell) == 307, 'first actual prime in tail')
    candidate = candidate_certificate()
    need(candidate == load(root/'results/candidate210.json'), 'candidate regenerated exactly')
    need(candidate['integer_comparison_pass'] and DELTA < F(candidate['c0_lower']),
         'strict comparison with exact candidate')
    worst_mean, worst_M, worst_p = max(mean_uppers)
    result = dict(
        status='PASS', normalization='Delta(B,N)=mu(B^N)/log(B)',
        all_residues_including_zero=True, norm='infinity', logarithms='natural',
        theorem_70p=dict(radical='70p', prime_min=11, prime_condition='p prime',
            all_positive_valuations=True, all_lengths_N_ge_1=True,
            uniform_upper=str(gamma), uniform_upper_decimal=out_decimal(gamma),
            strict_margin_to_2_over_85=str(DELTA-gamma),
            basic_uniform_upper=str(basic_max[0]), basic_uniform_upper_decimal=out_decimal(basic_max[0]),
            largest_bound_route=dict(route=final[1], B=final[2], N=1, p=final[3]),
            largest_bound_is_not_a_Delta_maximizer=True,
            marked_anchors=anchors, small_prime_routes=small_routes, infinite_regions=regions,
            frozen_finite_moduli=3518, analytic_moduli=2855, witness_moduli=663,
            small_prime_common_upper=str(small_upper), all_finite_witness_means_negative=True,
            worst_finite_witness_mean_upper=str(worst_mean),
            worst_finite_witness_mean_at=dict(M=worst_M, p=worst_p),
            historical_profile_caps_assumed=False),
        theorem_770ell=dict(radical='770*ell', scalar_cutoff=300, first_prime=307,
            all_positive_valuations=True, all_lengths_N_ge_1=True,
            exact_marked_anchor=str(F(q770['total_upper'])), relaxed_anchor='9/80',
            basic_upper=str(basic5), basic_upper_decimal=out_decimal(basic5),
            sharp_upper=str(sharp5), sharp_upper_decimal=out_decimal(sharp5),
            strict_margin_to_2_over_85=str(DELTA-sharp5),
            repair_radius=6, cutoff_changed_by_sharpening=False),
        reusable_results=dict(core='2*r*s for distinct odd primes r,s', repair_radius=4,
            periodic_C_upper='8', sharp_Phi_argument='5/(2p)',
            sharp_linear_remainder='15/(16p)', sharp_sqrtM_coefficient='298924/6615',
            fresh_anchor_still_required_for_other_cores=True),
        finite_states=receipt['states'], finite_representatives=receipt['representatives'],
        geometric_profiles=small['direct_profiles'],
        geometric_witnesses=small['geometric_witnesses_in_direct_profiles'],
        geometric_corner_evaluations=small['corner_evaluations_in_direct_profiles'],
        c0_lower=candidate['c0_lower'], c0_upper=candidate['c0_upper'],
        c0_interval_display=[candidate['c0_lower_decimal'], candidate['c0_upper_decimal']],
        exact_global_Cstar_proved=False, other_four_prime_supports_all_closed=False,
        all_five_or_more_prime_supports_closed=False, sign_classification_of_all_offsets=False,
        limits=dict(continuous_proof_in_documents=True, formal_assistant_proof=False,
            stored_aggregate_assembly_is_raw_replay=False,
            alternative_geometric_audit_exhaustive=False,
            producer_consumer_share_some_geometry=True,
            human_specialist_review_claimed=False))
    return result, sharpened

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result, rows = certificate()
    outputs = {'closure_certificate.json': json.dumps(result, indent=2)+'\n',
               'sharp_routes.csv': route_csv(rows)}
    for name, data in outputs.items():
        target = ROOT/'results'/name
        if args.check:
            need(target.read_bytes() == data.encode(), 'exact assembled output '+name)
        else:
            target.write_bytes(data.encode())
    print(json.dumps(dict(status='PASS', gamma70=result['theorem_70p']['uniform_upper_decimal'],
        gamma770_tail=result['theorem_770ell']['sharp_upper_decimal'],
        finite_moduli=3518, global_Cstar_identified=False), indent=2))

if __name__ == '__main__':
    main()
