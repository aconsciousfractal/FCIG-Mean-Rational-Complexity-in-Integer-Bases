"""Exact assembly of the two infinite regions and the finite N=1 complement.

Every geometric cap is a bound, not an observed value of the true mean.
This assembler never promotes an unsuccessful fine gate to a theorem.
"""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import argparse, csv, hashlib, json
from exact import need, log_lo, out_decimal
from candidate210 import certificate as candidate
from finite2310 import inventory
ROOT = Path(__file__).resolve().parent.parent
Q = 2310
PS = (2, 3, 5, 7, 11)
DELTA = F(2, 85)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def valuations(K):
    ee = []
    r = K
    for p in PS:
        e = 0
        while r % p == 0:
            r //= p
            e += 1
        ee.append(e)
    need(r == 1 and all((e >= 1 for e in ee)), 'full marked prime support of index')
    return ee

def certificate(root=ROOT):
    inv = inventory()
    need(inv == json.loads((root / 'results/finite_inventory.json').read_text()), 'entire N1 complement')
    domain = ''.join((f"{r['M']} {r['B']}\n" for r in inv['rows']))
    need(domain == (root / 'results/finite_domain.txt').read_text(), 'N1 physical domain')
    profiles = []
    for path in sorted((root / 'profiles/results').glob('*.json')):
        if path.name.endswith('.independent.json'):
            continue
        z = json.loads(path.read_text())
        if not isinstance(z, dict) or 'total_upper' not in z:
            continue
        need(z['q'] == Q, 'only declared full2310 anchors')
        K = z['transfer_index']
        ee = valuations(K)
        A = F(z['total_upper'])
        need(z['labels'] == 6912 and z['witnesses'] == z['cells'] * 6912 * (K // Q), 'marked branch and witness counts')
        if z.get('transfer_multiplier') is not None:
            need(z['k'] is None and z['transfer_multiplier'] == K // Q and (z['required_divisor_of_M'] == K) and z['valid_for_every_positive_multiple'], 'mixed index explicit domain')
        else:
            need(K == Q ** z['k'], 'legacy pure index domain')
        profiles.append(dict(file=str(path.relative_to(root)), sha256=sha(path), transfer_index=K, prime_valuations=ee, all_N_from=max(ee), mean_upper=str(A), mean_upper_decimal=out_decimal(A), normalized_upper=str(max(A, F(0)) / log_lo(Q)), cells=z['cells'], branches=6912 * (K // Q), angular_nodes=z['angular_nodes'], witnesses=z['witnesses'], witness_bytes=z['witness_bytes'], raw_sha256=z['witness_sha256']))
    need(profiles, 'new geometric anchors present')
    first = [r for r in profiles if r['transfer_index'] == Q]
    need(len(first) == 1, 'exactly one first-level2310 anchor')
    A1 = F(first[0]['mean_upper'])
    need(A1 < F(4, 15), 'coarse initial cap')
    high_base = max(A1, F(0)) / log_lo(100000)
    need(high_base < DELTA, 'all large bases at every N')
    eligible = [r for r in profiles if r['all_N_from'] <= 2]
    best2 = min(eligible, key=lambda r: F(r['mean_upper']))
    upper2 = max(F(best2['mean_upper']), F(0)) / log_lo(Q)
    fine2 = upper2 < DELTA
    rows = [{k: int(v) for k, v in z.items()} for z in csv.DictReader((root / 'results/finite_bounds.csv').open())]
    need([(r['M'], r['B']) for r in rows] == [(r['M'], r['B']) for r in inv['rows']], 'every finite residue table has its prescribed base')
    finite = []
    for r in rows:
        M, B = (r['M'], r['B'])
        num = 2 * r['sum_log_upper_ticks'] - M * r['log_M_lower_ticks']
        need(M == B and r['states'] == M and (r['representatives'] == M // 2 + 1), 'N1 all even states')
        need(num == r['mean_numerator'] and r['log_B_lower_ticks'] == r['log_M_lower_ticks'], 'centering and logB normalization')
        need(F(r['log_M_lower_ticks'], 10 ** 9) <= log_lo(M), 'independent modulus logarithm')
        u = F(max(num, 0), 2 * M * r['log_B_lower_ticks'])
        need(u < DELTA, 'every finite N1 bound below2/85')
        finite.append(dict(B=B, N=1, M=M, upper=str(u), mean_upper=str(F(num, 2 * M * 10 ** 9))))
    fin = max(finite, key=lambda r: F(r['upper']))
    fupper = F(fin['upper'])
    fcert = json.loads((root / 'results/finite_certificate.json').read_text())
    need(fcert['status'] == 'PASS' and fcert['moduli'] == 30 and (F(fcert['normalized_upper']) == fupper) and (fcert['states'] == inv['states']) and (fcert['representatives'] == inv['representatives']), 'consumed physical certificate reconciled')
    need(fcert['output_sha256'] == sha(root / 'results/finite_bounds.csv'), 'exact physical table custody')
    c = candidate()
    need(c == json.loads((root / 'results/candidate210.json').read_text()), 'candidate fresh exact minima')
    need(DELTA < F(c['c0_lower']), 'strict exact candidate separation')
    routes = [dict(region='N1_B_lt_100000', upper=str(fupper), worst_bound_base=fin['B']), dict(region='B_ge_100000_all_N', upper=str(high_base), mean_upper=str(A1), base_min=100000), dict(region='N_ge_2_all_B', upper=str(upper2), mean_upper=best2['mean_upper'], transfer_index=best2['transfer_index'], base_min=Q)]
    worst = max(routes, key=lambda r: F(r['upper']))
    Gamma = F(worst['upper'])
    need((Gamma < DELTA) == fine2, 'fine theorem promotion exactly follows all regional gates')
    new_wit = sum((r['witnesses'] for r in profiles))
    new_corner = sum((r['witnesses'] * r['angular_nodes'] * 4 for r in profiles))
    return dict(status='PASS', radical=Q, prime_support=list(PS), normalization='Delta(B,N)=mu(B^N)/log(B)', norm='infinity', logarithms='natural', all_residues_including_zero=True, all_positive_valuations=True, all_N_ge_1=True, uniform_upper=str(Gamma), uniform_upper_decimal=out_decimal(Gamma), full_fine_exclusion_below_2_over_85=fine2, full_exclusion_below_c0=Gamma < F(c['c0_lower']), strict_margin_to_2_over_85=str(DELTA - Gamma), worst_bound_route=worst, worst_bound_is_not_claimed_to_maximize_Delta=True, all_N1_excluded_below_2_over_85=True, N1_uniform_upper=str(max(high_base, fupper)), N1_uniform_upper_decimal=out_decimal(max(high_base, fupper)), mixed_transfer_lemma='rad(K)=q, every kernel K/q copies, exactlyK admitted branches per parent', transfer_domains='Every positive multiple of the declared index, including extra primes', regions=routes, profiles=profiles, finite_rows=finite, finite_moduli=len(rows), finite_states=inv['states'], finite_representatives=inv['representatives'], all_finite_witness_means_negative=all((F(r['mean_upper']) < 0 for r in finite)), geometric_profiles=len(profiles), geometric_witnesses_per_pass=new_wit, geometric_corner_evaluations_per_pass=new_corner, c0_lower=c['c0_lower'], c0_upper=c['c0_upper'], exact_global_Cstar_proved=False, all_five_prime_supports_excluded=False, sign_classification_of_all_offsets=False, continuous_proof_formalized=False, alternative_geometric_audit_exhaustive=False, physical_minimality_claimed=False, historical_numeric_profile_caps_assumed=False)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    z = certificate()
    s = json.dumps(z, indent=2) + '\n'
    path = ROOT / 'results/2310_certificate.json'
    if a.check:
        need(path.read_text() == s, 'frozen2310 theorem assembly')
    else:
        path.write_text(s, newline='\n')
    print(json.dumps({k: z[k] for k in ('status', 'uniform_upper_decimal', 'full_fine_exclusion_below_2_over_85', 'all_N1_excluded_below_2_over_85')}, indent=2))
