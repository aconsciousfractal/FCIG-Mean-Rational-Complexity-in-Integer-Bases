"""Pay 6006 and 154p, retain two 42pr gates, and assemble every prime-count sector."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import prod, isqrt
import argparse
import hashlib
import json
from exact import need, alternate_log, out_decimal
ROOT = Path(__file__).resolve().parents[1]

def prime(n):
    return n >= 2 and all((n % d for d in range(2, isqrt(n) + 1)))

def factors(n):
    result = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return result

def expression_interval(argument, constant, base):
    lo, hi = alternate_log(F(argument), terms=160)
    dl, dh = alternate_log(F(base), terms=160)
    return ((lo + F(constant)) / dh, (hi + F(constant)) / dl)

def independent_integral(row, periodic):
    """Different log series and finer square-root bounds for the final integral."""
    Q, k2, z = (F(row['base_min']), F(row['kappa_squared']), F(row['cut_z']))
    a, b, d = [F(periodic[k]) for k in ['quadratic_coefficient', 'linear_coefficient', 'constant_coefficient']]
    scale = 10 ** 40
    t = isqrt((k2 * scale * scale / Q).numerator // (k2 * scale * scale / Q).denominator)
    sl, sh = (F(t, scale), F(t + 1, scale))
    need(sl * sl <= k2 / Q < sh * sh, 'independent square-root enclosure')
    ll, lh = alternate_log(k2 * Q / z ** 2, terms=160)
    zl, zh = alternate_log(z, terms=160)
    ql, qh = alternate_log(Q, terms=160)
    E = a / 2 + b
    cl = a * k2 / (2 * z * z) + b * sl / z + d * ll / (2 * Q) - E / Q
    ch = a * k2 / (2 * z * z) + b * sh / z + d * lh / (2 * Q) - E / Q
    need(F(row['full_correction_lower']) < cl <= ch < F(row['full_correction_upper']), 'independent finite-endpoint correction enclosure')
    lo, hi = ((zl + cl) / qh, (zh + ch) / ql)
    need(F(row['lower']) < lo <= hi < F(row['upper']), 'independent normalized finite-endpoint integral enclosure')
    need(ll > 1 and b * b * k2 * Q > 4 * z * z * E * E, 'all-real-M monotonicity with both endpoint subtractions')
    need(row['monotone_in_M_for_every_real_M_ge_Q'] and row['tail_constant_is_not_integrated_to_infinity'] and (row['finite_endpoint_subtraction_E'] == str(E)), 'finite cutoff and subtraction are explicitly preserved')
    return dict(base_min=int(Q), cut_z=str(z), kappa_squared=str(k2), upper=row['upper'], upper_decimal=row['upper_decimal'], independent_log_series_terms=160, square_root_decimal_places=40, finite_endpoint_enclosure_checked=True, all_real_M_monotonicity_checked=True)

def certificate(root=ROOT):
    root = Path(root)
    inputs = []

    def bind(name):
        path = root / name
        inputs.append(dict(path=name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        return path

    def read(name):
        return json.loads(bind(name).read_text())
    for name in ['code/coverage.py', 'code/exact.py', 'code/assemble_family.py', 'global/verify_global.py', 'global/PROOF.md', 'repair42pr/PROOF.md', 'moments/PROOF.md', 'global/imports/TAIL_PROOF_4830.md', 'global/imports/H132_PROOF_4830.md', 'global/imports/PERIODIC_CDF_4290.md']:
        bind(name)
    old = read('global/imports/coverage_5610_102p.json')
    old_global = read('global/imports/global_5610_102p_certificate.json')
    source = read('global/imports/SOURCE_BINDINGS.json')
    source_manifest = read('global/imports/SOURCE_MANIFEST_5610_102P.json')
    need(hashlib.sha256((root / 'global/imports/SOURCE_MANIFEST_5610_102P.json').read_bytes()).hexdigest() == 'a0c939b13f16ada60e35d3e881c47da90689658ccb5756e1e2efcc1d76c3a48e', 'complete preceding frozen root manifest identity')
    new = read('results/closure_certificate.json')
    family = read('results/family154_certificate.json')
    g = read('global/results/global_certificate.json')
    replay = read('logs/global_replay.json')
    need(new['status'] == 'PASS' and new['radical'] == 6006 and (new['primes'] == [2, 3, 7, 11, 13]) and new['all_positive_valuations'] and new['all_N_ge_1'] and (F(new['Gamma6006']) < F(2, 85) < F(new['candidate']['c0_lower'])) and new['completed_full_geometric_replay']['completed_full_replay_verified'], 'new entire five-prime support has completed full primary geometry')
    need(family['status'] == 'PASS' and family['core'] == 154 and (family['core_primes'] == [2, 7, 11]) and (family['prime_minimum'] == 13) and family['all_positive_valuations'] and family['all_N_ge_1'] and (F(family['Gamma154']) < F(2, 85)) and family['completed_full_replay_verified'] and family['completed_independent_physical_consumption_verified'], 'new four-prime family has full geometry and independent physical consumption')
    need(old['new_5610_premise_paid'] and old['new_102p_family_premise_paid'] and (old['first_open_four_prime_radical'] == 2002) and (old['first_five_prime_radical_not_closed_vs_c0'] == 6006), 'historical premises and frontiers imported explicitly')
    historical = [(r['core'], r['new_prime_minimum']) for r in old['declared_four_prime_families']]
    need(historical == [(30, 7), (42, 11), (70, 11), (66, 13), (105, 11), (78, 17), (110, 13), (102, 19)], 'complete historical family inventory')
    families = historical + [(154, 13)]
    ps4 = [p for p in range(2, 2145 // 30 + 1) if prime(p)]
    four = sorted(((prod(t), list(t)) for t in combinations(ps4, 4) if prod(t) <= 2145))
    direct4 = [(n, factors(n)) for n in range(2, 2146) if len(factors(n)) == 4 and prod(factors(n)) == n]
    need(four == direct4 and len(four) == 53, 'two exhaustive four-prime support censuses')
    covered4 = []
    for q, support in four:
        matches = [c for c, pmin in families if q % c == 0 and prime(q // c) and (q // c >= pmin) and (q // c not in factors(c))]
        if q < 2145:
            need(bool(matches), 'every preceding four-prime support covered')
            covered4.append(dict(radical=q, primes=support, family=str(matches[0]) + 'p'))
        else:
            need(not matches and support == [3, 5, 11, 13], 'first residual four-prime support is odd')
    need([r for r in covered4 if r['radical'] < 2002] == old['covered_four_prime_radicals_below_first_open'], 'all historical four-prime coverage retained')
    need(any((r['radical'] == 2002 and r['family'] == '154p' for r in covered4)), 'new family pays the previous first residual support')
    closed5 = old['five_prime_radicals_closed_below_2_over_85'] + [6006]
    need(closed5 == [2310, 2730, 3570, 3990, 4290, 4830, 5610, 6006], 'eight complete five-prime closures')
    ps5 = [p for p in range(2, 6090 // 210 + 1) if prime(p)]
    five = sorted(((prod(t), list(t)) for t in combinations(ps5, 5) if prod(t) <= 6090))
    direct5 = [(n, factors(n)) for n in range(2, 6091) if len(factors(n)) == 5 and prod(factors(n)) == n]
    need(five == direct5 and [q for q, _ in five] == closed5 + [6090] and (five[-1][1] == [2, 3, 5, 7, 29]), 'two exhaustive five-prime censuses including the support without5')
    need(g['status'] == 'PASS_WITH_DECLARED_HISTORICAL_AND_6006_154P_PREMISES', 'explicit conditional component becomes paid only in the master')
    premises = [dict(kind='fixed_radical', radical=6006, all_positive_valuations=True, all_N_ge_1=True, bound='Delta(B,N)<=2/85'), dict(kind='prime_family', core=154, core_primes=[2, 7, 11], prime_minimum=13, all_positive_valuations=True, all_N_ge_1=True, bound='Delta(B,N)<=2/85')]
    need(g['required_new_premises'] == premises, 'both exact new global premises are paid')
    need(g['all_positive_valuations'] and g['all_N_ge_1'] and (g['normalization'] == new['normalization'] == family['normalization']), 'same physical all-length normalization in all components')
    need(replay['status'] == 'PASS' and replay['component_bytes_unchanged'] and (replay['certificate_sha256'] == hashlib.sha256((root / 'global/results/global_certificate.json').read_bytes()).hexdigest()) and (replay['exact_assertions'] == g['exact_assertions'] == 30417), 'executed arithmetic replay uses every field of this final certificate')
    front = g['coverage']
    need(front['declared_four_prime_families'] == [dict(core=c, new_prime_minimum=p) for c, p in families] and front['first_open_four_prime_radical'] == 2145 and (front['first_open_five_prime_radical_after_6006'] == 6090) and ([(r['radical'], r['prime_support']) for r in front['all_four_prime_radicals_through_frontier']] == four) and ([(r['radical'], r['prime_support']) for r in front['all_five_prime_radicals_through_frontier']] == five), 'global component and independent master censuses agree')
    expr = '(log(20/3)+1261/2500)/log(6090)'
    need(g['global_upper_expression'] == expr and g['simple_strict_upper'] == '7/25', 'new clean global expression')
    gl, gh = expression_interval(F(20, 3), F(1261, 2500), 6090)
    ulo, uhi = (F(g['global_lower_enclosure']), F(g['global_upper']))
    need(ulo < gl <= gh < uhi < F(7, 25) and g['global_upper_decimal'] == out_decimal(uhi), 'independent160-term enclosure of the clean expression')
    periodic = g['periodic_CDF']
    need(periodic == old_global['periodic_CDF'] and periodic['period'] == 30030 and (periodic['weights_checked'] == 30030) and periodic['two_weight_algorithms_agree'] and periodic['all_M_ge_1'] and periodic['all_real_x_ge_0'], 'periodic arithmetic weights retained with universal scope')
    five_sector = g['five_prime_sector']
    others = g['other_finite_sectors']
    tail = g['infinite_tail']
    need((five_sector['omega'], five_sector['base_min'], five_sector['H']) == (5, 6090, 10) and F(five_sector['kappa_squared']) == F(103, 3), 'five-prime residual with new base minimum')
    need([(r['omega'], r['base_min'], r['H']) for r in others] == [(4, 2145, 6), (6, 30030, 14), (7, 510510, 22), (8, 9699690, 26), (9, 223092870, 60), (10, 6469693230, 132)], 'all remaining finite prime-count sectors')
    primes = [p for p in range(2, 38) if prime(p)]
    for row in others:
        need(F(row['kappa_squared']) == F(row['H'] ** 2 + 3, 3), 'uniform repaired witness coefficient')
        if row['omega'] >= 6:
            need(row['base_min'] == prod(primes[:row['omega']]), 'exact minimum primorial')
    J = g['new_J132']
    need(J == old_global['new_J132'] and J['H'] == 132 and (J['residual_primes_at_most'] == 9) and J['all_forbidden_classes'] and J['all_integer_intervals'] and all((r['upper_covered'] < r['length'] for r in J['cases'])), 'inherited uniform covering retained')
    need(tail['from_omega'] == 11 and tail['base_min'] == prod(primes[:11]) == 200560490130 and (tail['primorial_primes'] == primes[:11]) and (tail['next_prime_lower'] == 37) and (tail['H_at_anchor'] == 1024) and (F(tail['kappa_squared']) == F(1048579, 3)) and (F(tail['cut_z']) == 669) and tail['ratio_kappa_over_cut_is_constant'] and tail['monotone_for_every_real_v_ge_0'] and (F(tail['derivative_v_upper']) < 0) and (not tail['finite_cardinality_scan_used_as_proof']), 'unbounded prime-count tail is analytic')
    _, h2 = alternate_log(2, terms=160)
    _, ph = alternate_log(tail['base_min'], terms=160)
    rl, _ = alternate_log(37, terms=160)
    need(h2 * ph - rl * F(tail['C_lower']) < 0, 'independent sign for all-real-v derivative')
    integrals = [independent_integral(r, periodic) for r in [five_sector, *others, tail]]
    stronger = F(five_sector['upper'])
    previous = F(old['unrounded_integration_upper'])
    need(stronger == F(g['stronger_bound_before_clean_rounding_upper']) and max((F(r['upper']) for r in [*others, tail])) < stronger < ulo and (stronger < previous), 'new exact integral dominates all sectors and improves preceding cap')
    need(F(new['candidate']['c0_upper']) < F(1, 40) < stronger and F(five_sector['full_correction_upper']) < F(1261, 2500), 'candidate and clean correction checked')
    ge6 = max([F(r['upper']) for r in others if r['omega'] >= 6] + [F(tail['upper'])])
    need(ge6 == F(tail['upper']) == F(g['uniform_omega_ge_6_upper']), 'uniform omega>=6 cap')
    far = old['inherited_far_tail']
    need(F(far['prime_count_upper']) < F(2, 85) and F(far['large_base_upper']) < F(2, 85), 'historical far-tail exclusions retained')
    joint = read('repair42pr/certificate.json')
    joint_replay = read('logs/repair42pr_replay.json')
    need(joint['status'] == 'PASS_UNIFORM_REPAIR_AND_TWO_INFINITE_GATES_WITH_DECLARED_HISTORICAL_ANCHORS' and joint['exact_assertions'] == joint_replay['exact_assertions'] == 118341 and (joint_replay['status'] == 'PASS') and (joint_replay['certificate_sha256'] == hashlib.sha256((root / 'repair42pr/certificate.json').read_bytes()).hexdigest()), 'new double-prime repair has a complete exact replay')
    gates = joint['infinite_gates']['gates']
    need([(r['p_and_r_minimum'], r['all_N_at_least'], r['anchor_divisor']) for r in gates] == [(229, 1, 42), (83, 2, 1764)] and all((r['all_positive_valuations'] for r in gates)) and all((F(r['normalized_upper']) < F(2, 85) for r in gates)) and (not joint['infinite_gates']['full_family_closed']) and (not joint['infinite_gates']['anchor_1764_is_negative']), 'exact scope of both five-prime infinite gates, without sign inflation')
    moments = read('moments/certificate.json')
    need(moments['status'] == 'PASS' and moments['all_M_negative_part_exact_from_two_moments'] and (not moments['global_cap_changed']) and (not moments['external_recent_theorems_used_as_premises']), 'moment identity is proved but no numerical global gain inferred from its finite examples')
    return dict(status='PASS_WITH_DECLARED_HISTORICAL_DEPENDENCIES', newly_closed_radical=6006, new_6006_premise_paid=True, new_154p_family_premise_paid=True, new_fixed_certificate_sha256=hashlib.sha256((root / 'results/closure_certificate.json').read_bytes()).hexdigest(), new_family_certificate_sha256=hashlib.sha256((root / 'results/family154_certificate.json').read_bytes()).hexdigest(), all_positive_valuations=True, all_N_ge_1=True, normalization=new['normalization'], first_open_four_prime_radical=2145, first_open_four_prime_support=[3, 5, 11, 13], declared_four_prime_families=[dict(core=c, new_prime_minimum=p) for c, p in families], historical_30p_includes_candidate_and_is_not_claimed_strictly_below_2_over_85=True, covered_four_prime_radicals_below_first_open=covered4, previous_first_open_five_prime_radical=6006, first_five_prime_radical_not_closed_vs_c0=6090, first_open_five_prime_support=[2, 3, 5, 7, 29], five_prime_radicals_closed_below_2_over_85=closed5, exact_five_prime_frontier=[dict(radical=q, prime_support=s, closed_vs_c0=q in closed5) for q, s in five], two_exhaustive_radical_inventories_agree=True, new_42pr_infinite_regions=gates, full_42pr_family_closed=False, moment_negative_part_identity_proved=True, moment_identity_has_not_yet_improved_the_global_numeric_cap=True, global_lower_candidate=new['candidate']['c0_lower'], global_upper_expression=expr, global_upper=str(uhi), global_upper_decimal=out_decimal(uhi), simple_strict_global_upper='7/25', global_numerical_bound_changed_this_task=True, previous_global_upper_decimal=old['global_upper_decimal'], previous_unrounded_integration_upper_decimal=old['unrounded_integration_upper_decimal'], unrounded_integration_upper=str(stronger), unrounded_integration_upper_decimal=out_decimal(stronger), independent_high_precision_logarithm_series_used=True, independent_finite_endpoint_integrals=integrals, reused_periodic_CDF=periodic, rechecked_uniform_covering=J, other_prime_count_sectors=others, infinite_prime_count_tail=tail, all_omega_ge_6_upper=str(ge6), all_omega_ge_6_upper_decimal=out_decimal(ge6), inherited_far_tail=far, previous_global_and_sector_proofs_are_explicit_imports=True, previous_large_geometry_or_covering_replayed_this_task=False, all_four_prime_supports_closed=False, all_five_prime_supports_closed=False, uniform_210p_family_closed=False, exact_Cstar_equals_c0_proved=False, recent_literature_theorems_used_as_external_numerical_premises=False, offset_sign_classification_claimed=False, next_four_prime_family='rad(B)=165p,p>=13 prime; proposed', next_five_prime_target='6090=2*3*5*7*29; no full profile produced in this task', remaining_domain_for_a_value_greater_than_c0=old['remaining_domain_for_a_value_greater_than_c0'], inputs=inputs)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(certificate(), indent=2) + '\n', newline='\n')
