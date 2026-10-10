"""Exact standalone assembly of the shaped-CDF global supplement.

The theory peer supplies frozen proof and rational-logarithm arithmetic.
This script reads peers, performs exact calculations and never writes to
them. The complete historical/profile replay obligations remain explicit.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, importlib.util, json
ROOT = Path(__file__).resolve().parent
CHECKS = 0

def need(c, m):
    global CHECKS
    CHECKS += 1
    if not c:
        raise AssertionError(m)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def certificate(theory):
    global CHECKS
    CHECKS = 0
    bindings = json.loads((ROOT / 'THEORY_BINDINGS.json').read_text())
    actual_inputs = []
    for item in bindings['files']:
        p = theory / item['theory_relative_path']
        need(sha(p) == item['sha256'], 'immutable theory source: ' + item['theory_relative_path'])
        actual_inputs.append(dict(path='../theory/' + item['theory_relative_path'], sha256=item['sha256']))
    manifest = json.loads((theory / 'MANIFEST.json').read_text())
    for item in manifest['files']:
        p = theory / item['path']
        need(p.stat().st_size == item['bytes'] and sha(p) == item['sha256'], 'entire frozen theory component')
    sp = importlib.util.spec_from_file_location('p72_shaped_cdf_exact_arithmetic', theory / 'verify_theory.py')
    T = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(T)
    prior = json.loads((theory / 'global_covered_certificate.json').read_text())
    cdf = json.loads((theory / 'theory_certificate.json').read_text())['cdf']
    need(prior['status'] == 'PASS_WITH_DECLARED_COMPONENT_DEPENDENCIES', 'previous global assembly status')
    need(cdf['status'] == 'PASS' and cdf['alpha'] == '3/2' and cdf['M1_included'], 'PCDF includes every seedM>=1')
    need(prior['five_prime_frontier'] == 4290 and prior['five_prime_sector']['H'] == 10, 'exact shaped-sector domain')
    need(prior['all_states_including_zero'] and prior['all_positive_valuations'] and prior['all_N_ge_1'], 'all previous model quantifiers')
    need(prior['normalization'] == 'mu(B^N)/log(B)', 'normalization stayslogB')
    need([a['q'] for a in prior['new_anchor_profiles']] == [2730, 3570, 3990], 'three paid anchor supports')
    nonlinear_a = F(29, 21)
    nonlinear_b = F(11, 12)
    need(nonlinear_a > 1 and nonlinear_b > 0, 'nonlinear CDF extends bytrivial bounds outside1<=Y<M')
    H = prior['five_prime_sector']['H']
    kappa2 = F(H * H + 3, 3)
    need(kappa2 == F(103, 3), 'full-mark W repair parameter')
    root_kappa_upper = 6
    root_modulus_lower = 65
    need(kappa2 < root_kappa_upper ** 2, 'strict kappa<6')
    need(root_modulus_lower ** 2 < prior['five_prime_frontier'], 'strictsqrtM>65 fromM>=4290')
    A = nonlinear_a * kappa2
    C = nonlinear_b * F(root_kappa_upper, root_modulus_lower)
    need(A == F(2987, 63) and C == F(11, 130), 'tail constants from exact rational algebra')
    cut = 7
    correction = A / (2 * cut * cut) + C / cut
    need(correction == F(2987, 6174) + F(11, 910) == F(99503, 200655), 'layercake integrals atz7')
    need(A / (cut * cut) + C / cut == F(393161, 401310) < 1, 'chosen cut lies beyondunit cap')
    l7, h7 = T.ln_series(cut)
    lb, hb = T.ln_series(prior['five_prime_frontier'])
    need(l7 > 0 and lb > 0, 'positive logarithms andnormalization')
    lower = (l7 + correction) / hb
    upper = (h7 + correction) / lb
    need(upper < F(73, 250), 'strict73/250 gate')
    need(upper.numerator * 250 < upper.denominator * 73, 'samegate by integercrossmultiplication')
    other = []
    for row in prior['other_sectors']:
        u = F(row['upper'])
        need(u < lower, 'shaped five-prime sector dominates prior sector')
        other.append(dict(omega=row['omega'], upper=row['upper'], upper_decimal=row['upper_decimal']))
    anchors = []
    for row in prior['new_anchor_profiles']:
        need(F(row['normalized_upper']) < lower, 'shapedsector dominates existing anchor')
        anchors.append(dict(q=row['q'], upper=row['normalized_upper'], source=row['input']))
    need(F(1, 40) < lower, 'shapedsector dominates imported c0 upper')
    need(upper < F(prior['global_upper']), 'newshape genuinely improves previousglobalupper')
    tail = prior['infinite_tail']
    need(F(tail['derivative_upper']) < 0 and tail['from_omega'] == 10, 'complete infinite-tail proof retained')
    need(F(prior['uniform_omega_ge_6']['upper']) < lower, 'everyomega>=6 dominated')
    return dict(status='PASS_WITH_DECLARED_BASELINE_DEPENDENCIES', theorem='c0<=Cstar<=Ushape=(log7+99503/200655)/log4290<73/250', exact_global_upper_expression='(log(7)+99503/200655)/log(4290)', global_lower_enclosure=str(lower), global_upper=str(upper), global_upper_decimal=T.upper_decimal(upper), simple_strict_upper='73/250', nonlinear_CDF=dict(theorem='For every M>=1,x>=0: Pr(lambda_s^2/M<=x)<=min(1,3*x/2,29*x/21+11*sqrt(x/M)/12)', all_integer_moduli_M_ge_1=True, all_real_thresholds_x_ge_0=True, middle_domain='1<=Y=sqrt(Mx)<M: formula6of frozenPCDFproof dividedbyM', small_domain='Y<1: eventempty', large_domain='Y>=M: x>=M>=1 and29x/21>1', source_formula_number='(6)', linear_CDF_sharpness_not_contradicted=True), improved_sector=dict(omega=5, base_min=4290, all_N_ge_1=True, modulus_relation='M=B^N>=B>=4290', H=H, kappa_squared=str(kappa2), tail_A=str(A), strict_tail_C_upper=str(C), cut_z=cut, rational_mean_correction=str(correction), strict_mean_bound_expression='log(7)+99503/200655', mean_upper_decimal=T.upper_decimal(h7 + correction)), normalization='mu(B^N)/log(B)', all_states_including_zero=True, all_positive_valuations=True, all_N_ge_1=True, other_sectors_preserved=other, anchor_bounds_preserved=anchors, candidate_upper_import='c0<1/40', infinite_tail_preserved=tail, integer_gate=dict(left=str(upper.numerator * 250), right=str(upper.denominator * 73), positive_difference=str(upper.denominator * 73 - upper.numerator * 250)), theory_inputs=actual_inputs, exact_assertions=CHECKS, prior_global_upper=str(prior['global_upper']), new_geometric_flows=0, required_baseline_dependencies=prior['required_component_validation'], this_supplement_performs_geometry_or_covering_replays=False, exact_Cstar_equals_c0_proved=False, formalized=False)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--theory', type=Path, default=ROOT.parent / 'theory')
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    z = certificate(a.theory)
    s = json.dumps(z, indent=2) + '\n'
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(s, newline='\n')
        print(json.dumps({k: z[k] for k in ['status', 'global_upper_decimal', 'simple_strict_upper', 'exact_assertions']}, indent=2))
    else:
        print(s, end='')
