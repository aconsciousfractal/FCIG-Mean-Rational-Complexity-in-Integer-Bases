"""Exact all-valuation/all-length assembly for radical3570 and3990.

The master pays every named component and the two new full geometric replays.
This arithmetic assembler binds actual files and never promotes a JSON alone
to a geometric or continuous proof.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
from exact import need, log_lo, out_decimal, prime_trial
from candidate210 import certificate as candidate
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def certificate():
    inputs = []

    def read(rel):
        path = ROOT / rel
        inputs.append(dict(path=rel, sha256=sha(path)))
        return json.loads(path.read_text())

    def bind(rel):
        value = sha(ROOT / rel)
        inputs.append(dict(path=rel, sha256=value))
        return value
    c = candidate()
    need(c == read('results/candidate210.json'), 'candidate reconstructed from its definition')
    delta = F(2, 85)
    c0lo = F(c['c0_lower'])
    need(delta < c0lo, 'strict candidate separation')
    theory = read('imports/previous2730/theory/theory_certificate.json')
    need(theory['status'] == 'PASS' and theory['cdf']['alpha'] == '3/2' and theory['cdf']['M1_included'], 'primitive CDF has every positive seed includingone')
    bind('imports/previous2730/theory/PRIMITIVE_CDF.md')
    bind('imports/previous2730/theory/imports/UNEQUAL_VALUATION_TRANSFER.md')
    plan = read('profiles/profiles_plan.json')
    need(len(plan) == 2 and sorted((r['p'] for r in plan)) == [17, 19], 'exactlytwo current full marks')
    inv = read('physical/results/finite_inventory.json')
    finite = read('physical/results/finite_certificate.json')
    need(inv['status'] == finite['status'] == 'PASS', 'physical inventory and consumer passed')
    need(inv['states'] == finite['states'] == 14700 and inv['representatives'] == finite['representatives'] == 7353, 'full finite state and representative masses')
    need(finite['moduli'] == 3 and finite['normalization'] == 'mu(B^N)/log(B)', 'finite scope and normalization')
    need([(r['M'], r['B'], r['radical'], r['N']) for r in inv['rows']] == [(3570, 3570, 3570, 1), (7140, 7140, 3570, 1), (3990, 3990, 3990, 1)], 'exact physical domain')
    sectors = []
    for q, cutoff, bases in [(3570, 10000, [3570, 7140]), (3990, 6000, [3990])]:
        p = q // 210
        primes = [2, 3, 5, 7, p]
        need(all((prime_trial(x) for x in primes)) and len(set(primes)) == 5, 'complete distinct prime support')
        means = []
        row = next((r for r in plan if r['p'] == p))
        need(row['multiplier'] == 2, 'new index is exactly2q')
        oldstem = f'q{q}_K{q}_coarse'
        for parent, stem, K in [('imports/previous2730/profiles', oldstem, q), ('profiles', row['stem'], 2 * q)]:
            base = read(f'{parent}/results/{stem}.json')
            alt = read(f'{parent}/results/{stem}.alpha3over2.json')
            need(base['status'] == 'PASS' and base['q'] == q and (base['transfer_index'] == K), 'full marked geometric identity')
            need(base['required_divisor_of_M'] == K and base['all_positive_multiples_including_extra_primes'] and base['all_residues_including_zero'], 'all-multiple, all-state scope')
            need(base['normalization'] == 'log B', 'profile normalization')
            need(alt['q'] == q and alt['transfer_index'] == K and (alt['alpha'] == '3/2'), 'exact alternate CDF scope')
            need(alt['status'] == 'ARITHMETIC_PASS_CONDITIONAL_ON_CDF_LEMMA' and alt['unchanged_geometric_contract'], 'explicit CDF premise and unchanged continuous cells')
            need(alt['base_total_upper'] == base['total_upper'], 'same original consumed profile')
            cellhash = bind(f'{parent}/results/{stem}.cells')
            need(alt['source_cells_sha256'] == cellhash, 'same complete continuous cells')
            A = F(alt['total_upper'])
            need(A > 0, 'positive explicit upper mean')
            means.append(dict(index=K, mean_upper=str(A), mean_upper_decimal=out_decimal(A), stem=stem, cell_sha256=cellhash))
        A1, A2 = [F(r['mean_upper']) for r in means]
        invq = next((r for r in inv['domains'] if r['radical'] == q))
        physq = next((r for r in finite['domains'] if r['radical'] == q))
        need(invq['prime_support'] == primes and invq['base_upper_exclusive'] == physq['base_upper_exclusive'] == cutoff and invq['all_positive_valuations'], 'finite and infinite domains meet')
        need([r['B'] for r in invq['rows']] == bases and physq['moduli'] == len(bases), 'complete finite complement')
        physical = F(physq['normalized_upper'])
        large = A1 / log_lo(cutoff)
        high = A2 / log_lo(q)
        need(physical < delta and large < delta and (high < delta), 'allthree exhaustive simple regions')
        need(q * q % (2 * q) == 0, '2q divides every B^N for N>=2')
        simple = max(physical, large, high)
        need(q % 4 == 2 and min(primes[1:]) == 3, 'exact odd-quotient minimum')
        physical_anchor = next((r for r in finite['rows'] if r['B'] == r['M'] == q and r['radical'] == q))
        anchor = F(physical_anchor['normalized_upper'])
        odd = A1 / log_lo(3 * q)
        gamma = max(anchor, high, odd)
        need(gamma <= simple and gamma < delta < c0lo, 'refined exhaustive valuation comparison')

        def case(label, u, **extra):
            return dict(case=label, upper=str(u), upper_decimal=out_decimal(u), **extra)
        sectors.append(dict(radical=q, prime_support=primes, all_positive_valuations=True, all_N_ge_1=True, mean_profiles=means, simple_cases=[case(f'N=1 and B<{cutoff}', physical, bases=bases), case(f'N=1 and B>={cutoff}', large), case('N>=2', high)], simple_Gamma=str(simple), simple_Gamma_upper_decimal=out_decimal(simple), refined_cases=[case(f'B={q}, N=1', anchor), case('N*v2(B)>=2', high), case(f'N=1,v2(B)=1,B!={q}; hence B>={3 * q}', odd)], refined_case_partition_exhaustive=True, geometric_domain='2q|M; rad(M)=q for identification with mu', Gamma=str(gamma), Gamma_upper_decimal=out_decimal(gamma), exact_margin_to_2_over_85=str(delta - gamma), exact_margin_to_c0_lower=str(c0lo - gamma), redundant_physical_moduli=[7140] if q == 3570 else []))
    return dict(status='PASS_WITH_DECLARED_COMPONENT_DEPENDENCIES', theorem='For rad(B) in {3570,3990}, Delta(B,N)<=Gamma_rad(B)<2/85<c0 for every N>=1', radicals=[3570, 3990], all_positive_valuations=True, all_N_ge_1=True, all_states_including_zero=True, norm='max(abs(a),d)', logarithms='natural', normalization='mu(B^N)/log(B)', simple_strict_upper='2/85', c0_lower=str(c0lo), sectors=sectors, inputs=inputs, required_components=['profiles', 'physical', 'theory', 'audit', 'inherited primitive CDF and first-index profiles'], this_assembler_performs_geometric_replay=False, valuation_split_is_not_monotonicity_in_q=True, offset_sign_classification_claimed=False, exact_Cstar_equals_c0_proved=False, formalized=False)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    z = certificate()
    s = json.dumps(z, indent=2) + '\n'
    if a.output:
        a.output.write_text(s, newline='\n')
        print(json.dumps({'status': z['status'], 'radicals': [{'q': r['radical'], 'Gamma': r['Gamma_upper_decimal']} for r in z['sectors']]}, indent=2))
    else:
        print(s, end='')
