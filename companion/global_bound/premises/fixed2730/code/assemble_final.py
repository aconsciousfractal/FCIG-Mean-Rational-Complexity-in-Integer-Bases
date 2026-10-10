"""Exact arithmetic assembly for all positive valuations of radical 2730.

Continuous and computational premises are paid by the master verifier and
the named components. This script never treats an unbound JSON cap as a
standalone geometric proof.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
from exact import need, log_lo, out_decimal
from candidate210 import certificate as candidate
ROOT = Path(__file__).resolve().parents[1]

def load(relative):
    return json.loads((ROOT / relative).read_text())

def digest(relative):
    return hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()

def certificate():
    inputs = []

    def read(p):
        inputs.append(dict(path=p, sha256=digest(p)))
        return load(p)
    c = candidate()
    need(c == read('results/candidate210.json'), 'definition-based candidate is reconstructed')
    delta = F(2, 85)
    c0lo = F(c['c0_lower'])
    need(delta < c0lo, 'strict separation from candidate')
    theory = read('theory/theory_certificate.json')
    need(theory['status'] == 'PASS' and theory['cdf']['alpha'] == '3/2' and theory['cdf']['M1_included'], 'new primitive CDF includes the physical seed1')
    profile_rows = []
    for stem, K, kind in [('q2730_K2730_coarse', 2730, 'all_M_divisible_by2730'), ('q2730_K16380_hybrid_AG8', 16380, 'all_M_divisible_by16380')]:
        base = read('profiles/results/' + stem + '.json')
        alt = read('profiles/results/' + stem + '.alpha3over2.json')
        need(base['status'] == 'PASS' and base['q'] == 2730 and (base['transfer_index'] == K), 'full marked profile identity')
        need(base['required_divisor_of_M'] == K and base['all_positive_multiples_including_extra_primes'] and base['all_residues_including_zero'], 'exact all-multiple transfer scope')
        need(base['normalization'] == 'log B', 'base normalization preserved')
        need(alt['q'] == 2730 and alt['transfer_index'] == K and (alt['alpha'] == '3/2'), 'alternative profile index and CDF')
        need(alt['base_total_upper'] == base['total_upper'], 'same original profile assembly')
        cells = 'profiles/results/' + stem + '.cells'
        need(alt['source_cells_sha256'] == digest(cells), 'same complete continuous cells')
        inputs.append(dict(path=cells, sha256=digest(cells)))
        a = F(alt['total_upper'])
        need(a > 0, 'positive explicit profile upper')
        profile_rows.append(dict(index=K, mean_upper=str(a), upper_decimal=out_decimal(a), domain=kind, cell_sha256=digest(cells)))
    A1, A6 = (F(r['mean_upper']) for r in profile_rows)
    cutoff = 20000
    large = A1 / log_lo(cutoff)
    high = A6 / log_lo(2730)
    need(large < delta and high < delta, 'both infinite base/length branches lie below2/85')
    need(2730 ** 2 % 16380 == 0, '6q divides B^N for every positive valuation whenN>=2')
    inventory = read('physical/results/finite_inventory.json')
    finite = read('physical/results/finite_certificate.json')
    need(inventory['radical'] == 2730 and inventory['base_upper_exclusive'] == cutoff and inventory['all_positive_valuations'], 'finite domain and large-base cutoff meet')
    expected = [2730, 5460, 8190, 10920, 13650, 16380, 19110]
    need([r['B'] for r in inventory['rows']] == expected and all((r['M'] == r['B'] and r['N'] == 1 for r in inventory['rows'])), 'exact finite N1 domain')
    need(finite['status'] == 'PASS' and finite['radical'] == 2730 and (finite['moduli'] == 7), 'complete finite certificate')
    need(finite['states'] == inventory['states'] == sum(expected) == 76440 and finite['representatives'] == inventory['representatives'] == 38227, 'all physical state mass')
    physical = F(finite['normalized_upper'])
    need(physical < delta, 'finite normalized upper')
    gamma = max(large, high, physical)
    need(gamma == large and gamma < delta < c0lo, 'closed all-valuation all-length comparison')
    return dict(status='PASS_WITH_DECLARED_COMPONENT_DEPENDENCIES', theorem='rad(B)=2730 implies Delta(B,N)<=Gamma2730<2/85<c0 for every N>=1', radical=2730, prime_support=[2, 3, 5, 7, 13], all_positive_valuations=True, all_N_ge_1=True, all_states_including_zero=True, norm='max(abs(a),d)', logarithms='natural', normalization='mu(B^N)/log(B)', profile_mean_uppers=profile_rows, cases=[dict(case='N=1 and B<20000', upper=str(physical), upper_decimal=out_decimal(physical), moduli=7), dict(case='N=1 and B>=20000', upper=str(large), upper_decimal=out_decimal(large)), dict(case='N>=2 and every admissible B', upper=str(high), upper_decimal=out_decimal(high))], Gamma2730=str(gamma), Gamma2730_upper_decimal=out_decimal(gamma), simple_strict_upper='2/85', c0_lower=str(c0lo), exact_margin_to_2_over_85=str(delta - gamma), exact_margin_to_c0_lower=str(c0lo - gamma), inputs=inputs, required_components=['profiles', 'theory', 'physical', 'audit'], this_assembler_performs_geometric_replay=False, classification_of_offset_signs_claimed=False, exact_Cstar_equals_c0_proved=False)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    z = certificate()
    s = json.dumps(z, indent=2) + '\n'
    if a.output:
        a.output.write_text(s, newline='\n')
        print(json.dumps({k: z[k] for k in ('status', 'Gamma2730', 'Gamma2730_upper_decimal', 'simple_strict_upper')}, indent=2))
    else:
        print(s, end='')
