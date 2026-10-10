"""Exact supplementary gates; no geometry or witness search replay."""
from pathlib import Path
from fractions import Fraction as F
from math import prod, factorial
import argparse, hashlib, json
HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(label)

def ceildiv(a, b):
    return (a + b - 1) // b

def rotated_sq(v, t):
    x, y = v
    return max((x + t * y) ** 2, (y - t * x) ** 2) / (1 + t * t)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'supplement_certificate.json')
    args = parser.parse_args()
    sources = json.loads((HERE / 'PROVENANCE.json').read_text())['sources']
    for source in sources:
        raw = (HERE / source['local_copy']).read_bytes()
        need(len(raw) == source['bytes'] and hashlib.sha256(raw).hexdigest() == source['sha256'], 'source bridge')
    odd = [ceildiv(76, p) for p in [3, 5, 7, 11, 13, 17, 19]]
    even = [ceildiv(38, p) for p in [3, 5, 7, 11, 13, 17]]
    need(odd == [26, 16, 11, 7, 6, 5, 4] and sum(odd) == 75 < 76, 'H76 odd case')
    need(even == [13, 8, 6, 4, 3, 3] and sum(even) == 37 < 38, 'H76 even case')
    eu = sum((F(1, factorial(n)) for n in range(9)), F()) + F(1, 8 * factorial(8))
    need(eu == F(876809, 322560) < F(49, 18), 'exponential upper')
    need(sum((F(1, factorial(n)) for n in range(4)), F()) == F(8, 3), 'exponential lower')
    need(F(29 * (76 ** 2 + 3), 42) == F(167591, 42), 'K8 coefficient')
    need(F(29 * (256 ** 2 + 3), 42) == F(1900631, 42), 'K9 coefficient')
    need(prod([2, 3, 5, 7, 11, 13, 17, 19, 23]) == 223092870, 'P9 product')
    need(F(4) ** 2 < F(8, 3) ** 3, 'ln4<3/2')
    need(F(223092870) < F(8, 3) ** 20, 'lnP9<20')
    need(F(29) > F(49, 18) ** 3, 'ln29>3')
    need(F(1900631, 42) > F(49, 18) ** 10, 'lnK9>10')
    derivative_upper = F(3, 2) * 20 - 3 * (1 + 10)
    need(derivative_upper == -3 < 0, 'tail derivative')
    rows = []
    for k, P, K in [(7, 510510, F(41963, 42)), (8, 9699690, F(167591, 42)), (9, 223092870, F(1900631, 42))]:
        left = (49 * K.numerator) ** 50
        right = P ** 31 * (18 * K.denominator) ** 50
        need(left < right, 'sector below31/100')
        rows.append(dict(k=k, P=P, K=str(K), integer_left=left, integer_right=right, strict_upper='31/100'))
    lower_left = (8 * 3161) ** 50
    lower_right = 30030 ** 31 * (3 * 14) ** 50
    need(lower_left > lower_right, 'U6 above31/100')
    upper_left = (49 * 3161) ** 8
    upper_right = 30030 ** 5 * (18 * 14) ** 8
    need(upper_left < upper_right, 'U6 below5/16')
    vectors = [(2, 1), (2, -1)]
    minima = []
    for t in [-F(1, 10), F(1, 10)]:
        values = [rotated_sq(v, t) for v in vectors]
        need(min(values) == F(361, 101), 'endpoint squared norm')
        minima.append([str(x) for x in values])
    need(min((rotated_sq(v, F(0)) for v in vectors)) == 4, 'midpoint squared norm')
    need(101 > 100, 'sqrt101>10')
    need(F(1, 19) / F(20, 19) == F(1, 20) > F(1, 100), 'log gap exceeds curvature budget')
    need(all(((x - 2 * y) % 4 == 0 and y % 2 for x, y in vectors)), 'two admissible marked lattice witnesses')
    old_vectors = [(3, 1), (3, -1)]
    old = [[rotated_sq(v, t) for t in [-F(1, 10), F(1, 10)]] for v in old_vectors]
    bank = [rotated_sq((3, 0), t) for t in [-F(1, 10), F(1, 10)]]
    need(bank == [F(900, 101)] * 2 and all((max(bank) < max(row) for row in old)), 'branchwise greedy endpoint improvement')
    need(all((old[0][i] * old[1][i] == F(899, 101) ** 2 < bank[i] ** 2 for i in [0, 1])), 'branchwise choices can realign logarithmic means')
    output = dict(status='PASS_WITH_DECLARED_IMPORTS', checks=CHECKS, higher_prime_tail=dict(theorem='omega(B)>=6 implies Delta<=U6<5/16', U6='(1+ln(3161/14))/(2ln30030)', all_valuations=True, all_N_ge1=True, odd_cover_counts=odd, even_cover_counts=even, sector_gates=rows, lower_U6=dict(integer_left=lower_left, integer_right=lower_right, strict_lower='31/100'), upper_U6=dict(integer_left=upper_left, integer_right=upper_right, strict_upper='5/16'), F9_derivative_strict_upper=str(derivative_upper), Kanold_only_k_ge9=True, geometry2310_required=False), angular_bank=dict(written_continuous_proof_required=True, implemented_profile_claimed=False, same_witness_both_endpoints=True, min_only_after_continuous_upper=True, negative_node_min_regression=dict(vectors=vectors, endpoint_norm_squares=minima, center_min_norm_square='4', log_gap_strict_lower='1/20', curvature_strict_upper='1/100'), greedy_aggregation_regression=dict(old_vectors=old_vectors, new_vector=[3, 0], old_exp_mean_square='899/101', bank_exp_mean_square='900/101', scope='actual vector norms, not a complete HNF transfer profile')), sources=sources, historical_imports=['SV/W/MEAN', 'Jacobsthal18,38', 'Kanold only k>=9'], historical_replays_performed=False, geometry_replayed=False, finite_checks_replace_infinite_proofs=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status=output['status'], checks=CHECKS, U6_upper='5/16', node_min_negative=True), indent=2))
if __name__ == '__main__':
    main()
