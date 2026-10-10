"""PRELIMINARY_NOT_A_FAMILY_CLOSURE: fresh scalar tables for core154."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, isqrt, prod
from hashlib import sha256
import argparse, json
Q = 154
PRIMES = (2, 7, 11)
PMIN = 13
CHECKS = 0

def need(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise ArithmeticError(label)

def certificate():
    global CHECKS
    CHECKS = 0
    need(prod(PRIMES) == Q and all((all((p % d for d in range(2, isqrt(p) + 1))) for p in PRIMES)), 'squarefree prime core')
    grouped = {}
    obstructions = []
    maxradius = 0
    for z in range(Q):
        if gcd(z, Q) == 1:
            i, j = (1, 0)
        else:
            choices = [i for i in range(1, 5) if gcd(i * z - 1, Q) == 1]
            need(bool(choices), 'whole-half radius4 candidate')
            i, j = (min(choices), -1)
        need(gcd(j + i * z, Q) == 1, 'positive half unit residue')
        need(max(abs(F(j)), abs(F(j) + F(i, 2))) <= 1, 'closed positive half, affine endpoint bound')
        need(gcd(-j - i * z, Q) == 1, 'reflected unit residue')
        need(max(abs(F(-j)), abs(F(-j) - F(i, 2))) <= 1, 'closed negative half, affine endpoint bound')
        grouped.setdefault((i, j), []).append(z)
        maxradius = max(maxradius, i)
        candidates = []
        for h in range(1, 4):
            for k in range(-4, 3):
                if abs(F(k) + F(h, 4)) <= 1:
                    candidates.append((h, k, gcd(k + h * z, Q)))
        need(len(candidates) == 6, 'all radius3 candidate pairs at tau1/4')
        if all((g > 1 for _, _, g in candidates)):
            obstructions.append(dict(z=z, tau='1/4', candidates=[dict(i=h, j=k, gcd=g) for h, k, g in candidates]))
    need(maxradius == 4 and [r['z'] for r in obstructions] == [11, 105], 'sharp radius within the specified ansatz')
    phi = sum((gcd(r, Q) == 1 for r in range(Q)))
    rho = F(phi, Q)
    A = F(5, 2) * rho
    need(phi == 60 and sum((r for r in range(Q) if gcd(r, Q) == 1)) == Q * phi // 2, 'unit count and symmetry')
    U = S = 0
    C = F(0)
    rows = []
    for j in range(Q):
        if j and gcd(j, Q) == 1:
            U += 1
            S += j
        eps = U - rho * j
        D = 2 * j * U + S - A * j * j
        endpoint = D / j if j else F(0)
        rows.append(dict(j=j, U=U, S=S, epsilon=str(eps), D=str(D), endpoint=str(endpoint)))
        C = max(C, 3 * eps, endpoint)
    need((rho, A, C) == (F(30, 77), F(75, 77), F(4274, 1001)), 'fresh derived density coefficients')
    maxima = []
    for r in rows:
        j, U, S = (r['j'], r['U'], r['S'])
        eps, D = (F(r['epsilon']), F(r['D']))
        need(3 * eps <= C and D <= C * j, 'all-block slope and intercept')
        need(F(5, 2) * Q * phi - A * Q * Q == 0, 'quadratic cancellation')
        need(2 * j * phi + 3 * Q * U - 2 * A * Q * j == 3 * Q * eps, 'linear coefficient identity')
        for n in (0, 1, 2, 100):
            T = Q * n + j
            Ut = phi * n + U
            St = Q * phi * n * n // 2 + Q * n * U + S
            need(2 * T * Ut + St - A * T * T == 3 * Q * eps * n + D, 'exact block polynomial')
        if F(r['endpoint']) == C:
            maxima.append(r)
    need(len(maxima) == 1 and maxima[0]['j'] == 117, 'unique endpoint maximizer')
    signed = [(prod((PRIMES[i] for i in range(3) if mask >> i & 1)), (-1) ** mask.bit_count()) for mask in range(8)]
    U = S = 0
    for T in range(1541):
        if T and gcd(T, Q) == 1:
            U += 1
            S += T
        iu = sum((sign * (T // d) for d, sign in signed))
        isu = sum((sign * d * (T // d) * (T // d + 1) // 2 for d, sign in signed))
        need((U, S) == (iu, isu), 'independent Mobius count and sum')
        need(2 * T * U + S <= A * T * T + C * T, 'integer counting regression')
    K, J = (F(46339, 250000), F(51, 10))
    need((2 * A, 4 * A * K, C * J) == (F(150, 77), F(139017, 192500), F(108987, 5005)), 'insertion coefficients')
    need(4 * A / PMIN <= F(4, 3) and F(12, PMIN) <= F(4, 3) and (4 < PMIN), 'all prime-domain and integral conditions')
    double_obstructions = []
    for q, z in [(462, 275), (546, 65)]:
        candidates = []
        for i in range(1, 5):
            for j in range(-4, 3):
                if abs(F(j) + F(i, 4)) <= 1:
                    common = gcd(j + i * z, q)
                    need(common > 1, 'radius4 cannot simply be reused for the second prime')
                    candidates.append(dict(i=i, j=j, gcd=common))
        need(len(candidates) == 9, 'all feasible pairs at tau1/4 through radius4')
        need(gcd(5 * z - 2, q) == 1 and abs(F(-2) + F(5, 4)) <= 1, 'radius5 repairs this single obstruction only')
        double_obstructions.append(dict(mark=q, z=z, tau='1/4', radius4_candidates=candidates, one_radius5_pair=dict(i=5, j=-2), scope='ansatz on arbitrary residue/shear; no claim of realizable minimum or global geometry optimality'))
    canonical = json.dumps(rows, sort_keys=True, separators=(',', ':')).encode()
    return dict(status='PRELIMINARY_NOT_A_FAMILY_CLOSURE', finite_algebra_checks='PASS', checks=CHECKS, core=Q, core_primes=list(PRIMES), new_prime_min=PMIN, density=dict(phi=phi, rho=str(rho), A=str(A), C=str(C), all_integer_T_proved_by_block_coefficients=True, residues=Q, maximum_row=maxima[0], maximum_slope=str(max((3 * F(r['epsilon']) for r in rows))), period_rows_sha256=sha256(canonical).hexdigest(), C_optimal_only_for_fixed_A_polynomial_count=True, actual_CDF_optimality_claimed=False), repair=dict(radius=4, positive_half_residues=Q, closed_half_intervals=2 * Q, all_tau_in_minus_half_to_half=True, grouped_table=[dict(i=i, j=j, z_values=zs) for (i, j), zs in sorted(grouped.items())], radius3_counterexamples=obstructions, radius_minimality_scope='only 1<=i<=R and |j+i*tau|<=1; no universal geometric optimality claim'), insertion=dict(Phi_definition='ln(1+x)+x*ln(1+1/x)', periodic_Phi_numerator=str(2 * A), periodic_linear_numerator=str(4 * A * K), periodic_sqrtM_coefficient=str(C * J), free_Phi_numerator='6', free_linear_numerator=str(12 * K), inherited_radius4_integral_caps=dict(K=str(K), J4=str(J)), p_divides_M_required=True, q_divides_M_required=False, p_squared_in_M_allowed=True), potential_profiles=dict(projective_labels=prod((p + 1 for p in PRIMES)), anchor_indices=[Q, 2 * Q, Q * Q], squared_anchor_branches=Q * prod((p + 1 for p in PRIMES)), squared_anchor_records_at2688cells=2688 * Q * prod((p + 1 for p in PRIMES)), squared_anchor_bytes=8 * 2688 * Q * prod((p + 1 for p in PRIMES)), cheap_mixed_indices=[dict(multiplier=m, index=Q * m, records=2688 * 288 * m, bytes=8 * 2688 * 288 * m, divides_all_B_to_N_for_N_ge2=True) for m in [2, 7, 11]], these_are_domain_sizes_not_completed_profiles=True), double_repair42pr_obstructions=double_obstructions, no_profiles_generated=True, no_anchor_caps_claimed=True, no_infinite_family_partition_certified=True, family154p_closed=False, radical2002_closed_here=False, global_Cstar_improved=False, double_insertion42pr_proved=False)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    z = certificate()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    print(json.dumps({k: z[k] for k in ['status', 'finite_algebra_checks', 'checks', 'core', 'density', 'insertion', 'potential_profiles']}, indent=2))
