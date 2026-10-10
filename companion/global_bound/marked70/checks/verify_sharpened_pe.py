#!/usr/bin/env python3
"""Rational constants for the coordinate-sensitive insertion improvement.

The pointwise and integration arguments are analytic, in SHARPENED_PE.md.
This certifies their rational constants and the same frozen gate corners.
It does not modify an inventory, a profile, or a previously accepted gate.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

from verify_repair70 import fs, need, log_bounds, phi_bounds, C70, PE70
from verify_gates70 import decimal_ceiling, decimal_floor

ROOT = Path(__file__).resolve().parent
SHARP70_C = F(971503, 44247)
SHARP770_C = F(1755677, 43120)


def atan_bounds(x, last_even_index=20):
    """Alternating atan series; an even final index gives an upper bound."""
    x = F(x)
    need(0 < x < 1 and last_even_index % 2 == 0, "atan series hypotheses")
    upper = sum(((-1)**j)*x**(2*j+1)/F(2*j+1)
                for j in range(last_even_index+1))
    lower = upper-x**(2*last_even_index+3)/F(2*last_even_index+3)
    return lower, upper


def sharp_pe70_upper(p, M):
    """Numerical sufficient-expression upper, also valid at scalar corners."""
    need(p > 0 and M >= 1, "positive sharp70 expression parameters")
    need(F(24, 7*p) <= F(4, 3), "integration range y<=4/3 for70")
    root = isqrt(M)
    need(root*root <= M < (root+1)*(root+1), "sharp70 integer root")
    _, phi_hi = phi_bounds(F(12, 7*p))
    return phi_hi+F(9, 14*p)+SHARP70_C/root


def sharp_pe770_upper(p, M):
    need(p > 0 and M >= 1, "positive sharp770 expression parameters")
    need(F(360, 77*p) <= F(4, 3), "integration range y<=4/3 for770")
    root = isqrt(M)
    need(root*root <= M < (root+1)*(root+1), "sharp770 integer root")
    _, phi_hi = phi_bounds(F(180, 77*p))
    return phi_hi+F(135, 154*p)+SHARP770_C/root


def main():
    al, ah = atan_bounds(F(3, 4))
    ll, lh = log_bounds(F(5, 2))
    Klo, Khi = al-lh/2, ah-ll/2
    need(0 < Klo <= Khi < F(3, 16), "exact upper for atan3/4 minus half ln5/2")
    l2, h2 = log_bounds(2)
    need(h2 < F(25, 36), "ln2 rational upper")
    need(l2 > F(2, 3), "ln2 rational lower for the alternative remainder check")
    need(F(99, 70)**2 > 2, "sqrt2 rational upper")
    need(F(26, 15)**2 > 3, "sqrt3 rational upper")
    need(F(49, 20)**2 > 6, "sqrt6 rational upper")
    need(F(4, 3)**2+1 == F(5, 3)**2, "pointwise transition square")
    need((1+F(4, 3)/2)**2 == 1+F(4, 3)**2, "continuous pointwise transition")
    need(F(4, 3)/F(5, 3) == F(4, 5) == 1/(1+F(1, 2)**2),
         "derivative comparison endpoint")
    need(1-F(1, 4)-F(1, 4)**2 == F(11, 16) > 0,
         "Cauchy norm polynomial on small t")
    need(F(2, 3) > 0, "minimum of 3t²-2t+1")

    # Exact checks over deliberately mixed scale points; not a universal proof.
    point_records = []
    for t in (F(0), F(1, 10), F(1, 4), F(1, 2), F(3, 4), F(1)):
        for v in (F(1, 100), F(1, 4), F(1), F(4, 3), F(3, 2), F(2), F(10), F(100)):
            f1 = 1+v*t/(1+t*t)
            f2 = t+v/(1+t*t)
            envelope_squared = (1+v/2)**2 if v <= F(4, 3) else 1+v*v
            need(max(f1*f1, f2*f2) <= envelope_squared, "pointwise mixed-scale regression")
            need(envelope_squared < (1+v)**2, "strict improvement over triangular envelope")
            point_records.append({"t": fs(t), "v": fs(v),
                                  "envelope_squared": fs(envelope_squared)})

    coeff4 = F(22, 7)*F(99, 70)+F(26, 15)*F(25, 36)
    coeff6 = F(22, 7)*F(26, 15)+F(3, 2)*F(99, 70)*F(25, 36)
    need(coeff4 == F(74731, 13230), "r4 universal sqrt coefficient")
    need(coeff6 == F(11627, 1680), "r6 universal sqrt coefficient")
    need(C70*coeff4 == SHARP70_C, "sharp70 exact sqrt coefficient")
    need(F(453, 77)*coeff6 == SHARP770_C, "sharp770 exact sqrt coefficient")
    need(F(24, 7)*F(3, 16) == F(9, 14), "sharp70 linear correction")
    need(F(360, 77)*F(3, 16) == F(135, 154), "sharp770 linear correction")
    need(8*coeff4 == F(298924, 6615), "uniform2rs sqrt coefficient")
    need(F(24, 35) < F(4, 3) and F(360, 539) < F(4, 3),
         "range condition follows from insertion hypotheses70and770")

    gates = []
    frozen = [
        ("70_p_ge_1000_all_N", 70, 1000, 70000, 70000, F(7, 50)),
        ("70_p_ge_73_N_ge_2", 70, 73, 5110, 5110**2, F(0)),
        ("70_p_ge_73_B_ge_3000000", 70, 73, 3000000, 3000000, F(7, 50)),
        ("770_ell_ge_300_all_N", 770, 300, 231000, 231000, F(9, 80)),
    ]
    for label, q, p, B, M, anchor in frozen:
        improved = sharp_pe70_upper(p, M) if q == 70 else sharp_pe770_upper(p, M)
        loglo, loghi = log_bounds(B)
        upper = (anchor+improved)/loglo
        need(upper < F(2, 85), "same frozen corner passes sharpened expression")
        oldarg = F(24, 7*p) if q == 70 else F(360, 77*p)
        oldC = PE70 if q == 70 else F(453, 10)
        oldphi_lo, _ = phi_bounds(oldarg)
        old_expression_lower = (anchor+oldphi_lo+oldC/isqrt(M))/loghi
        need(upper < old_expression_lower, "new bound improves old bound at the same corner")
        margin = F(2, 85)-upper
        gates.append({"route": label, "q": q, "p_scalar": p,
                      "B_min": B, "M_min": M, "anchor": fs(anchor),
                      "new_upper": fs(upper),
                      "new_upper_decimal_ceiling": decimal_ceiling(upper),
                      "old_expression_lower": fs(old_expression_lower),
                      "margin_from_delta": fs(margin),
                      "margin_decimal_floor": decimal_floor(margin),
                      "new_upper_is_value_of_Delta": False})

    output = {
        "status": "PASS", "analytic_proof": "SHARPENED_PE.md",
        "transition_v": fs(F(4, 3)), "transition_X_over_r": fs(F(3, 4)),
        "linear_remainder_constant": {"exact_expression": "atan(3/4)-(1/2)*ln(5/2)",
                                      "lower": fs(Klo), "upper": fs(Khi),
                                      "rational_cap": fs(F(3, 16))},
        "sqrt_constant_r4": fs(coeff4), "sqrt_constant_r6": fs(coeff6),
        "sharp70": {"Phi_argument": "12/(7p)", "linear_remainder": "9/(14p)",
                    "sqrtM_coefficient": fs(SHARP70_C)},
        "sharp770": {"Phi_argument": "180/(77*ell)", "linear_remainder": "135/(154*ell)",
                     "sqrtM_coefficient": fs(SHARP770_C)},
        "uniform2rs": {"Phi_argument_upper": "5/(2p)", "linear_remainder_upper": "15/(16p)",
                       "sqrtM_coefficient": fs(F(298924, 6615))},
        "pointwise_regression_pairs": len(point_records), "pointwise_records": point_records,
        "frozen_gates": gates,
        "dependency_scope": "same repair tables and marked-profile contracts; no new geometric profile",
        "inventory_or_cutoff_changed": False,
        "claim_scope": "analytic lemma with exact finite rational constants; finite regressions are separate"
    }
    path = ROOT/"sharpened_pe_certificate.json"
    path.write_text(json.dumps(output, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "PASS", "certificate": path.name,
                      "new_gate_uppers": {g["route"]: g["new_upper_decimal_ceiling"] for g in gates},
                      "linear_constant_upper": decimal_ceiling(Khi)}, indent=2))


if __name__ == "__main__":
    main()
