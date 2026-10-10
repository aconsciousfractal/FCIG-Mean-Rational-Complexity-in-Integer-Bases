#!/usr/bin/env python3
"""Exact PE70 gates, conditional on the separately certified profile interface.

This reads the two frozen profile summaries.  It checks their arithmetic
use and hashes; it does not rerun the geometric producer or consumer.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

from verify_repair70 import fs, log_bounds, need, phi_bounds, PE70

ROOT = Path(__file__).resolve().parent
DELTA = F(2, 85)


def decimal_ceiling(x, places=15):
    x = F(x)
    scale = 10**places
    value = -((-x.numerator*scale)//x.denominator)
    sign = "-" if value < 0 else ""
    value = abs(value)
    return f"{sign}{value//scale}.{value%scale:0{places}d}"


def decimal_floor(x, places=15):
    x = F(x)
    scale = 10**places
    value = x.numerator*scale//x.denominator
    sign = "-" if value < 0 else ""
    value = abs(value)
    return f"{sign}{value//scale}.{value%scale:0{places}d}"


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def import_anchor(k, upper):
    path = ROOT/"imports"/f"q70_k{k}.json"
    raw = path.read_bytes()
    source = json.loads(raw)
    need(source["q"] == 70 and source["k"] == k, "profile mark/index")
    need(source["transfer_index"] == 70**k, "profile divisibility domain")
    anchor = F(source["total_upper"])
    need(anchor < upper, "strict simplification of certified profile cap")
    return {"path": path.relative_to(ROOT).as_posix(),
            "sha256": sha256(raw).hexdigest(), "bytes": len(raw),
            "divisibility_required": 70**k,
            "total_upper": fs(anchor), "simple_cap": fs(upper),
            "geometry_replayed_here": False}


def main():
    anchors = [import_anchor(1, F(7, 50)), import_anchor(2, F(0))]
    records = []
    gates = [
        ("p_ge_1000_all_N", 1000, 70000, 70000, F(7, 50), 1),
        ("p_ge_73_N_ge_2", 73, 5110, 5110**2, F(0), 2),
        ("p_ge_73_B_ge_3000000", 73, 3000000, 3000000, F(7, 50), 1),
    ]
    for name, p, B, M, anchor, k in gates:
        root = isqrt(M)
        need(root*root <= M < (root+1)*(root+1), "integer root enclosure")
        argument = F(24, 7*p)
        phi_lo, phi_hi = phi_bounds(argument)
        logB_lo, logB_hi = log_bounds(B)
        need(logB_lo > 0 and phi_lo > 0, "positive denominator and numerator")
        numerator_hi = anchor+phi_hi+PE70/root
        need(numerator_hi > 0, "positive route numerator")
        upper = numerator_hi/logB_lo
        margin = DELTA-upper
        need(margin > 0, "strict gate "+name)
        records.append({
            "route": name, "p_scalar": p, "B_min": B, "M_min": M,
            "profile_k": k, "simple_anchor": fs(anchor),
            "sqrt_lower": root, "Phi_argument": fs(argument),
            "Phi_lower": fs(phi_lo), "Phi_upper": fs(phi_hi),
            "logB_lower": fs(logB_lo), "logB_upper": fs(logB_hi),
            "numerator_upper": fs(numerator_hi), "Delta_upper": fs(upper),
            "Delta_upper_decimal_ceiling": decimal_ceiling(upper),
            "margin_from_delta": fs(margin),
            "margin_decimal_floor": decimal_floor(margin),
        })

    p = 71
    B = 70*p
    phi_lo, _ = phi_bounds(F(24, 7*p))
    _, lnB_hi = log_bounds(B)
    negative_lower = (phi_lo+PE70/B)/lnB_hi
    need(negative_lower > DELTA, "limited failure of the zero-anchor rational bound at p71,N2")

    low_gap = [p for p in range(71, 74) if prime(p)]
    high_gap = [p for p in range(997, 1010) if prime(p)]
    need(low_gap == [71, 73] and high_gap == [997, 1009], "adjacent prime gates")
    finite_primes = [p for p in range(73, 1000) if prime(p)]
    small_primes = [p for p in range(11, 73) if prime(p)]
    need(len(finite_primes) == 148 and len(small_primes) == 16, "finite prime partition")
    max_route = max(records, key=lambda x: F(x["Delta_upper"]["numerator"],
                                           x["Delta_upper"]["denominator"]))
    output = {
        "status": "PASS", "q": 70, "delta": fs(DELTA),
        "normalization": "Delta(B,N)=mu(B^N)/ln(B)",
        "domain": "B=2^a*5^b*7^c*p^d; a,b,c,d,N>=1; p>=11 prime",
        "profile_imports": anchors, "gates": records,
        "largest_gate_upper_route": max_route["route"],
        "largest_gate_upper": max_route["Delta_upper"],
        "largest_gate_upper_is_true_Delta_value": False,
        "finite_complement": {"N": 1, "p_min": 73, "p_upper_exclusive": 1000,
                              "B_upper_exclusive": 3000000,
                              "prime_count": len(finite_primes), "primes": finite_primes},
        "separate_small_prime_list": small_primes,
        "adjacent_primes": {"lower": low_gap, "upper": high_gap},
        "limited_negative": {
            "p": 71, "N": 2, "B": B, "M": B*B, "simple_anchor": fs(0),
            "rational_sufficient_bound_lower": fs(negative_lower),
            "exceeds_delta": True,
            "interpretation": "failure of this chosen sufficient bound; not a counterexample for Delta",
        },
        "claim_scope": "exact gate arithmetic with separately certified marked-profile inputs",
        "finite_complement_closed_here": False,
    }
    path = ROOT/"gates70_certificate.json"
    path.write_text(json.dumps(output, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "PASS", "certificate": path.name,
                      "gate_uppers": {r["route"]: r["Delta_upper_decimal_ceiling"]
                                      for r in records},
                      "margin_floors": {r["route"]: r["margin_decimal_floor"]
                                        for r in records}}, indent=2))


if __name__ == "__main__":
    main()
