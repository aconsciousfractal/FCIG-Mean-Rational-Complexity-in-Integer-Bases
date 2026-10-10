#!/usr/bin/env python3
"""Exact scalar upper bounds for insertion into the mark 102.

The analytic statement concerns prime p >= 19, p | M.  The scalar
functions also accept integral p >= 19 that need not be prime, and M
that need not be a multiple of p: this is intentional for monotone
boundary evaluations.  validate_regime() checks actual theorem inputs.

All logarithms use a rational enclosure with an explicit positive tail.
The method is transparently adapted from the frozen P72_66P exact
support; this is not claimed as an independent arithmetic library.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd, isqrt

MARK = 102
RADIUS = 4
RHO = F(16, 51)
PERIODIC_A = F(40, 51)
PERIODIC_C = F(14744, 4029)
K_UPPER = F(46339, 250000)
J4_UPPER = F(51, 10)
PERIODIC_PHI_NUMERATOR = F(80, 51)
PERIODIC_LINEAR_NUMERATOR = F(92678, 159375)
PERIODIC_ROOT_COEFFICIENT = F(7372, 395)
FREE_PHI_NUMERATOR = F(6)
FREE_LINEAR_NUMERATOR = F(139017, 62500)


def _require(test, message):
    if not test:
        raise ValueError(message)


def _integer(n, name, minimum=1):
    _require(type(n) is int and n >= minimum,
             f"{name} must be an integer >= {minimum}")


def is_prime(n):
    return type(n) is int and n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def validate_regime(p, M, require_final_mark=False):
    """Check insertion hypotheses, optionally rad(M) = 102p.

    The insertion lemma itself does not require 102 | M.
    The condition rad(M) = 102p suffices for equality of denominator constraints and means;
    this validator does not characterize incidental equality of two numerical means.
    """
    _integer(p, "p", 19)
    _integer(M, "M")
    _require(is_prime(p) and gcd(p, MARK) == 1, "new prime outside mark 102")
    _require(M % p == 0, "new prime must divide M")
    if require_final_mark:
        _require(M % (MARK * p) == 0, "all four target primes must divide M")
        residual = M
        for prime in (2, 3, 17, p):
            while residual % prime == 0:
                residual //= prime
        _require(residual == 1, "extra prime prevents denominator-set identification with the physical mean")
    return True


@lru_cache(maxsize=None)
def log_bounds(x, terms=80):
    """Rational lo <= ln(x) <= hi, with binary reduction and tail."""
    x = F(x)
    _require(x > 0, "logarithm argument must be positive")
    _integer(terms, "terms")
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1

    def series(z):
        power, zz, total = z, z * z, F(0)
        for j in range(terms):
            total += 2 * power / (2 * j + 1)
            power *= zz
        remainder = 2 * power / ((2 * terms + 1) * (1 - zz))
        return total, total + remainder

    lo, hi = series((x - 1) / (x + 1))
    lo2, hi2 = series(F(1, 3))
    return (lo + k * lo2, hi + k * hi2) if k >= 0 else (lo + k * hi2, hi + k * lo2)


@lru_cache(maxsize=None)
def phi_upper(x):
    x = F(x)
    _require(x > 0, "Phi argument must be positive")
    return log_bounds(1 + x)[1] + x * log_bounds(1 + 1 / x)[1]


def sqrt_lower(M):
    _integer(M, "M")
    return isqrt(M)


def periodic_upper(p, M):
    """Upper for Phi(80/(51p)) + 92678/(159375p) + 7372/(395sqrt(M))."""
    _integer(p, "p", 19)
    _integer(M, "M")
    _require(4 * PERIODIC_A / p <= F(4, 3), "periodic integration range")
    return (phi_upper(PERIODIC_PHI_NUMERATOR / p)
            + PERIODIC_LINEAR_NUMERATOR / p
            + PERIODIC_ROOT_COEFFICIENT / sqrt_lower(M))


def free_upper(p):
    """Upper for Phi(6/p) + 139017/(62500p), independent of M."""
    _integer(p, "p", 19)
    _require(F(12, p) <= F(4, 3), "free integration range")
    return phi_upper(FREE_PHI_NUMERATOR / p) + FREE_LINEAR_NUMERATOR / p


def insertion_upper(p, M):
    return min(periodic_upper(p, M), free_upper(p))


def normalized_upper(p, B, M, anchor):
    """Uniform scalar upper: max(0, anchor + insertion cost) / ln(B).

    This performs scalar arithmetic only; the caller must prove the
    marked anchor on its domain and, for a physical conclusion, rad(M)
    = 102p and M = B**N.  The anchor may have either sign.  Clipping the
    complete numerator at zero preserves monotonicity in p, M, and B;
    clipping the anchor itself would unnecessarily weaken the bound.
    """
    _integer(B, "B", 2)
    _integer(M, "M", B)
    anchor = F(anchor)
    return max(F(0), anchor + insertion_upper(p, M)) / log_bounds(B)[0]


def fraction_json(x):
    x = F(x)
    return {"numerator": x.numerator, "denominator": x.denominator}


def decimal_ceiling(x, places=15):
    x = F(x)
    scale = 10 ** places
    n = -((-x.numerator * scale) // x.denominator)
    sign = "-" if n < 0 else ""
    n = abs(n)
    return f"{sign}{n // scale}.{n % scale:0{places}d}"
