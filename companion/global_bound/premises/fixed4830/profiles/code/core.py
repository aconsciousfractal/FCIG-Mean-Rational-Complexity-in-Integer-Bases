"""Explicit five-prime mark; no numeric profile is inherited."""
from math import prod
from exact import need
CORE_PRIMES=(2,3,5,7,23)
CORE=prod(CORE_PRIMES)
CORE_LABELS=prod(p+1 for p in CORE_PRIMES)
def factors(p):
 need(p==1,'fixed4830 mark selector')
 return CORE_PRIMES
def mark(p):
 ps=factors(p)
 return prod(ps),prod(q+1 for q in ps)
