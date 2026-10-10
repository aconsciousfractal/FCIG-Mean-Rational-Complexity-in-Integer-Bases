"""Projective label counts are derived from the declared squarefree core."""
from math import prod
from exact import need,prime_trial
CORE_PRIMES=(2,5,13)
CORE=prod(CORE_PRIMES)
CORE_LABELS=prod(p+1 for p in CORE_PRIMES)
MIN_NEW_PRIME=17
def factors(p):
 need(p==1 or (p>=MIN_NEW_PRIME and prime_trial(p) and p not in CORE_PRIMES),'new prime or sentinel1')
 return CORE_PRIMES if p==1 else CORE_PRIMES+(p,)
def mark(p):
 ps=factors(p)
 return prod(ps),prod(q+1 for q in ps)

CORE_MULTIPLIERS=(1,2,5,10,13,26,65,130)
