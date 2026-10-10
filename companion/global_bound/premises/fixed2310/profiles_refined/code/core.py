"""Projective label counts are derived from the declared squarefree core."""
from math import prod
from exact import need,prime_trial
CORE_PRIMES=(2,3,5,7)
CORE=prod(CORE_PRIMES)
CORE_LABELS=prod(p+1 for p in CORE_PRIMES)
MIN_NEW_PRIME=11
def factors(p):
 need(p==1 or (p>=MIN_NEW_PRIME and prime_trial(p) and p not in CORE_PRIMES),'new prime or sentinel1')
 return CORE_PRIMES if p==1 else CORE_PRIMES+(p,)
def mark(p):
 ps=factors(p)
 return prod(ps),prod(q+1 for q in ps)

def transfer_index(p,k=1,multiplier=None):
 ps=factors(p);q=prod(ps)
 if multiplier is None:
  need(k in (1,2),'supported legacy pure exponent')
  return q**k
 need(isinstance(multiplier,int) and multiplier>=1,'positive integer transfer multiplier')
 rest=multiplier
 for pr in ps:
  while rest%pr==0:rest//=pr
 need(rest==1,'transfer index has precisely the full marked prime support')
 return q*multiplier
