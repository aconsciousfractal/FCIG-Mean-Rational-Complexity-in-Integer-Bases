"""Exact assembly of the normalized comparison, not an exact value of C*."""
from pathlib import Path
from fractions import Fraction as F
import csv,json
from exact import need,log_lo,out_decimal,sieve
from inventory import DELTA,PMIN,PMAX,BMAX

def finite_certificate(path,rows):
 expected=[r for r in rows if r['route']=='witness'];got=list(csv.DictReader(Path(path).open()))
 need([int(r['M']) for r in got]==[r['M'] for r in expected],'full witness-domain match')
 values=[]
 for r,e in zip(got,expected):
  n={k:int(v) for k,v in r.items()};M=n['M'];need(n['p']==e['p'] and n['states']==M and n['representatives']==M//2+1,'finite state weights and support')
  need(n['mean_numerator']==2*n['sum_log_upper_ticks']-M*n['log_M_lower_ticks'],'finite exact mean assembly')
  need(n['mean_numerator']<0,'all new finite upper means are negative')
  need(n['stencil'] in (0,1) and 1<=n['max_height']<=M//2,'finite witness status')
  u=F(n['mean_numerator'],2*10**9*M);values.append((u,M,n['p']))
 worst=max(values)
 return dict(status='PASS',moduli=len(got),states=sum(int(r['M']) for r in got),representatives=sum(int(r['representatives']) for r in got),all_upper_means_negative=True,maximum_mean_upper=str(worst[0]),maximum_mean_upper_decimal=out_decimal(worst[0]),worst_majorant_modulus=worst[1],worst_majorant_prime=worst[2],refined_stencils=sum(int(r['stencil']) for r in got),dense_logarithms_verified=4000000,logarithm_secondary_seeds=1024,logarithm_secondary_recurrence_steps=3998976,witness_bytes_consumed=8*sum(int(r['representatives']) for r in got))

def assemble(profiles,rows,inventory,finite):
 ps=[p for p in sieve(210) if p>=37]
 need(len(ps)==35,'35 primes in the residual block')
 need([(z['p'],z['k']) for z in profiles]==[(p,1) for p in ps]+[(1,1),(1,2)],'complete profile domain')
 for z in profiles:
  u=F(z['total_upper']);cap=F(11,100) if z['p']>1 else (F(7,25) if z['k']==1 else F(61,500))
  need(u<cap,'strict profile cap')
 routes=[dict(kind='complete_mark_profile',p=z['p'],q=z['q'],ratio=z['normalized_upper']) for z in profiles if z['p']>1]
 routes += [dict(kind=z['name'],ratio=z['normalized_upper']) for z in inventory['infinite_regions']]
 for r in rows:
  if r['route']=='analytic':routes.append(dict(kind='finite_analytic_N1',M=r['M'],p=r['p'],ratio=r['normalized_upper']))
 routes.append(dict(kind='finite_witness_N1',ratio='0'))
 # The old seven-prime theorem is an explicitly scoped internal import.
 routes.append(dict(kind='S7_internal_import',primes=[11,13,17,19,23,29,31],ratio='233/10000'))
 gamma=max(F(r['ratio']) for r in routes)
 need(gamma<DELTA,'uniform strict separation of the supremum')
 need(finite['moduli']==inventory['witness_moduli'] and finite['states']==inventory['witness_states'],'finite/infinite partition join')
 low=[p for p in sieve(210) if p>=11]
 need(low==[11,13,17,19,23,29,31]+ps,'all primes below 211 covered, disjointly')
 out=dict(status='PASS',new_full_mark_primes=ps,new_full_mark_radicals=[30*p for p in ps],new_full_mark_supports=35,all_p_ge_37_closed_all_valuations_and_lengths=True,all_30p_p_ge_11_closed_with_S7_import=True,old_S7_import_replayed=False,radical_210_replayed=False,other_four_prime_supports_closed=False,global_Cstar_identified=False,mu_4_over_25_uniform_all_p_ge_37_proved=False,common_upper=str(gamma),common_upper_decimal=out_decimal(gamma),comparison_lower_candidate='2/85',margin_to_2_over_85=str(DELTA-gamma),worst_routes=[r for r in routes if F(r['ratio'])==gamma],geometric_profiles=len(profiles),geometric_cells=len(profiles)*2688,geometric_witnesses=sum(z['witnesses'] for z in profiles),geometric_corner_evaluations=sum(z['witnesses']*z['angular_nodes']*4 for z in profiles),finite_partition=inventory,finite_certificate=finite,profiles=profiles)
 return out
