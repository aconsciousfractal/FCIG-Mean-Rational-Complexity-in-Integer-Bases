#!/usr/bin/env python3
"""Exact finite-to-general assembly for every prime 11<=p<=71 and marked70 caps."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
from exact import need,sieve,log_lo,out_decimal
from assemble import profile
ROOT=Path(__file__).resolve().parents[1]
def assemble(root=ROOT):
 plan=json.loads((root/'profiles_plan.json').read_text())
 primes=[p for p in sieve(71) if p>=11]
 need(sorted(z['p'] for z in plan if z['p']>1)==primes,'all sixteen small primes, exactly once')
 need(sorted((z['p'],z['k']) for z in plan if z['p']==1)==[(1,1),(1,2)],'both marked70 divisibilities')
 rows=[];routes=[];marks=[]
 for z in plan:
  p,k,ag=z['p'],z['k'],z['ag'];need(z['stem']==(f'q{70*p}' if p>1 else f'q70_k{k}'),'unambiguous result stems')
  got=profile(p,root/'results'/(z['stem']+'.cells'),root,k,ag)
  stored=json.loads((root/'results'/(z['stem']+'.json')).read_text())
  need(set(stored)==set(got)|{'witness_bytes','witness_sha256'},'profile schema')
  need(all(stored[x]==v for x,v in got.items()),'exact reintegration of every stored field')
  need(stored['witness_bytes']==8*got['witnesses'] and len(stored['witness_sha256'])==64,'witness metadata only')
  cap=F(got['total_upper']);rows.append(got)
  if p==1:
   simple=F(z['simple_cap']);need(simple==({1:F(7,50),2:F(0)}[k]),'declared marked70 input')
   need(cap<simple,'strict marked70 simple cap')
   marks.append({'q':70,'k':k,'divisor':70**k,'exact_upper':str(cap),'upper_decimal':out_decimal(cap),'simple_cap':str(simple),'valid_for':'nu_70(M) for every positive M divisible by divisor'})
  else:
   need(k==1 and ag==4,'full mark one-step configuration');q=70*p
   norm=max(F(0),cap)/log_lo(q);need(norm<F(2,85),'strict all-exponents and all-lengths route')
   routes.append({'p':p,'q':q,'domain':'B=2^a 5^b 7^c p^d; a,b,c,d,N>=1','required_divisor_of_M':q,'base_min':q,'mean_upper':str(cap),'normalized_upper':str(norm),'normalized_upper_decimal':out_decimal(norm)})
 gamma=max(F(z['normalized_upper']) for z in routes)
 return {'status':'PASS','closed_primes':primes,'closed_radicals':[70*p for p in primes],'all_positive_prime_exponents':True,'all_lengths_N_ge_1':True,'marked70_caps':marks,'routes':routes,'common_upper':str(gamma),'common_upper_decimal':out_decimal(gamma),'comparison_threshold':'2/85','strict_margin':str(F(2,85)-gamma),'direct_profiles':len(rows),'geometric_witnesses_in_direct_profiles':sum(z['witnesses'] for z in rows),'corner_evaluations_in_direct_profiles':sum(4*z['angular_nodes']*z['witnesses'] for z in rows),'historical_profile_caps_assumed':False,'entire_70p_family_closed_by_this_component_alone':False,'remaining_tail_starts_at_prime':73,'global_Cstar_identified':False,'second_geometric_audit_exhaustive':False}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args();z=assemble();s=json.dumps(z,indent=2)+'\n';out=ROOT/'results/small_prime_closure.json'
 if a.check:need(out.read_text()==s,'exact closure report')
 else:out.write_text(s)
 print('PASS small-prime closure',z['common_upper_decimal'],flush=True)
if __name__=='__main__':main()
