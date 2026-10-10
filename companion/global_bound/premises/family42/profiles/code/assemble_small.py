#!/usr/bin/env python3
"""Exact assembly of every radical42p with11<=p<=103 and the marked42 inputs."""
from pathlib import Path
from fractions import Fraction as F
import json,argparse
from exact import need,log_lo,out_decimal,sieve
ROOT=Path(__file__).resolve().parents[1]
def assemble(root=ROOT):
 root=Path(root);plan=json.loads((root/'profiles_plan.json').read_text());rows=[];routes=[];marks=[]
 need(sorted(z['p'] for z in plan if z['p']>1)==[p for p in sieve(103) if p>=11],'exact declared prime domain')
 for z in plan:
  v=json.loads((root/z['folder']/'results'/(z['stem']+'.json')).read_text());need(v['p']==z['p'] and v['k']==z['k'],'profile identity');rows.append(v)
  cap=F(v['total_upper']);q=v['q']
  if z['p']==1:
   simple=F(z['simple_cap']);need(cap<simple,'marked42simplecap');marks.append({'q':42,'k':z['k'],'divisor':42**z['k'],'exact_upper':str(cap),'upper_decimal':out_decimal(cap),'simple_cap':str(simple),'valid_for':'nu_42(M) for every M divisible by divisor'});continue
  base=z.get('base_min',q);ratio=max(F(0),cap)/log_lo(base);need(ratio<F(2,85),'normalized direct route')
  routes.append({'p':z['p'],'q':q,'route':'N=1,B>=924' if z['p']==11 else 'all_N_all_exponents','base_min':base,'mean_upper':str(cap),'normalized_upper':str(ratio),'normalized_upper_decimal':out_decimal(ratio)})
 it=json.loads((root/'iterated/results/integrated.json').read_text());v=next(z for z in it if z['k']==2);cap=F(v['total_upper']);ratio=max(F(0),cap)/log_lo(462);need(ratio<F(2,85),'q462N>=2')
 routes.append({'p':11,'q':462,'route':'N>=2_all_exponents','base_min':462,'mean_upper':str(cap),'normalized_upper':str(ratio),'normalized_upper_decimal':out_decimal(ratio)})
 finite=json.loads((root/'iterated/results/finite_base462.json').read_text());need(finite['status']=='PASS' and finite['inventory']==[{'B':462,'exponents':[1,1,1,1]}] and finite['base_cutoff']==924,'finite residual completeness');ratio=F(finite['normalized_upper']);need(ratio<F(2,85),'finite residual gate')
 routes.append({'p':11,'q':462,'route':'N=1_B=462','base_min':462,'mean_upper':finite['mean_upper'],'normalized_upper':str(ratio),'normalized_upper_decimal':out_decimal(ratio)})
 gamma=max(F(z['normalized_upper']) for z in routes);need(gamma<F(2,85),'strict finite maximum of complete routes')
 return {'status':'PASS','closed_primes':[p for p in sieve(103) if p>=11],'closed_radicals':[42*p for p in sieve(103) if p>=11],'all_positive_prime_exponents':True,'all_lengths_N_ge_1':True,'historical_S113_profile_used':False,'historical_T3_closure_used':False,'marked42_caps':marks,'routes':routes,'common_upper':str(gamma),'common_upper_decimal':out_decimal(gamma),'comparison_threshold':'2/85','strict_margin':str(F(2,85)-gamma),'direct_profiles':len(rows),'selected_iterated_profiles':1,'iterated_auxiliary_horizon1_also_generated':True,'geometric_witnesses_in_direct_profiles':sum(z['witnesses'] for z in rows),'corner_evaluations_in_direct_profiles':sum(4*z['angular_nodes']*z['witnesses'] for z in rows),'initial_witnesses_iterated_profile':2688*1152,'entire_42p_family_closed_by_this_component_alone':False,'global_Cstar_identified':False,'second_geometric_audit_exhaustive':False}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'results/small_prime_closure.json');ap.add_argument('--check',action='store_true');a=ap.parse_args();data=(json.dumps(assemble(),indent=2)+'\n').encode()
 if a.check:need(a.output.read_bytes()==data,'stored exact small-prime closure')
 else:a.output.write_bytes(data)
 print(json.loads(data)['common_upper_decimal'],flush=True)
