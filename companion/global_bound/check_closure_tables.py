#!/usr/bin/env python3
"""Rebuild a family's all-valuation inventory and every rational scalar gate.

Retained physical caps are checked arithmetically. Their primary shortest-vector
records are independently regenerated and consumed only by the full route.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,importlib,json,sys
def need(condition,label):
 if not condition:raise ValueError(label)
def verify(component):
 root=Path(component).resolve();sources=list((root/'code').glob('inventory[0-9]*.py'))
 need(len(sources)==1,'unique family inventory source')
 sys.path.insert(0,str(root/'code'));module=importlib.import_module(sources[0].stem)
 rows,meta=module.make();expected=json.loads((root/'results/inventory_certificate.json').read_text())
 need(meta==expected,'complete regenerated inventory and infinite gates')
 need(len({r['M'] for r in rows})==len(rows),'unique finite moduli')
 with (root/'results/finite_inventory.csv').open(newline='') as stream:stored=list(csv.DictReader(stream))
 need(stored==[{k:str(v) for k,v in r.items()} for r in rows],'every valuation, route and scalar bound')
 witnesses=[r for r in rows if r['route']=='witness']
 need((root/'results/finite_domain.txt').read_text()==''.join(f"{r['M']} {r['p']}\n" for r in witnesses),'complete physical residual domain')
 core=getattr(module,'CORE',42);anchors=getattr(module,'ANCHORS',None)
 if anchors is None:anchors={core:module.A1,core**2:module.A2}
 profiles=[]
 for path in (root/'profiles').rglob('*.json'):
  if 'imports' in path.relative_to(root/'profiles').parts:continue
  try:value=json.loads(path.read_text())
  except (ValueError,UnicodeError):continue
  if isinstance(value,dict) and {'q','transfer_index','total_upper'}<=set(value):profiles.append((path,value))
 for K,cap in anchors.items():
  found=[value for path,value in profiles if value['q']==core and value['transfer_index']==K]
  need(found and min(F(v['total_upper']) for v in found)<cap,'signed geometric anchor obligation '+str(K))
 small=getattr(module,'SMALL',None)
 if small is None:small=[p for p in module.sieve(module.PMIN-1) if p>=({42:11,66:13,105:11}[core])]
 for p in small:
  if core==42 and p==11:
   sys.path.insert(0,str(root/'profiles/code'))
   from assemble_small import assemble
   result=assemble(root/'profiles')
   need(result==json.loads((root/'profiles/results/small_prime_closure.json').read_text()),'all three exceptional q462 regions')
   continue
  found=[value for path,value in profiles if value['q']==core*p and value['transfer_index']==core*p]
  need(found and max(F(0),min(F(v['total_upper']) for v in found))/module.log_lo(core*p)<F(2,85),'full small-prime mark '+str(p))
 with (root/'results/finite_bounds.csv').open(newline='') as stream:bounds=[{k:int(v) for k,v in r.items()} for r in csv.DictReader(stream)]
 need([(r['M'],r['p']) for r in bounds]==[(r['M'],r['p']) for r in witnesses],'all and only physical residual caps')
 for row in bounds:
  M=row['M'];need(row['states']==M and row['representatives']==M//2+1,'complete physical reflection mass')
  numerator=2*row['sum_log_upper_ticks']-M*row['log_M_lower_ticks']
  need(numerator==row['mean_numerator'],'exact physical logarithmic mean')
  need(row['log_M_lower_ticks']>0 and F(row['log_M_lower_ticks'],10**9)<=module.log_lo(M),'directed physical denominator')
  need(85*max(0,numerator)<4*M*row['log_M_lower_ticks'],'strict finite physical cap')
 for row in rows:
  if row['route']=='analytic':need(F(row['normalized_upper'])<F(2,85),'finite analytic cap')
 for gate in meta['infinite_regions']:need(F(gate['normalized_upper'])<F(2,85),'all-N infinite scalar region')
 return {'status':'PASS','core':core,'all_valuations_enumerated':True,'finite_rows':len(rows),'physical_caps':len(bounds),'analytic_rows':len(rows)-len(bounds),'infinite_regions':meta['infinite_regions'],'geometric_anchor_indices':list(anchors),'small_primes':small,'physical_primary_records_reconstructed_here':False}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--component',required=True,type=Path);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
 need(not args.output.exists(),'fresh output required');result=verify(args.component)
 args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print('PASS',result['core'],result['finite_rows'])
if __name__=='__main__':main()
