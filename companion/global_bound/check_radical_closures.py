#!/usr/bin/env python3
"""Pay eleven five-prime radical exclusions, including all valuations and N."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import argparse,csv,json,sys
from check import log_interval
ROOT=Path(__file__).resolve().parent
def need(condition,label):
 if not condition:raise ValueError(label)
def cap_profile(path,q,K,alpha=True):
 value=json.loads(Path(path).read_text())
 need(value['q']==q and value.get('transfer_index',value.get('K'))==K,'full marked profile scope')
 if alpha:need(value.get('alpha',value.get('cdf_alpha'))=='3/2','universal primitive CDF3/2')
 return F(value['total_upper'])
def normalized(A,B):return max(F(0),A)/log_interval(B)[0]
def finite(root,q,threshold):
 domain=json.loads((root/'results/finite_inventory.json').read_text())
 domains=domain.get('domains',[domain]);domain=next(d for d in domains if d['radical']==q)
 need(domain['base_upper_exclusive']==threshold,'finite all-valuation threshold')
 ps=domain['prime_support'] if 'prime_support' in domain else [p for p in (2,3,5,7,11,13,17,19) if q%p==0]
 need(prod(ps)==q,'exact five-prime radical')
 expected=[]
 for B in range(q,threshold,q):
  r=B
  for p in ps:
   while r%p==0:r//=p
  if r==1:expected.append(B)
 need([r['B'] for r in domain['rows']]==expected,'independent factorization of every N1 valuation')
 with (root/'results/finite_bounds.csv').open(newline='') as stream:rows=[{k:int(v) for k,v in z.items()} for z in csv.DictReader(stream)]
 rows=[r for r in rows if r['B'] in expected]
 need([r['B'] for r in rows]==expected,'each and only finite physical base')
 bounds=[]
 for r in rows:
  B=r['B'];M=r['M'];need(M==B and r['states']==M and r['representatives']==M//2+1,'all N1 physical states')
  n=2*r['sum_log_upper_ticks']-M*r['log_M_lower_ticks'];need(n==r['mean_numerator'],'finite centered mean')
  need(0<F(r['log_B_lower_ticks'],10**9)<=log_interval(B)[0],'outward physical normalization')
  ratio=F(max(n,0),2*M*r['log_B_lower_ticks']);need(ratio<F(2,85),'strict finite radical cap');bounds.append(ratio)
 return max(bounds,default=F(0)),len(rows)
def refined2310(root):
 refined=root/'profiles_refined';sys.path.insert(0,str(refined));sys.path.insert(0,str(refined/'code'))
 from assemble_refined import assemble,refined_grid
 value=assemble();den,old,grid,mapping=refined_grid()
 need(value==json.loads((refined/'results/profile.json').read_text()),'full retained hybrid assembly')
 v=[F(x,10**9) for x in value['monotone_profile_ticks']]
 radii=[F(old[b*21][2],den) for b in range(90)]+[F(grid[b*21][2],den) for b in range(148)]
 need(len(v)==len(radii)==238 and all(v[i]>=v[i+1] for i in range(len(v)-1)),'signed isotone radial envelope')
 alpha=F(3,2);rho=F(old[0][1],den)
 bulk=v[-1]+sum((v[i]-v[i+1])*min(F(1),alpha*radii[i]) for i in range(len(v)-1))
 cusp_at_rho=log_interval(2309)[1]+log_interval(13860)[1]/2-log_interval(rho)[0]/2+F(5,2)*rho
 cusp=alpha*rho*(max(cusp_at_rho-v[0],0)+F(1,2))
 return bulk+cusp
def verify(root=ROOT):
 premises=root/'premises';rows=[]
 A=refined2310(premises/'fixed2310')
 mixed=[(2310,13860,100000,premises/'fixed2310',A,premises/'fixed2310/profiles/results/q2310.json'),
        (2730,16380,20000,premises/'fixed2730/physical',cap_profile(premises/'fixed2730/profiles/results/q2730_K16380_hybrid_AG8.alpha3over2.json',2730,16380),premises/'fixed2730/profiles/results/q2730_K2730_coarse.alpha3over2.json'),
        (3570,7140,10000,premises/'fixed3570_3990/physical',cap_profile(premises/'fixed3570_3990/profiles/results/q3570_K7140_angular_tail_AG8.alpha3over2.json',3570,7140),premises/'fixed2730/profiles/results/q3570_K3570_coarse.alpha3over2.json'),
        (3990,7980,6000,premises/'fixed3570_3990/physical',cap_profile(premises/'fixed3570_3990/profiles/results/q3990_K7980_coarse.alpha3over2.json',3990,7980),premises/'fixed2730/profiles/results/q3990_K3990_coarse.alpha3over2.json')]
 for q,K,Bmin,physical,A,first in mixed:
  need(q*q%K==0,'all N>=2 valuations admit mixed transfer index')
  upper2=normalized(A,q);upper_large=normalized(cap_profile(first,q,q,False),Bmin);upper_finite,count=finite(physical,q,Bmin)
  upper=max(upper2,upper_large,upper_finite);need(upper<F(2,85),'entire mixed-index radical closure '+str(q)+' '+str([float(x) for x in (upper2,upper_large,upper_finite)]))
  rows.append({'q':q,'K_for_N_ge_2':K,'N1_base_threshold':Bmin,'finite_N1_caps':count,'upper_N_ge_2':str(upper2),'upper_large_base':str(upper_large),'upper_finite_N1':str(upper_finite),'upper':str(upper),'all_positive_valuations':True,'all_N_ge_1':True})
 for q,alias in [(4290,'fixed4290'),(4830,'fixed4830'),(5610,'fixed5610_family102'),(6006,'fixed6006_family154'),(6090,'fixed6090_family165'),(6270,'fixed6270_family130'),(6510,'fixed6510_family114')]:
  paths=list((premises/alias/'profiles/results').glob(f'q{q}_K{q}.json'));need(len(paths)==1,'unique direct full mark')
  upper=normalized(cap_profile(paths[0],q,q),q);need(upper<F(2,85),'direct full-radical gate')
  rows.append({'q':q,'K':q,'upper':str(upper),'all_positive_valuations':True,'all_N_ge_1':True})
 need([r['q'] for r in rows]==[2310,2730,3570,3990,4290,4830,5610,6006,6090,6270,6510],'all eleven exclusions exactly once')
 return {'schema':'eleven-radical-closures-v1','status':'PASS','rows':rows,'retained_geometric_caps_reintegrated_separately':True,'primary_geometry_reconstructed_here':False,'physical_primary_records_reconstructed_here':False}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();need(not args.output.exists(),'fresh output required')
 value=verify();args.output.write_text(json.dumps(value,indent=2)+'\n',encoding='utf8',newline='\n');print('PASS eleven radical closures')
if __name__=='__main__':main()
