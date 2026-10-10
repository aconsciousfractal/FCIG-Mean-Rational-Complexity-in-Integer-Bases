#!/usr/bin/env python3
"""Rebuild and consume every new certificate. Expected files are never overwritten."""
from pathlib import Path
import argparse,subprocess,tempfile,sys,os,json,time,hashlib,shutil,csv
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
from certify import *
from checks import run_checks
from negative import run_negative

def call(cmd,env=None):
 print('+',' '.join(map(str,cmd)),flush=True)
 subprocess.run(list(map(str,cmd)),check=True,env=env)
def same(a,b):need(Path(a).read_bytes()==Path(b).read_bytes(),'regenerated file differs: '+str(b))
def compare_json(path,data,work):
 p=Path(work)/Path(path).name;write_json(p,data);same(p,path)
def domain_join():
 from itertools import combinations
 ps=[n for n in range(2,48) if all(n%d for d in range(2,isqrt(n)+1))];need(len(ps)==15,'all primes up to 47');need(not any(all(n%d for d in range(2,isqrt(n)+1)) for n in range(48,53)),'no prime between 47 and 53')
 full={prod(t) for t in combinations(ps,3) if t[:2]!=(2,3)}
 old=set(map(int,(ROOT/'imports/closed_353_radicals.txt').read_text().split()))
 new={prod(t) for t in triples(ROOT)}
 need(len(old)==353 and len(new)==89 and not(old&new) and full==old|new,'complete disjoint 442-support partition')
 return dict(verified_set_identity=True,previous_supports=353,new_supports=89,total=442,
             G53_and_F6p_replayed=False,previous_analytic_proofs_reaudited=False)
def verify(sanitize=False,threads=2,keep=False):
 t0=time.monotonic();need(1<=threads<=8,'thread range');grid(ROOT);domain_join()
 work=Path(tempfile.mkdtemp(prefix='cstar_three_primes_'));print('WORK',work,flush=True)
 env={**os.environ,'CSTAR_THREADS':str(threads)}
 try:
  flags=['-O3','-std=c++17','-fopenmp']
  if sanitize:flags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
  for name in ['produce','consume','audit']:
   call(['g++',*flags,ROOT/'code'/f'{name}.cpp','-o',work/name])
  cf=['-O3','-std=c++17']+(['-fsanitize=undefined','-fno-sanitize-recover=all'] if sanitize else [])
  call(['g++',*cf,ROOT/'code/census.cpp','-o',work/'census'])
  profiles=work/'profiles';profiles.mkdir();digests={};allcases=cases(ROOT)
  for t,k in allcases:
   pre=profiles/stem(t,k);args=[*map(str,t),str(k),ROOT/'data',pre]
   call([work/'produce',*args],env);call([work/'consume',*args],env)
   same(str(pre)+'.checked.cells',str(pre)+'.cells');same(str(pre)+'.cells',ROOT/'results/profiles'/(stem(t,k)+'.cells'))
   call([work/'audit',*args],env);same(str(pre)+'.audit',ROOT/'results/profiles'/(stem(t,k)+'.audit'))
   z=Path(str(pre)+'.witness');h=hashlib.sha256()
   with open(z,'rb') as f:
    while b:=f.read(8*1024*1024):h.update(b)
   digests[stem(t,k)]=h.hexdigest()
   if t==(2,5,7) and k==1:
    neg=run_negative(work/'consume',work/'produce',pre,work/'negative',ROOT)
    compare_json(ROOT/'results/negative_checks.json',neg,work)
   z.unlink()
  compare_json(ROOT/'results/witness_digests.json',digests,work)
  finite=work/'finite';finite.mkdir();dom=finite/'domain.txt'
  rows=finite_domain(ROOT);dom.write_text(''.join(f"{r['M']} {r['q']}\n" for r in rows))
  call([work/'census','produce',dom,finite]);call([work/'census','consume',dom,finite])
  fc=finite_certificate(finite,ROOT);compare_json(ROOT/'results/finite.json',fc,work)
  cand=candidate();compare_json(ROOT/'results/candidate.json',cand,work)
  closure=assemble(profiles,fc,ROOT);compare_json(ROOT/'results/closure.json',closure,work)
  with open(work/'summary.csv','w') as f:
   wr=csv.writer(f);wr.writerow(['r','s','p','q','k','mu_upper','mu_upper_decimal','direct_radical_gate'])
   for row in closure['profiles']:wr.writerow(row['primes']+[row['q'],row['k'],row['total_upper'],row['upper_decimal'],row['direct_radical_gate']])
  same(work/'summary.csv',ROOT/'results/summary.csv')
  chk=run_checks(profiles,ROOT);compare_json(ROOT/'results/independent_checks.json',chk,work)
  report=dict(status='PASS',new_supports_closed=89,finite_profiles=93,witness_records=closure['witness_records'],
   finite_integer_states=fc['states'],expected_files_overwritten=False,python_optimized=not __debug__,
   undefined_behavior_sanitizer=sanitize,threads=threads,imported_theorems_replayed=False,
   all_three_prime_supports_excluded_with_declared_imports=True,global_Cstar_identified=False,
   seconds=time.monotonic()-t0,work_directory=str(work),set_join=domain_join())
  print(json.dumps(report,indent=2),flush=True);return report
 finally:
  if not keep:shutil.rmtree(work)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--sanitize',action='store_true');p.add_argument('--threads',type=int,default=2);p.add_argument('--keep-work',action='store_true');p.add_argument('--report');a=p.parse_args()
 report=verify(a.sanitize,a.threads,a.keep_work)
 if a.report:write_json(a.report,report)
