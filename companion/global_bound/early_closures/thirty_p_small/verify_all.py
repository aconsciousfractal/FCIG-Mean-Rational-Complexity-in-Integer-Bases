#!/usr/bin/env python3
"""Regenerate the mark-30 and selected mark-30p certificates.

Requires Python >=3.10 and GNU g++ C++17/OpenMP/__int128; no network or solver.
Default mode never overwrites expected results. --generate is explicit development mode,
not a comparison replay. Witness files are rebuilt in a temporary directory.
"""
from pathlib import Path
import argparse,subprocess,tempfile,sys,os,json,time,hashlib,shutil,csv
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
import arithmetic as ar, reference_family as ref, independent as ind, check_routes as cr
import certificate as cert,new_math as nm,negative
from exact import need

def data_bytes(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def run(cmd,env):
 print('+',' '.join(map(str,cmd)),flush=True);subprocess.run(list(map(str,cmd)),check=True,env=env)
def put(name,data,generate):
 p=ROOT/'results'/name
 if generate:p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
 else:need(p.is_file() and p.read_bytes()==data,'regeneration mismatch: '+name)
 print(('GENERATED ' if generate else 'PASS regenerated ')+name,flush=True)
def putj(name,data,g):put(name,data_bytes(data),g)
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()

def main(a):
 need(sys.version_info>=(3,10),'Python >=3.10');need(sys.byteorder=='little','little-endian initial witness dump')
 need(1<=a.threads<=4,'declared thread range');need(shutil.which('g++') is not None,'GNU C++ compiler required')
 t0=time.monotonic();cells,cellmeta=ar.make_cells(ROOT)
 need(ar.encode_cells(cells)==(ROOT/'data/cells.txt').read_text(),'whole-cell bounds regenerated')
 _,logmeta=ref.log_checks(ROOT/'data/log_table.txt')
 env={**os.environ,'CSTAR_THREADS':str(a.threads),'PYTHONDONTWRITEBYTECODE':'1'}
 flags=['g++','-O3','-std=c++17','-fopenmp']
 if a.sanitize:flags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
 w=Path(tempfile.mkdtemp(prefix='cstar_30p_'));print('WORK',w,flush=True)
 try:
  profiles=w/'profiles';profiles.mkdir();checks=[];common={};common_hash={};digests={};banks=[];refdir=w/'reference';refdir.mkdir()
  run(flags+[ROOT/'code/reference_min.cpp','-o',w/'reference_min'],env)
  for p,four,T,mk in cert.configurations():
   q=(30 if four else 6)*p;prefix=profiles/cert.stem(p,four)
   run(flags+[f'-DFOUR_PRIMES={int(four)}',f'-DPRIME={p}',f'-DDIR_T={T}','-DAUDIT_SPARSE=1',ROOT/'code/uniform.cpp','-o',w/'uniform'],env)
   run([w/'uniform',ROOT/'data/log_table.txt',ROOT/'data/partition.txt',ROOT/'data/shape_centres.tsv',prefix,mk,ROOT/'data/cells.txt'],env)
   exts=['.audit','.cellcosts','.points','.bank']+[f'.{s}{k}' for k in range(1,mk+1) for s in ['cells','radial']]
   for ext in exts:put('profiles/'+cert.stem(p,four)+ext,Path(str(prefix)+ext).read_bytes(),a.generate)
   ind.T=T;ind.J=4*T;cr.T=T;cr.J=4*T
   checks.append(ind.run_one(ROOT,p,prefix,four))
   if T not in common:
    common[T]=ind.common_costs(ROOT,prefix,cells);common_hash[T]=sha(str(prefix)+'.cellcosts')
   else:need(sha(str(prefix)+'.cellcosts')==common_hash[T],'same geometry costs across marks with equal direction count')
   digests[str(q)]=sha(str(prefix)+'.witnesses')
   for k in range(1,min(mk,2)+1):run([w/'reference_min',q**k,q,refdir/f'M{q**k}_q{q}.min'],env)
   if q==30:run([w/'reference_min',27000,30,refdir/'M27000_q30.min'],env)
   banks.append(cert.physical_bank(q,prefix,refdir,mk));Path(str(prefix)+'.witnesses').unlink()
   ar.logs.cache_clear()
  pr=cert.collect(profiles);rows,routes,cuts=cert.inventory(pr)
  finite=w/'finite';finite.mkdir();domain=finite/'domain.txt';domain.write_text(''.join(f"{z['M']} {z['q']}\n" for z in rows))
  run(flags+[ROOT/'code/census.cpp','-o',w/'census'],env)
  run([w/'census','produce',domain,finite],env);run([w/'census','consume',domain,finite],env)
  fc=cert.finite_certificate(finite,rows);cl=cert.closure(pr,fc)
  payloads={
   'uniform_profiles.json':pr,'finite_certificate.json':fc,'closure.json':cl,
   'repair_table.json':nm.repair_table(),'simple_repair.json':nm.simple_table(),'periodic_certificate.json':nm.periodic(),
   'prime_repair_banks.json':nm.physical_checks(),'candidate.json':nm.candidate(),
   'large_prime_gates.json':nm.large_prime_gates(),
   'independent_checks.json':dict(status='PASS',cell_inputs=cellmeta,log_table=logmeta,common_costs=common,per_mark=checks,physical_banks=banks),
   'witness_digests.json':digests,
   'negative_checks.json':negative.run(w/'census',finite,rows,w/'negative'),
  }
  for name,data in payloads.items():putj(name,data,a.generate)
  import io
  st=io.StringIO(newline='');cw=csv.writer(st);cw.writerow(['q','k','directions','mean_upper','upper_decimal','direct_radical_gate'])
  for z in pr:cw.writerow([z['q'],z['k'],z['directions'],z['total_upper'],z['upper_decimal'],z['radical_gate']])
  put('summary.csv',st.getvalue().encode(),a.generate)
  report=dict(status='PASS',generated_expected_files=a.generate,python_optimized=not __debug__,
              undefined_behavior_sanitizer=a.sanitize,closed_primes=list(cert.PS),finite_moduli=fc['moduli'],
              finite_states=fc['states'],uniform_profiles=len(pr),all_selected_sectors_closed=True,
              entire_30p_family_closed=False,global_Cstar_identified=False,
              historical_closure_theorems_used=False,threads=a.threads,seconds=time.monotonic()-t0)
  print(json.dumps(report,indent=2),flush=True)
  if a.report:Path(a.report).write_bytes(data_bytes(report))
 finally:
  if not a.keep_work:shutil.rmtree(w)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sanitize',action='store_true');p.add_argument('--threads',type=int,default=4)
 p.add_argument('--generate',action='store_true',help='development only: write expected outputs; not a comparison replay')
 p.add_argument('--keep-work',action='store_true');p.add_argument('--report');a=p.parse_args()
 try:main(a)
 except (Exception,KeyboardInterrupt) as e:print('FAIL '+type(e).__name__+': '+str(e),file=sys.stderr,flush=True);raise SystemExit(1)
