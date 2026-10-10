"""Reconstruct every profile witness from the preserved algorithms.

Every stage has a 600-second timeout and keeps stdout/stderr. Producer and
consumer cover all records in all eighteen profiles. Independent algebra,
stratified geometric audits and complete small physical banks are separate.
"""
from pathlib import Path
import os,sys,subprocess,time,json,hashlib,argparse,uuid
ROOT=Path(__file__).resolve().parent
RUN=None
def run(name,args):
 start=time.monotonic(); log=RUN/'logs'/f'{name}.log'
 with log.open('w') as f:
  f.write('COMMAND '+repr([str(x) for x in args])+'\n');f.flush()
  try:
   z=subprocess.run([str(x) for x in args],stdout=f,stderr=subprocess.STDOUT,timeout=600,env={**os.environ,'CSTAR_THREADS':'4','PYTHONDONTWRITEBYTECODE':'1'})
   result=dict(stage=name,returncode=z.returncode,seconds=time.monotonic()-start,timeout_seconds=600,command=[str(x) for x in args],log=log.relative_to(RUN).as_posix())
  except subprocess.TimeoutExpired:
   result=dict(stage=name,returncode=None,status='TIMEOUT',seconds=time.monotonic()-start,timeout_seconds=600,command=[str(x) for x in args],log=log.relative_to(RUN).as_posix())
  stages.append(result);(RUN/'logs/reconstruction_stages.json').write_text(json.dumps(stages,indent=2)+'\n')
 if result['returncode']!=0:raise RuntimeError(result)
 print('PASS',name,result['seconds'],flush=True)
stages=[]
def main():
 global RUN
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--work-dir',type=Path);a=ap.parse_args()
 RUN=(a.work_dir or ROOT/'runs'/('geometry-'+uuid.uuid4().hex)).resolve()
 if RUN.exists():raise ValueError('work directory must be new')
 missing=[];parent=RUN
 while not parent.exists():missing.append(parent);parent=parent.parent
 for p in reversed(missing):print(p.resolve(),flush=True);p.mkdir()
 for name in ('bin','work','logs'):
  p=RUN/name;print(p.resolve(),flush=True);p.mkdir()
 code=ROOT/'profiles/code';data=ROOT/'profiles/data'
 run('compile_algebra',['g++','-O3','-std=c++17','-fopenmp',code/'algebra.cpp','-o',RUN/'bin/algebra'])
 for k,ag in ((1,4),(1,16),(2,16)):
  for mode in ('produce','consume'):
   run(f'compile_{mode}_{k}_{ag}',['g++','-O3','-std=c++17','-fopenmp',f'-DDEGREE_EXP={k}',f'-DANGULAR={ag}',*(['-DBUILD_CONSUMER=1'] if mode=='consume' else []),code/'angle.cpp','-o',RUN/f'bin/{mode}_{k}_{ag}'])
 plan=json.loads((ROOT/'profiles/profiles_plan.json').read_text())
 for row in plan:
  p,k,ag,stem=row['p'],row['k'],row['ag'],row['stem']
  pre=RUN/'work'/stem
  for mode in ('produce','consume'):run(f'{stem}_{mode}',[RUN/f'bin/{mode}_{k}_{ag}',mode,p,data,pre])
  for ext in ('.cells','.branches','.audit'):
   if Path(str(pre)+ext).read_bytes()!=Path(str(pre)+'.checked'+ext).read_bytes():raise ArithmeticError('producer consumer mismatch '+stem+ext)
  # Verify exact byte identity of rebuilt cell caps against authentic cells.
  if Path(str(pre)+'.cells').read_bytes()!=(ROOT/'profiles/results'/f'{stem}.cells').read_bytes():raise ArithmeticError('authentic cell mismatch '+stem)
  for ext in ('.branches','.audit'):
   if Path(str(pre)+ext).read_bytes()!=(ROOT/'profiles/results'/(stem+ext)).read_bytes():raise ArithmeticError('compact output mismatch '+stem+ext)
  expected=json.loads((ROOT/'profiles/results'/(stem+'.json')).read_text());w=Path(str(pre)+'.witness')
  with w.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
  if expected['witness_bytes']!=w.stat().st_size or expected['witness_sha256']!=h:raise ArithmeticError('full witness mismatch '+stem)
  tensor=RUN/'work'/(stem+'.tensor')
  run(f'{stem}_tensor',[RUN/'bin/algebra',p,k,data,tensor])
  if tensor.read_bytes()!=(ROOT/'profiles/results'/(stem+'.tensor')).read_bytes():raise ArithmeticError('tensor output mismatch '+stem)
  audit=RUN/'work'/(stem+'.independent.json')
  run(f'{stem}_audit',[sys.executable,'-B',*(['-O'] if not __debug__ else []),code/'audit_one.py',p,ag,'--k',k,'--prefix',ROOT/'profiles/results'/stem,'--output',audit])
  if audit.read_bytes()!=(ROOT/'profiles/results'/(stem+'.independent.json')).read_bytes():raise ArithmeticError('independent audit output mismatch '+stem)
 print('RECONSTRUCTION COMPLETE',flush=True)
if __name__=='__main__':main()
