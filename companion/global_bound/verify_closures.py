#!/usr/bin/env python3
"""Reintegrate retained finite profiles and validate numerical closure premises.

This mode does not reconstruct primary geometric witnesses. See
reconstruct_closures.py for the complete producer/consumer route.
"""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time,shutil,signal
ROOT=Path(__file__).resolve().parent
def need(condition,label):
 if not condition:raise ValueError(label)
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def guarded(command,log,env,limit):
 process=subprocess.Popen(list(map(str,command)),stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=(os.name=='posix'))
 try:return process.wait(timeout=limit)
 except subprocess.TimeoutExpired:
  if os.name=='posix':
   try:os.killpg(process.pid,signal.SIGKILL)
   except ProcessLookupError:pass
  else:process.kill()
  process.wait();print('STAGE TIMEOUT; process group terminated',file=log,flush=True);return 124
def profile_components(root=ROOT):
 paths=[]
 for path in (root/'premises').rglob('verify_profiles.py'):
  relative=path.relative_to(root/'premises')
  if any(part in ('imports','baseline','originals') for part in relative.parts):continue
  paths.append(path)
 paths.extend(path for path in (root/'premises').rglob('verify_anchor.py') if 'imports' not in path.relative_to(root/'premises').parts)
 paths.extend(path for path in (root/'premises').rglob('verify_refined.py') if 'imports' not in path.relative_to(root/'premises').parts)
 return sorted(paths)
def expected_components(root=ROOT,full=False):
 result=[p.parent.relative_to(root).as_posix() for p in profile_components(root)]
 result.extend(p.parent.parent.relative_to(root).as_posix()+'/all_valuations' for p in sorted((root/'premises').rglob('inventory[0-9]*.py')) if p.parent.name=='code' and 'imports' not in p.relative_to(root/'premises').parts)
 result.extend(p.relative_to(root).as_posix() for p in sorted((root/'premises').rglob('verify_insertion[0-9]*.py')) if 'imports' not in p.relative_to(root/'premises').parts)
 for alias,scripts in [('family105',['verify_theory105.py','verify_gates105.py']),('family66',['verify_theory66.py','verify_gates66.py']),('family42',['verify_repair42.py','verify_repair462.py'])]:result.extend(alias+'/'+script for script in scripts)
 result.extend(['eleven_radical_closures','early_prime_count_and30p_closures','marked70_all_valuations'])
 if full:
  physical=[p for p in (root/'premises').rglob('run_finite[0-9]*.py') if 'imports' not in p.relative_to(root/'premises').parts]
  physical.extend([root/'premises/fixed2310/code/finite2310.py',root/'premises/fixed2730/physical/code/finite2730.py',root/'premises/fixed3570_3990/physical/code/finite.py'])
  for path in sorted(physical):
   result.append(path.relative_to(root).as_posix())
   if (path.parent/'independent_physical.py').exists():result.append((path.parent/'independent_physical.py').relative_to(root).as_posix())
 return result
def verify(output,root=ROOT,full=False,threads=4):
 output=Path(output).resolve();need(not output.exists(),'fresh output directory required')
 print('VERIFIED ABSOLUTE OUTPUT',output,flush=True);output.mkdir(parents=True)
 temporary=output/'temporary';print('VERIFIED ABSOLUTE TEMPORARY WORK',temporary,flush=True);temporary.mkdir()
 rows=[];env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','TMPDIR':str(temporary),'TMP':str(temporary),'TEMP':str(temporary)}
 if not __debug__:env['PYTHONOPTIMIZE']='1'
 catalog_path=root/'closure_input_catalog.json';catalog=json.loads(catalog_path.read_text())
 need(len({r['path'] for r in catalog['files']})==len(catalog['files']),'unique catalog inventory')
 for row in catalog['files']:
  path=root/row['path'];need(path.resolve().is_relative_to(root.resolve()),'contained catalog path')
  need(path.stat().st_size==row['bytes'] and digest(path)==row['sha256'],'reader input binding '+row['path'])
 for index,path in enumerate(profile_components(root)):
  alias=path.parent.relative_to(root).as_posix();name=f'profile_{index:02d}'
  command=[sys.executable,'-B',str(path),'--rebuild' if full else '--quick','--report',str(output/(name+'.json'))]
  if path.name=='verify_anchor.py' and not full:command.remove('--quick')
  if full:command.extend(['--threads',str(threads),'--sanitize','--keep-work'])
  if path.name=='verify_anchor.py' and full:
   command.remove('--sanitize')
  if path.name=='verify_refined.py':command.extend(['--baseline-root',str(path.parent.parent/'profiles')])
  print('CHECK',alias,flush=True);start=time.monotonic()
  with (output/(name+'.log')).open('w',encoding='utf8',newline='\n') as log:
   returncode=guarded(command,log,env,3600 if full else 600)
  row={'component':alias,'returncode':returncode,'seconds':round(time.monotonic()-start,3),'log_sha256':digest(output/(name+'.log'))}
  rows.append(row);(output/'progress.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
  if returncode:print((output/(name+'.log')).read_text()[-12000:],flush=True)
 need(all(row['returncode']==0 for row in rows),'all retained profile numerical checks')
 def run_step(alias,command):
  name=f'premise_{len(rows):02d}';start=time.monotonic();print('CHECK',alias,flush=True)
  with (output/(name+'.log')).open('w',encoding='utf8',newline='\n') as log:
   returncode=guarded(command,log,env,3600 if full else 600)
  row={'component':alias,'returncode':returncode,'seconds':round(time.monotonic()-start,3),'log_sha256':digest(output/(name+'.log'))};rows.append(row)
  (output/'progress.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
  if returncode:print((output/(name+'.log')).read_text()[-12000:],flush=True)
  need(returncode==0,'numerical premise '+alias)
 # The all-valuation inventories independently factor every eligible multiple.
 inventories=[p for p in (root/'premises').rglob('inventory[0-9]*.py') if p.parent.name=='code' and 'imports' not in p.relative_to(root/'premises').parts]
 for index,path in enumerate(sorted(inventories)):
  component=path.parent.parent
  run_step(component.relative_to(root).as_posix()+'/all_valuations',[sys.executable,'-B',root/'check_closure_tables.py','--component',component,'--output',output/f'family_{index:02d}.json'])
 # These finite predicates verify residue repair, density and integral constants.
 for index,path in enumerate(sorted((root/'premises').rglob('verify_insertion[0-9]*.py'))):
  if 'imports' in path.relative_to(root/'premises').parts:continue
  run_step(path.relative_to(root).as_posix(),[sys.executable,'-B',path,'--output',output/f'insertion_{index:02d}.json'])
 for alias,scripts in [('family105',['verify_theory105.py','verify_gates105.py']),('family66',['verify_theory66.py','verify_gates66.py']),('family42',['verify_repair42.py','verify_repair462.py'])]:
  copied=output/(alias+'_checks');print('VERIFIED ABSOLUTE CHECK WORK',copied.resolve(),flush=True);shutil.copytree(root/'premises'/alias/'checks',copied)
  sources=[{'local_copy':p.relative_to(copied).as_posix(),'sha256':digest(p),'bytes':p.stat().st_size} for p in (copied/'imports').rglob('*') if p.is_file()]
  (copied/'PROVENANCE.json').write_text(json.dumps({'schema':'reader-input-bindings-v1','sources':sources},indent=2)+'\n',encoding='utf8',newline='\n')
  for script in scripts:run_step(alias+'/'+script,[sys.executable,'-B',copied/script])
 run_step('eleven_radical_closures',[sys.executable,'-B',root/'check_radical_closures.py','--output',output/'radical_closures.json'])
 early=root/'early_closures/verify_early.py';marked=root/'marked70/verify.py'
 need(early.exists() and marked.exists(),'complete early and70-family companion inputs')
 script="import importlib.util,json,sys;from pathlib import Path;p=Path(sys.argv[1]);s=importlib.util.spec_from_file_location('early_closure_reader',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);z=m.verify(root=p.parent);Path(sys.argv[2]).write_text(json.dumps(z,indent=2)+'\\n',encoding='utf8',newline='\\n')"
 run_step('early_prime_count_and30p_closures',[sys.executable,'-B','-c',script,early,output/'early_closures.json'])
 run_step('marked70_all_valuations',[sys.executable,'-B',marked,'--rebuild' if full else '--quick','--output-dir',output/'marked70_work','--report',output/'marked70.json'])
 if full:
  # Primary physical producer and independent consumer use the exact residual domains.
  physical=[p for p in (root/'premises').rglob('run_finite[0-9]*.py') if 'imports' not in p.relative_to(root/'premises').parts]
  physical.extend([root/'premises/fixed2310/code/finite2310.py',root/'premises/fixed2730/physical/code/finite2730.py',root/'premises/fixed3570_3990/physical/code/finite.py'])
  for index,path in enumerate(sorted(physical)):
   work=output/f'physical_{index:02d}'
   run_step(path.relative_to(root).as_posix(),[sys.executable,'-B',path,'--work-dir',work,'--threads',threads,'--compare'])
   independent=path.parent/'independent_physical.py'
   if independent.exists():run_step(independent.relative_to(root).as_posix(),[sys.executable,'-B',independent,'--root',path.parent.parent,'--raw-dir',work,'--output',output/f'independent_physical_{index:02d}.json'])
 need([r['component'] for r in rows]==expected_components(root,full),'complete ordered numerical component inventory')
 report={'status':'PASS','mode':'full_primary_geometric_reconstruction' if full else 'retained_caps_exact_integration_and_independent_HNF','primary_witnesses_reconstructed':full,'primary_witness_scope':'later family/radical components and marked70','early_primary_witnesses_reconstructed':False,'early_retained_arithmetic_verified':True,'geometric_sampling_is_not_a_full_check':True,'python_optimized':not __debug__,'input_catalog_sha256':digest(catalog_path),'catalog_files':len(catalog['files']),'all_consumed_scalar_closures_reconstructed':True,'components':rows}
 (output/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
 return report
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
 verify(args.output)
if __name__=='__main__':main()
