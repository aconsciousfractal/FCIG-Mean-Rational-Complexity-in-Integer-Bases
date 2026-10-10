"""Marked-profile checker. Exact stored reintegration; full replay is explicit."""
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
from checks import input_checks,algebra_check,geometric_check
from assemble import grid
from assemble_small import assemble
from exact import need
def verify():
 start=time.monotonic();inp=input_checks(ROOT);z=assemble(ROOT)
 need(z==json.loads((ROOT/'results/small_prime_closure.json').read_text()),'small closure exact equality')
 den,gr=grid(ROOT);audits=[]
 for row in json.loads((ROOT/'profiles_plan.json').read_text()):
  p,k,ag=row['p'],row['k'],row['ag'];pre=ROOT/'results'/row['stem']
  br,algebra=algebra_check(p,k,pre)
  geo=geometric_check(p,k,ag,pre,gr,den,br)
  audits.append(dict(stem=row['stem'],algebra=algebra,stratified_geometry=geo))
 return dict(status='PASS',mode='exact_stored_integrals_and_independent_algebra',direct_profiles_reintegrated=z['direct_profiles'],common_upper=z['common_upper'],input_checks=inp,audits=audits,geometric_witnesses_replayed_here=False,geometric_second_audit_exhaustive=False,seconds=time.monotonic()-start)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');ap.add_argument('--report',type=Path);a=ap.parse_args();z=verify()
 if a.report:a.report.write_text(json.dumps(z,indent=2)+'\n')
 print(json.dumps({k:v for k,v in z.items() if k!='audits'},indent=2))
if __name__=='__main__':main()
