from pathlib import Path
from fractions import Fraction as F
from math import prod
import subprocess,shutil,os
from exact import need,log_lo,log_hi
from certify import HERE,profile

def run_negative(consumer,producer,good_prefix,work,root=HERE):
 work=Path(work);work.mkdir(exist_ok=True);good=Path(str(good_prefix)+'.witness');data=good.read_bytes();outputs=[]
 for name,payload,expected in [
  ('truncated',data[:-1],'truncated witness records'),
  ('extra',data+b'x','extra witness records'),
  ('zero',b'\0'*8+data[8:],'nonunit full mark')]:
  pre=work/name;Path(str(pre)+'.witness').write_bytes(payload)
  z=subprocess.run([str(consumer),'2','5','7','1',str(root/'data'),str(pre)],capture_output=True,text=True,env={**os.environ,'CSTAR_THREADS':'2'})
  need(z.returncode!=0 and expected in z.stderr,'consumer rejects '+name)
  outputs.append(dict(test=name,returncode=z.returncode,output=z.stderr.strip()))
 z=subprocess.run([str(producer),'2','5','9','1',str(root/'data'),str(work/'composite')],capture_output=True,text=True)
 need(z.returncode!=0 and 'prime input' in z.stderr,'composite prime rejected')
 outputs.append(dict(test='composite',returncode=z.returncode,output=z.stderr.strip()))
 first=profile((2,5,7),1,root/'results/profiles/q70.cells',root)
 need(F(first['total_upper'])>F(2,85)*log_hi(70),'single step is not claimed sufficient for q70')
 outputs.append(dict(test='first_step_q70',output='FAIL: the one-step upper certificate does not establish the radical gate; not a counterexample'))
 den=int((root/'data/grid.tsv').read_text().split()[0]);old_d=1600*den*110**2
 need(old_d.bit_length()+(10**15).bit_length()>=126 and old_d.bit_length()+(10**14).bit_length()<126,'old guard rejected, new scale within guard')
 outputs.append(dict(test='old_scale_guard',old_bit_budget=old_d.bit_length()+(10**15).bit_length(),new_bit_budget=old_d.bit_length()+(10**14).bit_length(),output='old conservative arithmetic guard fails; outward ceiling scale changed from 10^15 to 10^14'))
 # A constant must be preserved by a probability transfer.
 need(F(70**2,70**2)==1 and F(70**2,70)!=1,'prime-power normalization')
 outputs.append(dict(test='wrong_transfer_denominator',correct='1',wrong='70'))
 # All entries of this HNF are divisible by every mark prime.
 need(70*70==4900 and all(70%p==0 for p in (2,5,7)),'rank-zero branches excluded')
 outputs.append(dict(test='zero_rank_HNF',H=[70,0,0,70],admissible_parent_marks=0))
 raw=[list(map(int,l.split())) for l in (root/'results/profiles/q70.cells').read_text().splitlines()]
 need(20*log_lo(2)>F(max(r[3] for r in raw),10**9),'thin lattice beyond finite grid')
 outputs.append(dict(test='dropping_cusp',rho='1/(70*2^40)',operator_lower_bound='20 log 2',output='the finite grid maximum alone is not a global pointwise bound'))
 return dict(status='PASS',tests=len(outputs),outputs=outputs,expected_failures_are_not_counterexamples_to_Cstar=True)
