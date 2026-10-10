"""Deliberate invalid inputs and bounded failures; no failure is a counterexample to C*."""
from pathlib import Path
from fractions import Fraction as F
import subprocess,struct,shutil,os
from exact import need,phi_hi,log_hi


def run(root,tmp,finite,consumer,valid_geometric):
    root=Path(root);tmp=Path(tmp);tmp.mkdir(exist_ok=True,parents=True)
    data=root/'data';records=[];env={**os.environ,'CSTAR_THREADS':'1'}
    def reject(name,cmd,expected):
        p=subprocess.run(list(map(str,cmd)),stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=env)
        text=p.stdout+p.stderr
        need(p.returncode!=0 and expected in text,'negative input not rejected: '+name)
        print('EXPECTED FAILURE',name,repr(text.strip()),flush=True)
        records.append(dict(name=name,returncode=p.returncode,required_output=expected))
    # A finite angle grid without the curvature allowance underestimates this fixed vector.
    values=[F(max((4+3*n)**2,(-n+12)**2),16+n*n) for n in range(-4,5)]
    discrete=max(values);need(discrete<10,'angular-grid counterexample')
    records.append(dict(name='angular_grid_without_error',grid_squared_max=str(discrete),continuous_squared_max='10',witness=[1,3],missing_angle_tangent='-1/3',status='bounded shortcut rejected'))
    # The old scalar prime-erasure estimate does not establish the proposed raw .16 gate.
    upper=F(7,25)+phi_hi(F(8,111))+F(7304,399*33)
    need(upper>F(4,25) and upper>F(2,85)*log_hi(1110),'bounded erasure failure')
    records.append(dict(name='excision_bound_at_p37_M1110',upper=str(upper),status='does not establish 4/25 or the normalized target; not a counterexample'))
    # Geometric consumer: valid record then three distinct corruptions.
    good=Path(valid_geometric).read_bytes();pre=tmp/'geo'
    for name,blob,err in [('geometry_truncated',good[:-1],'truncated witness file'),('geometry_extra',good+b'X','extra witness data'),('geometry_nonunit',bytes(8)+good[8:],'full mark unit')]:
        Path(str(pre)+'.witness').write_bytes(blob)
        reject(name,[consumer,'consume',1,data,pre],err)
    Path(str(pre)+'.witness').unlink()
    reject('geometry_composite_prime',[consumer,'consume',9,data,pre],'prime input')
    # An independent finite stream is small enough to inspect completely in every mutation.
    pre=tmp/'finite';cmd=[finite,'produce',6330,211,data,pre,0]
    p=subprocess.run(list(map(str,cmd)),capture_output=True,text=True)
    need(p.returncode==0,'valid finite mutation source')
    good=Path(str(pre)+'.bin').read_bytes()
    p=subprocess.run(list(map(str,[finite,'consume',6330,211,data,pre,0])),capture_output=True,text=True)
    need(p.returncode==0 and 'PASS consumed' in p.stdout,'valid finite source consumption')
    records.append(dict(name='valid_finite_mutation_control',status='PASS',states=6330))
    bad_unit=struct.pack('<ii',0,2)+good[8:]
    bad_cong=good[:8]+struct.pack('<ii',0,1)+good[16:]
    for name,blob,err in [('finite_truncated',good[:-1],'truncated witnesses'),('finite_extra',good+b'X','extra witnesses'),('finite_nonunit',bad_unit,'full denominator unit'),('finite_congruence',bad_cong,'state congruence')]:
        Path(str(pre)+'.bin').write_bytes(blob)
        reject(name,[finite,'consume',6330,211,data,pre,0],err)
    Path(str(pre)+'.bin').unlink()
    reject('finite_composite_prime',[finite,'consume',6390,213,data,pre,0],'prime input')
    reject('finite_fifth_prime',[finite,'consume',44310,211,data,pre,0],'exact radical 30p')
    dd=tmp/'duplicate.txt';dd.write_text('6330 211\n6330 211\n')
    reject('finite_duplicate_modulus',[finite,'run',dd,data,tmp,tmp/'bad.csv',1],'duplicate modulus')
    return dict(status='PASS',tests=records,negative_test_count=len(records),global_counterexample_claimed=False)
