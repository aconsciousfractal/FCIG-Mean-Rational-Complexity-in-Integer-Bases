"""Re-run deliberately invalid inputs. A rejected majorant is not a counterexample."""
from pathlib import Path
import subprocess,struct,shutil
from fractions import Fraction as F
from exact import need
import new_math

def run(exe,finite_folder,rows,work):
    work=Path(work);work.mkdir();r=rows[0];M=r['M'];q=r['q'];orig=(Path(finite_folder)/f'M{M}.bin').read_bytes();out=[]
    variants={
        'nonunit_denominator':struct.pack('<qq',0,2)+orig[16:],
        'false_congruence':struct.pack('<qq',1,1)+orig[16:],
        'truncated_record':orig[:-1],
        'extra_record':orig+struct.pack('<qq',0,1),
    }
    for name,data in variants.items():
        d=work/name;d.mkdir();(d/f'M{M}.bin').write_bytes(data);domain=d/'domain.txt';domain.write_text(f'{M} {q}\n')
        p=subprocess.run([str(exe),'consume',str(domain),str(d)],capture_output=True,text=True)
        need(p.returncode!=0,'corrupt stream unexpectedly accepted');out.append(dict(test=name,returncode=p.returncode,stderr=p.stderr.strip()))
    for name,domaintext in [('duplicate_modulus',f'{M} {q}\n{M} {q}\n'),('wrong_radical',f'{M} {q+1}\n')]:
        d=work/name;d.mkdir();(d/f'M{M}.bin').write_bytes(orig);domain=d/'domain.txt';domain.write_text(domaintext)
        p=subprocess.run([str(exe),'consume',str(domain),str(d)],capture_output=True,text=True)
        need(p.returncode!=0,'invalid domain unexpectedly accepted');out.append(dict(test=name,returncode=p.returncode,stderr=p.stderr.strip()))
    for name,action in [
        ('false_relaxed_equality',lambda:need(new_math.direct(330,30)[0][30]==new_math.direct(330,330)[0][30],'false equality after deleting p')),
        ('illegal_ambient_division',lambda:need((0-30*1)%330==0,'division by ambient gcd leaves lattice')),
        ('prime_already_in_mark',lambda:new_math.repair(330,30,5,0,11)),
        ('composite_prime',lambda:new_math.repair(330,30,9,0,11)),
    ]:
        try:action()
        except ArithmeticError as e:out.append(dict(test=name,output='ArithmeticError: '+str(e)))
        else:raise ArithmeticError('negative fixture unexpectedly accepted: '+name)
    return dict(status='PASS',rejected_inputs=len(out),tests=out,
                coefficient_obstruction=new_math.repair_table()['bounded_negative'],
                insufficient_gate=new_math.large_prime_gates()['negative_p97'])
