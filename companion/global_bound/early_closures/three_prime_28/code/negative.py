"""Executed rejection tests. A failed upper-bound gate is never a counterexample."""
from pathlib import Path
from fractions import Fraction as F
import shutil,subprocess,struct
import arithmetic as ar
import finite

def run(executable,work):
    out={}
    try:ar.need(F(1,8)<F(2*ar.logs(30)[0],85*ar.LS),'k=2 uniform cap alone does not close B=30; finite case required')
    except ArithmeticError as e:out['insufficient_uniform_gate']='ArithmeticError: '+str(e)
    else:raise ArithmeticError('negative gate unexpectedly passed')
    try:ar.need(finite.primes_to_113()[:-1]==finite.primes_to_113(),'prime block missing 113')
    except ArithmeticError as e:out['incomplete_domain']='ArithmeticError: '+str(e)
    else:raise ArithmeticError('domain deletion accepted')
    root=Path(work);root.mkdir(exist_ok=True);dom=root/'domain.txt';dom.write_text('30 30\n')
    good=Path(root/'M30.bin');valid=b''.join(struct.pack('<qq',s,1) for s in range(16))
    tests={'nonunit':struct.pack('<qq',0,5)+valid[16:], 'wrong_state':struct.pack('<qq',1,1)+valid[16:], 'truncation':valid[:-1], 'extra_record':valid+struct.pack('<qq',0,1)}
    expected={'nonunit':'full-modulus denominator unit','wrong_state':'state congruence','truncation':'truncated witness','extra_record':'extra witness records'}
    for label,data in tests.items():
        good.write_bytes(data)
        cp=subprocess.run([str(executable),'consume',str(dom),str(root)],capture_output=True,text=True)
        ar.need(cp.returncode!=0 and expected[label] in cp.stderr,'mutated witness stream rejected: '+label)
        out[label]=cp.stderr.strip()
    good.unlink()
    return out
