"""Actual mutation tests of the witness consumer and mathematical domain guards."""
import json,struct,subprocess,tempfile
from pathlib import Path
from central import *
from inventory import primitive_base,write_json

def run(binary,work,out):
    binary=Path(binary).resolve();work=Path(work).resolve();out=Path(out)
    M=430;original=(work/f'M{M}.bin').read_bytes()
    need(len(original)==16*(M//2+1),'mutation fixture complete')
    records=[]
    mutations=[
      ('nonunit_denominator',struct.pack('<qq',0,43)+original[16:],'FAIL full-modulus denominator unit'),
      ('zero_denominator',struct.pack('<qq',0,0)+original[16:],'FAIL full-modulus denominator unit'),
      ('negative_denominator',struct.pack('<qq',0,-1)+original[16:],'FAIL full-modulus denominator unit'),
      ('false_congruence',struct.pack('<qq',1,1)+original[16:],'FAIL state congruence'),
      ('truncated_record',original[:-1],'FAIL truncated witness'),
      ('extra_byte',original+b'\0','FAIL extra witness records'),
      ('extra_record',original+struct.pack('<qq',0,1),'FAIL extra witness records'),
    ]
    with tempfile.TemporaryDirectory(prefix='cstar353-negative-') as tmp:
        t=Path(tmp)
        def call(label,expect):
            c=subprocess.run([str(binary),'consume',str(t/'domain.txt'),str(t)],capture_output=True,text=True,timeout=10)
            need(c.returncode==1 and c.stderr.strip()==expect,'mutation not rejected correctly: '+label)
            records.append(dict(name=label,exit_code=c.returncode,stderr=c.stderr.strip()))
        for label,data,expect in mutations:
            (t/'domain.txt').write_text('430 430\n');(t/'M430.bin').write_bytes(data);call(label,expect)
        (t/'M430.bin').write_bytes(original)
        for label,domain,expect in [
         ('duplicate_modulus','430 430\n430 430\n','FAIL duplicate modulus'),
         ('wrong_radical','430 10\n','FAIL exact radical support'),
         ('empty_domain','','FAIL empty domain'),
         ('missing_field','430\n','FAIL domain parse'),
         ('extra_field','430 430 1\n','FAIL extra domain fields'),
         ('modulus_over_limit','200000000 10\n','FAIL modulus domain')]:
            (t/'domain.txt').write_text(domain);call(label,expect)
    def reject(label,fun):
        try:fun()
        except ArithmeticError as e:records.append(dict(name=label,stderr=str(e)))
        else:raise ArithmeticError('expected failure did not occur: '+label)
    reject('fourth_prime_in_modulus',lambda:true_upper(70*3,2,5,7))
    reject('composite_prime_parameter',lambda:central_constants(2,5,9))
    reject('missing_support_prime',lambda:true_upper(10,2,5,7))
    # Artificial cutoffs may be off support. They are scalar profile evaluations only.
    value=improved(1110,3,5,37)
    reject('off_support_cutoff_used_as_physical_modulus',lambda:true_upper(1110,3,5,37))
    # Explicit logical falsifier of division by log(M) at length two.
    need(F(3,200)<DELTA<F(3,100),'false normalization comparison')
    records.append(dict(name='log_M_is_not_log_B',formal_model='M=70^2; U=(3/100)*log(70)',wrong_ratio='3/200',correct_ratio='3/100',threshold='2/85',actual_mean_counterexample=False))
    # Signs: negative bounds cannot simply be divided by the smaller base.
    need(F(-1,2)>F(-1),'negative ratio order')
    records.append(dict(name='negative_ratio_direction',left='-1/2',right='-1',inequality='left > right; clip negative upper bounds at zero for a common upper ratio'))
    # It is essential to retain the third prime in the central count.
    def P(q,T):return sum(1 for d in range(1,T+1) if gcd(d,q)==1 for a in range(1,T+1) if gcd(a,d)==1)
    a=P(10,43);b=P(430,43);need(a-b==42,'third-prime removal witness')
    records.append(dict(name='marked_pairs_are_not_full_pairs',q=10,Q=430,T=43,Pq=a,PQ=b,difference=a-b))
    # A declared nonclosed support stays open, not a false positive.
    reject('radical70_leading_gate',lambda:need(coeffs(2,5,7)[0]<DELTA*log_lo(70),'leading majorant at radical 70 is insufficient; not a counterexample'))
    obj=dict(status='PASS',tests=len(records),controls=records,all_errors_actually_executed=True)
    write_json(out/'negative_checks.json',obj);print(json.dumps(obj,indent=2),flush=True);return obj
if __name__=='__main__':
    import sys
    run(sys.argv[1],sys.argv[2],sys.argv[3])
