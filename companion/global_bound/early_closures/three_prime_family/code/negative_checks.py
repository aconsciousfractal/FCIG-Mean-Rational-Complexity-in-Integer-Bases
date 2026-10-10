"""Mutation tests of the executable consumer. Each failure is actually invoked."""
from pathlib import Path
import json, struct, subprocess, tempfile
from exact import need


def run(binary, work):
    binary=Path(binary).resolve();work=Path(work).resolve()
    original=(work/'M762.bin').read_bytes()
    need(len(original)==16*(762//2+1),'mutation fixture size')
    out=[]
    def call(root,label,expected):
        p=subprocess.run([str(binary),'consume',str(root/'domain.txt'),str(root)],
                         text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=10)
        need(p.returncode==1,label+': a malformed input was not rejected')
        need(p.stderr.strip()==expected,label+': unexpected failure')
        out.append({'name':label,'exit_code':p.returncode,'stderr':p.stderr.strip()})
    mutations=[
        ('nonunit_denominator',struct.pack('<qq',0,127)+original[16:],'FAIL full-modulus denominator unit'),
        ('zero_denominator',struct.pack('<qq',0,0)+original[16:],'FAIL full-modulus denominator unit'),
        ('negative_denominator',struct.pack('<qq',0,-1)+original[16:],'FAIL full-modulus denominator unit'),
        ('wrong_congruence',struct.pack('<qq',1,1)+original[16:],'FAIL state congruence'),
        ('truncated_byte',original[:-1],'FAIL truncated witness'),
        ('extra_byte',original+b'\x00','FAIL extra witness records'),
        ('extra_record',original+struct.pack('<qq',0,1),'FAIL extra witness records'),
        ('signed_min_numerator',struct.pack('<qq',-(1<<63),1)+original[16:],'FAIL state congruence'),
    ]
    with tempfile.TemporaryDirectory(prefix='cstar-negative-') as folder:
        root=Path(folder)
        for label,data,expected in mutations:
            (root/'domain.txt').write_text('762 762\n')
            (root/'M762.bin').write_bytes(data)
            call(root,label,expected)
        (root/'M762.bin').unlink()
        call(root,'missing_witness','FAIL missing witness')
        (root/'M762.bin').write_bytes(original)
        domains=[
            ('duplicate_modulus','762 762\n762 762\n','FAIL duplicate modulus'),
            ('missing_domain_field','762\n','FAIL domain parse'),
            ('extra_domain_field','762 762 1\n','FAIL extra domain fields'),
            ('wrong_radical','762 6\n','FAIL exact radical support'),
            ('empty_domain','','FAIL empty domain'),
        ]
        for label,data,expected in domains:
            (root/'domain.txt').write_text(data)
            call(root,label,expected)
    return {'status':'PASS','tests':len(out),'controls':out,
            'domain':'one complete valid modulus M=762, mutated independently',
            'a_rejection_is_not_a_negative_mathematical_result':True}

if __name__=='__main__':
    import sys
    result=run(sys.argv[1],sys.argv[2])
    Path(sys.argv[3]).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2),flush=True)
