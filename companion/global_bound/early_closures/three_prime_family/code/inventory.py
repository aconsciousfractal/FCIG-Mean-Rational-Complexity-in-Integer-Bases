"""Exhaustive finite complement of the analytic theorem, two independent inventories."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import csv,json
from exact import need,analytic_upper,out_decimal
CUTOFF=150000

def sieve(n):
    flags=bytearray(b'\1')*(n+1)
    flags[:2]=b'\0\0'
    for p in range(2,isqrt(n)+1):
        if flags[p]:flags[p*p:n+1:p]=b'\0'*((n-p*p)//p+1)
    return [i for i in range(2,n+1) if flags[i]]

def inventories():
    ps=[p for p in sieve((CUTOFF-1)//6) if p>=127]
    d={};a=2
    while a*3*127<CUTOFF:
        b=3
        while a*b*127<CUTOFF:
            for p in ps:
                z=a*b*p
                while z<CUTOFF:
                    need(z not in d,'duplicate exponential inventory')
                    d[z]=p;z*=p
            b*=3
        a*=2
    e={}
    for M in range(762,CUTOFF,6):
        r=M
        while r%2==0:r//=2
        while r%3==0:r//=3
        if r<127:continue
        p=next((t for t in range(5,isqrt(r)+1,2) if r%t==0),r)
        if p<127:continue
        while r%p==0:r//=p
        if r==1:e[M]=p
    need(d==e,'independent complete inventories disagree')
    return d,ps

def produce(out,work):
    out=Path(out);work=Path(work);out.mkdir(parents=True,exist_ok=True);work.mkdir(parents=True,exist_ok=True)
    d,ps=inventories();analytic=[];finite=[]
    for M,p in sorted(d.items()):
        hi=analytic_upper(M,p)
        row={'M':M,'p':p,'upper':str(hi),'upper_decimal':out_decimal(hi)}
        if hi<F(3,20):analytic.append(row)
        else:finite.append(row)
    for name,rows in [('analytic_small',analytic),('finite_domain',finite)]:
        with (out/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['M','p','upper','upper_decimal']);w.writeheader();w.writerows(rows)
    (work/'domain.txt').write_text(''.join(f"{r['M']} {6*r['p']}\n" for r in finite))
    result={'status':'PASS','strict_modulus_cutoff':CUTOFF,'min_prime':127,'inventory_methods':['prime sieve and positive exponent products','all multiples of 6, trial division of stripped cofactor'],'prime_parameters':len(ps),'all_moduli':len(d),'all_state_mass':sum(d),'analytic_small_moduli':len(analytic),'finite_moduli':len(finite),'finite_state_mass':sum(r['M'] for r in finite),'finite_representatives':sum(r['M']//2+1 for r in finite),'largest_finite_modulus':max(r['M'] for r in finite),'largest_finite_prime':max(r['p'] for r in finite),'largest_analytic_upper':str(max(F(r['upper']) for r in analytic))}
    (out/'inventory.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result

if __name__=='__main__':
    import sys
    print(json.dumps(produce(sys.argv[1],sys.argv[2]),indent=2),flush=True)
