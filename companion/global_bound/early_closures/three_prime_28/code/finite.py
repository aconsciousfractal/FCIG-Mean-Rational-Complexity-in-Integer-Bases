from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json,hashlib
import arithmetic as ar

def primes_to_113():
    return [n for n in range(5,114) if all(n%d for d in range(2,isqrt(n)+1))]
def smooth(p,limit):
    out=set();a=2
    while a*3*p<limit:
        b=3
        while a*b*p<limit:
            c=p
            while a*b*c<limit:out.add(a*b*c);c*=p
            b*=3
        a*=2
    return sorted(out)
def domains():
    cases=[];inventory=[]
    for p,cut in ((5,200000),(7,10000),(11,250),(13,90)):
        xs=smooth(p,cut);ys=[]
        for n in range(6*p,cut,6*p):
            r=n
            for t in (2,3,p):
                while r%t==0:r//=t
            if r==1:ys.append(n)
        ar.need(xs==ys,'two full finite base censuses')
        inventory.append({'p':p,'N':1,'B_strict_cutoff':cut,'bases':xs})
        cases.extend({'p':p,'B':B,'N':1,'M':B} for B in xs)
    xs=smooth(5,210)
    ar.need(xs==[30,60,90,120,150,180],'complete N=2 residue')
    inventory.append({'p':5,'N':2,'B_strict_cutoff':210,'bases':xs})
    cases.extend({'p':5,'B':B,'N':2,'M':B*B} for B in xs)
    mods={x['M']:6*x['p'] for x in cases}
    return cases,mods,inventory

def file_hash(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        while chunk:=f.read(1<<20):h.update(chunk)
    return h.hexdigest()
def upper_decimal(x,places=12):
    z=ar.ceildiv(x.numerator*10**places,x.denominator)
    sign='-' if z<0 else '';z=abs(z)
    return sign+str(z//10**places)+'.'+str(z%10**places).zfill(places)

def census_certificate(folder):
    folder=Path(folder);cases,mods,inventory=domains();rows=[];means={};histout={}
    for M,q in sorted(mods.items()):
        hist=[tuple(map(int,l.split())) for l in (folder/f'M{M}.hist').read_text().splitlines()]
        ar.need(all(1<=h<=M//2 and n>0 for h,n in hist),'valid histogram')
        ar.need(all(hist[j][0]<hist[j+1][0] for j in range(len(hist)-1)),'sorted distinct heights')
        ar.need(sum(n for h,n in hist)==M,'full histogram mass')
        upper=sum(n*ar.logs(h)[1] for h,n in hist)
        mu=F(upper,M*ar.LS)-F(ar.logs(M)[0],2*ar.LS)
        means[M]=mu;histout[str(M)]=hist
        rows.append({'M':M,'q':q,'representatives':M//2+1,'mass':M,'histogram_bins':len(hist),'witness_stream_sha256':file_hash(folder/f'M{M}.bin'),'mu_upper':str(mu),'mu_upper_decimal':upper_decimal(mu)})
    tests=[]
    for x in cases:
        mu=means[x['M']];gate=F(2*ar.logs(x['B'])[0],85*ar.LS)
        ar.need(mu<gate,f"finite mean gate B={x['B']} N={x['N']}")
        delta_hi=mu/F(ar.logs(x['B'])[0],ar.LS) if mu>=0 else mu/F(ar.logs(x['B'])[1],ar.LS)
        tests.append(dict(x,mu_upper=str(mu),delta_upper=str(delta_hi),positive_margin=str(gate-mu)))
    return {'status':'PASS','inventory':inventory,'moduli':rows,'case_comparisons':tests,'counts':{'cases':len(cases),'distinct_moduli':len(mods),'full_state_mass':sum(mods),'representatives':sum(M//2+1 for M in mods),'histogram_bins':sum(r['histogram_bins'] for r in rows)},'worst_delta_upper':str(max(F(x['delta_upper']) for x in tests))},histout

if __name__=='__main__':
    import sys
    cert,hist=census_certificate(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n')
    Path(sys.argv[3]).write_text(json.dumps(hist,sort_keys=True,separators=(',',':'))+'\n')
    print(cert['counts'],'worst_delta',upper_decimal(F(cert['worst_delta_upper'])))
