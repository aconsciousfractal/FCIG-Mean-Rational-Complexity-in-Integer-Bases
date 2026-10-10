#!/usr/bin/env python3
"""Regenerate all decisive finite results. Python>=3.10 and g++ C++17/OpenMP.
No network, third-party Python packages, floating-point gates or assertions.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,hashlib,json,os,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
import arithmetic as ar
import checks

def run(cmd):
    print('RUN '+' '.join(map(str,cmd)),flush=True)
    subprocess.run(list(map(str,cmd)),check=True,cwd=ROOT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS','4')))

def manifest():
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS' and p.suffix not in ['.pyc','.pyo']}
    seen=set()
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        h,name=line.split('  ',1);p=(ROOT/name).resolve()
        ar.need(p.is_relative_to(ROOT) and name not in seen,'manifest path or duplicate');seen.add(name)
        ar.need(p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==h,'manifest '+name)
    ar.need(seen==actual,'manifest complete inventory')
    print('PASS manifest',len(seen),flush=True)

def perform(tmp:Path,sanitize=False):
    cells,cellmeta=ar.make_cells(ROOT)
    ar.need(ar.encode_cells(cells)==(ROOT/'data/cells.txt').read_text(),'all condition/cell data regenerated')
    ms=ar.moduli();ar.need(''.join(str(m)+'\n' for m in ms)==(ROOT/'data/moduli.txt').read_text(),'both finite-domain enumerations')
    flags=['g++','-O3','-std=c++17','-fopenmp']
    if sanitize:flags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
    for name in ['uniform','produce_witnesses','verify_witnesses']:
        run(flags+[ROOT/'code'/f'{name}.cpp','-o',tmp/name])
    prefix=tmp/'uniform.csv'
    run([tmp/'uniform',ROOT/'data/log_table.txt',ROOT/'data/partition.txt',ROOT/'data/shape_centres.tsv',prefix,3,ROOT/'data/cells.txt'])
    uniform=ar.uniform_certificate(prefix)
    independent=checks.run(ROOT,prefix,cells);independent.update(cellmeta)
    w=tmp/'w';w.mkdir()
    run([tmp/'produce_witnesses',ROOT/'data/moduli.txt',w])
    hist=tmp/'histograms.csv';run([tmp/'verify_witnesses',ROOT/'data/moduli.txt',w,hist])
    finite=ar.finite_certificate(hist,w)
    # A real mutation must be rejected by the independent witness consumer.
    target=w/'w210.bin';old=target.read_bytes();mut=bytearray(old);mut[4:8]=(2).to_bytes(4,'little')
    try:
        target.write_bytes(mut)
        r=subprocess.run([str(tmp/'verify_witnesses'),str(ROOT/'data/moduli.txt'),str(w),str(tmp/'mutated.csv')],capture_output=True,text=True)
        ar.need(r.returncode!=0 and 'invalid witness: congruence/unit/range' in r.stderr,'mutation rejected for semantic reason')
    finally:target.write_bytes(old)
    independent['intentional_nonunit_mutation']={'exit_code':r.returncode,'output':r.stderr.strip()}
    ar.need(F(uniform['rows'][2]['total_upper'])>F(3,25),'tighter 0.12 gate is not certified')
    independent['tighter_3_over_25_gate']='FAIL: the certified majorant exceeds 3/25; no such uniform claim is made'
    conclusion=ar.conclusion(uniform,finite)
    return {'uniform_certificate.json':uniform,'finite_certificate.json':finite,'independent_checks.json':independent,'conclusion.json':conclusion},prefix,hist

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sanitize',action='store_true',help='all GNU undefined-behaviour checks, including signed overflow');args=ap.parse_args()
    ar.need(sys.version_info>=(3,10),'Python>=3.10');manifest();start=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='cstar_radical210_') as work:
        certs,prefix,hist=perform(Path(work),args.sanitize)
        for name,obj in certs.items():
            ar.need(obj==json.loads((ROOT/'results'/name).read_text()),'certificate byte-equivalent JSON '+name)
            print('PASS regenerated '+name,flush=True)
        ar.need(hist.read_bytes()==(ROOT/'results/histograms.csv').read_bytes(),'all witness histograms regenerated')
        with (ROOT/'results/finite_bounds.csv').open(newline='') as f:
            table=list(csv.DictReader(f))
        expected=[{'B_at_N1':str(r['M']),'admissible_witness_mu_upper':r['mu_upper'],
                   'admissible_witness_Delta_upper':r['Delta_upper'],
                   'below_2_over_85':str(r['below_2_over_85'])}
                  for r in certs['finite_certificate.json']['rows']]
        ar.need(table==expected,'all finite summary-table rows agree with regenerated certificate')
        for k in [1,2,3]:
            for term in ['cells','radial']:
                ar.need(Path(str(prefix)+f'.{term}{k}').read_bytes()==(ROOT/'results'/f'{term}{k}.tsv').read_bytes(),'full cell/radial table '+str(k))
    print(json.dumps({'status':'PASS','closed_radical':210,'unique_maximizing_pair':[210,1],'Cstar_global_value_determined':False,'sanitizer':args.sanitize,'seconds':round(time.monotonic()-start,3)},indent=2),flush=True)

if __name__=='__main__':
    try:main()
    except (Exception,KeyboardInterrupt) as e:
        print('FAIL '+type(e).__name__+': '+str(e),file=sys.stderr,flush=True);raise SystemExit(1)
