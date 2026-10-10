#!/usr/bin/env python3
"""Exact replay of the 28 complete marked-radical families and the combined frontier.

Python >= 3.10; g++ with C++17, OpenMP, __int128. No network or solver.
Normal mode never rewrites expected results. --generate is explicit developer mode.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile, time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
import arithmetic as ar
import finite, certify, independent, banks, combined, negative
import checks_mod6, routes_mod6, excision

def jb(x): return (json.dumps(x,sort_keys=True,ensure_ascii=False,indent=2)+'\n').encode()
def run(cmd):
    print('RUN '+' '.join(map(str,cmd)),flush=True)
    subprocess.run(list(map(str,cmd)),check=True,cwd=ROOT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',CSTAR_THREADS=os.environ.get('CSTAR_THREADS','4'),OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS','2')))
def manifest():
    seen=set()
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        h,n=line.split('  ',1);p=(ROOT/n).resolve()
        ar.need(p.is_relative_to(ROOT) and n not in seen,'manifest path/duplicate');seen.add(n)
        ar.need(p.is_file() and finite.file_hash(p)==h,'manifest '+n)
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS' and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo')}
    ar.need(actual==seen,'complete manifest inventory')
    print('PASS manifest '+str(len(seen)),flush=True)
def put(n,data,generate):
    p=ROOT/'results'/n
    if generate:p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
    else:ar.need(p.is_file() and p.read_bytes()==data,'exact regenerated match '+n)
    print(('GENERATED ' if generate else 'PASS regenerated ')+n,flush=True)

def inherited(w,flags,cells,meta,generate):
    prefix=w/'mod6';run(flags+[ROOT/'imports/mod6/uniform.cpp','-o',w/'uniform_mod6'])
    run([w/'uniform_mod6',ROOT/'data/log_table.txt',ROOT/'data/partition.txt',ROOT/'data/shape_centres.tsv',prefix,5,ROOT/'data/cells.txt'])
    u=ar.uniform_certificate(prefix)
    tests=checks_mod6.run(ROOT,prefix,cells);tests.update(meta)
    tests['all_geometric_routes']=routes_mod6.check(ROOT,prefix)
    tests['complete_physical_banks']=excision.bank_checks(prefix)
    try:ar.need(F(u['rows'][4]['total_upper'])<F(-39,1000),'stronger uniform cap -39/1000 is not certified')
    except ArithmeticError as e:tests['negative_uniform_gate']='ArithmeticError: '+str(e)
    else:raise ArithmeticError('stronger-cap negative test did not fail')
    pe=excision.finite_checks();cand=excision.candidate();gates=excision.family_gates(u)
    for n,x in [('uniform_mod6.json',u),('independent_mod6.json',tests),('prime_excision_checks.json',pe),('candidate.json',cand),('large_prime_gates.json',gates)]:put(n,jb(x),generate)
    for k in range(1,6):
        for kind in ('cells','radial'):put(f'import_mod6/{kind}{k}.tsv',Path(str(prefix)+f'.{kind}{k}').read_bytes(),generate)
    return gates

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--sanitize',action='store_true');ap.add_argument('--generate',action='store_true',help='Explicit developer mode: write expected results; not a comparison replay')
    args=ap.parse_args();ar.need(sys.version_info>=(3,10),'Python >= 3.10');ar.need(shutil.which('g++') is not None,'g++ required')
    if not args.generate:manifest()
    start=time.monotonic()
    flags=['g++','-O3','-std=c++17','-fopenmp']
    if args.sanitize:flags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
    cells,meta=ar.make_cells(ROOT);ar.need(ar.encode_cells(cells)==(ROOT/'data/cells.txt').read_text(),'all cell-condition data regenerated')
    with tempfile.TemporaryDirectory(prefix='cstar_28_complete_') as td:
        w=Path(td);rows=[];checks=[];physical=[];common_hash=None
        for p,maxk in certify.configurations():
            pd=w/f'p{p:03}';pd.mkdir();prefix=pd/'engine';exe=pd/'uniform'
            run(flags+[f'-DPRIME={p}','-DDIR_T=16','-DAUDIT_SPARSE=1',ROOT/'code/uniform.cpp','-o',exe])
            run([exe,ROOT/'data/log_table.txt',ROOT/'data/partition.txt',ROOT/'data/shape_centres.tsv',prefix,maxk,ROOT/'data/cells.txt'])
            ch=finite.file_hash(Path(str(prefix)+'.cellcosts'))
            if common_hash is None:
                common_hash=ch;common=dict(meta,**independent.common_costs(ROOT,prefix,cells));put('common_independent.json',jb(common),args.generate)
            else:ar.need(ch==common_hash,'identical common geometric costs for all p')
            for k,cap in enumerate(certify.rational_caps(p),1):
                ar.need(k<=maxk,'selected horizon');rows.append(certify.profile(p,k,prefix,cap))
                for kind in ('cells','radial'):put(f'profiles/p{p:03}.{kind}{k}.tsv',Path(str(prefix)+f'.{kind}{k}').read_bytes(),args.generate)
            checks.append(independent.run_one(ROOT,p,prefix));physical.append(banks.run(p,prefix))
            print(f'PASS p={p}: all declared horizons, all initial witnesses, independent routes and complete physical banks',flush=True)
            shutil.rmtree(pd);ar.logs.cache_clear()
        ar.need(len(rows)==33,'exact profile count')
        put('uniform_certificates.json',jb({'status':'PASS','direction_count':64,'complete_prime_list':finite.primes_to_113(),'profiles':rows}),args.generate)
        put('independent_checks.json',jb({'status':'PASS','per_prime':checks}),args.generate)
        put('physical_banks.json',jb({'status':'PASS','per_prime':physical}),args.generate)
        cases,mods,inventory=finite.domains();cw=w/'census';cw.mkdir();domain=cw/'domain.txt';domain.write_text(''.join(f'{M} {q}\n' for M,q in sorted(mods.items())))
        run(flags+[ROOT/'code/census.cpp','-o',w/'census_engine'])
        run([w/'census_engine','produce',domain,cw]);run([w/'census_engine','consume',domain,cw])
        census,hist=finite.census_certificate(cw)
        put('census_certificate.json',jb(census),args.generate)
        put('witness_histograms.json',(json.dumps(hist,sort_keys=True,separators=(',',':'))+'\n').encode(),args.generate)
        put('negative_tests.json',jb(negative.run(w/'census_engine',w/'mutations')),args.generate)
        cover=certify.covering_gates(rows,census);put('closure_certificate.json',jb(cover),args.generate)
        old=inherited(w,flags,cells,meta,args.generate)
        put('combined_frontier.json',jb(combined.certificate(cover,old)),args.generate)
    print(json.dumps({'status':'PASS','sanitizer':args.sanitize,'generated_expected_files':args.generate,'all_28_radicals_closed_all_lengths':True,'all_radical_6p_pairs_N_ge_5_excluded':True,'entire_6p_family_closed':False,'global_Cstar_identified':False,'seconds':round(time.monotonic()-start,3)},indent=2),flush=True)
if __name__=='__main__':
    try:main()
    except (Exception,KeyboardInterrupt) as e:
        print('FAIL '+type(e).__name__+': '+str(e),file=sys.stderr,flush=True);raise SystemExit(1)
