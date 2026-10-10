#!/usr/bin/env python3
"""Regenerate all new evidence, then compare with immutable expected results.

S7 (the seven previously closed small-prime sectors) is a scoped internal import,
not a claim of a fresh replay. No global value of C* is identified by this script.
"""
from pathlib import Path
import argparse,csv,hashlib,io,json,os,shutil,subprocess,sys,tempfile,time,platform
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
from exact import need
from assemble import profile
from checks import cases,stem,run as independent
from inventory import make
from closure import finite_certificate,assemble
from repair import repair_table,simple_table,periodic,physical_checks,candidate
from negative import run as negative

def jsonbytes(obj):return (json.dumps(obj,sort_keys=True,indent=2)+'\n').encode()
def csvbytes(rows):
    s=io.StringIO(newline='');w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows);return s.getvalue().encode()
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def summarybytes(profiles):
    rows=[{k:z[k] for k in ('p','q','k','transfer_index','angular_nodes','total_upper','upper_decimal','witnesses')} for z in profiles]
    return csvbytes(rows)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--threads',type=int,choices=range(1,5),default=4)
    a.add_argument('--sanitize',action='store_true',help='GNU undefined-behaviour checks; fail on first finding')
    a.add_argument('--generate',action='store_true',help='Explicit development mode: create expected files, NOT a comparison replay')
    a.add_argument('--keep-work',action='store_true')
    a.add_argument('--work-dir',type=Path,help='New, empty temporary work directory')
    a.add_argument('--report',type=Path,help='Write a metadata report outside the expected-results directory')
    args=a.parse_args();need(sys.version_info>=(3,10),'Python 3.10 or later');need(sys.byteorder=='little','little-endian native witness records')
    need(shutil.which('g++') is not None,'GNU C++ compiler is required')
    if args.report:need(ROOT/'results' not in args.report.resolve().parents,'do not write a report over expected results')
    temp=args.work_dir or Path(tempfile.mkdtemp(prefix='cstar30p-'))
    if args.work_dir:temp.mkdir(parents=True,exist_ok=False)
    out=temp/'profiles';out.mkdir();bdir=temp/'bin';bdir.mkdir();fin=temp/'finite';fin.mkdir()
    env={**os.environ,'CSTAR_THREADS':str(args.threads)};checked=[];started=time.monotonic()
    def compare(rel,data):
        p=ROOT/rel
        if args.generate:p.parent.mkdir(exist_ok=True,parents=True);p.write_bytes(data)
        else:need(p.is_file() and p.read_bytes()==data,'regenerated file mismatch: '+rel)
        checked.append(rel)
    def put(name,obj):compare('results/'+name,jsonbytes(obj));print('PASS file',name,flush=True)
    def command(cmd):
        print('RUN',' '.join(map(str,cmd)),flush=True)
        subprocess.run(list(map(str,cmd)),check=True,env=env)
    flags=['-O3','-std=c++17','-fopenmp']
    if args.sanitize:flags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
    def compile_(name,source,extra=()):
        dst=bdir/name;command(['g++',*flags,*extra,ROOT/'code'/source,'-o',dst]);return dst
    report={}
    try:
        bins={}
        for ag,k in ((4,1),(16,1),(16,2)):
            for consumer in (False,True):
                extra=[f'-DANGULAR={ag}',f'-DDEGREE_EXP={k}']+(['-DBUILD_CONSUMER=1'] if consumer else [])
                bins[(ag,k,consumer)]=compile_(f'{"consume" if consumer else "produce"}_{ag}_{k}','angle.cpp',extra)
        algebra=compile_('algebra','algebra.cpp');finite=compile_('finite','finite.cpp')
        ds=[];profiles=[];validgeo=temp/'negative_model.witness'
        for p,k,ag in cases():
            name=stem(p,k);prefix=out/name
            command([bins[(ag,k,False)],'produce',p,ROOT/'data',prefix])
            command([bins[(ag,k,True)],'consume',p,ROOT/'data',prefix])
            for ext in ('.cells','.branches','.audit'):
                v=Path(str(prefix)+ext).read_bytes()
                need(v==Path(str(prefix)+'.checked'+ext).read_bytes(),'producer/consumer mismatch '+name+ext)
                compare('results/profiles/'+name+ext,v)
            command([algebra,p,k,ROOT/'data',str(prefix)+'.tensor'])
            compare('results/profiles/'+name+'.tensor',Path(str(prefix)+'.tensor').read_bytes())
            wf=Path(str(prefix)+'.witness');ds.append(dict(p=p,q=30*p if p>1 else 30,k=k,bytes=wf.stat().st_size,sha256=digest(wf)))
            if p==1 and k==1:shutil.copyfile(wf,validgeo)
            wf.unlink()
            profiles.append(profile(p,Path(str(prefix)+'.cells'),ROOT,k,ag))
        put('uniform_profiles.json',profiles);put('witness_digests.json',ds)
        compare('results/summary.csv',summarybytes(profiles))
        rows,inv=make();put('inventory_certificate.json',inv)
        compare('results/finite_inventory.csv',csvbytes(rows))
        wr=[r for r in rows if r['route']=='witness'];compare('results/witness_inventory.csv',csvbytes(wr))
        domain=''.join(f"{r['M']} {r['p']}\n" for r in wr).encode();compare('data/finite_domain.txt',domain)
        domainfile=temp/'domain.txt';domainfile.write_bytes(domain)
        fcsv=temp/'finite_bounds.csv';command([finite,'run',domainfile,ROOT/'data',fin,fcsv,args.threads])
        compare('results/finite_bounds.csv',fcsv.read_bytes());fc=finite_certificate(fcsv,rows);put('finite_certificate.json',fc)
        put('independent_checks.json',independent(out,ROOT))
        for name,fn in [('repair_table.json',repair_table),('simple_repair.json',simple_table),('periodic_certificate.json',periodic),('prime_repair_banks.json',physical_checks),('candidate.json',candidate)]:put(name,fn())
        put('negative_checks.json',negative(ROOT,temp/'negative',finite,bins[(16,1,True)],validgeo))
        cl=assemble(profiles,rows,inv,fc);put('closure.json',cl)
        # The seven small-prime imports are joined by verify_early.py;
        # this case reconstructs only its own declared numerical inputs.
        # All expected results must be covered by regeneration; no stale unverified extras.
        actual={p.relative_to(ROOT).as_posix() for p in (ROOT/'results').rglob('*') if p.is_file()}
        need(actual<={x for x in checked if x.startswith('results/')},'all expected results regenerated')
        report=dict(status='GENERATED' if args.generate else 'PASS',generated_expected_files=args.generate,
           python_optimized=not __debug__,undefined_behavior_sanitizer=args.sanitize,threads=args.threads,
           regenerated_compared_files=len(checked),expected_result_files=len(actual),
           geometric_profiles=37,new_full_mark_supports=35,
           finite_moduli=fc['moduli'],finite_states=fc['states'],finite_representatives=fc['representatives'],
           all_p_ge_37_closed_all_valuations_and_lengths=True,entire_30p_family_closed_with_S7_import=True,
           old_S7_import_replayed=False,global_Cstar_identified=False,
           strong_mu_4_over_25_claimed=False,common_upper=cl['common_upper'],
           runtime_seconds=round(time.monotonic()-started,3),python=sys.version.split()[0],platform=platform.platform())
        if args.report:args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_bytes(jsonbytes(report))
        print(json.dumps(report,indent=2),flush=True)
    finally:
        if args.keep_work:print('Retained work directory:',temp,flush=True)
        else:shutil.rmtree(temp,ignore_errors=True)

if __name__=='__main__':
    try:main()
    except (ArithmeticError,RuntimeError,subprocess.CalledProcessError,OSError) as e:
        print(type(e).__name__+': '+str(e),file=sys.stderr,flush=True);sys.exit(1)
