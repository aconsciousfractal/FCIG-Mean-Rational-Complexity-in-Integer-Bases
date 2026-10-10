#!/usr/bin/env python3
"""Offline complete regeneration. Expected results are never overwritten."""
import argparse,hashlib,json,pathlib,shutil,subprocess,sys,tempfile,time

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sanitize',action='store_true');ap.add_argument('--report');ap.add_argument('--keep-work',action='store_true');a=ap.parse_args()
    root=pathlib.Path(__file__).resolve().parent
    need(sys.version_info>=(3,10),'Python >=3.10 required')
    cc=shutil.which('g++');need(cc is not None,'GNU g++ with C++17 and __int128 required')
    manifest=root/'SHA256SUMS'
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            digest,name=line.split('  ',1);p=root/name
            need(p.is_file() and sha(p)==digest,'manifest mismatch: '+name)
    tmp=pathlib.Path(tempfile.mkdtemp(prefix='cstar353-replay-'))
    work=tmp/'work';out=tmp/'results';logs=tmp/'logs'
    for p in (work,out,logs):p.mkdir()
    flags=['-O','-B'] if sys.flags.optimize else ['-B'];start=time.monotonic();runs=[]
    def run(label,cmd):
        t=time.monotonic();print('RUN',label,flush=True)
        with (logs/(label+'.log')).open('w') as f:
            p=subprocess.run([str(x) for x in cmd],cwd=root,stdout=f,stderr=subprocess.STDOUT)
        if p.returncode:
            print((logs/(label+'.log')).read_text()[-10000:],file=sys.stderr)
            raise RuntimeError(label+' failed; retained work '+str(tmp))
        runs.append(dict(stage=label,exit_code=0,seconds=round(time.monotonic()-t,3)))
    try:
        run('inventory',[sys.executable,*flags,root/'code/inventory.py',out,work])
        run('cost_audit',[sys.executable,*flags,root/'code/cost_audit.py',out])
        cflags=['-O3','-std=c++17']
        if a.sanitize:cflags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
        run('compile',[cc,*cflags,root/'code/census.cpp','-o',work/'census'])
        run('producer',[work/'census','produce',work/'domain.txt',work])
        run('consumer',[work/'census','consume',work/'domain.txt',work])
        run('certify',[sys.executable,*flags,root/'code/certify.py',work,out])
        run('independent_checks',[sys.executable,*flags,root/'code/checks.py',out])
        run('negative_checks',[sys.executable,*flags,root/'code/negative_checks.py',work/'census',work,out])
        expected=sorted(p.name for p in (root/'results').iterdir() if p.is_file())
        generated=sorted(p.name for p in out.iterdir() if p.is_file())
        need(expected==generated,'expected/generated file sets differ')
        for name in expected:
            need((root/'results'/name).read_bytes()==(out/name).read_bytes(),'regeneration mismatch: '+name)
        # Independent consistency check of the descriptive worst-route locator.
        import csv
        cl=json.loads((out/'closure.json').read_text());label=cl['worst_majorant_route']
        if label.startswith('analytic M='):
            M=int(label.split('=')[1]);rr=[x for x in csv.DictReader((out/'analytic_moduli.csv').open()) if int(x['M'])==M]
            need(len(rr)==1 and int(rr[0]['radical'])==cl['worst_route_radical'],'worst-route metadata mismatch')
        report=dict(status='PASS',expected_results_overwritten=False,all_expected_results_regenerated=len(expected),python_optimized=bool(sys.flags.optimize),undefined_behavior_sanitizer=a.sanitize,runs=runs,seconds=round(time.monotonic()-start,3),selected_supports_closed=cl['selected_supports_closed'],remaining_three_prime_supports_after_declared_imports=89,G53_large_prime_closure_replayed=False,F6p_replayed=False,global_Cstar_identified=False)
        if a.report:pathlib.Path(a.report).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
        print(json.dumps(report,indent=2),flush=True)
        if a.keep_work:print('Retained replay directory:',tmp,flush=True)
        else:shutil.rmtree(tmp)
    except BaseException:
        print('Replay directory retained:',tmp,file=sys.stderr);raise
if __name__=='__main__':main()
