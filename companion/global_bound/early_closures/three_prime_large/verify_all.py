#!/usr/bin/env python3
"""Offline full regeneration. Never overwrite the distributed expected results."""
import argparse,hashlib,json,os,pathlib,shutil,subprocess,sys,tempfile,time

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def need(b,msg):
    if not b:raise RuntimeError(msg)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sanitize',action='store_true');ap.add_argument('--report');ap.add_argument('--keep-work',action='store_true');args=ap.parse_args()
    root=pathlib.Path(__file__).resolve().parent
    need(sys.version_info>=(3,10),'Python >=3.10 is required')
    cc=shutil.which('g++');need(cc is not None,'GNU g++ with C++17 and __int128 is required')
    manifest=root/'SHA256SUMS'
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            digest,name=line.split('  ',1);p=root/name
            need(p.is_file() and sha(p)==digest,'manifest mismatch: '+name)
    temp=pathlib.Path(tempfile.mkdtemp(prefix='cstar-g53-replay-'));work=temp/'work';out=temp/'results';logs=temp/'logs'
    for d in (work,out,logs):d.mkdir()
    start=time.monotonic();runs=[]
    flags=['-O','-B'] if sys.flags.optimize else ['-B']
    def run(label,cmd):
        t=time.monotonic();print('RUN',label,flush=True)
        with (logs/(label+'.log')).open('w') as f:
            p=subprocess.run([str(x) for x in cmd],cwd=root,stdout=f,stderr=subprocess.STDOUT)
        if p.returncode:
            print((logs/(label+'.log')).read_text()[-8000:],file=sys.stderr)
            raise RuntimeError(label+' failed with exit '+str(p.returncode)+'; retained work '+str(temp))
        runs.append({'stage':label,'exit_code':0,'seconds':round(time.monotonic()-t,3)})
    try:
        run('inventory',[sys.executable,*flags,root/'code/inventory.py',out,work])
        cflags=['-O3','-std=c++17']
        if args.sanitize:cflags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
        run('compile',[cc,*cflags,root/'code/census.cpp','-o',work/'census'])
        run('producer',[work/'census','produce',work/'domain.txt',work])
        run('consumer',[work/'census','consume',work/'domain.txt',work])
        run('certify',[sys.executable,*flags,root/'code/certify.py',work,out])
        run('independent_checks',[sys.executable,*flags,root/'code/checks.py',out,out/'independent_checks.json'])
        run('negative_checks',[sys.executable,*flags,root/'code/negative_checks.py',work/'census',work,out/'consumer_negative_checks.json'])
        expected=sorted(p.name for p in (root/'results').iterdir() if p.is_file())
        generated=sorted(p.name for p in out.iterdir() if p.is_file())
        need(expected==generated,'expected/generated file sets differ')
        for name in expected:need((root/'results'/name).read_bytes()==(out/name).read_bytes(),'regenerated mismatch: '+name)
        report={'status':'PASS','expected_results_overwritten':False,'all_expected_results_regenerated':len(expected),'python_optimized':bool(sys.flags.optimize),'undefined_behavior_sanitizer':args.sanitize,'runs':runs,'seconds':round(time.monotonic()-start,3),'all_three_prime_supports_with_largest_prime_ge_53_excluded_after_declared_F6p_import':True,'remaining_supports':442,'F6p_import_replayed':False,'global_Cstar_identified':False}
        if args.report:pathlib.Path(args.report).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
        print(json.dumps(report,indent=2),flush=True)
        if args.keep_work:print('Retained replay directory:',temp)
        else:shutil.rmtree(temp)
    except BaseException:
        print('Replay directory retained:',temp,file=sys.stderr);raise
if __name__=='__main__':main()
