#!/usr/bin/env python3
"""Rebuild and compare every new mathematical result without rewriting the expected files.
Requires Python >=3.10 and GNU C++17 with __int128, no third-party Python modules.
The imported S113 theorem is an explicit dependency ONLY of the combined corollary.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
from exact import need
import inventory,certify,checks,negative_checks,family

EXPECTED=[
 'analytic_small.csv','candidate.json','constants.json','finite_certificate.csv',
 'finite_domain.csv','finite_summary.json','independent_checks.json','inventory.json',
 'logarithm_checks.json','consumer_negative_checks.json','family_splice.json',
]

def write_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def command(args,log):
    print('RUN',' '.join(map(str,args)),flush=True)
    with log.open('w') as f:
        p=subprocess.run(list(map(str,args)),stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:
        print(log.read_text()[-16000:],flush=True)
        raise ArithmeticError('subprocess failed: '+str(p.returncode))
    print('PASS',log.name,flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sanitize',action='store_true',help='GNU undefined-behavior checks, no recovery')
    ap.add_argument('--compiler',default=os.environ.get('CXX','g++'))
    ap.add_argument('--keep-work',type=Path,help='use a new directory and retain the large regenerated streams')
    args=ap.parse_args()
    need(sys.version_info>=(3,10),'Python 3.10 or later required')
    need(shutil.which(args.compiler) is not None,'GNU C++ compiler not found')
    start=time.monotonic()
    if args.keep_work:
        work=args.keep_work.resolve();work.mkdir(parents=True,exist_ok=False);tmp=None
    else:
        tmp=tempfile.TemporaryDirectory(prefix='cstar-family-6p-');work=Path(tmp.name)
    try:
        result=work/'results';result.mkdir()
        print('Temporary work:',str(work),flush=True)
        # Inventory contains the mathematical finite-domain generation, checked two ways.
        inv=inventory.produce(result,work)
        print('PASS inventory',json.dumps(inv,sort_keys=True),flush=True)
        binary=work/'census'
        flags=['-O2','-std=c++17','-Wall','-Wextra']
        if args.sanitize:flags+=['-fsanitize=undefined','-fno-sanitize-recover=all']
        command([args.compiler,*flags,ROOT/'code/census.cpp','-o',binary],work/'compile.log')
        command([binary,'produce',work/'domain.txt',work],work/'produce.log')
        command([binary,'consume',work/'domain.txt',work],work/'consume.log')
        cert=certify.run(work,result)
        print('PASS exact means and all logarithm inputs',flush=True)
        write_json(result/'independent_checks.json',checks.run())
        print('PASS independent finite banks',flush=True)
        write_json(result/'consumer_negative_checks.json',negative_checks.run(binary,work))
        print('PASS consumer mutation controls',flush=True)
        write_json(result/'family_splice.json',family.run(ROOT))
        print('PASS declared-import interface',flush=True)
        actual=sorted(p.name for p in result.iterdir())
        need(actual==sorted(EXPECTED),'new result file inventory mismatch')
        expected=sorted(p.name for p in (ROOT/'results').iterdir())
        need(expected==sorted(EXPECTED),'expected result file inventory mismatch')
        for name in EXPECTED:
            a=(result/name).read_bytes();b=(ROOT/'results'/name).read_bytes()
            need(a==b,'regenerated result differs: '+name)
            print('MATCH',name,hashlib.sha256(a).hexdigest(),flush=True)
        print(json.dumps({'status':'PASS','generated_expected_files':False,
            'new_result_files_compared':len(EXPECTED),
            'new_G127_finite_complement_states':cert['all_state_mass'],
            'new_G127_finite_complement_moduli':cert['moduli'],
            'all_new_certificates_regenerated':True,
            'consumer_mutations_replayed':True,
            'S113_full_proof_replayed_this_task':False,
            'entire_6p_family_closed_with_imported_S113':True,
            'global_Cstar_identified':False,'GNU_sanitizers_enabled':args.sanitize,
            'python_assertions_enabled':__debug__,
            'elapsed_seconds':round(time.monotonic()-start,3)},indent=2),flush=True)
    finally:
        if tmp is not None:tmp.cleanup()

if __name__=='__main__':
    try:main()
    except (ArithmeticError,OSError,subprocess.SubprocessError) as e:
        print(type(e).__name__+': '+str(e),file=sys.stderr,flush=True)
        sys.exit(1)
