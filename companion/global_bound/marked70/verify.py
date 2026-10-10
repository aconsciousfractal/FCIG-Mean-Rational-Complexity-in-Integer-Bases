"""Check compact marked-70 certificates; optionally reconstruct every raw witness."""
from pathlib import Path, PurePosixPath
import argparse, hashlib, json, subprocess, sys, time, uuid

ROOT=Path(__file__).resolve().parent
RUNTIME={'runs','bin','work','logs','__pycache__'}

def ensure(path):
    missing=[]
    while not path.exists():
        missing.append(path);path=path.parent
    for p in reversed(missing):
        print(p.resolve(),flush=True);p.mkdir()

def pairs(items):
    result={}
    for k,v in items:
        if k in result:raise ValueError('duplicate manifest key')
        result[k]=v
    return result

def manifest_check():
    expected=json.loads((ROOT/'MANIFEST_SHA256.json').read_bytes(),object_pairs_hook=pairs)
    if not isinstance(expected,dict) or not expected:raise ValueError('nonempty manifest required')
    for name,h in expected.items():
        pp=PurePosixPath(name)
        if pp.is_absolute() or pp.as_posix()!=name or any(v in ('.','..') for v in pp.parts) or '\\' in name or ':' in name:
            raise ValueError('noncanonical member')
        p=ROOT.joinpath(*pp.parts)
        if not p.resolve().is_relative_to(ROOT) or not p.is_file() or any(q.is_symlink() for q in [p,*p.parents] if q.is_relative_to(ROOT)):raise ValueError('member containment')
        if not isinstance(h,str) or len(h)!=64 or any(c not in '0123456789abcdef' for c in h):raise ValueError('canonical sha256')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise ValueError('member digest '+name)
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
            and p.relative_to(ROOT).parts[0] not in RUNTIME and '__pycache__' not in p.parts
            and p.suffix not in ('.pyc','.pyo') and p.relative_to(ROOT).as_posix()!='MANIFEST_SHA256.json'}
    if actual!=set(expected):raise ValueError('exact public payload inventory')
    return len(expected)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group();mode.add_argument('--quick',action='store_true');mode.add_argument('--rebuild',action='store_true')
    ap.add_argument('--output-dir',type=Path);ap.add_argument('--report',type=Path)
    a=ap.parse_args();count=manifest_check()
    out=(a.output_dir or ROOT/'runs'/uuid.uuid4().hex).resolve()
    if out.exists():raise ValueError('output directory must be new')
    if out.is_relative_to(ROOT) and out.relative_to(ROOT).parts[0] not in RUNTIME:raise ValueError('output must not overlap scientific payload')
    ensure(out);rows=[];flag=['-O'] if not __debug__ else []
    stages=[('profiles',ROOT/'profiles/verify_profiles.py',['--quick','--report',str(out/'profiles.json')]),
            ('repair70',ROOT/'checks/verify_repair70.py',[]),('repair770',ROOT/'checks/verify_repair770.py',[]),
            ('gates70',ROOT/'checks/verify_gates70.py',[]),('sharpened',ROOT/'checks/verify_sharpened_pe.py',[]),
            ('closure',ROOT/'review/check_final_closure.py',[])]
    if a.rebuild:
        stages[:0]=[('geometric_full',ROOT/'reconstruct_profiles.py',['--work-dir',str(out/'geometric')]),
                    ('finite_full',ROOT/'code/run_finite70.py',['--work-dir',str(out/'finite'),'--threads','4','--compare']),
                    ('negative_fresh',ROOT/'profiles/code/negative_profiles.py',['--consumer',str(out/'geometric/bin/consume_1_4'),'--output',str(out/'negative.json'),'--log',str(out/'negative.log')])]
    for name,path,args in stages:
        start=time.monotonic();limit=3600 if name=='geometric_full' else 600
        with (out/(name+'.log')).open('w',encoding='utf8') as log:
            try:
                z=subprocess.run([sys.executable,'-B',*flag,str(path),*args],stdout=log,stderr=subprocess.STDOUT,timeout=limit)
                row=dict(stage=name,returncode=z.returncode,seconds=time.monotonic()-start,timeout_seconds=limit)
            except subprocess.TimeoutExpired:
                row=dict(stage=name,returncode=None,status='TIMEOUT',seconds=time.monotonic()-start,timeout_seconds=limit)
        rows.append(row);(out/'STAGES.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8')
        if row['returncode']!=0:raise RuntimeError('stage failed: '+name+'; output retained at '+str(out))
        print('PASS',name,flush=True)
    if manifest_check()!=count:raise ValueError('payload changed during checking')
    report=dict(status='PASS',payload_files=count,stages=rows,
                full_geometry_rebuilt_in_this_command=a.rebuild,finite_streams_replayed_in_this_command=a.rebuild,
                quick_geometric_scope='96 stratified witness samples per each of 18 profiles; complete HNF algebra and stored reintegration',
                quick_geometric_audit_exhaustive=False,raw_streams_in_public_payload=False)
    (out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    if a.report:
        p=a.report.resolve()
        if p.is_relative_to(ROOT) and p.relative_to(ROOT).as_posix() in json.loads((ROOT/'MANIFEST_SHA256.json').read_text()):raise ValueError('report overlaps input')
        ensure(p.parent);p.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
