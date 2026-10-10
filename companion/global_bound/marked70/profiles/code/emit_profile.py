#!/usr/bin/env python3
from pathlib import Path
import sys,json,hashlib,argparse
from assemble import profile
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('p',type=int);ap.add_argument('ag',type=int);ap.add_argument('k',type=int);ap.add_argument('prefix',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args()
v=profile(a.p,Path(str(a.prefix)+'.cells'),ROOT,a.k,a.ag)
w=Path(str(a.prefix)+'.witness');h=hashlib.sha256()
with w.open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
v['witness_bytes']=w.stat().st_size;v['witness_sha256']=h.hexdigest()
a.output.write_text(json.dumps(v,indent=2)+'\n')
