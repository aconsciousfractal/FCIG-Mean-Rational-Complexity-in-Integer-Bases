#!/usr/bin/env python3
"""Reintegrate stored cell caps. This does not replay the geometric witnesses."""
from pathlib import Path
import argparse,json
from assemble import profile
from exact import need
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('p',type=int);ap.add_argument('ag',type=int);ap.add_argument('k',type=int);ap.add_argument('prefix',type=Path);a=ap.parse_args()
got=profile(a.p,Path(str(a.prefix)+'.cells'),ROOT,a.k,a.ag);expected=json.loads(Path(str(a.prefix)+'.json').read_text())
need(set(expected)==set(got)|{'witness_bytes','witness_sha256'},'exact stored profile schema')
for key,value in got.items():need(expected[key]==value,'stored integral field '+key)
need(expected['witness_bytes']==8*got['witnesses'] and len(expected['witness_sha256'])==64,'witness metadata schema only')
print('PASS stored integral',got['q'],a.k,got['upper_decimal'],flush=True)
