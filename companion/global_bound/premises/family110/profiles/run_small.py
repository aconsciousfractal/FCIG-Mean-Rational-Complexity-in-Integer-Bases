#!/usr/bin/env python3
"""Sequentially generate the selected complete small-prime prefix."""
from pathlib import Path
import subprocess,sys,json
ROOT=Path(__file__).resolve().parent
for z in json.loads((ROOT/'profiles_plan.json').read_text()):
 if z['p']==1:continue
 stem=z['stem']
 with (ROOT/'logs'/f'{stem}_generation.log').open('w') as f:
  subprocess.run([sys.executable,'-B',ROOT/'generate.py','--stem',stem,'--threads','4'],stdout=f,stderr=subprocess.STDOUT,check=True)
 print('Completed',stem,flush=True)
