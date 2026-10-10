#!/usr/bin/env python3
"""Sequentially generate the remaining five small-prime profiles."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parent
for p in (19,23,29,31,37):
 stem=f'q{78*p}'
 with (ROOT/'logs'/f'{stem}_generation.log').open('w') as f:
  subprocess.run([sys.executable,'-B',ROOT/'generate.py','--stem',stem,'--threads','3'],stdout=f,stderr=subprocess.STDOUT,check=True)
 print('Completed',stem,flush=True)
