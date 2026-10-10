"""Exact assembly checker for the current marked-70 closure."""
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from assemble_closure import certificate,route_csv
from inventory70 import make
from exact import need
def verify():
 rows,meta=make();need(meta==json.loads((ROOT/'results/inventory_certificate.json').read_text()),'inventory regenerated exactly')
 z,rr=certificate(ROOT)
 need(z==json.loads((ROOT/'results/closure_certificate.json').read_text()),'closure exact regeneration')
 need(route_csv(rr).encode()==(ROOT/'results/sharp_routes.csv').read_bytes(),'sharp routes exact regeneration')
 return dict(status='PASS',checker_provenance='exact numerical assembly',finite_moduli=len(rows),gamma70=z['theorem_70p']['uniform_upper_decimal'],gamma770_tail=z['theorem_770ell']['sharp_upper_decimal'],finite_streams_replayed_here=False)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
