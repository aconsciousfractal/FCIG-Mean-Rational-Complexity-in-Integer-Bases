"""Reconstruct search-space costs only: no large-state enumeration is claimed here."""
from pathlib import Path
import csv
from central import *
from inventory import power_tuples,primitive_base,write_json

def legacy(M,r,s,p):
    d=parameters(r,s);u=F(d['u']);A=isqrt(M);l=log_hi(A)
    return F(d['T0_tick'],CS)+phi_hi(5*u/p)+u*(1+ELL_HI)/p+(F(d['C0_tick'],CS)+F(d['C1_tick'],CS)*l)/A+(F(d['C2_tick'],CS)+F(d['C3_tick'],CS)*l)/(A*A)

def run(root,out):
    root=Path(root);out=Path(out);rows=list(csv.DictReader((root/'imports/frontier_442.csv').open()))
    naive_mass=normalized_mass=naive_n=normalized_n=0;old_tails=0
    for row in rows:
        if row['tail_from_M']=='not certified':continue
        prs=tuple(int(row[k]) for k in ('r','s','p'));Q=prod(prs);K=int(row['tail_from_M']);target=DELTA*log_lo(Q)
        need(legacy(K,*prs)<target,'legacy scalar cutoff recomputed')
        for M,es in power_tuples(prs,K):
            u=legacy(M,*prs);B,g,_=primitive_base(M,prs)
            if u>=target:naive_n+=1;naive_mass+=M
            if u>=DELTA*log_lo(B):normalized_n+=1;normalized_mass+=M
        old_tails+=1
    inv=__import__('json').loads((out/'inventory.json').read_text())
    result=dict(status='PASS',old_tail_supports=old_tails,legacy_radical_normalized_pending_moduli=naive_n,legacy_radical_normalized_state_mass=naive_mass,legacy_primitive_base_pending_moduli=normalized_n,legacy_primitive_base_state_mass=normalized_mass,new_supports=inv['selected_supports'],new_pending_moduli=inv['witness_moduli'],new_state_mass=inv['witness_state_mass'],large_legacy_states_enumerated=False)
    need(naive_mass==605095848801 and normalized_mass==3511413042,'cost audit expected identity')
    write_json(out/'cost_audit.json',result);return result
if __name__=='__main__':
    import sys,json
    print(json.dumps(run(Path(__file__).resolve().parents[1],sys.argv[1]),indent=2))
