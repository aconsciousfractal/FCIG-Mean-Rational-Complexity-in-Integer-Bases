"""Reintegrate all fixed4290 profiles without assuming an exclusion passes."""
from pathlib import Path
from fractions import Fraction as F
import argparse, json
from exact import need, log_lo, out_decimal
from assemble import profile
from core import CORE, CORE_PRIMES, CORE_LABELS
ROOT = Path(__file__).resolve().parents[1]

def assemble(root=ROOT):
    plan = json.loads((root / 'profiles_plan.json').read_text())
    contract = json.loads((root / 'profile_contract.json').read_text())
    need(CORE == 4290 and CORE_PRIMES == (2, 3, 5, 11, 13) and (CORE_LABELS == 12096), 'fixed five-prime scientific core')
    need(contract['core'] == CORE and tuple(contract['core_primes']) == CORE_PRIMES and (contract['cdf_alpha'] == '3/2'), 'declared core and CDF')
    need(len({r['stem'] for r in plan}) == len(plan) and plan, 'unique nonempty profile plan')
    need(set((r['multiplier'] for r in plan)) <= set(contract['allowed_multipliers']) and all((r['p'] == 1 for r in plan)), 'fixed mark mixed indices')
    rows = []
    marked = []
    for z in plan:
        m, ag = (z['multiplier'], z['ag'])
        need(ag in (4, 8, 16) and z['stem'] == f'q4290_K{CORE * m}', 'audited angle and unambiguous stem')
        got = profile(1, root / 'results' / (z['stem'] + '.cells'), root, m, ag)
        stored = json.loads((root / 'results' / (z['stem'] + '.json')).read_text())
        need(set(stored) == set(got) | {'witness_bytes', 'witness_sha256'} and all((stored[k] == v for k, v in got.items())), 'exact complete integral schema and fields')
        need(stored['witness_bytes'] == 8 * got['witnesses'] and len(stored['witness_sha256']) == 64, 'record count and digest schema')
        cap = F(got['total_upper'])
        norm = max(F(0), cap) / log_lo(CORE)
        rows.append(got)
        marked.append({'q': CORE, 'required_divisor': CORE * m, 'exact_upper': str(cap), 'upper_decimal': out_decimal(cap), 'normalized_upper': str(norm), 'normalized_upper_decimal': out_decimal(norm), 'strict_gate_2_over85': norm < F(2, 85), 'strict_margin': str(F(2, 85) - norm), 'scope': 'all positive M divisible by required_divisor, extra primes allowed; marked mean nu4290, equality with physical mean under rad(M)=4290'})
    direct = any((z['required_divisor'] == CORE and z['strict_gate_2_over85'] for z in marked))
    return {'status': 'PASS', 'cdf_alpha': '3/2', 'marked_caps': marked, 'profiles': len(rows), 'geometric_witnesses': sum((r['witnesses'] for r in rows)), 'geometric_corner_evaluations': sum((4 * r['angular_nodes'] * r['witnesses'] for r in rows)), 'arbitrary_prime_valuations': True, 'whole_radical_directly_closed_here': direct, 'all_N_positive_if_direct_gate': direct, 'geometric_witnesses_replayed_here': False, 'global_Cstar_identified': False}
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    z = assemble()
    s = json.dumps(z, indent=2) + '\n'
    p = ROOT / 'results/profile_certificate.json'
    if a.check:
        need(p.read_text() == s, 'exact profile certificate')
    else:
        p.write_text(s, newline='\n')
    print('PASS profiles assembly', z['profiles'], 'profiles')
