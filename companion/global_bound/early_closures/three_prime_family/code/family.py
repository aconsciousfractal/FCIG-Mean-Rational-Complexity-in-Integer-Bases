"""Exact interface to the PRIOR internal S113 theorem; not a replay of its proof."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json
from exact import need,log_lo,out_decimal


def run(root):
    root=Path(root)
    parent=json.loads((root/'imports/S113_closure_certificate.json').read_text())
    primes=[p for p in range(5,127) if all(p%d for d in range(2,isqrt(p)+1))]
    need(parent['primes']==primes and len(primes)==28,'S113 prime interface')
    need(parent['all_lengths_and_positive_valuations'] is True,'S113 quantifiers')
    need(parent['clean_rational_cap']=='117/5000','S113 bound interface')
    old=F(parent['uniform_delta_upper']);new=F(3,20)/log_lo(762)
    need(old<F(117,5000) and new<F(117,5000)<F(2,85),'common family cap')
    composites={str(n):next(d for d in range(2,isqrt(n)+1) if n%d==0) for n in range(114,127)}
    need(all(1<d<int(n) and int(n)%d==0 for n,d in composites.items()),'prime gap 113--127')
    imports={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted((root/'imports').glob('*')) if p.is_file()}
    return {'status':'PASS_WITH_DECLARED_INTERNAL_THEOREM_IMPORT',
            'S113_internal_theorem_reused':True,'S113_full_proof_replayed_this_task':False,
            'new_G127_depends_on_S113':False,'S113_primes':primes,
            'composite_witnesses_114_to_126':composites,
            'new_p_ge127_all_valuations_all_lengths_excluded':True,
            'new_delta_upper':str(new),'new_delta_upper_decimal':out_decimal(new),
            'S113_imported_delta_upper':str(old),'combined_uniform_cap':'117/5000',
            'combined_below_candidate_lower':'117/5000 < 2/85 < c0',
            'entire_6p_family_closed_with_S113':True,
            'global_Cstar_identified':False,'four_prime_radicals_covered':False,
            'import_file_sha256':imports,
            'warning':'These checks validate the theorem interface and the rational splice, not the imported geometric proof.'}

if __name__=='__main__':
    import sys
    obj=run(sys.argv[1]);Path(sys.argv[2]).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    print(json.dumps(obj,indent=2),flush=True)
