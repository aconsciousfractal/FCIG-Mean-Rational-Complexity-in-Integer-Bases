"""Independent normalization and full-replay binding of the6510 theorem."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import argparse
import json
from common import need,read,sha,canonical
from exact import log_lo,out_decimal
from candidate210 import certificate as candidate
from custody import replay

ROOT=Path(__file__).resolve().parents[1]
Q=6510
PRIMES=[2,3,5,7,31]


def certificate(root=ROOT):
    root=Path(root)
    p=root/'profiles'
    contract=read(p/'profile_contract.json');plan=read(p/'profiles_plan.json')
    summary=read(p/'results/profile_certificate.json');profile=read(p/'results/q6510_K6510.json')
    need(prod(PRIMES)==Q and contract['core']==Q and contract['core_primes']==PRIMES,
         'actual complete five-prime support6510')
    need(plan==[dict(p=1,multiplier=1,ag=8,stem='q6510_K6510')],
         'exact one-profile plan without additional index')
    need(summary['status']=='PASS' and summary['profiles']==1
         and summary['arbitrary_prime_valuations'] and summary['whole_radical_directly_closed_here']
         and summary['all_N_positive_if_direct_gate'], 'whole support scope from transfer')
    need(profile['q']==profile['transfer_index']==Q and profile['alpha']=='3/2'
         and profile['labels']==18432 and profile['cells']==2688 and profile['angular_nodes']==17,
         'complete projective transfer geometry')
    records=prod(p+1 for p in PRIMES)*2688
    corners=records*4*17
    need(profile['witnesses']==summary['geometric_witnesses']==records==49545216
         and summary['geometric_corner_evaluations']==corners==3369074688,
         'complete distinct records and continuous corner checks')
    A=F(profile['total_upper'])
    need(A==F(profile['bulk_upper'])+F(profile['cusp_upper'])
         ==F(24921625818115921811023309563933918363563,
              154742504910672534362390528000000000000000),
         'new mean is assembled from signed bulk and cusp')
    gamma=A/log_lo(Q)
    need(0<A<F(2,85)*log_lo(Q) and gamma<F(2,85),
         'actual sufficient gate; no unnecessary A<4/25 hypothesis')
    cap=summary['marked_caps'][0]
    need(len(summary['marked_caps'])==1 and cap['q']==cap['required_divisor']==Q
         and F(cap['exact_upper'])==A and F(cap['normalized_upper'])==gamma
         and cap['strict_gate_2_over85'] and F(profile['normalized_upper'])==gamma,
         'profile certificate agrees with independent lnB normalization')
    c=candidate()
    need(c==read(root/'results/candidate210.json') and gamma<F(2,85)<F(c['c0_lower']),
         'exact candidate reconstructed from all210 minima')
    custody=replay(root,'profiles')
    e=custody['expectation']
    need((e['profile_count'],e['geometric_records'],e['corner_evaluations'],e['complete_profile_comparisons'])==(1,49545216,3369074688,10),
         'exact new6510 full replay counts')
    names=['code/assemble_fixed.py','code/custody.py','code/run_final_replay.py',
           'code/common.py','code/exact.py','PROOF_6510.md',
           'profiles/PROOF.md','profiles/results/profile_certificate.json',
           'profiles/results/q6510_K6510.json','results/candidate210.json']
    return dict(status='PASS',radical=Q,primes=PRIMES,
                theorem='rad(B)=6510 implies Delta(B,N)<=Gamma6510<2/85<c0',
                all_positive_valuations=True,all_N_ge_1=True,
                normalization='Delta(B,N)=mu(B^N)/ln(B)',
                mean_upper=str(A),mean_upper_decimal=out_decimal(A),
                Gamma6510=str(gamma),Gamma6510_upper_decimal=out_decimal(gamma),
                margin_to_2_over_85=str(F(2,85)-gamma),candidate=c,
                transfer=dict(mark=Q,actual_index=Q,seed='L=M/6510',
                              no_seed_coprimality_assumption=True,all_positive_multiples=True,
                              physical_identification='rad(M)=6510',normalizing_base_minimum=Q),
                geometry=dict(records=records,corner_evaluations=corners,
                              raw_bytes=profile['witness_bytes'],raw_sha256=profile['witness_sha256']),
                completed_full_replay_verified=True,completed_replay=custody,
                second_profile_required=False,finite_complement_required=False,
                exact_Cstar_equals_c0_proved=False,offset_sign_classification_claimed=False,
                inputs=[dict(path=n,sha256=sha(root/n)) for n in names])


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.write_bytes(canonical(certificate()))
