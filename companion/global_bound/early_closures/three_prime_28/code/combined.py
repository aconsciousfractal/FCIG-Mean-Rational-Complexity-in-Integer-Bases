from fractions import Fraction as F
from math import isqrt
import arithmetic as ar
import excision
from finite import primes_to_113

def certificate(new,old):
    ps=primes_to_113();ar.need(new['primes']==ps==old['N_ge_5_unexcluded_primes_by_this_filter'],'complete finite/large-prime interface')
    cuts=[]
    for r in old['rows']:
        k=r['k'];count=r['remaining_prime_count_for_N_ge_k']-len(ps);ar.need(count>=0,'nonnegative residual count')
        if k<=4:
            small=[n for n in range(127,r['p_threshold']) if all(n%d for d in range(2,isqrt(n)+1))] if k>=2 else None
            if small is not None:ar.need(len(small)==count,'explicit residual prime range')
            cuts.append({'N':k,'p_min':127,'p_strict_upper':r['p_threshold'],'prime_count':count,'largest_prime':r['largest_unexcluded_prime_for_N_ge_k'],'min_base_2_3_valuation_at_most':4//k,'primes':small})
        else:ar.need(count==0,'all N>=5 closed')
    upper=F(393,1000)+excision.psi_bounds(127)[1]
    bound=10**11;lo=F(ar.logs(bound)[0],ar.LS);margin=F(2,85)*lo-upper
    ar.need(margin>0,'explicit residual base cutoff')
    return {'status':'PASS','all_radical_6p_pairs_N_ge_5_excluded':True,'whole_6p_family_closed':False,'remaining_necessary_conditions':{'p_min':127,'N_min':1,'N_max':4,'B_strict_upper':bound,'min_v2_v3_of_B_to_N_at_most':4},'base_cutoff_mu_upper':str(upper),'base_cutoff_margin':str(margin),'residual_rows':cuts,'source_of_large_prime_filter':'replayed mod6 profiles plus elementary prime-excision theorem','global_Cstar_identified':False,'no_actual_competitor_claimed':True}
