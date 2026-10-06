"""Check every reconstructed semiprime histogram and ledger with distinct log bounds."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from archive import data, entry
from verify_all_central import strict_load
from verify_retained import hist_check, log32

S = 2**64
COUNT = 4596
fine_log = lru_cache(maxsize=None, typed=True)(log32)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def verify(path, pin):
    fine_log.cache_clear()
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == pin, 'full reconstruction byte pin')
    payload = strict_load(raw)
    rows = payload['rows']
    old = data(28, 'SEMIPRIME_BRIDGE.json')
    domain_record = entry(28, 'SEMIPRIME_DOMAIN.json')
    need(old['domain_sha256'] == domain_record['sha256'], 'historical bridge/domain byte binding')
    domain = data(28, 'SEMIPRIME_DOMAIN.json')['rows']
    expected = [[r[k] for k in ('M','p','r')] for r in domain]
    actual = [[r[k] for k in ('M','p','r')] for r in rows]
    need(payload['scope'] == 'full' and payload['status'] == 'PORTABLE_E032_RECONSTRUCTION_PASS', 'completed full reconstruction')
    historic_profiles = [[r[k] for k in ('M','p','r')] for r in old['rows']]
    need(all(type(v) is int for collection in (actual,expected,historic_profiles) for profile in collection for v in profile), 'exact integer original profiles')
    need(len(rows) == len(domain) == len(old['rows']) == COUNT and actual == expected == historic_profiles
         and len({tuple(v) for v in actual}) == COUNT, 'all original ordered obligations')
    histograms = consumed_ticks = mass = 0
    for row, historic in zip(rows, old['rows']):
        M = row['M']
        hist = row['upper_histogram']
        hist_check(hist, M)
        digest = hashlib.sha256(json.dumps(hist,separators=(',',':')).encode()).hexdigest()
        need(digest == row['upper_histogram_sha256'] == historic['upper_histogram_sha256'], 'original histogram identity')
        ticks = row['log_tick_ledger']
        need(all(type(t) is list and len(t)==3 and all(type(v) is int for v in t) for t in ticks), 'exact log ledger')
        ledger = {n:(L,U) for n,L,U in ticks}
        need(len(ledger)==len(ticks) and list(ledger)==sorted(ledger)
             and set(ledger)=={M}|{h for h,c in hist}, 'complete sorted consumed log ledger')
        for n,(L,U) in ledger.items():
            lo,hi = fine_log(n)
            need(L <= lo <= hi <= U, 'distinct K32 enclosure contained in consumed ticks')
        q = 2*sum(c*ledger[h][1] for h,c in hist)-M*ledger[M][0]
        denominator = 2*M*S
        need(type(row['upper_centered_log_numerator']) is int and type(historic['upper_centered_log_numerator']) is int and q == row['upper_centered_log_numerator']
             == historic['upper_centered_log_numerator'] < 0, 'every reconstructed outward sign')
        need(type(row['upper_centered_log_denominator']) is int and type(historic['upper_centered_log_denominator']) is int
             and denominator == row['upper_centered_log_denominator'] == historic['upper_centered_log_denominator'], 'exact normalization')
        histograms += 1
        consumed_ticks += len(ticks)
        mass += M
    return {'status':'COMPLETE_SEMIPRIME_HISTOGRAM_AND_DISTINCT_LOG_PASS',
            'input_sha256':pin,'histograms_checked':histograms,'consumed_tick_entries':consumed_ticks,
            'distinct_log_integers':fine_log.cache_info().currsize,'upper_envelope_entries':mass,
            'every_outward_sum_negative':True,'pointwise_dominance_producer_and_proof_required':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('sha256')
    args=parser.parse_args()
    print(json.dumps(verify(args.input,args.sha256),sort_keys=True,separators=(',',':')))


if __name__=='__main__':main()
