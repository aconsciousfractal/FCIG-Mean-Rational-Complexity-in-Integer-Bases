"""Reconstruct the pointwise semiprime histograms and consumed log ledgers.

Default fixed acceptance sample: M=6,10,77,221. --all selects all 4596.
--output writes a fresh streamed artifact and prints only its compact identity.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
from archive import data, entry
import sixth_length_math as core

FULL_COUNT = 4596
SAMPLE = (6, 10, 77, 221)
STATUS = 'PORTABLE_E032_RECONSTRUCTION_PASS'
_uncached_log_ticks = core.log_ticks


@lru_cache(maxsize=None, typed=True)
def _cached_log_ticks(n):
    return _uncached_log_ticks(n)


def log_ticks(n):
    """Memoize only exact positive integer inputs; preserve the original recipe."""
    if type(n) is not int or n < 1:
        raise ValueError('log input must be a positive exact integer')
    return _cached_log_ticks(n)


# The certificate and the exported ledger must consume the same exact recipe.
core.log_ticks = log_ticks


def complete_domain(bridge, domain, domain_sha256):
    """Compare complete ordered obligations with the distinct pinned domain."""
    if bridge.get('domain_sha256') != domain_sha256:
        raise ValueError('semiprime domain byte binding')
    profiles = []
    for payload in (bridge, domain):
        rows = payload.get('rows')
        if type(rows) is not list or len(rows) != FULL_COUNT:
            raise ValueError('complete semiprime domain count')
        current = []
        for row in rows:
            if type(row) is not dict or not all(key in row for key in ('M', 'p', 'r')):
                raise ValueError('semiprime domain profile fields')
            profile = tuple(row[key] for key in ('M', 'p', 'r'))
            M, p, q = profile
            if (not all(type(n) is int for n in profile) or not 2 <= p < q
                    or M != p*q or not (q < 6613 if p == 2 else M < 19600)):
                raise ValueError('exact semiprime domain profile')
            current.append(profile)
        if (len(set(current)) != len(current)
                or [M for M,p,q in current] != sorted({M for M,p,q in current})):
            raise ValueError('unique ordered semiprime domain')
        profiles.append(current)
    if profiles[0] != profiles[1]:
        raise ValueError('complete ordered semiprime domain identity')
    return bridge['rows']


def compare_record(record, expected):
    for key in ('upper_histogram_sha256', 'upper_centered_log_numerator',
                'upper_centered_log_denominator'):
        if type(record[key]) is not type(expected[key]) or record[key] != expected[key]:
            raise ValueError('E032 reconstruction mismatch '+str(record['M'])+'/'+key)


def reconstruct(M,p,q):
    T=isqrt(2*M)+1
    if q<=T:
        exact,h,evidence=core.central(M,M)
        upper=[v or min(q,max(1,min(s,M-s))) for s,v in enumerate(exact)]
    else:
        upper,h,exact,evidence=core.envelope(M,M)
        upper=[min(v,q) for v in upper]
    hist=Counter(upper);cert=core.certify(hist,M)
    histogram=cert['histogram'];digest=hashlib.sha256(json.dumps(histogram,separators=(',',':')).encode()).hexdigest()
    ledger=[[n,*core.log_ticks(n)] for n in sorted({M}|set(hist))]
    return dict(M=M,p=p,r=q,upper_histogram=histogram,log_tick_ledger=ledger,
                upper_histogram_sha256=digest,upper_centered_log_numerator=cert['upper_centered_log_numerator'],
                upper_centered_log_denominator=cert['upper_centered_log_denominator'])


def reconstructed_rows(selected, progress=False):
    total = len(selected)
    for index, expected in enumerate(selected, 1):
        row = reconstruct(expected['M'], expected['p'], expected['r'])
        compare_record(row, expected)
        yield row
        if progress and (index % 250 == 0 or index == total):
            print('PROGRESS '+str(index)+'/'+str(total), file=sys.stderr, flush=True)


def result(rows, scope):
    return dict(status=STATUS, scope=scope, rows=rows,
                historical_job_replayed=False, new_portable_reconstruction=True)


def export_file(path, rows, scope, total):
    """Exclusive streamed output; a failed run leaves clearly incomplete JSON."""
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError('reconstruction output must be fresh')
    if not path.parent.is_dir() or path.parent.is_symlink():
        raise ValueError('reconstruction output requires an existing nonsymlink parent')
    digest = hashlib.sha256()
    output_bytes = 0
    count = 0
    with path.open('xb') as stream:
        def write(text):
            nonlocal output_bytes
            payload = text.encode('utf-8')
            stream.write(payload)
            digest.update(payload)
            output_bytes += len(payload)
        write('{"historical_job_replayed":false,"new_portable_reconstruction":true,"rows":[')
        for row in rows:
            if count >= total:
                raise ValueError('reconstruction output has excess rows')
            if count:
                write(',')
            write(json.dumps(row, sort_keys=True, separators=(',', ':')))
            count += 1
            if count % 250 == 0:
                stream.flush()
        if count != total:
            raise ValueError('reconstruction output has missing rows')
        write('],"scope":'+json.dumps(scope)+',"status":'+json.dumps(STATUS)+'}\n')
    return dict(status=STATUS, scope=scope, reconstructed_rows=count,
                output_sha256=digest.hexdigest(), output_bytes=output_bytes,
                historical_job_replayed=False, new_portable_reconstruction=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--output', type=Path, help='fresh JSON artifact; stdout contains its compact identity')
    a=ap.parse_args()
    bridge=data(28,'SEMIPRIME_BRIDGE.json')
    domain=data(28,'SEMIPRIME_DOMAIN.json')
    old=complete_domain(bridge, domain, entry(28,'SEMIPRIME_DOMAIN.json')['sha256'])
    selected=old if a.all else [r for r in old if r['M'] in SAMPLE]
    if not a.all and tuple(r['M'] for r in selected) != SAMPLE:
        raise ValueError('fixed reconstruction sample identity')
    scope='full' if a.all else 'fixed four-modulus acceptance sample'
    rows=reconstructed_rows(selected, progress=a.all or a.output is not None)
    payload=export_file(a.output, rows, scope, len(selected)) if a.output is not None else result(list(rows), scope)
    print(json.dumps(payload,sort_keys=True,separators=(',',':')))


if __name__=='__main__':main()
