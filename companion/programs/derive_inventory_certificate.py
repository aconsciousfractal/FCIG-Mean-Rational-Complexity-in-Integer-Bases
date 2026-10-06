"""Derive the two-method inventory certificate from complete pinned metadata.

No height, CDF, logarithmic sign, inventory producer or scientific oracle
is imported. Summary arithmetic is recomputed from the five full lists.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path

from archive import entry, read_object

FIELDS = ('all_rows', 'retained_rows', 'optional_r27_qualifying_retained_rows',
          'supplemental_ratio_qualifying_retained_rows', 'active_rows')
COUNTS = dict(zip(FIELDS, (110922, 19239, 5361, 7857, 11382)))
BRANCHES = ('p2', 'p3', 'p_ge5')
REASONS = ('closed_radical', 'common_valuation_factor', 'squarefree', 'branch_cutoff')
SCHEMA = 'p72-r36-source-only-frontier-v1'
STATUS = 'R36_TWO_METHOD_DOMAIN_PASS'
ROW_KEYS = {'M', 'p', 'q', 'a', 'b', 'branch', 'cutoff', 'reasons'}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def strict_load(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            need(key not in result, 'duplicate inventory JSON key')
            result[key] = value
        return result
    def constant(value):
        raise ValueError('nonfinite inventory JSON constant: ' + value)
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)


def pinned(path, pin):
    need(type(pin) is str and len(pin) == 64 and all(c in '0123456789abcdef' for c in pin),
         'exact source artifact digest')
    raw = Path(path).read_bytes()
    need(sha(raw) == pin, 'source artifact byte pin')
    return strict_load(raw)


@lru_cache(maxsize=None)
def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def profile(rows):
    return dict(count=len(rows), mass=sum(r['M'] for r in rows),
                min_M=min((r['M'] for r in rows), default=None),
                max_M=max((r['M'] for r in rows), default=None),
                radical_count=len({(r['p'], r['q']) for r in rows}))


def ratio_qualifies(row, optional=False):
    coefficient = 1500 if optional else {'p2':225, 'p3':1300, 'p_ge5':400}[row['branch']]
    return row['b'] == 1 and row['M'] >= 19600 and row['q'] >= coefficient * row['p']**row['a']


def validate_rows(rows, field):
    need(type(rows) is list and len(rows) == COUNTS[field], 'complete inventory list count: ' + field)
    previous = 0
    for row in rows:
        need(type(row) is dict and set(row) == ROW_KEYS, 'complete inventory row fields')
        need(all(type(row[key]) is int for key in ('M','p','q','a','b','cutoff')),
             'exact integer inventory profile')
        M,p,q,a,b = (row[key] for key in ('M','p','q','a','b'))
        need(previous < M < 360000 and 2 <= p < q <= M and 1 <= a <= 18 and 1 <= b <= 18,
             'unique strictly sorted bounded inventory profiles')
        need(prime(p) and prime(q) and M == p**a * q**b, 'exact two-prime inventory factorization')
        branch,cutoff = ('p2',147456) if p == 2 else ('p3',262144) if p == 3 else ('p_ge5',360000)
        reasons = []
        if p*q in (6,10,14,15):reasons.append('closed_radical')
        if gcd(a,b) >= 2:reasons.append('common_valuation_factor')
        if a == b == 1:reasons.append('squarefree')
        if M >= cutoff:reasons.append('branch_cutoff')
        need(type(row['branch']) is str and row['branch'] == branch and row['cutoff'] == cutoff
             and type(row['reasons']) is list and row['reasons'] == reasons,
             'original inventory classification metadata')
        previous = M


def validate_inventory(document, algorithm):
    need(type(document) is dict and document.get('schema') == SCHEMA
         and document.get('status') == 'SOURCE_ONLY' and document.get('algorithm') == algorithm,
         'complete producer status/schema/algorithm')
    need(document.get('reason_order') == list(REASONS), 'inventory reason order')
    domain = document.get('domain')
    expected_domain = ({'M_lower_inclusive':2,'M_upper_exclusive':360000,'omega':2,
                        'prime_order':'p<q','valuations':'a,b>=1'}
                       if algorithm == 'ordered_primes_and_powers'
                       else {'M_min_inclusive':2,'M_max_exclusive':360000})
    need(type(domain) is dict and domain == expected_domain
         and all(type(v) is int for v in domain.values() if not isinstance(v,str)),
         'complete fixed inventory domain')
    for field in FIELDS:
        need(field in document, 'missing full inventory list: ' + field)
        validate_rows(document[field], field)
    all_rows = document['all_rows']
    retained = [r for r in all_rows if not r['reasons']]
    optional = [r for r in retained if ratio_qualifies(r, optional=True)]
    supplemental = [r for r in retained if ratio_qualifies(r)]
    active = [r for r in retained if not ratio_qualifies(r)]
    need([retained, optional, supplemental, active] == [document[field] for field in FIELDS[1:]],
         'complete inventory subset/partition identities')


def derive(generator, factor, generator_pin, factor_pin, template):
    validate_inventory(generator, 'ordered_primes_and_powers')
    validate_inventory(factor, 'integer-traversal-smallest-prime-factor-sieve')
    for field in FIELDS:
        need(generator[field] == factor[field], 'two-method full-array disagreement: ' + field)
    need(type(template) is dict and template.get('status') == STATUS
         and template.get('normalized_full_list_fields') == list(FIELDS)
         and type(template.get('supplemental_import')) is str and template['supplemental_import'],
         'pinned baseline attribution template')
    active = generator['active_rows']
    return dict(active_by_branch=[dict(branch=branch, **profile([r for r in active if r['branch'] == branch]))
                                  for branch in BRANCHES],
                factor_output='inventory_factor_normal1.json', factor_sha256=factor_pin,
                general_N1_closed=False,
                generator_output='inventory_generate_normal1.json', generator_sha256=generator_pin,
                normalized_full_list_fields=list(FIELDS),
                scientific_heights_CDFs_or_signs_evaluated=False, scientific_means_evaluated=False,
                status=STATUS, summaries={field:profile(generator[field]) for field in FIELDS},
                supplemental_import=template['supplemental_import'])


def fresh_output(path, document):
    target = Path(path).absolute()
    need(not target.exists() and not target.is_symlink(), 'fresh inventory certificate output required')
    need(target.parent.is_dir() and all(not p.is_symlink() for p in target.parents),
         'existing nonsymlink output parent chain required')
    raw = (json.dumps(document, sort_keys=True, indent=2) + '\n').encode()
    with target.open('xb') as stream:
        stream.write(raw)
    return dict(status=STATUS, certificate_sha256=sha(raw), certificate_bytes=len(raw),
                generator_sha256=document['generator_sha256'], factor_sha256=document['factor_sha256'],
                complete_full_list_fields=list(FIELDS))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('generator', type=Path)
    parser.add_argument('generator_sha256')
    parser.add_argument('factor', type=Path)
    parser.add_argument('factor_sha256')
    parser.add_argument('output', type=Path)
    parser.add_argument('--baseline-template', type=Path)
    parser.add_argument('--baseline-template-sha256')
    args = parser.parse_args()
    need(bool(args.baseline_template) == bool(args.baseline_template_sha256),
         'baseline template path and digest supplied together')
    need(not args.output.exists() and not args.output.is_symlink(), 'fresh inventory certificate output required')
    generator = pinned(args.generator, args.generator_sha256)
    factor = pinned(args.factor, args.factor_sha256)
    if args.baseline_template is None:
        template = strict_load(read_object(entry(36, 'FRONTIER_CERTIFICATE.json', 'results')))
    else:
        template = pinned(args.baseline_template, args.baseline_template_sha256)
    result = derive(generator, factor, args.generator_sha256, args.factor_sha256, template)
    print(json.dumps(fresh_output(args.output, result), sort_keys=True, separators=(',', ':')))


if __name__ == '__main__':main()
