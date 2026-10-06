"""Distinct exact central enumeration; no producer or interval-oracle import.

Every signed legal pair in the central square is visited, including pairs
with nontrivial content. Only central minima are computed; zero in the vector
means that the state's true height exceeds floor(sqrt(M)).
"""
import argparse
import ast
from array import array
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import sys

from archive import PACKAGE, entry, read_object

BATCHES = {
    29: ('RADICAL6_AGGREGATE.json', 'RADICAL6_DOMAIN.json', 143),
    30: ('RADICAL10_AGGREGATE.json', 'RADICAL10_DOMAIN.json', 97),
    31: ('PARAMETRIC_AGGREGATE.json', 'PARAMETRIC_DOMAIN.json', 132),
    32: ('LENGTH_FOUR_AGGREGATE.json', 'LENGTH_FOUR_DOMAIN.json', 3),
    33: ('LENGTH_THREE_AGGREGATE.json', 'LENGTH_THREE_DOMAIN.json', 35),
    35: ('LENGTH_TWO_AGGREGATE.json', 'LENGTH_TWO_DOMAIN.json', 99),
    39: ('BRANCH_AGGREGATE.json', 'BRANCH_EXECUTION_DOMAIN.json', 1087),
    40: ('BRANCH_AGGREGATE.json', 'BRANCH_DOMAIN.json', 1422),
    41: ('BRANCH_AGGREGATE.json', 'BRANCH_DOMAIN.json', 1901),
}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def exact(value):
    return type(value) is int


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def strict_load(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            need(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def constant(value):
        raise ValueError('nonfinite JSON constant: ' + value)
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)


def factorization(M):
    need(exact(M) and M >= 2, 'positive exact modulus >=2')
    factors = []
    value, prime = M, 2
    while prime * prime <= value:
        if value % prime == 0:
            exponent = 0
            while value % prime == 0:
                value //= prime
                exponent += 1
            factors.append((prime, exponent))
        prime += 1
    if value > 1:
        factors.append((value, 1))
    return factors


def proper_divisors(M):
    values = [1]
    for prime, exponent in factorization(M):
        values = [d * prime**j for d in values for j in range(exponent + 1)]
    return sorted(d for d in values if d < M)


def enumerate_central(M):
    """Return exact CDF arrays for every proper gcd stratum and every T<=k.

    No primitivity test, Mobius sum, collision formula, Farey interval or
    scientific source import is used. A nonprimitive legal pair is a valid
    witness and must not be confused with a separate random state.
    """
    need(exact(M) and 2 <= M < 2**32, 'exact uint32-sized modulus >=2')
    k = isqrt(M)
    heights = array('I', [0]) * M
    need(heights.itemsize == 4, 'uint32 array representation')
    heights[0] = 1
    denominator_units = pair_visits = signed_visits = improvements = 0
    for d in range(1, k + 1):
        if gcd(d, M) != 1:
            continue
        denominator_units += 1
        inverse = pow(d, -1, M)
        residue = 0
        for a in range(1, k + 1):
            residue = (residue + inverse) % M
            height = max(a, d)
            pair_visits += 1
            # At a self-complementary state both visits update the same atom.
            for state in (residue, (-residue) % M):
                signed_visits += 1
                if heights[state] == 0 or height < heights[state]:
                    heights[state] = height
                    improvements += 1
    divisors = proper_divisors(M)
    bins = {g: array('I', [0]) * (k + 1) for g in divisors}
    state_visits = covered = 0
    for state in range(1, M):
        state_visits += 1
        height = heights[state]
        if height:
            g = gcd(state, M)
            need(g <= k and 1 <= height <= k, 'covered nonzero central state')
            bins[g][height] += 1
            covered += 1
    cdf = array('I', [0]) * (k + 1)
    for T in range(1, k + 1):
        total = 1  # The zero state has height one, and no nonzero gcd layer.
        for g in divisors:
            bins[g][T] += bins[g][T - 1]
            total += bins[g][T]
        cdf[T] = total
    need(pair_visits == k * denominator_units and signed_visits == 2 * pair_visits,
         'complete legal signed-square visit accounting')
    need(cdf[k] == covered + 1 and pair_visits <= k*k <= M,
         'central coverage and exact pair budget')
    vector_bytes = heights.tobytes()
    if sys.byteorder != 'little':
        copy = array('I', heights)
        copy.byteswap()
        vector_bytes = copy.tobytes()
    result = {
        'M': M, 'k': k, 'zero_count': 1, 'cdf': list(cdf),
        'strata': [[g, list(bins[g])] for g in divisors],
        'central_vector_sha256': digest(vector_bytes),
        'central_vector_encoding': 'uint32-little-endian; zero means height>k',
        'counters': {'denominator_units': denominator_units,
                     'positive_pair_visits': pair_visits,
                     'signed_witness_visits': signed_visits,
                     'height_improvements': improvements,
                     'nonzero_state_scan': state_visits,
                     'central_thresholds': k,
                     'proper_stratum_thresholds': k * len(divisors),
                     'forced_zero_stratum_thresholds': k * sum(g > k for g in divisors)},
        'true_far_minima_computed': False,
    }
    return result


def row_profile(batch, row):
    p = 2 if batch in (29, 30) else row.get('p')
    q = (3 if batch == 29 else 5) if batch in (29, 30) else row.get('r')
    values = [row.get('M'), p, q, row.get('a'), row.get('b')]
    need(all(exact(v) for v in values), 'exact complete profile')
    M, p, q, a, b = values
    need(1 < p < q and a >= 1 and b >= 1 and M == p**a * q**b,
         'original ordered two-prime profile')
    need(factorization(M) == [(p, a), (q, b)], 'complete two-prime factorization')
    return values


def bind_domain(batch, payload, domain):
    need(batch in BATCHES, 'supported aggregate batch')
    count = BATCHES[batch][2]
    profiles = domain.get('profiles' if batch >= 39 else 'rows')
    rows = payload.get('rows')
    need(type(profiles) is list and type(rows) is list
         and len(profiles) == len(rows) == count, 'complete batch domain cardinality')
    need(all(type(p) is list and len(p) == 5 and all(exact(v) for v in p)
             for p in profiles), 'exact domain profiles')
    need(len({tuple(p) for p in profiles}) == count, 'duplicate domain profile')
    need([row_profile(batch, row) for row in rows] == profiles,
         'ordered original domain/aggregate identity')
    return profiles


def bind_validation(batch, payload, source_bytes=None):
    rows = payload.get('validation', [])
    need(type(rows) is list, 'validation row list')
    if batch >= 39:
        need(not rows, 'final branch has no aggregate validation extras')
        return []
    if source_bytes is None:
        source_bytes = read_object(entry(batch, 'verify_aggregate.py', 'results'))
    tree = ast.parse(source_bytes)
    values = [node.value for node in tree.body if isinstance(node, ast.Assign)
              and any(isinstance(t, ast.Name) and t.id == 'VALIDATION' for t in node.targets)]
    need(len(values) == 1, 'one literal producer VALIDATION declaration')
    constants = ast.literal_eval(values[0])
    need(type(constants) in (tuple, list), 'literal fixed producer validation domain')
    expected = [[M, 2, 3 if batch == 29 else 5] for M in constants] if batch in (29,30) else [list(v) for v in constants]
    need(all(len(v) == 3 and all(exact(x) for x in v) for v in expected)
         and len({tuple(v) for v in expected}) == len(expected), 'unique exact validation obligations')
    profiles = []
    for row in rows:
        need(type(row) is dict and exact(row.get('M')), 'exact validation modulus')
        fs = factorization(row['M'])
        need(len(fs) == 2, 'two-prime validation model')
        (p,a),(q,b) = fs
        profile = [row['M'],p,q,a,b]
        for key, value in zip(('M','p','r','a','b'), profile):
            if key in row:
                need(exact(row[key]) and row[key] == value, 'validation original profile field')
        profiles.append(profile)
    need([p[:3] for p in profiles] == expected, 'complete ordered fixed validation domain')
    return profiles


def compare_row(profile, row, actual, producer_arrays):
    M, p, q, a, b = profile
    k = isqrt(M)
    need(actual['M'] == M and actual['k'] == k and row.get('central_k') == k,
         'original central modulus/cutoff')
    all_strata = {g: values for g, values in actual['strata']}
    headers = [g for g in all_strata if g <= k]
    need(type(producer_arrays) is list and len(producer_arrays) == len(headers),
         'complete producer central stratum list')
    need(all(type(pair) is list and len(pair) == 2 and exact(pair[0])
             and type(pair[1]) is list for pair in producer_arrays), 'producer array shape')
    need([g for g, values in producer_arrays] == headers,
         'ordered producer stratum headers')
    for g, values in producer_arrays:
        need(len(values) == k + 1 and all(exact(v) for v in values),
             'complete exact producer prefix array')
        need(values == all_strata[g], 'every central threshold/stratum: M=%d g=%d' % (M, g))
    need(all(not any(values) for g, values in all_strata.items() if g > k),
         'all omitted producer strata forced zero')
    need(actual['cdf'][k] == row.get('exact_central_states'), 'retained central total')
    if 'central_stratum_totals' in row:
        need(row['central_stratum_totals'] == [[g, v[k]] for g, v in all_strata.items()],
             'all retained proper stratum totals')
    histogram = row.get('upper_histogram')
    need(type(histogram) is list and all(type(pair) is list and len(pair) == 2
         and all(exact(v) for v in pair) and pair[0] >= 1 and pair[1] > 0
         for pair in histogram), 'retained histogram exact shape')
    need([h for h, c in histogram] == sorted({h for h, c in histogram})
         and sum(c for h, c in histogram) == M, 'retained histogram order/mass')
    central_bins = {h: c for h, c in histogram if h <= k}
    prefix = 0
    for T in range(1, k + 1):
        prefix += central_bins.get(T, 0)
        need(prefix == actual['cdf'][T], 'every retained global central CDF threshold')
    return {'profile': profile, 'central_vector_sha256': actual['central_vector_sha256'],
            'central_arrays_sha256': digest(canonical(actual['strata'])),
            'central_cdf_sha256': digest(canonical(actual['cdf'])),
            'exact_central_states': actual['cdf'][k],
            'all_proper_strata_checked': len(all_strata),
            'producer_array_strata_checked': len(headers),
            'all_thresholds_checked': k, 'counters': actual['counters'],
            'status': 'ALL_CENTRAL_THRESHOLDS_AND_STRATA_PASS'}


def load_batch(batch, certificate=None, certificate_pin=None):
    aggregate_name, domain_name, count = BATCHES[batch]
    aggregate_record = entry(batch, aggregate_name, 'results')
    domain_record = entry(batch, domain_name, 'results')
    if certificate is None:
        payload = strict_load(read_object(aggregate_record))
        actual_pin = aggregate_record['sha256']
    else:
        raw = Path(certificate).read_bytes()
        need(digest(raw) == certificate_pin, 'fresh aggregate certificate byte pin')
        payload = strict_load(raw)
        actual_pin = digest(raw)
    domain = strict_load(read_object(domain_record))
    need(payload.get('domain_sha256') == domain_record['sha256'], 'aggregate/domain byte pin')
    profiles = bind_domain(batch, payload, domain)
    pins = {'aggregate_sha256': aggregate_record['sha256'],
            'consumed_certificate_sha256': actual_pin,
            'domain_sha256': domain_record['sha256'],
            'provenance_sha256': digest((PACKAGE / 'provenance.json').read_bytes()),
            'checker_sha256': digest(Path(__file__).read_bytes())}
    return payload, profiles, pins


def load_producer_arrays(batch, path, pin, profiles):
    raw = Path(path).read_bytes()
    need(digest(raw) == pin, 'fresh producer-array export byte pin')
    data = strict_load(raw)
    need(data.get('status') == 'PRODUCER_CENTRAL_ARRAY_CAPTURE_COMPLETE'
         and data.get('batch') == batch, 'producer-array capture status/batch')
    record = entry(batch, 'verify_aggregate.py', 'results')
    read_object(record)
    source_map = strict_load((PACKAGE / 'portable_source_map.json').read_bytes())
    adapters = [r for r in source_map if r.get('locator') == record['locator']]
    need(len(adapters) == 1 and adapters[0]['public_source_sha256'] == record['sha256'],
         'captured producer source/adapter identity')
    adapter_pin = adapters[0]['portable_sha256']
    need(digest((PACKAGE / 'src' / adapter_pin).read_bytes()) == adapter_pin
         and data.get('producer_adapter_sha256') == adapter_pin, 'captured producer adapter bytes')
    rows = data.get('rows')
    need(type(rows) is list, 'producer capture row list')
    by_key = {}
    for row in rows:
        need(type(row) is dict and all(exact(row.get(key)) for key in ('M','p','r','central_k')),
             'exact producer capture profile')
        key = (row['M'],row['p'],row['r'])
        need(key not in by_key and row['central_k'] == isqrt(row['M']),
             'unique producer capture/cutoff')
        by_key[key] = row
    # profiles is the deduplicated union of primary and consumed validation
    # obligations. Unrelated extras cannot count as successful coverage.
    need(set(by_key) == {tuple(p[:3]) for p in profiles},
         'all and only required primary/validation producer capture profiles')
    return [by_key[tuple(p[:3])] for p in profiles], digest(raw)


def fresh_output(path):
    target = Path(path).absolute()
    need(not target.is_symlink() and not target.exists(), 'output must be fresh')
    need(target.parent.is_dir(), 'output parent must already exist')
    for parent in target.parents:
        need(not parent.is_symlink(), 'symlink output parent')
    return target.open('x', encoding='utf-8', newline='\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', required=True, type=int, choices=sorted(BATCHES))
    parser.add_argument('--start', type=int, default=0, help='zero-based inclusive index')
    parser.add_argument('--stop', type=int, help='exclusive index; default complete batch')
    parser.add_argument('--producer-arrays')
    parser.add_argument('--producer-arrays-sha256')
    parser.add_argument('--certificate', help='fresh source-produced aggregate JSON')
    parser.add_argument('--certificate-sha256')
    parser.add_argument('--output', required=True, help='fresh JSON file; existing parent required')
    args = parser.parse_args()
    need(bool(args.producer_arrays) == bool(args.producer_arrays_sha256),
         'producer-array path and digest supplied together')
    need(bool(args.certificate) == bool(args.certificate_sha256),
         'certificate path and digest supplied together')
    payload, profiles, pins = load_batch(args.batch, args.certificate, args.certificate_sha256)
    validations = bind_validation(args.batch, payload)
    pins['validation_domain_sha256'] = digest(canonical(validations))
    capture_profiles = list(profiles)
    seen_profiles = {tuple(p) for p in profiles}
    for profile in validations:
        if tuple(profile) not in seen_profiles:
            capture_profiles.append(profile)
            seen_profiles.add(tuple(profile))
    stop = len(profiles) if args.stop is None else args.stop
    need(0 <= args.start < stop <= len(profiles), 'nonempty bounded chunk indices')
    external = None
    if args.producer_arrays:
        need(args.batch < 39, 'final branches already retain central arrays in their certificates')
        external, export_pin = load_producer_arrays(args.batch, args.producer_arrays,
            args.producer_arrays_sha256, capture_profiles)
        external = {(r['M'],r['p'],r['r']): r for r in external}
        pins['producer_arrays_sha256'] = export_pin
    else:
        need(args.batch >= 39, 'earlier batch needs freshly exported source-produced central arrays')
        pins['producer_arrays_sha256'] = None
    report = {'schema': 'P72_ALL_CENTRAL_CHUNK_V1', 'batch': args.batch,
              'start': args.start, 'stop': stop, 'complete_batch_rows': len(profiles),
              'pins': pins, 'central_scope': 'every T=1..floor(sqrt(M)); all proper gcd strata; zero once',
              'algorithm': 'joint nonprimitive-inclusive signed legal-pair enumeration',
              'producer_imported': False, 'interval_oracle_imported': False,
              'true_far_minima_computed': False, 'rows': [],
              'validation_profiles': validations, 'validation_rows': []}
    cached_controls = {}
    validation_keys = {tuple(p[:3]) for p in validations}
    with fresh_output(args.output) as output:
        try:
            for index in range(args.start, stop):
                row = payload['rows'][index]
                arrays = external[tuple(profiles[index][:3])]['central_strata'] if external is not None else row.get('central_strata')
                actual = enumerate_central(profiles[index][0])
                if tuple(profiles[index][:3]) in validation_keys:
                    cached_controls[tuple(profiles[index][:3])] = actual
                checked = compare_row(profiles[index], row, actual, arrays)
                checked['domain_index'] = index
                report['rows'].append(checked)
                print('CENTRAL %d %d/%d M=%d' % (args.batch, index+1, stop, profiles[index][0]),
                      file=sys.stderr, flush=True)
            for index, profile in enumerate(validations):
                key = tuple(profile[:3])
                actual = cached_controls.get(key)
                if actual is None:
                    actual = enumerate_central(profile[0])
                arrays = external[key]['central_strata']
                checked = compare_row(profile, payload['validation'][index], actual, arrays)
                checked['validation_index'] = index
                report['validation_rows'].append(checked)
            report['status'] = 'COMPLETE_REQUESTED_CHUNK_PASS'
            report['coverage_indices'] = list(range(args.start, stop))
            report['complete_batch_verified'] = args.start == 0 and stop == len(profiles)
        except Exception as error:
            report['status'] = 'FAILED_PARTIAL_CHUNK'
            report['complete_batch_verified'] = False
            report['error'] = {'type': type(error).__name__, 'message': str(error)}
            output.write(canonical(report).decode() + '\n')
            raise
        output.write(canonical(report).decode() + '\n')
    print(json.dumps({k: report[k] for k in ('status','batch','start','stop','complete_batch_verified')},
                     sort_keys=True, separators=(',', ':')))


if __name__ == '__main__':
    main()
