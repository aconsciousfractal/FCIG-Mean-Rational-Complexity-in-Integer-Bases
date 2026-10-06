"""Check distributed complete evidence; recount commands are separate."""
import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'companion/programs'))
from archive import PACKAGE, entry, read_object
from reconstruction_chain import REQUIRED_JOB_IDS, relative
from verify_all_central import BATCHES, load_batch, bind_validation, load_producer_arrays, strict_load
from verify_full_semiprime import verify as verify_semiprime

ROOT = PACKAGE.parent
BATCH_IDS = (29,30,31,32,33,35,39,40,41)
COUNTS = (143,97,132,3,35,99,1087,1422,1901)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def verify_observations(index, plan):
    jobs = plan['jobs']
    need({j['id'] for j in jobs} == REQUIRED_JOB_IDS and len(jobs) == 37,
         'complete unique 37-job plan')
    need(index['normal_optimized_scientific_bytes_identical'] is True,
         'reported normal/optimized identity')
    need(set(index['observed_runs']) == {'normal','optimized'}, 'both observed modes')
    for mode, rows in index['observed_runs'].items():
        need(len(rows) == 37 and [r['id'] for r in rows] == [j['id'] for j in jobs],
             'complete ordered observed job coverage')
        for row, job in zip(rows, jobs):
            record = index['records']['chain/' + job['id']]
            need(row['status'] == 'PASS' and type(row['exit_code']) is int
                 and row['exit_code'] == 0 and row['retained_math_equal'] is True
                 and row['raw_retained_equal'] is True, 'observed job success')
            need(record['sha256'] == row['consumed_output_sha256']
                 == job['output']['retained_sha256'], 'generated/consumed/original byte identity')


def check(index_path=None):
    path = ROOT/'companion/reconstruction/index.json' if index_path is None else Path(index_path)
    index = strict_load(path.read_bytes())
    need(index['schema'] == 'COMPLETE_PROOF_EVIDENCE_V1'
         and index['mathematical_results_changed'] is False, 'complete evidence schema')
    need(index['version'] == (ROOT/'VERSION').read_text().strip(), 'evidence version')
    for name, pin in index['source_pins'].items():
        name = relative(name)
        source = ROOT/name
        need(source.resolve().is_relative_to(ROOT) and not source.is_symlink(), 'source path')
        need(sha(source.read_bytes()) == pin, 'executed scientific source identity: ' + name)
    need(sha((PACKAGE/'provenance.json').read_bytes()) == index['provenance_sha256']
         and sha((PACKAGE/'portable_source_map.json').read_bytes()) == index['source_map_sha256'],
         'original public scientific graph identity')
    plan = strict_load((PACKAGE/'reconstruction_plan.json').read_bytes())
    expected_labels = {'chain/'+j['id'] for j in plan['jobs']}
    expected_labels |= {'semiprime/full','semiprime/logcheck','inventory/frontier'}
    expected_labels |= {'central/'+str(b) for b in BATCH_IDS}
    expected_labels |= {'capture/'+str(b) for b in BATCH_IDS[:6]}
    expected_labels |= {'pointwise/'+k for k in ('semiprime','base6','sixth')}
    need(set(index['records']) == expected_labels, 'all and only complete evidence records')
    verify_observations(index, plan)
    for label, record in index['records'].items():
        need(relative(record['storage']).startswith('objects/'), 'relative object store')
        read_object(record)
    main = controls = thresholds = cells = control_thresholds = control_cells = 0
    main_moduli = set()
    control_moduli = set()
    for batch, count in zip(BATCH_IDS, COUNTS):
        report = strict_load(read_object(index['records']['central/'+str(batch)]))
        producer, profiles, pins = load_batch(batch)
        validations = bind_validation(batch, producer)
        need(report['status'] == 'COMPLETE_REQUESTED_CHUNK_PASS'
             and report['complete_batch_verified'] is True and report['batch'] == batch
             and report['start'] == 0 and report['stop'] == report['complete_batch_rows'] == count,
             'complete aggregate report')
        need(report['coverage_indices'] == list(range(count))
             and [r['domain_index'] for r in report['rows']] == list(range(count))
             and [r['profile'] for r in report['rows']] == profiles, 'exact central primary domain')
        need(report['validation_profiles'] == validations
             and [r['validation_index'] for r in report['validation_rows']] == list(range(len(validations)))
             and [r['profile'] for r in report['validation_rows']] == validations,
             'separate complete validation domain')
        for k,v in pins.items():
            need(report['pins'][k] == v, 'source and consumed certificate binding')
        if batch < 39:
            record = index['records']['capture/'+str(batch)]
            with tempfile.TemporaryDirectory(prefix='central_capture_') as directory:
                capture = Path(directory)/'capture.json'
                capture.write_bytes(read_object(record))
                combined = list(profiles)
                for p in validations:
                    if p not in combined: combined.append(p)
                load_producer_arrays(batch,capture,record['sha256'],combined)
            need(report['pins']['producer_arrays_sha256'] == record['sha256'], 'capture/report binding')
        for validation, rows in ((False,report['rows']),(True,report['validation_rows'])):
            for row in rows:
                need(row['status'] == 'ALL_CENTRAL_THRESHOLDS_AND_STRATA_PASS'
                     and type(row['all_thresholds_checked']) is int
                     and row['all_thresholds_checked'] == isqrt(row['profile'][0]),
                     'complete central threshold record')
                if validation:
                    controls += 1
                    control_thresholds += row['all_thresholds_checked']
                    control_cells += row['counters']['proper_stratum_thresholds']
                    control_moduli.add(row['profile'][0])
                else:
                    main += 1
                    thresholds += row['all_thresholds_checked']
                    cells += row['counters']['proper_stratum_thresholds']
                    main_moduli.add(row['profile'][0])
    census = dict(batch_count=9,main_profile_obligations=main,distinct_main_moduli=len(main_moduli),
                  validation_comparisons=controls,distinct_validation_moduli=len(control_moduli),
                  main_thresholds=thresholds,main_proper_stratum_threshold_cells=cells,
                  validation_thresholds=control_thresholds,
                  validation_proper_stratum_threshold_cells=control_cells)
    need(census == index['central_counts'], 'complete central census')
    pointwise = {}
    for kind,count,source in (('semiprime',4596,'semiprime/full'),('base6',8,'chain/base6'),('sixth',3,'chain/sixth')):
        record = index['records']['pointwise/'+kind]
        report = strict_load(read_object(record))
        premise = strict_load(read_object(index['records'][source]))
        need(report['status'] == 'COMPLETE_POINTWISE_CENTRAL_CHUNK_PASS'
             and report['complete_domain_verified'] is True and report['start'] == 0
             and report['stop'] == report['full_domain_count'] == count
             and report['kind'] == kind and len(premise['rows']) == count,
             'complete pointwise report')
        need(report['input_sha256'] == index['records'][source]['sha256']
             and [r['domain_index'] for r in report['rows']] == list(range(count))
             and [r['M'] for r in report['rows']] == [r['M'] for r in premise['rows']],
             'pointwise exact domain and consumed input')
        need(all(r['status']=='ALL_POINTWISE_CENTRAL_THRESHOLDS_PASS'
                 and r['k']==isqrt(r['M']) for r in report['rows'])
             and report['all_thresholds_checked']==sum(r['k'] for r in report['rows']),
             'all pointwise central threshold records')
        pointwise[kind] = report['all_thresholds_checked']
    record = index['records']['semiprime/full']
    with tempfile.TemporaryDirectory(prefix='semiprime_evidence_') as directory:
        full = Path(directory)/'full.json'
        full.write_bytes(read_object(record))
        logcheck = verify_semiprime(full,record['sha256'])
    need(logcheck == strict_load(read_object(index['records']['semiprime/logcheck'])),
         'distinct log replay agrees with recorded full check')
    need(read_object(index['records']['inventory/frontier'])
         == read_object(entry(36,'FRONTIER_CERTIFICATE.json')), 'derived inventory certificate identity')
    run = subprocess.run([sys.executable,'-B',*(['-O'] if sys.flags.optimize else []),
                          str(PACKAGE/'programs/verify_retained.py')],
                         check=True,capture_output=True)
    retained = strict_load(run.stdout)
    need(retained['status'] == 'PORTABLE_RETAINED_ARITHMETIC_PASS'
         and retained['final_bridge_rows'] == 4410 and retained['two_prime_remaining'] == 0,
         'retained arithmetic and domain closure')
    return dict(status='COMPLETE_DISTRIBUTED_EVIDENCE_AND_DISTINCT_LOG_PASS',
                evidence_records=len(expected_labels),observed_reconstruction_jobs_per_mode=37,
                central_counts=census,pointwise_thresholds=pointwise,
                semiprime=logcheck,central_recount_executed_by_this_command=False,
                full_reconstruction_executed_by_this_command=False,
                retained_arithmetic_rechecked=True,formal_proof_or_originality_certified=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(check(),sort_keys=True,separators=(',',':')))


if __name__=='__main__': main()
