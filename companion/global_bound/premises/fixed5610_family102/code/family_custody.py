"""Bind the final full 102p replay, including a separate complete raw consumer."""
from pathlib import Path
import json
from exact import need
from custody import sha


def family_replay(root):
    root = Path(root)
    component, logs = root/'family102p', root/'logs'
    manifest_path = component/'MANIFEST_SHA256.json'
    manifest = json.loads(manifest_path.read_text())
    present = {str(p.relative_to(component)) for p in component.rglob('*')
               if p.is_file() and p != manifest_path}
    need(set(manifest) == present, 'complete immutable family payload inventory')
    for name,digest in manifest.items():
        path = component/name
        need(path.resolve().is_relative_to(component.resolve()) and sha(path) == digest,
             'unchanged final family payload '+name)

    def read(suffix):
        return json.loads((logs/('family_full_replay'+suffix)).read_text())

    result = read('.json')
    receipt = read('_RECEIPT.json')
    contract = read('_CONTRACT.json')
    independent = read('_independent_physical.json')
    physical = read('_physical_RECEIPT.json')
    physical_contract = read('_physical_CONTRACT.json')
    source_manifest = read('_SOURCE_MANIFEST.json')
    root_contract = json.loads((logs/'family_full_root_CONTRACT.json').read_text())
    root_receipt = json.loads((logs/'family_full_root_RECEIPT.json').read_text())
    need(source_manifest == manifest and
         (logs/'family_full_root_SOURCE_MANIFEST.json').read_bytes() == manifest_path.read_bytes(),
         'executed family inventory is exactly the final inventory')
    need(root_receipt['status'] == 'PASS' and root_receipt['process_returncode'] == 0 and
         root_receipt['contract_sha256'] == sha(logs/'family_full_root_CONTRACT.json'),
         'root completed full family process receipt')
    need(root_receipt['family_manifest_sha256'] == root_contract['family_manifest_sha256'] ==
         result['family_manifest_sha256'] == contract['family_manifest_sha256'] == sha(manifest_path),
         'root, family precontract, executed report and final manifest identity')
    need(root_receipt['executed_payloads_verified_unchanged_at_completion'] and
         root_receipt['executed_payload_files'] == contract['source_file_count'] ==
         result['manifest_files'] == len(manifest), 'all family payloads verified after execution')
    for name,digest in root_receipt['evidence_sha256'].items():
        need(sha(logs/name) == digest, 'complete external execution evidence '+name)
    need(result['status'] == receipt['status'] == 'PASS' and
         result['mode'] == receipt['mode'] == contract['mode'] == 'full', 'actual complete family mode')
    need(result['source_inputs_unchanged'] and receipt['all_source_inputs_unchanged'] and
         receipt['all_commands_succeeded'], 'full reconstruction preserved its inputs')
    for key,suffix in [('report_sha256','.json'),('contract_sha256','_CONTRACT.json'),
                       ('source_manifest_sha256','_SOURCE_MANIFEST.json'),('transcript_sha256','.log'),
                       ('independent_physical_report_sha256','_independent_physical.json')]:
        need(receipt[key] == sha(logs/('family_full_replay'+suffix)),
             'family receipt file identity '+suffix)
    need(result['closure_certificate_sha256'] == sha(component/'results/closure_certificate.json') and
         result['insertion_certificate_sha256'] == sha(component/'insertion/results/insertion102_certificate.json'),
         'executed family closure and insertion are final bytes')
    need(result['core'] == 102 and result['prime_min'] == 19 and
         result['all_positive_valuations'] and result['all_N_ge_1'], 'complete family scope in replay')

    plan = json.loads((component/'profiles/profiles_plan.json').read_text())
    pc = json.loads((component/'profiles/profile_contract.json').read_text())
    need(pc['core'] == 102 and pc['core_primes'] == [2,3,17], 'correct complete three-prime anchor mark')
    need(sorted((r['p'],r['multiplier'],r['ag']) for r in plan) ==
         sorted([(1,1,16),(1,2,16),(1,102,8),(19,1,4),(23,1,4),(29,1,4),(31,1,4)]),
         'complete predeclared seven-profile family plan')
    records = sum(2688*216*(r['p']+1 if r['p'] > 1 else r['multiplier']) for r in plan)
    corners = sum(2688*216*(r['p']+1 if r['p'] > 1 else r['multiplier'])*4*(2*r['ag']+1) for r in plan)
    pr = result['profiles']
    need(pr['status'] == 'PASS' and pr['mode'] == 'full_witness_rebuild' and
         pr['primary_geometric_witnesses_replayed'] and pr['negative_tests_reexecuted'] and
         pr['sanitize'] and pr['threads'] == contract['threads'] == root_contract['threads'] == 4,
         'actual full UBSan family geometry and fixture replay')
    need(pr['profiles'] == root_contract['profiles'] == 7 and
         pr['geometric_witnesses'] == root_contract['geometric_records'] == records == 122508288 and
         pr['geometric_corner_evaluations'] == root_contract['corner_evaluations'] == corners == 6472617984,
         'all seven profile records and corner evaluations counted once')
    need(pr['compared_files'] == root_contract['required_geometric_comparisons'] == 9*len(plan)+1 == 64,
         'all family geometric comparisons completed')
    need(pr['commands'] and all(r['returncode'] == 0 for r in pr['commands']),
         'every geometric subprocess succeeded')
    expected = sorted((r['p'],r['multiplier'],r['ag']) for r in plan)
    for role in ['produce','consume']:
        actual = []
        for row in pr['commands']:
            argv = row['argv']
            if len(argv) >= 4 and argv[1] == role:
                binary = Path(argv[0]).name
                need(binary.startswith(role), 'actual geometric binary role')
                actual.append((int(argv[2]),int(argv[3]),int(binary[len(role):])))
        need(sorted(actual) == expected, 'actual '+role+' invocation for every family profile')
    need(result['commands'] and all(row['returncode'] == 0 for row in result['commands']),
         'all family verifier subprocesses succeeded')
    scripts = [Path(row['argv'][2]).name for row in result['commands']]
    need(scripts == ['verify_profiles.py','verify_insertion102.py','run_finite102.py','independent_physical.py'],
         'complete ordered reconstruction and separate raw-consumer invocations')

    need(result['physical_primary_records_reconstructed'] and result['physical_independent_records_consumed'],
         'physical generation and independent raw consumption actually completed')
    need(result['physical_replay'] == physical and physical['status'] == 'PASS' and
         physical['returncode'] == 0 and physical['expected_byte_match'] and
         physical['all_stream_sha256_match'] and physical['inputs_unchanged'],
         'all physical streams reproduced and consumed against frozen outputs')
    need(physical['contract_sha256'] == sha(logs/'family_full_replay_physical_CONTRACT.json') and
         physical['output_sha256'] == sha(component/'results/finite_bounds.csv'),
         'completed physical process input and output identity')
    need('-fsanitize=undefined' in physical_contract['compile'] and
         '-fno-sanitize-recover=all' in physical_contract['compile'] and
         physical_contract['compare_expected'], 'actual UBSan physical compilation and frozen comparison')
    for name,digest in physical_contract['inputs'].items():
        need(sha(component/name) == digest, 'physical executed source identity '+name)
    need((physical['moduli'],physical['states'],physical['representatives'],physical['stream_bytes']) ==
         (117,4055316,2027775,16222200), 'complete physical residual mass')
    need(independent['status'] == 'PASS' and independent['independent_full_raw_consumer'] and
         independent['consumer_has_no_producer_imports'] and not independent['assumes_minimality'] and
         not independent['uses_lattice_reduction'], 'separate complete constructive physical consumer')
    need((independent['moduli'],independent['states'],independent['fully_consumed_records'],independent['raw_bytes']) ==
         (117,4055316,2027775,16222200), 'independent consumer covers every physical record')
    need(independent['finite_scope_N'] == 1 and not independent['all_N_claimed'] and
         independent['negative_upper_means'] == 117, 'truthful finite scope and all negative upper means')
    need(result['independent_physical_report_sha256'] == sha(logs/'family_full_replay_independent_physical.json') and
         result['independent_physical_source_sha256'] == contract['independent_physical_source_sha256'] ==
         sha(component/'code/independent_physical.py'), 'actual independent consumer and report identity')
    for name,digest in independent['inputs'].items():
        need(sha(component/name) == digest, 'independent physical source/table identity '+name)
    streams = json.loads((component/'results/finite_stream_manifest.json').read_text())
    projected = [dict(M=r['M'],p=r['p'],records=r['records'],bytes=8*r['records'],sha256=r['sha256'])
                 for r in independent['results']]
    need(projected == streams, 'independent full reader consumed exactly the frozen streams')
    need(result['independent_physical_checks'] == independent['checks'], 'independent check total identity')
    return dict(status='PASS',completed_full_replay_verified=True,
        completed_independent_physical_consumption_verified=True,
        family_manifest_sha256=sha(manifest_path),executed_payload_files=len(manifest),
        report_sha256=sha(logs/'family_full_replay.json'),
        root_receipt_sha256=sha(logs/'family_full_root_RECEIPT.json'),
        independent_physical_report_sha256=sha(logs/'family_full_replay_independent_physical.json'),
        profiles=7,geometric_records=records,corner_evaluations=corners,geometric_comparisons=64,
        physical_moduli=117,physical_states=4055316,physical_records=2027775,
        independent_physical_checks=independent['checks'],runtime_seconds=result['runtime_seconds'])
