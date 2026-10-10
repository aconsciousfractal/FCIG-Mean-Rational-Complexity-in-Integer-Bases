"""Bind the final full154p replay and its transparent count-only recovery.

Adapted from P72_5610_102P_20261009/code/family_custody.py. The original
outer driver is preserved as failed; the complete scientific subprocess
and all immutable records are validated independently of that counter.
"""
from pathlib import Path
import json,re,difflib
from exact import need
from custody import sha


def family_replay(root):
    root = Path(root)
    component, logs = root/'family154p', root/'logs'
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
    # The scientific subprocess completed PASS; its outer driver failed only
    # because it expected91 comparisons and the new mixed-index check makes92.
    # Preserve the original failed source/precontract and validate the transparent
    # count-only finalization. No scientific rerun or family mutation pays this gate.
    recovery_path=logs/'family_full_root_RECOVERY.json'
    recovery=json.loads(recovery_path.read_text())
    need(root_receipt['status']=='PASS_SCIENTIFIC_REPLAY_RECOVERED' and
         root_receipt['scientific_replay_status']=='PASS' and root_receipt['recovery_status']=='PASS' and
         root_receipt['outer_driver_status']=='FAILED_COUNT_ASSERTION' and
         root_receipt['scientific_command_returncode']==0 and root_receipt['outer_driver_returncode']==1 and
         root_receipt['contract_sha256']==sha(logs/'family_full_root_CONTRACT.json'),
         'completed scientific replay distinguished from failed outer count assertion')
    need(recovery['status']=='COUNT_EXPECTATION_CORRECTED' and
         recovery['scientific_replay_status']=='PASS' and recovery['outer_driver_status']=='FAILED_COUNT_ASSERTION' and
         recovery['scientific_command_returncode']==0 and recovery['outer_driver_returncode']==1 and
         recovery['original_precontract_unchanged'] and not recovery['scientific_changes'] and
         not recovery['scientific_rerun'] and not root_receipt['scientific_replay_rerun'],
         'transparent arithmetic-only counter recovery without changing or rerunning science')
    need(root_receipt['recovery_note_sha256']==sha(recovery_path) and recovery['all190_payload_files_unchanged'],
         'source-bound count recovery report')
    original=logs/recovery['original_driver_copy'];operational=root/recovery['operational_driver_path']
    need(sha(original)==root_contract['driver_sha256']==recovery['precommitted_driver_sha256']==
         root_receipt['executed_driver_sha256'], 'byte-identical actually executed precommitted driver')
    old=original.read_text();need(len(re.findall(r'\b91\b',old))==3,'exact three original count occurrences')
    corrected=re.sub(r'\b91\b','92',old)
    need(operational.read_text()==corrected and sha(operational)==recovery['corrected_driver_sha256']==
         root_receipt['corrected_operational_driver_sha256'], 'only91-to92 operational source correction')
    diff=''.join(difflib.unified_diff(old.splitlines(keepends=True),corrected.splitlines(keepends=True),
         fromfile='executed_precommitted_driver91.py',tofile='operational_driver92.py'))
    need((logs/recovery['driver_diff']).read_text()==diff and
         sha(logs/recovery['driver_diff'])==recovery['driver_diff_sha256'], 'exact preserved source diff')
    need(sha(logs/recovery['recovery_script_copy'])==recovery['recovery_script_sha256'] and
         recovery['original_driver_archived_before_patch'] and recovery['operational_patch_after_driver_termination'],
         'finalizer source and completed-process correction order preserved')
    need(recovery['observed_outer_failure_log']=='family_full_root_OUTER_FAILURE.log',
         'actual outer-driver failure transcript retained')
    failure=(logs/recovery['observed_outer_failure_log']).read_bytes()
    need(b'AssertionError' in failure and b"['compared_files']==91" in failure,
         'observed failure was precisely the outer comparison-count equality')
    need(recovery['precommitted_expected_comparisons']==root_contract['required_geometric_comparisons']==
         root_receipt['precommitted_expected_geometric_comparisons']==91 and
         recovery['actual_completed_comparisons']==root_receipt['geometric_comparisons']==92 and
         recovery['comparison_decomposition']==dict(profile_comparisons=90,mixed_index_certificate=1,
                                                   negative_fixture_certificate=1,total=92),
         'ninety geometry comparisons plus one mixed-index and one negative-suite comparison')
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
         result['insertion_certificate_sha256'] == sha(component/'insertion/results/insertion154_certificate.json'),
         'executed family closure and insertion are final bytes')
    need(result['core'] == 154 and result['prime_min'] == 13 and
         result['all_positive_valuations'] and result['all_N_ge_1'], 'complete family scope in replay')

    plan = json.loads((component/'profiles/profiles_plan.json').read_text())
    pc = json.loads((component/'profiles/profile_contract.json').read_text())
    need(pc['core'] == 154 and pc['core_primes'] == [2,7,11], 'correct complete three-prime anchor mark')
    need(sorted((r['p'],r['multiplier'],r['ag']) for r in plan) ==
         sorted([(1,1,16),(1,2,16),(1,7,16),(1,11,16),(13,1,4),(17,1,4),(19,1,4),(23,1,4),(29,1,4),(31,1,4)]),
         'complete predeclared ten-profile family plan')
    records = sum(2688*288*(r['p']+1 if r['p'] > 1 else r['multiplier']) for r in plan)
    corners = sum(2688*288*(r['p']+1 if r['p'] > 1 else r['multiplier'])*4*(2*r['ag']+1) for r in plan)
    pr = result['profiles']
    need(pr['status'] == 'PASS' and pr['mode'] == 'full_witness_rebuild' and
         pr['primary_geometric_witnesses_replayed'] and pr['negative_tests_reexecuted'] and
         pr['sanitize'] and pr['threads'] == contract['threads'] == root_contract['threads'] == 3,
         'actual full UBSan family geometry and fixture replay')
    need(pr['profiles'] == root_contract['profiles'] == 10 and
         pr['geometric_witnesses'] == root_contract['geometric_records'] == records == 123088896 and
         pr['geometric_corner_evaluations'] == root_contract['corner_evaluations'] == corners == 5991874560,
         'all ten profile records and corner evaluations counted once')
    need(pr['compared_files']==root_receipt['geometric_comparisons']==9*len(plan)+2==92 and
         pr['compared_files']>root_contract['required_geometric_comparisons'] and
         pr['mixed_index_exact_checks']==61939,
         'all geometric comparisons plus the additional complete mixed-index comparison completed')
    need(recovery['scientific_report_sha256']==root_receipt['report_sha256']==sha(logs/'family_full_replay.json') and
         recovery['unchanged_family_manifest_sha256']==sha(manifest_path), 'recovery binds the actual complete report and payload')
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
    need(scripts == ['verify_profiles.py','verify_insertion154.py','run_finite154.py','independent_physical.py'],
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
         (46,1079232,539662,4317296), 'complete physical residual mass')
    need(independent['status'] == 'PASS' and independent['independent_full_raw_consumer'] and
         independent['consumer_has_no_producer_imports'] and not independent['assumes_minimality'] and
         not independent['uses_lattice_reduction'], 'separate complete constructive physical consumer')
    need((independent['moduli'],independent['states'],independent['fully_consumed_records'],independent['raw_bytes']) ==
         (46,1079232,539662,4317296), 'independent consumer covers every physical record')
    need(independent['finite_scope_N'] == 1 and not independent['all_N_claimed'] and
         independent['negative_upper_means'] == 46, 'truthful finite scope and all negative upper means')
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
        profiles=10,geometric_records=records,corner_evaluations=corners,geometric_comparisons=92,
        physical_moduli=46,physical_states=1079232,physical_records=539662,
        independent_physical_checks=independent['checks'],runtime_seconds=result['runtime_seconds'],
        scientific_command_returncode=0,outer_driver_returncode=1,
        outer_driver_status='FAILED_COUNT_ASSERTION',count_only_recovery_verified=True,
        original_expected_comparisons=91,actual_completed_comparisons=92,
        recovery_note_sha256=sha(recovery_path),mixed_index_exact_checks=61939)
