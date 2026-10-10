from pathlib import Path
import hashlib, json, os, resource, subprocess, time
R = Path(__file__).resolve().parent
assert not (R / 'work/profile.witness').exists()
files = [R / 'CONTRACT.json', *sorted((R / 'code').glob('*')), *sorted((R / 'data').glob('*')), *sorted((R / 'imports').glob('*'))]

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
manifest = {p.relative_to(R).as_posix(): digest(p) for p in files if p.is_file()}
(R / 'logs/generation_input_sha256.json').write_text(json.dumps(manifest, indent=2) + '\n', newline='\n')
env = {**os.environ, 'CSTAR_THREADS': '4', 'PYTHONDONTWRITEBYTECODE': '1'}
start = time.monotonic()
stages = []

def run(args):
    print('RUN', ' '.join(map(str, args)), flush=True)
    s = time.monotonic()
    subprocess.run(list(map(str, args)), check=True, env=env)
    stages.append(dict(command=list(map(str, args)), seconds=round(time.monotonic() - s, 3)))
for mode in ['produce', 'consume']:
    binary = R / 'work' / mode
    run(['g++', '-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all', *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), R / 'code/refine.cpp', '-o', binary])
    run([binary, mode, R / 'data', R / 'work/profile'])
    if mode == 'produce':
        got = digest(R / 'work/profile.witness')
        h = hashlib.sha256()
        source = R.parent / 'work_mixed/q2310_K13860_AG4.witness'
        with source.open('rb') as f:
            f.seek(90 * 21 * 41472 * 8)
            remaining = 777 * 41472 * 8
            while remaining:
                b = f.read(min(8 * 1024 * 1024, remaining))
                assert b
                h.update(b)
                remaining -= len(b)
        assert got == h.hexdigest(), 'complete identity with selected old witness span'
        (R / 'results/witness_metadata.json').write_text(json.dumps(dict(bytes=777 * 41472 * 8, records=777 * 41472, sha256=got, complete_byte_identity_with_baseline_span=True, baseline_cell_span=[1890, 2666]), indent=2) + '\n', newline='\n')
for name, d in manifest.items():
    assert digest(R / name) == d, 'frozen input mutated ' + name
rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024
assert rss < 2 * 1024 ** 3
report = dict(status='PASS', stages=stages, runtime_seconds=round(time.monotonic() - start, 3), peak_child_rss_bytes=rss, threads=4, sanitize=True, frozen_input_files=len(manifest), consumer_exhaustive=True, producer_computes_geometric_caps=False)
(R / 'logs/generation_report.json').write_text(json.dumps(report, indent=2) + '\n', newline='\n')
print(json.dumps(report, indent=2), flush=True)
