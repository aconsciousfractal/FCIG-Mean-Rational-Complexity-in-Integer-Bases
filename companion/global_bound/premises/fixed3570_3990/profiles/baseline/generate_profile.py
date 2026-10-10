"""Freeze a declared contract before a new K=2q geometric reconstruction."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, os, resource, subprocess, sys, time
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'code'))
from core import mark
from exact import need, log_lo

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while (b := f.read(8 * 1024 * 1024)):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int, choices=(17, 19))
    ap.add_argument('--threads', type=int, choices=(1, 2, 3), default=3)
    ap.add_argument('--prepare', action='store_true')
    ap.add_argument('--work-dir', type=Path)
    a = ap.parse_args()
    need(sys.version_info >= (3, 10) and __debug__, 'Python3.10+ with assertions enabled')
    q, nl = mark(a.p)
    m = 2
    K = m * q
    style = 'coarse'
    ag = 8
    stem = f'q{q}_K{K}_{style}'
    plan = json.loads((R / 'profiles_plan.json').read_text())
    need(dict(stem=stem, p=a.p, multiplier=m, style=style, angular=ag) in plan, 'declared profile plan')
    paths = [R / 'generate_profile.py', R / 'profiles_plan.json', R / 'REUSED_SOURCES.json', *sorted((R / 'code').glob('*')), *sorted((R / 'data').glob('*')), *sorted((R / 'imports').glob('*'))]
    inputs = {p.relative_to(R).as_posix(): sha(p) for p in paths if p.is_file()}
    contract = dict(p=a.p, q=q, multiplier=m, transfer_index=K, style=style, source_cells=2667, output_cells=2667, labels=nl, branches=nl * m, geometric_witness_records=2667 * nl * m, witness_bytes=2667 * nl * m * 8, geometric_corner_evaluations=2667 * nl * m * 9 * 4, coarse_band_range=[0, 126], refined_band_range=None, radial_parts=1, fine_angular=None, coarse_angular=4, gauss_crop=5200308914369309, threads=a.threads, sanitize=True, producer='only fixed coefficient pairs', consumer='all coefficients and all nodes of all continuous rectangles', normalization='log B', full_mark=list((2, 3, 5, 7, a.p)), universal_module_domain='every positive multiple of K', target_mu_lower=str(F(2, 85) * log_lo(q)), target_CDF_coefficient='3/2', CDF_proof_source='imports/PRIMITIVE_CDF.md', inputs=inputs)
    cp = R / 'contracts' / (stem + '.json')
    if a.prepare:
        need(not cp.exists(), 'contract must not overwrite an earlier record')
        cp.write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
        print(json.dumps(dict(status='CONTRACT_FROZEN', stem=stem, contract_sha256=sha(cp), records=contract['geometric_witness_records'], target_mu_lower=contract['target_mu_lower']), indent=2), flush=True)
        return
    need(cp.exists() and json.loads(cp.read_text()) == contract, 'pre-existing unchanged run contract')
    need(a.work_dir is not None, 'external fresh work directory required')
    w = a.work_dir.resolve()
    need(not w.is_relative_to(R) and (not w.exists()), 'fresh external generation workspace')
    w.mkdir(parents=True)
    pre = w / 'profile'
    out = R / 'results' / stem
    cs = sha(cp)
    env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    start = time.monotonic()
    stages = []
    flags = ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all', f'-DINDEX_MULTIPLIER={m}', f'-DFINE_ANGULAR={ag}']

    def run(cmd):
        cmd = list(map(str, cmd))
        print('RUN', ' '.join(cmd), flush=True)
        t = time.monotonic()
        subprocess.run(cmd, check=True, env=env)
        stages.append(dict(command=cmd, seconds=round(time.monotonic() - t, 3)))
    run(['g++', *flags, R / 'code/algebra.cpp', '-o', w / 'algebra'])
    run([w / 'algebra', a.p, 1, R / 'data', str(pre) + '.tensor'])
    for mode in ('produce', 'consume'):
        binary = w / mode
        run(['g++', *flags, *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), R / 'code/profile.cpp', '-o', binary])
        run([binary, mode, a.p, style, R / 'data', pre])
    run([sys.executable, '-B', R / 'code/integrate.py', a.p, m, ag, style, pre, str(pre) + '.json'])
    run([sys.executable, '-B', R / 'code/integrate_alpha.py', a.p, m, style, '--angular', ag, '--prefix', pre, '--output', str(pre) + '.alpha3over2.json'])
    for name, d in inputs.items():
        need(sha(R / name) == d, 'input unchanged ' + name)
    need(sha(cp) == cs, 'contract unchanged after generation')
    meta = dict(bytes=Path(str(pre) + '.witness').stat().st_size, sha256=sha(str(pre) + '.witness'))
    need(meta['bytes'] == contract['witness_bytes'], 'exact complete witness extent')
    for ext in ('.cells', '.audit', '.branches', '.tensor', '.json', '.alpha3over2.json'):
        Path(str(out) + ext).write_bytes(Path(str(pre) + ext).read_bytes())
    Path(str(out) + '.witness_metadata.json').write_text(json.dumps(meta, indent=2) + '\n', newline='\n')
    report = dict(status='PASS', stem=stem, contract_sha256=cs, contract_prepared_before_run=True, inputs_unchanged=True, stages=stages, runtime_seconds=round(time.monotonic() - start, 3), peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024, witness=meta)
    (R / 'logs' / (stem + '.generation.json')).write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report, indent=2), flush=True)
if __name__ == '__main__':
    main()
