"""Create the run contract before compiling, then freeze all numerical inputs."""
from pathlib import Path
import argparse, hashlib, json, os, resource, subprocess, time
R = Path(__file__).resolve().parent

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        while (b := f.read(8 * 1024 * 1024)):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('multiplier', type=int, choices=(1, 6))
    ap.add_argument('--angular', type=int, default=8, choices=(8, 16))
    ap.add_argument('--threads', type=int, default=2, choices=(1, 2, 3, 4))
    a = ap.parse_args()
    assert a.p in (13, 17, 19)
    style = 'coarse' if a.multiplier == 1 else 'hybrid'
    q = 210 * a.p
    stem = f'q{q}_K{q * a.multiplier}_{style}' + ('' if style == 'coarse' else f'_AG{a.angular}')
    w = R / 'work' / stem
    w.mkdir()
    pre = w / 'profile'
    out = R / 'results' / stem
    paths = [R / 'generate_profile.py', R / 'code/integrate.py', R / 'code/assemble.py', R / 'code/core.py', R / 'code/exact.py', R / 'code/kernel.hpp', R / 'code/search.hpp', R / 'code/profile.cpp', R / 'code/algebra.cpp', *sorted((R / 'data').glob('*')), *sorted((R / 'imports').glob('*'))]
    inputs = {p.relative_to(R).as_posix(): sha(p) for p in paths if p.is_file()}
    contract = dict(p=a.p, q=q, multiplier=a.multiplier, transfer_index=q * a.multiplier, style=style, source_cells=2667, coarse_band_range=[0, 89] if style == 'hybrid' else [0, 126], refined_band_range=[90, 126] if style == 'hybrid' else None, radial_parts=4 if style == 'hybrid' else 1, fine_angular=a.angular if style == 'hybrid' else None, coarse_angular=4, gauss_crop=5200308914369309, threads=a.threads, sanitize=True, producer='only fixed coefficient pairs', consumer='all coefficients and all nodes of all continuous rectangles', inputs=inputs)
    cp = R / 'contracts' / (stem + '.json')
    assert not cp.exists()
    cp.write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
    cs = sha(cp)
    env = {**os.environ, 'CSTAR_THREADS': str(a.threads), 'PYTHONDONTWRITEBYTECODE': '1'}
    start = time.monotonic()
    stages = []
    flags = ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all', f'-DINDEX_MULTIPLIER={a.multiplier}', f'-DFINE_ANGULAR={a.angular}']

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
    run(['python3', '-B', R / 'code/integrate.py', a.p, a.multiplier, a.angular, style, pre, str(pre) + '.json'])
    for name, d in inputs.items():
        assert sha(R / name) == d, 'input changed ' + name
    assert sha(cp) == cs, 'contract changed'
    meta = dict(bytes=Path(str(pre) + '.witness').stat().st_size, sha256=sha(str(pre) + '.witness'))
    for ext in ('.cells', '.audit', '.branches', '.tensor', '.json'):
        Path(str(out) + ext).write_bytes(Path(str(pre) + ext).read_bytes())
    Path(str(out) + '.witness_metadata.json').write_text(json.dumps(meta, indent=2) + '\n', newline='\n')
    report = dict(status='PASS', stem=stem, contract_sha256=cs, inputs_unchanged=True, stages=stages, runtime_seconds=round(time.monotonic() - start, 3), peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024, witness=meta)
    (R / 'logs' / (stem + '.generation.json')).write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report, indent=2), flush=True)
if __name__ == '__main__':
    main()
