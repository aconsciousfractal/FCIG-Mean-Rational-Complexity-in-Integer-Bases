"""Precontract and compile separate UBSan producer/consumer in a fresh /tmp directory."""
from pathlib import Path
import hashlib, json, subprocess, time, tempfile, sys
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'code'))
from exact import need

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    src = {name: sha(R / name) for name in ('code/angle.cpp', 'code/algebra.cpp', 'code/kernel.hpp', 'code/search.hpp', 'build.py', 'profiles_plan.json')}
    flags = ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all']
    binary_dir = Path(tempfile.mkdtemp(prefix='p72-6510-initial-bin-'))
    plan = json.loads((R / 'profiles_plan.json').read_text())
    ags = sorted({z['ag'] for z in plan})
    commands = [['g++', *flags, str(R / 'code/algebra.cpp'), '-o', str(binary_dir / 'algebra')]]
    for ag in ags:
        for mode in ('produce', 'consume'):
            commands.append(['g++', *flags, f'-DANGULAR={ag}', *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), str(R / 'code/angle.cpp'), '-o', str(binary_dir / f'{mode}{ag}')])
    cp = R / 'logs/build_CONTRACT.json'
    need(not cp.exists(), 'fresh build precontract')
    cp.write_text(json.dumps({'source_sha256': src, 'flags': flags, 'commands': commands, 'consumer_search_code_included': False, 'binary_directory': str(binary_dir), 'binary_directory_is_temporary': True}, indent=2) + '\n', newline='\n')
    rows = []
    for cmd in commands:
        tick = time.monotonic()
        z = subprocess.run(cmd, text=True, capture_output=True)
        rows.append({'argv': cmd, 'exit_code': z.returncode, 'seconds': round(time.monotonic() - tick, 3), 'transcript': z.stdout + z.stderr})
        print(json.dumps(rows[-1]), flush=True)
        need(z.returncode == 0, 'every compilation succeeded')
        Path(cmd[-1]).chmod(448)
    need(all((sha(R / k) == v for k, v in src.items())), 'all compiled inputs unchanged')
    binaries = {p.name: sha(p) for p in binary_dir.iterdir() if p.is_file()}
    need(set(binaries) == {'algebra', *(f'{mode}{ag}' for ag in ags for mode in ('produce', 'consume'))}, 'exact selected executable inventory')
    (R / 'logs/build_RECEIPT.json').write_text(json.dumps({'status': 'PASS', 'contract_sha256': sha(cp), 'source_sha256': src, 'binaries': binaries, 'binary_directory': str(binary_dir), 'commands': rows}, indent=2) + '\n', newline='\n')
if __name__ == '__main__':
    main()
