"""Compile separate producer and consumer after freezing their source contract."""
from pathlib import Path
import hashlib, json, subprocess, time
R = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    src = {p.relative_to(R).as_posix(): sha(p) for p in (R / 'code/angle.cpp', R / 'code/algebra.cpp', R / 'code/kernel.hpp', R / 'code/search.hpp')}
    flags = ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all']
    commands = [['g++', *flags, str(R / 'code/algebra.cpp'), '-o', str(R / 'bin/algebra')]]
    for ag in (8,):
        for mode in ('produce', 'consume'):
            commands.append(['g++', *flags, f'-DANGULAR={ag}', *(['-DBUILD_CONSUMER=1'] if mode == 'consume' else []), str(R / 'code/angle.cpp'), '-o', str(R / 'bin' / f'{mode}{ag}')])
    cp = R / 'logs/build_CONTRACT.json'
    assert not cp.exists()
    cp.write_text(json.dumps({'source_sha256': src, 'flags': flags, 'commands': commands, 'consumer_search_code_included': False}, indent=2) + '\n', newline='\n')
    rows = []
    for cmd in commands:
        t = time.monotonic()
        z = subprocess.run(cmd, text=True, capture_output=True)
        rows.append({'argv': cmd, 'exit_code': z.returncode, 'seconds': round(time.monotonic() - t, 3), 'transcript': z.stdout + z.stderr})
        print(rows[-1], flush=True)
        assert z.returncode == 0
    assert all((sha(R / k) == v for k, v in src.items()))
    (R / 'logs/build_RECEIPT.json').write_text(json.dumps({'status': 'PASS', 'contract_sha256': sha(cp), 'source_sha256': src, 'binaries': {p.name: sha(p) for p in (R / 'bin').iterdir() if p.is_file()}, 'commands': rows}, indent=2) + '\n', newline='\n')
if __name__ == '__main__':
    main()
