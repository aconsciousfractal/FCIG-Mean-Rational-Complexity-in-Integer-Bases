"""Compile explicit110 geometry with a pre-run source contract and UBSan."""
from pathlib import Path
import hashlib, json, subprocess, time
ROOT = Path(__file__).resolve().parent
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
sources = {str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT / 'code').iterdir()) if p.suffix in ('.hpp', '.cpp')}
contract = {'status': 'PRECOMMITTED', 'source_sha256': sources, 'compiler': 'g++', 'flags': ['-O3', '-std=c++17', '-fopenmp', '-fsanitize=undefined', '-fno-sanitize-recover=all'], 'targets': ['algebra', 'produce4', 'consume4', 'produce8', 'consume8', 'produce16', 'consume16']}
(ROOT / 'logs/build_CONTRACT.json').write_text(json.dumps(contract, indent=2) + '\n', newline='\n')
commands = []
start = time.monotonic()
for name in contract['targets']:
    mode = 'consume' if name.startswith('consume') else 'produce'
    ag = 4 if name.endswith('4') else 8 if name.endswith('8') else 16
    src = 'algebra.cpp' if name == 'algebra' else 'angle.cpp'
    extra = [] if name == 'algebra' else [f'-DANGULAR={ag}'] + (['-DBUILD_CONSUMER=1'] if mode == 'consume' else [])
    cmd = ['g++', *contract['flags'], *extra, str(ROOT / 'code' / src), '-o', str(ROOT / 'bin' / name)]
    t = time.monotonic()
    r = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    commands.append({'argv': cmd, 'exit_code': r.returncode, 'seconds': round(time.monotonic() - t, 3), 'output': r.stdout})
    print(name, r.returncode, flush=True)
    if r.returncode:
        raise RuntimeError(r.stdout)
assert sources == {s: sha(ROOT / s) for s in sources}
receipt = {'status': 'PASS', 'source_sha256': sources, 'contract_sha256': sha(ROOT / 'logs/build_CONTRACT.json'), 'binaries': {name: sha(ROOT / 'bin' / name) for name in contract['targets']}, 'compiler_version': subprocess.check_output(['g++', '--version'], text=True).splitlines()[0], 'commands': commands, 'seconds': round(time.monotonic() - start, 3)}
(ROOT / 'logs/build_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n', newline='\n')
