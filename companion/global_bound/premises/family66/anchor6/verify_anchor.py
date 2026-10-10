"""Read-only exact q6 assembly; optional raw geometric rebuild in a temporary copy."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, os, shutil, struct, subprocess, sys, tempfile, time
from contextlib import contextmanager
from uuid import uuid4
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from exact import need
from assemble import profile

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

@contextmanager
def retained_workspace(keep):
    directory=Path(tempfile.gettempdir())/('marked6-'+uuid4().hex)
    print('VERIFIED ABSOLUTE ANCHOR WORK',directory.resolve(),flush=True)
    directory.mkdir();failed=True
    try:
        yield str(directory)
        failed=False
    finally:
        if keep or failed:print('Retained anchor work',directory.resolve(),flush=True)
        else:shutil.rmtree(directory)

def certificate(root=ROOT):
    p = root / 'results/q6.json'
    stored = json.loads(p.read_text())
    got = profile(1, root / 'results/q6.cells', root, 1, 16)
    need(all((stored[k] == v for k, v in got.items())), 'reintegrated marked-6 profile')
    need(F(stored['total_upper']) < F(79, 200), 'fresh anchor cap')
    witness = root / 'results/q6.witness'
    raw = witness.read_bytes()
    need(len(raw) == stored['witness_bytes'] == 258048 and sha(witness) == stored['witness_sha256'], 'included full q6 witness hash and length')
    branches = [list(map(int, l.split())) for l in (root / 'results/q6.branches').read_text().splitlines()]
    need(len(branches) == 12, 'twelve marked-6 branches')
    for n, (ci, cj) in enumerate(struct.iter_unpack('<ii', raw)):
        _, a, b, d, x, y = branches[n % 12]
        need((ci * x + cj * y) % 2 and (ci * x + cj * y) % 3, 'every raw record is unit-marked')
        u, v = (ci * a + cj * b, cj * d)
        need((u or v) and abs(u) < 1000000 and (abs(v) < 1000000), 'every raw record in range')
    audit = json.loads((root / 'results/q6.independent.json').read_text())
    need(audit['status'] == 'PASS' and audit['q'] == 6, 'independent marked-6 audit')
    need(audit['algebra']['branches'] == 12 and audit['algebra']['signed_direct_outputs'] == 24, 'all q6 projective labels and signed sums')
    return dict(status='PASS', q=6, all_multiples_of_6=True, upper=stored['total_upper'], relaxed_cap='79/200', raw_witnesses=32256, raw_bytes=len(raw), raw_sha256=sha(witness), all_raw_mark_checks=True, continuous_coverage_proved_in='proofs/PROOF.md; article appendix E1-E7', alternative_geometric_audit_exhaustive=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rebuild', action='store_true')
    ap.add_argument('--keep-work', action='store_true')
    ap.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=1)
    ap.add_argument('--report', type=Path)
    args = ap.parse_args()
    if args.report:
        need(not args.report.resolve().is_relative_to(ROOT), 'report outside frozen component')
    start = time.monotonic()
    result = certificate()
    with retained_workspace(args.keep_work) as td:
        work = Path(td)
        copy = work / 'anchor6'
        shutil.copytree(ROOT, copy)
        print('VERIFIED ABSOLUTE ANCHOR LOGS',(copy/'logs').resolve(),flush=True)
        (copy/'logs').mkdir(exist_ok=True)
        env = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'}
        log = work / 'audit.log'
        with log.open('w') as stream:
            subprocess.run([sys.executable, '-B', str(copy / 'code/audit_one.py'), '1', '16', '--k', '1', '--prefix', str(copy / 'results/q6')], stdout=stream, stderr=subprocess.STDOUT, check=True, env=env)
        need((copy / 'results/q6.independent.json').read_bytes() == (ROOT / 'results/q6.independent.json').read_bytes(), 'independent audit reproduced')
        if args.rebuild:
            with (work / 'rebuild.log').open('w') as stream:
                subprocess.run([sys.executable, '-B', str(copy / 'run_anchor.py'), '--compare', '--threads', str(args.threads)], stdout=stream, stderr=subprocess.STDOUT, check=True, env=env)
            replay = json.loads((copy / 'logs/replay_receipt.json').read_text())
            need(replay['status'] == 'PASS' and replay['sanitize'], 'geometric replay with UBSan')
            rebuilt = Path(replay['raw_work_directory']) / 'q6.witness'
            need(rebuilt.read_bytes() == (ROOT / 'results/q6.witness').read_bytes(), 'regenerated raw q6 byte identity')
            raw = rebuilt.read_bytes()
            negative = []
            for name, p, data, message in [('truncated', 1, raw[:-1], 'truncated witness file'), ('extra', 1, raw + b'\x00', 'extra witness data'), ('zero_nonunit', 1, b'\x00' * 8 + raw[8:], 'full mark unit'), ('wrong_mark', 13, raw, 'only the marked-6 anchor')]:
                prefix = work / ('negative_' + name)
                Path(str(prefix) + '.witness').write_bytes(data)
                got = subprocess.run([str(rebuilt.parent / 'consume'), 'consume', str(p), str(copy / 'data'), str(prefix)], capture_output=True, text=True, env=env)
                need(got.returncode != 0 and message in got.stderr, 'specific negative rejection ' + name)
                negative.append(name)
            result['negative_consumer_checks'] = negative
            shutil.rmtree(rebuilt.parent)
    result.update(independent_audit_recomputed=True, geometric_raw_replayed=args.rebuild, runtime_seconds=round(time.monotonic() - start, 3))
    if args.report:
        args.report.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(result, indent=2), flush=True)
if __name__ == '__main__':
    main()
