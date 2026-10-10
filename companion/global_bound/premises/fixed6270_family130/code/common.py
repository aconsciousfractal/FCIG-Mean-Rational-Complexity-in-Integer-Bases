"""Small transparent helpers for the 6270/130p checkpoint."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys


def need(test, message):
    if not test:
        raise AssertionError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def canonical(obj):
    return (json.dumps(obj, sort_keys=True, indent=2) + '\n').encode()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, Path(path))
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def manifest(root):
    root = Path(root).resolve()
    path = root/'MANIFEST_SHA256.json'
    obj = read(path)
    rows = obj.get('files', obj)
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*')
              if p.is_file() and p != path}
    need(actual == set(rows), 'complete immutable inventory '+str(root))
    for rel, meta in rows.items():
        p = root/rel
        need(p.resolve().is_relative_to(root), 'contained manifest path')
        if isinstance(meta, str):
            need(sha(p) == meta, 'immutable bytes '+rel)
        else:
            need(p.stat().st_size == meta['bytes'] and sha(p) == meta['sha256'],
                 'immutable bytes and size '+rel)
    return len(rows)


def check_imports(root):
    root = Path(root)
    bindings = read(root/'imports/SOURCE_BINDINGS.json')
    source = root/'imports/SOURCE_MANIFEST_6090_165P.json'
    expected = '4c6313f3d6c29b90cca7dd739f9df03a5e4da08c08b81086ad106b35bb1b026a'
    need(sha(source) == expected == bindings['source_manifest_sha256'],
         'complete preceding source manifest identity')
    need(bindings['source_checkpoint'] == 'P72_6090_165P_20261009'
         and bindings['status'] == 'PASS' and not bindings['old_large_raw_replayed_this_task'],
         'historical sources have declared scope')
    mf = read(source)['files']
    rows = bindings['files']
    need(len({r['copy'] for r in rows}) == len(rows), 'unique import destinations')
    for row in rows:
        p = root/row['copy']
        need(p.resolve().is_relative_to(root), 'contained import destination')
        meta = dict(bytes=p.stat().st_size,sha256=sha(p))
        need(meta == dict(bytes=row['bytes'],sha256=row['sha256']), 'import bytes preserved')
        if row['source'] != 'MANIFEST_SHA256.json':
            need(mf[row['source']] == meta, 'direct authentication by preceding root manifest')
        else:
            need(meta['sha256'] == expected, 'manifest import identity')
    return rows
