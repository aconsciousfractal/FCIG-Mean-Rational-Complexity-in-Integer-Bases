"""Explicit optional scientific replay in a fresh package-owned directory.

Usage: python programs/replay.py BATCH SCRIPT --work NAME [--optimized]
Read the reproduction guide before choosing a potentially expensive job.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from archive import PACKAGE, entry, manifest, read_object


def workspace_path(name):
    if not name or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in name):
        raise ValueError('work name must be a simple package-local directory name')
    parent = PACKAGE / 'replay_work'
    if parent.is_symlink():
        raise ValueError('replay work parent must not be a symlink')
    path = parent / name
    if not path.resolve().is_relative_to(PACKAGE.resolve()):
        raise ValueError('replay work directory escapes the package')
    if path.exists() or path.is_symlink():
        raise FileExistsError('replay requires a fresh work directory; already exists: ' + name)
    return path


def safe_locator(name):
    if type(name) is not str:
        raise ValueError('locator must be a relative POSIX string')
    rel = PurePosixPath(name)
    if (not name or rel.is_absolute() or '..' in rel.parts or '\\' in name
            or ':' in name or str(rel) != name or name == '.'):
        raise ValueError('unsafe or nonnormalized locator: ' + name)
    return rel


def digest(value):
    if type(value) is not str or re.fullmatch('[0-9a-f]{64}', value) is None:
        raise ValueError('exact SHA256 hexadecimal digest required')
    return value


def preflight(workspace, records, adapters):
    """Validate the complete source map before any work-tree write."""
    if type(records) is not list or type(adapters) is not list:
        raise ValueError('provenance and adapter lists required')
    root = workspace.resolve()
    sources = {}
    locators = set()
    destinations = set()
    objects = []
    for record in records:
        locator = str(safe_locator(record['locator']))
        if locator in locators:
            raise ValueError('duplicate provenance locator')
        locators.add(locator)
        for key in ('sha256', 'original_sha256', 'stored_sha256'):
            digest(record[key])
        storage = safe_locator(record['storage'])
        if not storage.parts or storage.parts[0] != 'objects':
            raise ValueError('public object storage required')
        target = (workspace / locator).resolve()
        if not target.is_relative_to(root) or target in destinations:
            raise ValueError('locator escapes replay workspace')
        destinations.add(target)
        objects.append((target, record))
        if record['kind'] == 'source':
            sources[locator] = record
    portable = []
    mapped = set()
    for adapter in adapters:
        locator = str(safe_locator(adapter['locator']))
        if locator in mapped or locator not in sources:
            raise ValueError('adapter map must be bijective with provenance sources')
        mapped.add(locator)
        for key in ('portable_sha256', 'public_source_sha256', 'original_sha256'):
            digest(adapter[key])
        source = sources[locator]
        if (adapter['public_source_sha256'] != source['sha256']
                or adapter['original_sha256'] != source['original_sha256']):
            raise ValueError('adapter source identity mismatch')
        original = workspace / locator
        target = original.with_name('_portable_' + original.name).resolve()
        if not target.is_relative_to(root) or target in destinations:
            raise ValueError('unsafe or colliding adapter destination')
        destinations.add(target)
        payload = (PACKAGE / 'src' / adapter['portable_sha256']).read_bytes()
        if hashlib.sha256(payload).hexdigest() != adapter['portable_sha256']:
            raise ValueError('portable source byte identity mismatch')
        portable.append((target, payload))
    if mapped != set(sources):
        raise ValueError('adapter map must be bijective with provenance sources')
    return objects, portable


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('batch', type=int)
    ap.add_argument('script')
    ap.add_argument('--optimized', action='store_true')
    ap.add_argument('--work', required=True)
    a = ap.parse_args()
    workspace = workspace_path(a.work)  # Reject reuse before any file write or scientific launch.
    source = entry(a.batch, a.script, 'reports' if a.batch == 24 else 'results')
    if source['kind'] != 'source':
        raise ValueError('scientific source job required')
    if os.name == 'nt':
        workspace = Path('\\\\?\\' + str(workspace.resolve()))
    records = manifest()['entries']
    adapters = json.loads((PACKAGE / 'portable_source_map.json').read_bytes())
    objects, portable = preflight(workspace, records, adapters)
    workspace.mkdir(parents=True, exist_ok=False)
    for target, e in objects:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(read_object(e))
    for target, payload in portable:
        target.write_bytes(payload)
    original = workspace / source['locator']
    job = original.with_name('_portable_' + original.name)
    argv = [sys.executable, '-B'] + (['-O'] if a.optimized else []) + [str(job)]
    if a.batch in (39, 40, 41) and a.script in ('run_branch.py', 'check_branch.py', 'verify_local_controls.py'):
        name = 'BRANCH_EXECUTION_DOMAIN.json' if a.batch == 39 else 'BRANCH_DOMAIN.json'
        domain = entry(a.batch, name)
        argv += [str(original.parent / name), domain['sha256']]
        if a.script == 'check_branch.py':
            certificate = entry(a.batch, 'BRANCH_AGGREGATE.json')
            argv += [str(original.parent / 'BRANCH_AGGREGATE.json'), certificate['sha256']]
        elif a.script == 'verify_local_controls.py':
            core = entry(a.batch, 'aggregate_core.py')
            argv += [str(original.parent / 'aggregate_core.py'), core['sha256']]
    elif a.batch == 39 and a.script == 'qualify_oracle.py':
        core = entry(39, 'compressed_oracle.py')
        argv += [str(original.parent / 'compressed_oracle.py'), core['sha256']]
    raise SystemExit(subprocess.run(argv, cwd=workspace, check=False).returncode)


if __name__ == '__main__':
    main()
