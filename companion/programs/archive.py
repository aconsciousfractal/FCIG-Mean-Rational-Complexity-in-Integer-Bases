"""Verified relative access to the public scientific objects."""
import gzip
import hashlib
import json
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]


def manifest():
    return json.loads((PACKAGE / 'provenance.json').read_bytes())


def entry(batch, name, area=None):
    matches = [e for e in manifest()['entries']
               if e['batch'] == batch and e['locator'].endswith('/' + name)
               and (area is None or e['locator'].startswith(area + '/'))]
    if len(matches) != 1:
        raise ValueError('ambiguous or missing public object: ' + str((batch, name)))
    return matches[0]


def read_object(record):
    path = (PACKAGE / record['storage']).resolve()
    if not path.is_relative_to((PACKAGE / 'objects').resolve()):
        raise ValueError('object path escapes public store')
    stored = path.read_bytes()
    if len(stored) != record['stored_bytes'] or hashlib.sha256(stored).hexdigest() != record['stored_sha256']:
        raise ValueError('stored object byte identity mismatch')
    payload = gzip.decompress(stored) if path.suffix == '.gz' else stored
    if len(payload) != record['bytes'] or hashlib.sha256(payload).hexdigest() != record['sha256']:
        raise ValueError('uncompressed scientific byte identity mismatch')
    return payload


def raw(batch, name, area=None):
    return read_object(entry(batch, name, area))


def data(batch, name, area=None):
    return json.loads(raw(batch, name, area))
