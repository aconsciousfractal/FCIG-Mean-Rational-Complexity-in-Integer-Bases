"""Package integrity and metadata controls; no mathematical source imports."""
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('package_manifest_checker', ROOT / 'scripts/check_manifest.py')
manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(manifest)
replay_spec = importlib.util.spec_from_file_location('replay_output_checker', ROOT / 'scripts/check_replay_output.py')
replay_output = importlib.util.module_from_spec(replay_spec)
replay_spec.loader.exec_module(replay_output)


class ReplayJSONTests(unittest.TestCase):
    def test_valid_json_preserves_canonical_identity(self):
        value = {'status': 'PASS', 'rows': [{'M': 6, 'count': 3}]}
        expected = hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                             ensure_ascii=True, allow_nan=False).encode('utf8')).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'output.json'
            path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf8')
            self.assertEqual(replay_output.canonical_digest(path), expected)

    def test_ambiguous_and_nonfinite_outputs_fail(self):
        corrupt = ('{"status":"FAIL","status":"PASS"}',
                   '{"rows":[{"M":900,"M":6}]}',
                   '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}',
                   '{"x":', '{} {}')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'output.json'
            for content in corrupt:
                path.write_text(content, encoding='utf8')
                with self.subTest(content=content), self.assertRaises(ValueError):
                    replay_output.canonical_digest(path)

    def test_expected_identity_file_rejects_duplicate_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'companion').mkdir()
            (root / 'companion/expected_replay.json').write_text(
                '{"outputs":{},"outputs":{"smoke":{"canonical_json_sha256":"ignored"}}}',
                encoding='utf8')
            with patch.object(replay_output, 'ROOT', root), patch.object(
                    replay_output.argparse.ArgumentParser, 'parse_args',
                    return_value=type('Args', (), {'name': 'smoke', 'outputs': []})()):
                with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
                    replay_output.main()


class ManifestTests(unittest.TestCase):
    def fixture(self, root):
        (root / 'reader.txt').write_bytes(b'reader')
        row = hashlib.sha256(b'reader').hexdigest() + '  reader.txt\n'
        (root / manifest.MANIFEST).write_text(row, encoding='utf-8')
        return row

    def test_valid_tree_with_local_generated_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            for name in ('.git/config', 'build/replay.json', 'companion/replay_work/job/result.json', '__pycache__/cache.pyc'):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b'local generated content')
            self.assertEqual(manifest.verify(root)['status'], 'MANIFEST_PASS')

    def test_empty_and_malformed_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = self.fixture(root)
            for content in ('', 'bad\n', row.replace('  reader', ' reader'), row.upper()):
                with self.subTest(content=content):
                    (root / manifest.MANIFEST).write_text(content, encoding='utf-8')
                    with self.assertRaises(ValueError):
                        manifest.verify(root)

    def test_duplicate_row(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = self.fixture(root)
            (root / manifest.MANIFEST).write_text(row + row, encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                manifest.verify(root)

    def test_unsafe_and_self_manifest_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            for name in ('../reader.txt', '/reader.txt', 'C:/reader.txt', 'folder\\reader.txt',
                         './reader.txt', 'folder//reader.txt', manifest.MANIFEST):
                with self.subTest(name=name):
                    (root / manifest.MANIFEST).write_text('0' * 64 + '  ' + name + '\n', encoding='utf-8')
                    with self.assertRaisesRegex(ValueError, 'unsafe'):
                        manifest.verify(root)

    def test_generated_paths_cannot_be_manifest_entries(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            for name in ('.git/config', 'build/result.json', 'companion/replay_work/job/result.json', '.venv/source.py'):
                with self.subTest(name=name):
                    (root / manifest.MANIFEST).write_text('0' * 64 + '  ' + name + '\n', encoding='utf-8')
                    with self.assertRaisesRegex(ValueError, 'excluded'):
                        manifest.verify(root)

    def test_missing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            (root / 'reader.txt').unlink()
            with self.assertRaisesRegex(ValueError, 'file-set mismatch.*missing'):
                manifest.verify(root)

    def test_tampered_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            (root / 'reader.txt').write_bytes(b'tampered')
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                manifest.verify(root)

    def test_unexpected_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            (root / 'extra.txt').write_bytes(b'extra')
            with self.assertRaisesRegex(ValueError, 'file-set mismatch.*unexpected'):
                manifest.verify(root)

    def test_symlink_file_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'reader_tree'
            root.mkdir()
            self.fixture(root)
            outside = Path(directory) / 'outside.txt'
            outside.write_bytes(b'reader')
            reader = root / 'reader.txt'
            reader.unlink()
            try:
                reader.symlink_to(outside)
            except (OSError, NotImplementedError) as error:
                self.skipTest('file symlink unsupported in this environment: ' + str(error))
            with self.assertRaisesRegex(ValueError, 'symlink'):
                manifest.verify(root)
            self.assertEqual(outside.read_bytes(), b'reader')


class PackageIdentityTests(unittest.TestCase):
    def manifest_rows(self):
        rows = {}
        for line in (ROOT / manifest.MANIFEST).read_text(encoding='utf-8').splitlines():
            digest, name = line.split('  ', 1)
            self.assertNotIn(name, rows)
            rows[name] = digest
        return rows

    def test_version_citation_and_changelog_agree(self):
        version = (ROOT / 'VERSION').read_text(encoding='utf-8').strip()
        self.assertRegex(version, r'^\d+\.\d+\.\d+$')
        citation = (ROOT / 'CITATION.cff').read_text(encoding='utf-8')
        versions = re.findall(r'^version:\s*([^\r\n]+)$', citation, flags=re.MULTILINE)
        self.assertEqual(versions, [version])
        self.assertIn('companion version ' + version, citation)
        self.assertIn('## ' + version + '\n', (ROOT / 'CHANGELOG.md').read_text(encoding='utf-8'))

    def test_license_scopes_and_citation_license_agree(self):
        license_text = (ROOT / 'LICENSE').read_text(encoding='utf-8')
        self.assertTrue(license_text.startswith('MIT License\n'))
        manuscript = (ROOT / 'LICENSE_MANUSCRIPT.md').read_text(encoding='utf-8')
        self.assertIn('SPDX-License-Identifier: CC-BY-4.0', manuscript)
        self.assertIn('https://creativecommons.org/licenses/by/4.0/', manuscript)
        scope = (ROOT / 'LICENSE_SCOPE.md').read_text(encoding='utf-8')
        self.assertIn('[LICENSE](LICENSE)', scope)
        self.assertIn('[LICENSE_MANUSCRIPT.md](LICENSE_MANUSCRIPT.md)', scope)
        self.assertIn('paper/', scope)
        self.assertIn('certificates/data', scope)
        citation = (ROOT / 'CITATION.cff').read_text(encoding='utf-8')
        self.assertEqual(re.findall(r'^license:\s*([^\r\n]+)$', citation, flags=re.MULTILINE), ['MIT'])
        self.assertIn('Original manuscript CC BY 4.0', citation)

    def test_companion_objects_match_package_manifest_identities(self):
        rows = self.manifest_rows()
        provenance = json.loads((ROOT / 'companion/provenance.json').read_bytes())
        self.assertEqual(provenance['schema'], 'integer-base-public-provenance-v1')
        locators = set()
        for record in provenance['entries']:
            locator = record['locator']
            self.assertNotIn(locator, locators)
            locators.add(locator)
            for name in (locator, record['origin_locator'], record['storage']):
                path = PurePosixPath(name)
                self.assertFalse(path.is_absolute())
                self.assertNotIn('..', path.parts)
                self.assertNotIn('\\', name)
                self.assertNotIn(':', name)
            storage = record['storage']
            self.assertTrue(storage.startswith('objects/'))
            self.assertEqual(PurePosixPath(storage).stem, record['sha256'])
            self.assertRegex(record['original_sha256'], r'^[0-9a-f]{64}$')
            self.assertRegex(record['sha256'], r'^[0-9a-f]{64}$')
            self.assertEqual(rows['companion/' + storage], record['stored_sha256'])
            self.assertEqual((ROOT / 'companion' / storage).stat().st_size, record['stored_bytes'])
            if not storage.endswith('.gz'):
                self.assertEqual(record['sha256'], record['stored_sha256'])
                self.assertEqual(record['bytes'], record['stored_bytes'])
        self.assertTrue(locators)
        complete = json.loads((ROOT / 'companion/reconstruction/index.json').read_bytes())
        records = provenance['entries'] + list(complete['records'].values())
        for record in complete['records'].values():
            self.assertEqual(rows['companion/' + record['storage']], record['stored_sha256'])
        self.assertEqual({name for name in rows if name.startswith('companion/objects/')},
                         {'companion/' + record['storage'] for record in records})

    def test_portable_map_source_and_adapter_bindings(self):
        rows = self.manifest_rows()
        provenance = json.loads((ROOT / 'companion/provenance.json').read_bytes())['entries']
        sources = {record['locator']: record for record in provenance if record['kind'] == 'source'}
        adapters = json.loads((ROOT / 'companion/portable_source_map.json').read_bytes())
        self.assertEqual(set(sources), {record['locator'] for record in adapters})
        self.assertEqual(len(adapters), len(sources))
        for adapter in adapters:
            source = sources[adapter['locator']]
            self.assertEqual(adapter['original_sha256'], source['original_sha256'])
            self.assertEqual(adapter['public_source_sha256'], source['sha256'])
            self.assertRegex(adapter['portable_sha256'], r'^[0-9a-f]{64}$')
            name = 'companion/src/' + adapter['portable_sha256']
            self.assertEqual(rows[name], adapter['portable_sha256'])
            self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(), adapter['portable_sha256'])
        self.assertEqual({name for name in rows if name.startswith('companion/src/')},
                         {'companion/src/' + adapter['portable_sha256'] for adapter in adapters})


if __name__ == '__main__':
    unittest.main()
