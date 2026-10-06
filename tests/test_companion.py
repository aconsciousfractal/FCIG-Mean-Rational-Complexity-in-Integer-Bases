"""Integrity and replay isolation checks without historical scientific sweeps."""
import ast
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
PROGRAMS = ROOT / 'companion/programs'
sys.path.insert(0, str(PROGRAMS))
import archive
import replay
import verify_retained as retained


class CompanionTests(unittest.TestCase):
    def test_required_complete_domains(self):
        expected = [(23, 'FINITE_BRIDGE.json'), (23, 'GENERAL_PRIME_BRIDGE.json'),
                    (28, 'SEMIPRIME_DOMAIN.json'), (28, 'SEMIPRIME_BRIDGE.json'),
                    (34, 'ANALYTIC_N2.json'), (37, 'ANALYTIC_N1.json')]
        expected += [(n, 'BRANCH_AGGREGATE.json') for n in (39, 40, 41)]
        for n, name in expected:
            self.assertEqual(archive.entry(n, name)['kind'], 'retained_data')

    def test_public_maps_and_adapters(self):
        entries = archive.manifest()['entries']
        self.assertEqual(len({e['locator'] for e in entries}), len(entries))
        maps = json.loads((ROOT / 'companion/portable_source_map.json').read_bytes())
        self.assertEqual({e['locator'] for e in entries if e['kind'] == 'source'}, {e['locator'] for e in maps})
        for e in maps:
            raw = (ROOT / 'companion/src' / e['portable_sha256']).read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(), e['portable_sha256'])
            ast.parse(raw)

    def test_replay_rejects_existing_directory_without_write(self):
        with tempfile.TemporaryDirectory() as directory:
            saved = replay.PACKAGE
            replay.PACKAGE = Path(directory)
            try:
                existing = replay.PACKAGE / 'replay_work/existing'
                existing.mkdir(parents=True)
                sentinel = existing / 'sentinel.txt'
                sentinel.write_bytes(b'unchanged')
                with self.assertRaises(FileExistsError):
                    replay.workspace_path('existing')
                self.assertEqual(sentinel.read_bytes(), b'unchanged')
                self.assertEqual(list(existing.iterdir()), [sentinel])
            finally:
                replay.PACKAGE = saved

    def test_replay_name_guard(self):
        for name in ('', '../outside', 'with/slash', 'with space', 'C:drive'):
            with self.assertRaises(ValueError):
                replay.workspace_path(name)

    def test_replay_rejects_symlink_parent_without_write(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / 'package'
            outside = Path(directory) / 'outside'
            package.mkdir()
            outside.mkdir()
            link = package / 'replay_work'
            try:
                link.symlink_to(outside, target_is_directory=True)
            except (OSError, NotImplementedError) as error:
                self.skipTest('directory symlink unsupported in this environment: ' + str(error))
            saved = replay.PACKAGE
            replay.PACKAGE = package
            try:
                with self.assertRaises(ValueError):
                    replay.workspace_path('fresh')
                self.assertEqual(list(outside.iterdir()), [])
                self.assertFalse((outside / 'fresh').exists())
            finally:
                replay.PACKAGE = saved

    def test_semiprime_complete_command_is_explicit(self):
        tree = ast.parse((PROGRAMS / 'reconstruct_semiprime.py').read_text())
        values = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant)]
        self.assertIn('--all', values)
        self.assertIn(4596, values)
        self.assertIn('upper_histogram_sha256', values)

    def test_prime_power_obligations_reject_duplicate_omission_and_reordering(self):
        ternary = [{'N': n, 'M': 3**n} for n in range(1, 9)]
        cases = ((5,3),(5,4),(5,5),(5,6),(7,3),(7,4),(7,5),(11,3),(11,4),
                 (13,3),(17,3),(19,3),(23,3))
        general = [{'p': p, 'N': n, 'M': p**n} for p, n in cases]
        # Histograms are deliberately absent: invalid obligations must fail first.
        for source in ('FINITE_BRIDGE.json', 'GENERAL_PRIME_BRIDGE.json'):
            for mutation in ('duplicate', 'omission', 'reorder', 'wrong_modulus', 'float'):
                payloads = {'FINITE_BRIDGE.json': copy.deepcopy(ternary),
                            'GENERAL_PRIME_BRIDGE.json': copy.deepcopy(general)}
                rows = payloads[source]
                if mutation == 'duplicate': rows[-1] = copy.deepcopy(rows[0])
                elif mutation == 'omission': rows.pop()
                elif mutation == 'reorder': rows[0], rows[1] = rows[1], rows[0]
                elif mutation == 'wrong_modulus': rows[0]['M'] += 1
                else: rows[0]['N'] = float(rows[0]['N'])
                with self.subTest(source=source, mutation=mutation), patch.object(
                        retained, 'data', side_effect=lambda batch, name: {'rows': payloads[name]}):
                    with self.assertRaisesRegex(ValueError, 'domain'):
                        retained.prime_power()

    def test_early_obligations_reject_duplicate_omission_and_reordering(self):
        fixtures = {26: [{'N': n, 'M': 6**n} for n in range(1, 9)],
                    27: [{'base': b, 'N': 6, 'M': b**6, 'radical': rad}
                         for b, rad in ((10,10),(12,6),(14,14))]}
        for batch, original in fixtures.items():
            for mutation in ('duplicate', 'omission', 'reorder'):
                rows = copy.deepcopy(original)
                if mutation == 'duplicate': rows[-1] = copy.deepcopy(rows[0])
                elif mutation == 'omission': rows.pop()
                else: rows[0], rows[1] = rows[1], rows[0]
                with self.subTest(batch=batch, mutation=mutation), patch.object(
                        retained, 'data', return_value={'rows': rows}):
                    with self.assertRaisesRegex(ValueError, 'domain'):
                        retained.early(batch, 'fixture.json', len(original))

    def test_aggregate_obligations_bind_to_distinct_domain(self):
        for batch in (29,30,31,32,33,35,39,40,41):
            q = 5 if batch == 30 else 3
            profiles = [[2**a * q, 2, q, a, 1] for a in (1, 2)]
            original = [dict(zip(('M','p','r','a','b'), values)) for values in profiles]
            if batch in (29,30):
                original = [{k: v for k, v in row.items() if k not in ('p','r')} for row in original]
            for mutation in ('duplicate', 'omission', 'reorder', 'wrong_profile'):
                rows = copy.deepcopy(original)
                if mutation == 'duplicate': rows[-1] = copy.deepcopy(rows[0])
                elif mutation == 'omission': rows.pop()
                elif mutation == 'reorder': rows[0], rows[1] = rows[1], rows[0]
                else: rows[0]['a'] += 1
                payload = {'rows': rows, 'domain_sha256': 'frozen-domain'}
                domain = {'profiles' if batch >= 39 else 'rows': profiles}
                def supplied(b, name):
                    return payload if name == 'fixture.json' else domain
                with self.subTest(batch=batch, mutation=mutation), patch.object(
                        retained, 'data', side_effect=supplied), patch.object(
                        retained, 'entry', return_value={'sha256': 'frozen-domain'}):
                    with self.assertRaisesRegex(ValueError, 'domain|aggregate count'):
                        retained.aggregate(batch, 'fixture.json', 2)

    def test_invalid_adapter_maps_never_write_or_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / 'package'
            package.mkdir()
            outside = package / '_portable_outside.py'
            outside.write_bytes(b'sentinel')
            payload = b'raise RuntimeError("must not launch")\n'
            pin = hashlib.sha256(payload).hexdigest()
            (package / 'src').mkdir()
            (package / 'src' / pin).write_bytes(payload)
            (package / 'src' / ('0' * 64)).write_bytes(b'tampered adapter')
            record = {'locator': 'results/batch_22/job.py', 'kind': 'source',
                      'storage': 'objects/' + pin, 'sha256': pin,
                      'original_sha256': pin, 'stored_sha256': pin}
            valid = {'locator': record['locator'], 'portable_sha256': pin,
                     'public_source_sha256': pin, 'original_sha256': pin}
            variants = [[], [valid, valid]]
            for name in ('../../outside.py', '/outside.py', 'C:/outside.py',
                         'results\\job.py', './results/batch_22/job.py',
                         'results//batch_22/job.py', 'results/batch_22/unknown.py'):
                variants.append([{**valid, 'locator': name}])
            for key in ('public_source_sha256', 'original_sha256', 'portable_sha256'):
                variants.append([{**valid, key: '0' * 64}])
                variants.append([{**valid, key: '../bad'}])
            for index, adapters in enumerate(variants):
                (package / 'portable_source_map.json').write_text(json.dumps(adapters), encoding='utf8')
                with self.subTest(index=index), patch.object(replay, 'PACKAGE', package), patch.object(
                        replay, 'entry', return_value=record), patch.object(
                        replay, 'manifest', return_value={'entries': [record]}), patch.object(
                        replay, 'read_object') as read, patch.object(
                        replay.subprocess, 'run') as launch, patch.object(
                        sys, 'argv', ['replay.py', '22', 'job.py', '--work', 'fresh']):
                    with self.assertRaises(ValueError):
                        replay.main()
                    read.assert_not_called()
                    launch.assert_not_called()
                    self.assertEqual(outside.read_bytes(), b'sentinel')
                    self.assertFalse((package / 'replay_work').exists())

    def test_valid_adapter_preflight_is_read_only(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            payload = b'pass\n'
            pin = hashlib.sha256(payload).hexdigest()
            (package / 'src').mkdir()
            (package / 'src' / pin).write_bytes(payload)
            record = dict(locator='results/batch_22/job.py', kind='source',
                          storage='objects/' + pin, sha256=pin, original_sha256=pin, stored_sha256=pin)
            adapter = dict(locator=record['locator'], portable_sha256=pin,
                           public_source_sha256=pin, original_sha256=pin)
            work = package / 'replay_work/fresh'
            with patch.object(replay, 'PACKAGE', package):
                objects, portable = replay.preflight(work, [record], [adapter])
            self.assertEqual(objects, [(work / record['locator'], record)])
            self.assertEqual(portable, [(work / 'results/batch_22/_portable_job.py', payload)])
            self.assertFalse(work.parent.exists())


if __name__ == '__main__':
    unittest.main()
