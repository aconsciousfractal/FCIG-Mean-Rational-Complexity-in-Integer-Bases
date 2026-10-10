"""Input guard regressions from the adversarial review of the new companion."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'companion/global_bound'
spec=importlib.util.spec_from_file_location('early_closure_guard',BASE/'early_closures/verify_early.py')
early=importlib.util.module_from_spec(spec)
spec.loader.exec_module(early)

class EarlyInputGuards(unittest.TestCase):
    def test_json_duplicate_and_nonfinite_values_rejected(self):
        for value in ('{"x":1,"x":2}','{"x":NaN}','{"x":Infinity}','{"x":-Infinity}'):
            with self.subTest(value=value),self.assertRaises(ValueError):
                early.parse(value)

    def test_case_escape_rejected_before_reading_a_sibling(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'root';root.mkdir()
            outside=Path(td)/'outside';outside.mkdir()
            data=outside/'note.txt';data.write_bytes(b'reader')
            ref=dict(case='../outside',path='note.txt',sha256=hashlib.sha256(b'reader').hexdigest())
            with self.assertRaisesRegex(ValueError,'declared case'):
                early.source(root,ref)
            self.assertEqual(data.read_bytes(),b'reader')

    def test_noncanonical_and_parent_member_paths_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);(root/'reader.txt').write_bytes(b'reader')
            for name in ('../reader.txt','/reader.txt','C:/reader.txt','a\\reader.txt',
                         './reader.txt','a//reader.txt'):
                with self.subTest(name=name),self.assertRaises(ValueError):
                    early.member(root,name)

    def test_complete_manifest_rejects_duplicate_missing_extra_and_tampered_inputs(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);data=root/'reader.txt';data.write_bytes(b'reader')
            row=hashlib.sha256(b'reader').hexdigest()+'  reader.txt\n'
            manifest=root/'SHA256SUMS';manifest.write_text(row,encoding='utf8')
            self.assertEqual(early.inventory(root),1)
            manifest.write_text(row+row,encoding='utf8')
            with self.assertRaisesRegex(ValueError,'duplicate'):early.inventory(root)
            manifest.write_text(row,encoding='utf8');data.unlink()
            with self.assertRaises(ValueError):early.inventory(root)
            data.write_bytes(b'reader');extra=root/'extra';extra.write_bytes(b'x')
            with self.assertRaisesRegex(ValueError,'inventory'):early.inventory(root)
            extra.unlink();data.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'digest'):early.inventory(root)

    def test_valid_retained_contract_with_exact_domains(self):
        result=early.verify()
        self.assertEqual(result['three_prime_partition'],[353,89])
        self.assertEqual(result['three_prime_finite_moduli'],11911)
        self.assertFalse(result['full_geometric_reconstruction_in_this_run'])

if __name__=='__main__':unittest.main()

