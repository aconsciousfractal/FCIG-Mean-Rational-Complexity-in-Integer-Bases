"""Bounded reconstruction/cache/domain/output controls, not the full sweep."""
from collections import Counter
import copy
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'companion/programs'))
import reconstruct_semiprime as semiprime


class FullSemiprimeTests(unittest.TestCase):
    def test_exact_log_cache_matches_original_and_reuses_values(self):
        semiprime._cached_log_ticks.cache_clear()
        original = semiprime._uncached_log_ticks
        with patch.object(semiprime, '_uncached_log_ticks', wraps=original) as raw:
            for n in (1, 2, 3, 6, 36, 1296, 19599):
                expected = original(n)
                self.assertEqual(semiprime.log_ticks(n), expected)
                self.assertEqual(semiprime.log_ticks(n), expected)
            self.assertEqual(raw.call_count, 7)
        semiprime._cached_log_ticks.cache_clear()

    def test_log_cache_rejects_inexact_and_invalid_inputs_before_lookup(self):
        with patch.object(semiprime, '_cached_log_ticks') as cache:
            for value in (0, -1, True, 1.0, Fraction(1), '1'):
                with self.subTest(value=value), self.assertRaises(ValueError):
                    semiprime.log_ticks(value)
            cache.assert_not_called()

    def test_small_capped_and_far_envelopes_dominate_defining_minima(self):
        for M,p,q in ((6,2,3), (22,2,11)):
            with self.subTest(M=M):
                row = semiprime.reconstruct(M,p,q)
                self.assertEqual(sum(c for h,c in row['upper_histogram']), M)
                ledger = {n:(L,U) for n,L,U in row['log_tick_ledger']}
                self.assertEqual(set(ledger), {M} | {h for h,c in row['upper_histogram']})
                self.assertEqual(2*sum(c*ledger[h][1] for h,c in row['upper_histogram'])-M*ledger[M][0],
                                 row['upper_centered_log_numerator'])
                self.assertLess(row['upper_centered_log_numerator'], 0)
                minimum = [semiprime.core.defining_minimum(M,s,M) for s in range(M)]
                if q <= semiprime.isqrt(2*M)+1:
                    exact,_,_ = semiprime.core.central(M,M)
                    upper = [h or min(q,max(1,min(s,M-s))) for s,h in enumerate(exact)]
                else:
                    upper,_,_,_ = semiprime.core.envelope(M,M)
                    upper = [min(h,q) for h in upper]
                self.assertTrue(all(h >= r for h,r in zip(upper, minimum)))
                self.assertEqual(row['upper_histogram'], [[h,c] for h,c in sorted(Counter(upper).items())])

    def test_distinct_complete_domain_rejects_count_order_duplicates_and_types(self):
        rows = [{'M':6,'p':2,'r':3}, {'M':10,'p':2,'r':5}, {'M':15,'p':3,'r':5}]
        bridge = {'domain_sha256':'bound-domain', 'rows':copy.deepcopy(rows)}
        domain = {'rows':copy.deepcopy(rows)}
        with patch.object(semiprime, 'FULL_COUNT', 3):
            self.assertEqual(semiprime.complete_domain(bridge,domain,'bound-domain'), rows)
            for target in ('bridge','domain'):
                for mutation in ('missing','duplicate','reorder','bool','float','wrong_product','changed_prime'):
                    b,d = copy.deepcopy(bridge),copy.deepcopy(domain)
                    altered = (b if target=='bridge' else d)['rows']
                    if mutation=='missing':altered.pop()
                    elif mutation=='duplicate':altered[-1]=copy.deepcopy(altered[0])
                    elif mutation=='reorder':altered.reverse()
                    elif mutation=='bool':altered[0]['p']=True
                    elif mutation=='float':altered[0]['M']=6.0
                    elif mutation=='wrong_product':altered[0]['M']=7
                    else:altered[-1]={'M':14,'p':2,'r':7}
                    with self.subTest(target=target,mutation=mutation), self.assertRaises(ValueError):
                        semiprime.complete_domain(b,d,'bound-domain')
            with self.assertRaisesRegex(ValueError,'byte binding'):
                semiprime.complete_domain(bridge,domain,'different-domain')

    def test_each_historical_digest_and_Q_field_is_checked(self):
        row={'M':6,'upper_histogram_sha256':'digest','upper_centered_log_numerator':-2,
             'upper_centered_log_denominator':12}
        semiprime.compare_record(row,copy.deepcopy(row))
        for key in ('upper_histogram_sha256','upper_centered_log_numerator','upper_centered_log_denominator'):
            expected=copy.deepcopy(row)
            expected[key]='other' if key=='upper_histogram_sha256' else expected[key]+1
            with self.subTest(key=key), self.assertRaisesRegex(ValueError,'reconstruction mismatch'):
                semiprime.compare_record(row,expected)

    def test_streamed_artifact_matches_legacy_json_and_compact_identity(self):
        rows=[{'M':6,'p':2,'r':3,'upper_histogram':[[1,5],[2,1]]}, {'M':10,'p':2,'r':5}]
        scope='fixture'
        expected=(json.dumps(semiprime.result(rows,scope),sort_keys=True,separators=(',',':'))+'\n').encode()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'complete.json'
            summary=semiprime.export_file(path,iter(rows),scope,len(rows))
            self.assertEqual(path.read_bytes(),expected)
            self.assertEqual(summary['output_sha256'],hashlib.sha256(expected).hexdigest())
            self.assertEqual(summary['output_bytes'],len(expected))
            self.assertEqual(summary['reconstructed_rows'],len(rows))
            self.assertNotIn('rows',summary)

    def test_existing_output_is_rejected_before_consuming_rows(self):
        def forbidden_rows():
            raise AssertionError('rows must not be consumed')
            yield None
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'sentinel.json'
            path.write_bytes(b'unchanged')
            with self.assertRaises(FileExistsError):
                semiprime.export_file(path,forbidden_rows(),'full',1)
            self.assertEqual(path.read_bytes(),b'unchanged')

    def test_failed_stream_retains_incomplete_artifact_without_success(self):
        def failed_rows():
            yield {'M':6}
            raise ValueError('fixture reconstruction failure')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'partial.json'
            with self.assertRaisesRegex(ValueError,'fixture reconstruction failure'):
                semiprime.export_file(path,failed_rows(),'full',2)
            with self.assertRaises(json.JSONDecodeError):
                json.loads(path.read_bytes())
            self.assertNotIn(b'"status"',path.read_bytes())
            with self.assertRaisesRegex(ValueError,'missing rows'):
                semiprime.export_file(Path(directory)/'missing.json',iter([{'M':6}]),'full',2)
            with self.assertRaisesRegex(ValueError,'excess rows'):
                semiprime.export_file(Path(directory)/'excess.json',iter([{'M':6},{'M':10}]),'full',1)

    def test_progress_is_on_stderr_after_historical_comparison(self):
        expected=[{'M':6,'p':2,'r':3,'upper_histogram_sha256':'d',
                   'upper_centered_log_numerator':-2,'upper_centered_log_denominator':12}]
        stderr=io.StringIO()
        with patch.object(semiprime,'reconstruct',return_value=copy.deepcopy(expected[0])), patch.object(sys,'stderr',stderr):
            self.assertEqual(list(semiprime.reconstructed_rows(expected,progress=True)),expected)
        self.assertEqual(stderr.getvalue(),'PROGRESS 1/1\n')


if __name__=='__main__':unittest.main()
