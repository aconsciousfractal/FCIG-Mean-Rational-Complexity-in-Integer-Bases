"""Bounded defining/interval qualification and rejection tests for central audit."""
import copy
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

PROGRAMS = Path(__file__).resolve().parents[1] / 'companion/programs'
sys.path.insert(0, str(PROGRAMS))
import verify_all_central as checker
import capture_central as capture
from central_oracle import central_counts


def defining_height(M, s):
    best = max(1, min(s, M-s))
    for d in range(1, best+1):
        if gcd(d, M) == 1:
            a = (s*d) % M
            best = min(best, max(d, min(a, M-a)))
    return best


def fixture(M):
    actual = checker.enumerate_central(M)
    fs = checker.factorization(M)
    (p, a), (q, b) = fs
    k = isqrt(M)
    histogram = []
    previous = 0
    for T in range(1, k+1):
        count = actual['cdf'][T]
        if count > previous:
            histogram.append([T, count-previous])
        previous = count
    if previous < M:
        histogram.append([k+1, M-previous])
    row = {'central_k': k, 'exact_central_states': actual['cdf'][k],
           'upper_histogram': histogram,
           'central_stratum_totals': [[g, values[k]] for g, values in actual['strata']]}
    arrays = [[g, values] for g, values in actual['strata'] if g <= k]
    return [M,p,q,a,b], row, actual, arrays


class AllCentralTests(unittest.TestCase):
    def test_definition_all_thresholds_and_strata_fixed_panel(self):
        # Prime/power/zero/small midpoint guards plus both two-prime parities.
        for M in (2,3,4,6,9,10,12,15,21,25,35,44,63,77,100,175,221,352,900):
            actual = checker.enumerate_central(M)
            heights = [defining_height(M,s) for s in range(M)]
            for T in range(1,isqrt(M)+1):
                self.assertEqual(actual['cdf'][T], sum(h<=T for h in heights), (M,T))
                for g, values in actual['strata']:
                    self.assertEqual(values[T], sum(h<=T and gcd(s,M)==g
                        for s,h in enumerate(heights) if s != 0), (M,T,g))
            self.assertEqual(actual['zero_count'],1)
            self.assertEqual(actual['cdf'][0],0)
            self.assertEqual(actual['counters']['signed_witness_visits'],
                             2*actual['counters']['positive_pair_visits'])

    def test_distinct_interval_oracle_all_thresholds_fixed_panel(self):
        for M in (6,10,12,15,18,21,35,44,63,77,100,175,221,352):
            (p,a),(q,b) = checker.factorization(M)
            actual = checker.enumerate_central(M)
            for T in range(1,isqrt(M)+1):
                interval = central_counts(M,p,q,a,b,T)
                self.assertEqual(actual['cdf'][T],interval['cdf'],(M,T))
                self.assertEqual([[g,values[T]] for g,values in actual['strata']],
                                 interval['strata'],(M,T))

    def test_square_endpoint_and_nonprimitive_witness(self):
        actual = checker.enumerate_central(100)
        # (a,d)=(2,2) is illegal here, but (6,3) is legal and nonprimitive.
        # State2 has the legal shorter witness (2,1); content is not a state.
        self.assertEqual(actual['cdf'][10],sum(defining_height(100,s)<=10 for s in range(100)))
        self.assertEqual(checker.enumerate_central(35)['counters']['positive_pair_visits'],20)
        self.assertEqual(checker.enumerate_central(2)['cdf'],[0,2])

    def test_selected_large_interval_boundaries(self):
        # A separate fixed qualifier domain, not a sample promoted to full audit.
        for M in (194481,359375,1000000):
            (p,a),(q,b) = checker.factorization(M)
            actual=checker.enumerate_central(M)
            for T in sorted({isqrt(M//2),isqrt(M)}):
                interval=central_counts(M,p,q,a,b,T)
                self.assertEqual(actual['cdf'][T],interval['cdf'],(M,T))
                self.assertEqual([[g,v[T]] for g,v in actual['strata']],interval['strata'],(M,T))

    def test_bad_modulus_rejected(self):
        for M in (True,False,0,1,-2,6.0,2**32):
            with self.subTest(M=M),self.assertRaises(ValueError):
                checker.enumerate_central(M)

    def test_consumed_arrays_and_histogram_checked_at_every_threshold(self):
        profile,row,actual,arrays = fixture(352)
        self.assertEqual(checker.compare_row(profile,row,actual,arrays)['status'],
                         'ALL_CENTRAL_THRESHOLDS_AND_STRATA_PASS')
        for mutation in ('interior','missing','reordered','float','global_histogram'):
            r = copy.deepcopy(row);v=copy.deepcopy(arrays)
            if mutation=='interior':v[0][1][2]+=1
            elif mutation=='missing':v.pop()
            elif mutation=='reordered':v.reverse()
            elif mutation=='float':v[0][1][1]=float(v[0][1][1])
            else:
                r['upper_histogram'][0][1]+=1
                r['upper_histogram'][-1][1]-=1
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                checker.compare_row(profile,r,actual,v)

    def test_omitted_strata_are_explicitly_forced_zero(self):
        profile,row,actual,arrays = fixture(352)
        self.assertTrue(any(g>actual['k'] for g,v in actual['strata']))
        bad=copy.deepcopy(actual)
        next(v for g,v in bad['strata'] if g>bad['k'])[1]=1
        with self.assertRaises(ValueError):
            checker.compare_row(profile,row,bad,arrays)

    def test_fresh_output_never_overwrites(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'existing.json';path.write_bytes(b'sentinel')
            with self.assertRaises(ValueError):checker.fresh_output(path)
            self.assertEqual(path.read_bytes(),b'sentinel')
            with checker.fresh_output(Path(directory)/'fresh.json') as stream:
                stream.write('fixture')
            self.assertEqual((Path(directory)/'fresh.json').read_text(),'fixture')

    def test_ambiguous_json_rejected(self):
        for raw in (b'{"x":1,"x":2}',b'{"x":{"a":1,"a":2}}',b'{"x":NaN}'):
            with self.assertRaises(ValueError):checker.strict_load(raw)

    def test_capture_binding_allows_controls_but_not_missing_or_duplicate_main(self):
        with tempfile.TemporaryDirectory() as directory:
            package=Path(directory);(package/'src').mkdir()
            source=b'producer identity fixture';adapter=b'adapter identity fixture'
            source_pin=hashlib.sha256(source).hexdigest();adapter_pin=hashlib.sha256(adapter).hexdigest()
            record={'locator':'results/batch_29/verify_aggregate.py','sha256':source_pin}
            (package/'src'/adapter_pin).write_bytes(adapter)
            (package/'portable_source_map.json').write_text(json.dumps([
                {'locator':record['locator'],'public_source_sha256':source_pin,'portable_sha256':adapter_pin}]))
            required={'M':6,'p':2,'r':3,'central_k':2,'central_strata':[]}
            control={'M':12,'p':2,'r':3,'central_k':3,'central_strata':[]}
            for mutation in ('valid','missing','duplicate','wrong_adapter','unrelated_extra'):
                document={'status':'PRODUCER_CENTRAL_ARRAY_CAPTURE_COMPLETE','batch':29,
                          'producer_adapter_sha256':adapter_pin,'rows':[required,control]}
                if mutation=='missing':document['rows']=[control]
                elif mutation=='duplicate':document['rows']=[required,required]
                elif mutation=='wrong_adapter':document['producer_adapter_sha256']='0'*64
                elif mutation=='unrelated_extra':document['rows'].append(
                    {'M':18,'p':2,'r':3,'central_k':4,'central_strata':[]})
                raw=json.dumps(document).encode();path=package/'capture.json';path.write_bytes(raw)
                with self.subTest(mutation=mutation),patch.object(checker,'PACKAGE',package),\
                     patch.object(checker,'entry',return_value=record),\
                     patch.object(checker,'read_object',return_value=source):
                    if mutation=='valid':
                        rows,pin=checker.load_producer_arrays(29,path,hashlib.sha256(raw).hexdigest(),[[6,2,3,1,1],[12,2,3,2,1]])
                        self.assertEqual(rows,[required,control])
                    else:
                        with self.assertRaises(ValueError):
                            checker.load_producer_arrays(29,path,hashlib.sha256(raw).hexdigest(),[[6,2,3,1,1],[12,2,3,2,1]])

    def test_consumed_validation_domain_is_source_literal_and_exact(self):
        source=b'VALIDATION = (12,18)\n'
        rows=[{'M':12},{'M':18}]
        self.assertEqual(checker.bind_validation(29,{'validation':rows},source),
                         [[12,2,3,2,1],[18,2,3,1,2]])
        for bad in ([rows[0]],rows[::-1],[rows[0],rows[0]],[{'M':12.0},rows[1]]):
            with self.assertRaises(ValueError):
                checker.bind_validation(29,{'validation':bad},source)
        with self.assertRaises(ValueError):
            checker.bind_validation(39,{'validation':[{'M':12}]})

    def test_capture_wrapper_preserves_return_objects_without_scientific_source(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'fixture.py';source.write_bytes(b'fixture identity only')
            output=Path(directory)/'capture.json'
            returned=({1:[0,1,2]}, {'sentinel':'original statistics'})
            module=SimpleNamespace(central=lambda *args:returned)
            seen=[]
            module.main=lambda:seen.append(module.central(6,object()))
            spec=SimpleNamespace(loader=SimpleNamespace(exec_module=lambda mod:None))
            with patch.object(capture.importlib.util,'spec_from_file_location',return_value=spec),\
                 patch.object(capture.importlib.util,'module_from_spec',return_value=module),\
                 patch.object(sys,'argv',['capture_central.py',str(source),
                     hashlib.sha256(source.read_bytes()).hexdigest(),'29',str(output)]):
                capture.main()
            self.assertIs(seen[0],returned)
            self.assertEqual(returned,({1:[0,1,2]},{'sentinel':'original statistics'}))
            self.assertEqual(json.loads(output.read_text())['rows'][0]['central_strata'],[[1,[0,1,2]]])

    def test_complete_domain_identity_rejects_omission_duplicate_and_reordering(self):
        # Metadata-only fixture: no central enumeration of these moduli.
        profiles=[[2**a*3**b,2,3,a,b] for a in range(1,12) for b in range(1,14)]
        original=[{'M':M,'a':a,'b':b} for M,p,q,a,b in profiles]
        self.assertEqual(checker.bind_domain(29,{'rows':original},{'rows':profiles}),profiles)
        for mutation in ('omission','duplicate','reorder','float'):
            rows=copy.deepcopy(original)
            if mutation=='omission':rows.pop()
            elif mutation=='duplicate':rows[-1]=copy.deepcopy(rows[0])
            elif mutation=='reorder':rows.reverse()
            else:rows[0]['a']=float(rows[0]['a'])
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                checker.bind_domain(29,{'rows':rows},{'rows':profiles})


if __name__=='__main__':
    unittest.main()
