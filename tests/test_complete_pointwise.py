"""Small metadata/negative guard fixtures, not large mathematical replay.

The synthetic histograms and mocked log/CDF values test consumer plumbing.
They are not claimed admissible height envelopes or numerical certificates.
"""
import copy
from functools import lru_cache
import hashlib
import io
import json
from math import isqrt
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

PROGRAMS = Path(__file__).resolve().parents[1] / 'companion/programs'
sys.path.insert(0,str(PROGRAMS))
import verify_full_semiprime as full
import verify_pointwise_central as pointwise


def raw_digest(value):
    raw=json.dumps(value,sort_keys=True,separators=(',',':')).encode()
    return raw,hashlib.sha256(raw).hexdigest()


def log_fixture():
    @lru_cache(maxsize=None,typed=True)
    def ticks(n):
        return (0,0) if n==1 else (1,1)
    return ticks


def semiprime_fixture():
    domain=[{'M':6,'p':2,'r':3},{'M':10,'p':2,'r':5},{'M':15,'p':3,'r':5}]
    rows=[]
    for profile in domain:
        M=profile['M'];hist=[[1,M]]
        rows.append({**profile,'upper_histogram':hist,
                     'upper_histogram_sha256':hashlib.sha256(json.dumps(hist,separators=(',',':')).encode()).hexdigest(),
                     'log_tick_ledger':[[1,0,0],[M,1,1]],
                     'upper_centered_log_numerator':-M,
                     'upper_centered_log_denominator':2*M*full.S})
    historic={'domain_sha256':'domain-pin','rows':copy.deepcopy(rows)}
    return {'scope':'full','status':'PORTABLE_E032_RECONSTRUCTION_PASS','rows':rows},historic,domain


def point_fixture():
    rows=[]
    for n in range(1,9):
        M=6**n;k=isqrt(M)
        hist=[[1,1],[3,1],[k+1,M-2]] if k>=3 else [[1,1],[k+1,M-1]]
        rows.append({'M':M,'N':n,'histogram':hist})
    return {'status':'E029_PASS','rows':rows}


def synthetic_central(M):
    k=isqrt(M)
    return {'k':k,'cdf':[0]+[1 if T<3 else 2 for T in range(1,k+1)],
            'central_vector_sha256':'a'*64,'counters':{}}


class CompletePointwiseTests(unittest.TestCase):
    def check_full(self,payload,historic,domain,pin_override=None):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'fixture.json';raw,pin=raw_digest(payload);path.write_bytes(raw)
            def supplied(batch,name):
                return historic if name=='SEMIPRIME_BRIDGE.json' else {'rows':domain}
            with patch.object(full,'COUNT',3),patch.object(full,'data',side_effect=supplied),\
                 patch.object(full,'entry',return_value={'sha256':'domain-pin'}),\
                 patch.object(full,'fine_log',log_fixture()):
                return full.verify(path,pin if pin_override is None else pin_override)

    def test_full_fixture_complete_consumed_ledger(self):
        payload,historic,domain=semiprime_fixture()
        result=self.check_full(payload,historic,domain)
        self.assertEqual(result['histograms_checked'],3)
        self.assertEqual(result['consumed_tick_entries'],6)
        self.assertTrue(result['every_outward_sum_negative'])
        self.assertTrue(result['pointwise_dominance_producer_and_proof_required'])

    def test_full_domain_binding_types_order_and_scope(self):
        for mutation in ('float_profile','duplicate','reorder','missing','scope','domain_pin','historic_order'):
            payload,historic,domain=semiprime_fixture()
            if mutation=='float_profile':payload['rows'][0]['p']=2.0
            elif mutation=='duplicate':payload['rows'][-1]=copy.deepcopy(payload['rows'][0])
            elif mutation=='reorder':payload['rows'].reverse()
            elif mutation=='missing':payload['rows'].pop()
            elif mutation=='scope':payload['scope']='sample'
            elif mutation=='domain_pin':historic['domain_sha256']='wrong'
            else:historic['rows'].reverse()
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                self.check_full(payload,historic,domain)

    def test_full_histogram_tick_sign_and_normalization_guards(self):
        for mutation in ('histogram','ledger_missing','ledger_extra','ledger_duplicate','ledger_float',
                         'tick_direction','positive_Q','wrong_Q','float_denominator','historic_float_denominator',
                         'historic_float_Q'):
            payload,historic,domain=semiprime_fixture();r=payload['rows'][0]
            if mutation=='histogram':r['upper_histogram'][0][1]-=1
            elif mutation=='ledger_missing':r['log_tick_ledger'].pop()
            elif mutation=='ledger_extra':r['log_tick_ledger'].insert(1,[2,1,1])
            elif mutation=='ledger_duplicate':r['log_tick_ledger'].append([6,1,1])
            elif mutation=='ledger_float':r['log_tick_ledger'][0][1]=0.0
            elif mutation=='tick_direction':r['log_tick_ledger'][-1][1:]=[2,2]
            elif mutation=='positive_Q':r['upper_centered_log_numerator']=1
            elif mutation=='wrong_Q':r['upper_centered_log_numerator']-=1
            elif mutation=='float_denominator':r['upper_centered_log_denominator']=float(r['upper_centered_log_denominator'])
            elif mutation=='historic_float_denominator':historic['rows'][0]['upper_centered_log_denominator']=float(historic['rows'][0]['upper_centered_log_denominator'])
            else:historic['rows'][0]['upper_centered_log_numerator']=-6.0
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                self.check_full(payload,historic,domain)

    def test_full_input_byte_pin(self):
        with self.assertRaises(ValueError):
            self.check_full(*semiprime_fixture(),pin_override='0'*64)

    def test_full_histogram_identity_rejected_with_correct_mass(self):
        payload,historic,domain=semiprime_fixture()
        payload['rows'][0]['upper_histogram']=[[1,5],[2,1]]
        with self.assertRaisesRegex(ValueError,'histogram identity'):
            self.check_full(payload,historic,domain)

    def run_pointwise(self,payload,start=0,stop=None,existing=False,pin_override=None):
        with tempfile.TemporaryDirectory() as directory:
            input_path=Path(directory)/'input.json';output=Path(directory)/'output.json'
            raw,pin=raw_digest(payload);input_path.write_bytes(raw)
            if existing:output.write_bytes(b'sentinel')
            argv=['verify_pointwise_central.py','--kind','base6','--input',str(input_path),
                  '--input-sha256',pin if pin_override is None else pin_override,
                  '--output',str(output),'--start',str(start)]
            if stop is not None:argv+=['--stop',str(stop)]
            with patch.object(sys,'argv',argv),patch.object(pointwise,'enumerate_central',side_effect=synthetic_central) as calls,\
                 patch('sys.stdout',new_callable=io.StringIO),patch('sys.stderr',new_callable=io.StringIO):
                try:pointwise.main()
                except Exception:
                    if existing:
                        self.assertEqual(output.read_bytes(),b'sentinel')
                        calls.assert_not_called()
                    raise
            return json.loads(output.read_text()),calls.call_count

    def test_pointwise_complete_domain_and_chunk_scope(self):
        report,calls=self.run_pointwise(point_fixture())
        self.assertTrue(report['complete_domain_verified']);self.assertEqual(calls,8)
        self.assertEqual([r['domain_index'] for r in report['rows']],list(range(8)))
        report,calls=self.run_pointwise(point_fixture(),1,2)
        self.assertFalse(report['complete_domain_verified']);self.assertEqual(calls,1)
        self.assertEqual(report['rows'][0]['domain_index'],1)

    def test_pointwise_interior_prefix_mutation_preserves_endpoint_but_fails(self):
        payload=point_fixture()
        # At M=36 change height3 to height2: total mass and CDF(k) unchanged.
        payload['rows'][1]['histogram'][1][0]=2
        with self.assertRaisesRegex(ValueError,'T=2'):
            self.run_pointwise(payload,1,2)

    def test_pointwise_histogram_domain_pin_and_fresh_output(self):
        for mutation in ('duplicate','reorder','hist_float','hist_duplicate','hist_mass','status','pin','existing'):
            payload=point_fixture()
            if mutation=='duplicate':payload['rows'][-1]=copy.deepcopy(payload['rows'][0])
            elif mutation=='reorder':payload['rows'].reverse()
            elif mutation=='hist_float':payload['rows'][0]['histogram'][0][1]=1.0
            elif mutation=='hist_duplicate':payload['rows'][0]['histogram'].append([1,1])
            elif mutation=='hist_mass':payload['rows'][0]['histogram'][0][1]+=1
            elif mutation=='status':payload['status']='FAILED_PARTIAL'
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                self.run_pointwise(payload,existing=mutation=='existing',
                                   pin_override='0'*64 if mutation=='pin' else None)

    def test_pointwise_unselected_nonexact_domain_row_rejected(self):
        payload=point_fixture()
        payload['rows'][-1]['M']=float(payload['rows'][-1]['M'])
        with self.assertRaises(ValueError):
            self.run_pointwise(payload,0,1)


if __name__=='__main__':unittest.main()
