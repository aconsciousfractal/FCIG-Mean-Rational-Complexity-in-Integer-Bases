"""Negative controls for complete observed-domain bindings."""
import copy
import importlib.util
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('complete_evidence_control',ROOT/'scripts/check_complete_evidence.py')
check=importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


class EvidenceControls(unittest.TestCase):
    def fixture(self):
        jobs=[{'id':name,'output':{'retained_sha256':'a'*64}}
              for name in sorted(check.REQUIRED_JOB_IDS)]
        rows=[{'id':j['id'],'status':'PASS','exit_code':0,'retained_math_equal':True,
               'raw_retained_equal':True,'consumed_output_sha256':'a'*64} for j in jobs]
        index={'normal_optimized_scientific_bytes_identical':True,
               'observed_runs':{'normal':copy.deepcopy(rows),'optimized':copy.deepcopy(rows)},
               'records':{'chain/'+j['id']:{'sha256':'a'*64} for j in jobs}}
        return index,{'jobs':jobs}

    def test_complete_domain(self):
        i,p=self.fixture()
        check.verify_observations(i,p)

    def test_missing_duplicate_reordered_jobs(self):
        for mutation in ('missing','duplicate','reorder'):
            i,p=self.fixture()
            rows=i['observed_runs']['normal']
            if mutation=='missing':rows.pop()
            elif mutation=='duplicate':rows[1]=rows[0]
            else:rows.reverse()
            with self.assertRaises(ValueError):check.verify_observations(i,p)

    def test_wrong_consumed_pin(self):
        i,p=self.fixture()
        i['observed_runs']['optimized'][0]['consumed_output_sha256']='b'*64
        with self.assertRaisesRegex(ValueError,'byte identity'):check.verify_observations(i,p)

    def test_failed_job_or_bool_exitcode(self):
        for key,value in (('status','FAILED'),('exit_code',False),('retained_math_equal',False)):
            i,p=self.fixture()
            i['observed_runs']['normal'][0][key]=value
            with self.assertRaisesRegex(ValueError,'success'):check.verify_observations(i,p)

    def test_missing_mode_and_extra_plan_job(self):
        i,p=self.fixture()
        del i['observed_runs']['optimized']
        with self.assertRaises(ValueError):check.verify_observations(i,p)
        i,p=self.fixture()
        p['jobs'].append(copy.deepcopy(p['jobs'][0]))
        with self.assertRaisesRegex(ValueError,'37-job'):check.verify_observations(i,p)


if __name__=='__main__':unittest.main()
