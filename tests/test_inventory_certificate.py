"""Small metadata fixtures for complete two-method inventory derivation."""
import copy
import hashlib
import json
from math import gcd
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'companion/programs'))
import derive_inventory_certificate as inventory


def fixture():
    rows=[]
    for p,q,a,b in ((2,3,1,1),(2,11,2,1),(3,7,2,1),(5,7,2,1),
                    (2,8191,2,1),(2,8191,3,1),(3,8191,2,1)):
        M=p**a*q**b
        branch,cutoff=('p2',147456) if p==2 else ('p3',262144) if p==3 else ('p_ge5',360000)
        reasons=[]
        if p*q in (6,10,14,15):reasons.append('closed_radical')
        if gcd(a,b)>=2:reasons.append('common_valuation_factor')
        if a==b==1:reasons.append('squarefree')
        if M>=cutoff:reasons.append('branch_cutoff')
        rows.append(dict(M=M,p=p,q=q,a=a,b=b,branch=branch,cutoff=cutoff,reasons=reasons))
    retained=[r for r in rows if not r['reasons']]
    optional=[r for r in retained if inventory.ratio_qualifies(r,optional=True)]
    supplemental=[r for r in retained if inventory.ratio_qualifies(r)]
    active=[r for r in retained if not inventory.ratio_qualifies(r)]
    arrays=dict(zip(inventory.FIELDS,(rows,retained,optional,supplemental,active)))
    generator=dict(schema=inventory.SCHEMA,status='SOURCE_ONLY',algorithm='ordered_primes_and_powers',
                   reason_order=list(inventory.REASONS),
                   domain={'M_lower_inclusive':2,'M_upper_exclusive':360000,'omega':2,
                           'prime_order':'p<q','valuations':'a,b>=1'},**copy.deepcopy(arrays))
    factor=dict(schema=inventory.SCHEMA,status='SOURCE_ONLY',algorithm='integer-traversal-smallest-prime-factor-sieve',
                reason_order=list(inventory.REASONS),domain={'M_min_inclusive':2,'M_max_exclusive':360000},
                **copy.deepcopy(arrays))
    template={'status':inventory.STATUS,'normalized_full_list_fields':list(inventory.FIELDS),
              'supplemental_import':'fixed attributed mathematical premise',
              'summaries':{'deliberately':'false, never copied'},'source_freeze_sha256':'omitted'}
    return generator,factor,template,{field:len(rows) for field,rows in arrays.items()}


class InventoryCertificateTests(unittest.TestCase):
    def test_summaries_are_derived_and_control_metadata_is_not_copied(self):
        generator,factor,template,counts=fixture()
        with patch.object(inventory,'COUNTS',counts):
            result=inventory.derive(generator,factor,'g','f',template)
        self.assertEqual(result['summaries']['all_rows'],dict(count=7,mass=172299,min_M=6,max_M=73719,radical_count=6))
        self.assertEqual(result['summaries']['active_rows']['mass'],74001)
        self.assertEqual([(r['branch'],r['count'],r['mass']) for r in result['active_by_branch']],
                         [('p2',1,44),('p3',2,73782),('p_ge5',1,175)])
        self.assertEqual(result['supplemental_import'],template['supplemental_import'])
        self.assertNotIn('source_freeze_sha256',result)
        self.assertFalse(result['general_N1_closed'])

    def test_each_full_field_is_compared_not_only_summaries(self):
        generator,factor,template,counts=fixture()
        for field in inventory.FIELDS:
            changed=copy.deepcopy(factor)
            changed[field][0]['M']+=1
            with self.subTest(field=field),patch.object(inventory,'COUNTS',counts),self.assertRaises(ValueError):
                inventory.derive(generator,changed,'g','f',template)

    def test_jointly_truncated_inputs_cannot_claim_complete(self):
        generator,factor,template,counts=fixture()
        for document in (generator,factor):document['all_rows'].pop()
        with patch.object(inventory,'COUNTS',counts),self.assertRaisesRegex(ValueError,'complete inventory list count'):
            inventory.derive(generator,factor,'g','f',template)

    def test_status_schema_and_algorithm_are_required(self):
        generator,factor,template,counts=fixture()
        for key in ('status','schema','algorithm'):
            for which in (0,1):
                documents=[copy.deepcopy(generator),copy.deepcopy(factor)]
                documents[which][key]='incomplete'
                with self.subTest(key=key,which=which),patch.object(inventory,'COUNTS',counts),self.assertRaises(ValueError):
                    inventory.derive(*documents,'g','f',template)

    def test_missing_duplicate_reordered_or_inexact_profiles_rejected(self):
        generator,factor,template,counts=fixture()
        for mutation in ('missing','duplicate','reorder','float','bool','composite','reason','partition'):
            changed=copy.deepcopy(generator)
            rows=changed['all_rows']
            if mutation=='missing':del changed['active_rows']
            elif mutation=='duplicate':rows[-1]=copy.deepcopy(rows[-2])
            elif mutation=='reorder':rows.reverse()
            elif mutation=='float':rows[0]['M']=6.0
            elif mutation=='bool':rows[0]['p']=True
            elif mutation=='composite':rows[1]['q']=9
            elif mutation=='reason':rows[0]['reasons']=[]
            else:changed['active_rows'][0]=copy.deepcopy(changed['retained_rows'][-1])
            with self.subTest(mutation=mutation),patch.object(inventory,'COUNTS',counts),self.assertRaises(ValueError):
                inventory.derive(changed,factor,'g','f',template)

    def test_ambiguous_JSON_and_wrong_artifact_pin_rejected(self):
        for raw in (b'{"x":1,"x":2}',b'{"x":{"a":1,"a":2}}',b'{"x":NaN}'):
            with self.assertRaises(ValueError):inventory.strict_load(raw)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.json';path.write_bytes(b'{"x":1}')
            self.assertEqual(inventory.pinned(path,hashlib.sha256(path.read_bytes()).hexdigest()),{'x':1})
            with self.assertRaisesRegex(ValueError,'byte pin'):inventory.pinned(path,'0'*64)

    def test_output_is_fresh_and_summary_binds_actual_bytes(self):
        generator,factor,template,counts=fixture()
        with patch.object(inventory,'COUNTS',counts):result=inventory.derive(generator,factor,'g','f',template)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'certificate.json'
            summary=inventory.fresh_output(path,result)
            raw=path.read_bytes()
            self.assertEqual(json.loads(raw),result)
            self.assertEqual(summary['certificate_sha256'],hashlib.sha256(raw).hexdigest())
            self.assertEqual(summary['certificate_bytes'],len(raw))
            with self.assertRaises(ValueError):inventory.fresh_output(path,result)
            self.assertEqual(path.read_bytes(),raw)


if __name__=='__main__':unittest.main()
