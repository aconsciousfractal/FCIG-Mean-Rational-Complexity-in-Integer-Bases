"""Verify retained early closure inputs, joins, rational caps and finite domains.

This does not re-execute the geometric witnesses. Full producer/consumer routes
are provided in the eight case drivers; their recorded reconstructions are
separately described in RECONSTRUCTION_HISTORY.json.
"""
from pathlib import Path, PurePosixPath
from itertools import combinations
from fractions import Fraction
from math import isqrt, gcd, prod
import argparse, csv, hashlib, io, json, re

ROOT=Path(__file__).resolve().parent
CASES=('thirty_p_large','thirty_p_small','radical210','three_prime_28',
       'three_prime_family','three_prime_large','three_prime_353','three_prime_89')

def need(x,message):
    if not x: raise ValueError(message)

def pairs(items):
    d={}
    for k,v in items:
        need(k not in d,'duplicate JSON key: '+k);d[k]=v
    return d

def parse(data):
    return json.loads(data,object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON: '+x)))

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while b:=f.read(1<<20):h.update(b)
    return h.hexdigest()

def member(home,name):
    need(type(name) is str and '\\' not in name,'POSIX member path')
    pp=PurePosixPath(name)
    need(not pp.is_absolute() and all(p not in ('','.','..') for p in pp.parts)
         and pp.as_posix()==name and ':' not in name,'safe canonical member')
    path=home.joinpath(*pp.parts)
    need(path.resolve().is_relative_to(home.resolve()) and path.is_file(),'member containment/existence')
    need(not any(p.is_symlink() for p in [path,*path.parents] if p.is_relative_to(home)),'no symbolic inputs')
    return path

def inventory(home):
    j=home/'MANIFEST_SHA256.json';t=home/'SHA256SUMS'
    need(j.exists()!=t.exists(),'one case manifest')
    if j.exists():expected=parse(j.read_bytes());manifest=j.name
    else:
        expected={}
        for line in t.read_text(encoding='utf8').splitlines():
            h,name=line.split('  ',1);need(name not in expected,'duplicate manifest member');expected[name]=h
        manifest=t.name
    need(type(expected) is dict and expected,'nonempty manifest')
    for name,h in expected.items():
        need(type(h) is str and re.fullmatch('[0-9a-f]{64}',h),'canonical sha256')
        need(sha(member(home,name))==h,'case input digest: '+name)
    actual={p.relative_to(home).as_posix() for p in home.rglob('*') if p.is_file()
            and p.relative_to(home).as_posix()!=manifest and '__pycache__' not in p.parts
            and p.suffix not in ('.pyc','.pyo')}
    need(actual==set(expected),'complete exact case inventory')
    return len(expected)

def source(root,ref):
    need(type(ref.get('case')) is str and ref['case'] in CASES,'declared case')
    home=root/ref['case']
    need(home.resolve().is_relative_to(root.resolve()),'case containment')
    path=member(home,ref['path'])
    need(sha(path)==ref['sha256'],'contract source digest')
    return path.read_bytes()

def pointer(d,p):
    need(p=='' or p.startswith('/'),'JSON pointer')
    for key in p.split('/')[1:]:
        key=key.replace('~1','/').replace('~0','~')
        d=d[int(key)] if isinstance(d,list) else d[key]
    return d

def selected(root,ref):
    data=source(root,ref);sel=ref['selector']
    need(sel['format'] in ('whitespace_integer_lines','csv'),'selector format')
    need(sel['transform'] in ('integer_set','ordered_prime_triples_to_set'),'selector transform')
    need(sel['format']!='whitespace_integer_lines' or sel['transform']=='integer_set','text selector transform')
    if sel['format']=='whitespace_integer_lines':values=[int(t) for t in data.split()]
    else:
        rows=list(csv.DictReader(io.StringIO(data.decode())));cols=sel['columns']
        values=[int(t[cols[0]]) if sel['transform']=='integer_set'
                else tuple(int(t[k]) for k in cols) for t in rows]
        if sel['transform']=='ordered_prime_triples_to_set':
            for v,t in zip(values,rows):
                need(v[0]<v[1]<v[2],'ordered prime triple')
                need(all(n>=2 and all(n%d for d in range(2,isqrt(n)+1)) for n in v),'prime triple')
                if 'radical' in t:need(int(t['radical'])==v[0]*v[1]*v[2],'support product')
    need(len(values)==len(set(values)),'duplicate support selector')
    return set(values)

def collect_M(o):
    if isinstance(o,dict):
        if 'M' in o and type(o['M']) is int:yield o['M']
        else:
            for v in o.values():yield from collect_M(v)
    elif isinstance(o,list):
        for v in o:yield from collect_M(v)

def verify(root=ROOT):
    root=Path(root).resolve()
    need(tuple(parse((root/'CASES.json').read_bytes()))==CASES,'eight-case catalog order/domain')
    inputs=sum(inventory(root/c) for c in CASES)
    spec=parse((root/'CONTRACT.json').read_bytes());checks=0
    need(spec['schema']=='early-closure-contract-1','contract schema')
    for c in spec['checks']:
        b=source(root,c['source']);op=c['operation']
        if 'json_pointer' in c:
            need(op in ('equals','array_length_equals'),'JSON check operation')
            x=pointer(parse(b),c['json_pointer'])
            if op=='array_length_equals':x=len(x)
        else:
            need(op in ('row_count_equals','integer_sum_equals','sum_floor_half_plus_one_equals'),'CSV check operation')
            rows=list(csv.DictReader(io.StringIO(b.decode())));nums=[int(t[c['column']]) for t in rows]
            x={'row_count_equals':lambda:len(rows),'integer_sum_equals':lambda:sum(nums),
               'sum_floor_half_plus_one_equals':lambda:sum(t//2+1 for t in nums)}[op]()
        need(type(x) is type(c['expected']) and x==c['expected'],'retained result contract');checks+=1
    for e in spec['imports']:
        need(e['property'] in ('bytes_equal','sets_equal'),'import property')
        if e['property']=='bytes_equal':need(source(root,e['child'])==source(root,e['parent']),'numeric import byte join')
        else:
            got=selected(root,e['child']);need(got==selected(root,e['parent']) and len(got)==e['count'],'numeric import set join')
    part=spec['partition'];universe=selected(root,part['universe']);parts=[selected(root,r) for r in part['parts']]
    primes=[n for n in range(2,48) if all(n%d for d in range(2,isqrt(n)+1))]
    defined={t for t in combinations(primes,3) if t[:2]!=(2,3)}
    need(universe==defined==parts[0]|parts[1] and not(parts[0]&parts[1])
         and [len(v) for v in parts]==[353,89],'complete three-prime partition')
    census=parse((root/'three_prime_28/results/census_certificate.json').read_bytes())
    mods=[int(t['M']) for t in census['moduli']]
    for c in ('three_prime_family','three_prime_large','three_prime_353'):
        mods += [int(t['M']) for t in csv.DictReader((root/c/'results/finite_domain.csv').open())]
    last=sorted(set(collect_M(parse((root/'three_prime_89/results/finite.json').read_bytes()))))
    need(len(last)==35,'last finite domain');mods+=last
    need(len(mods)==len(set(mods))==11911 and sum(mods)==707038477
         and sum(m//2+1 for m in mods)==353528303,'exact three-prime finite-domain aggregate')
    # Numeric parent import is joined independently of an imported-replay flag.
    small=parse((root/'thirty_p_small/results/closure.json').read_bytes())
    large=parse((root/'thirty_p_large/results/closure.json').read_bytes())
    need(parse((root/'thirty_p_large/imports/S7_SOURCE_CLOSURE.json').read_bytes())==small,'seven-prime numeric parent join')
    need(small['closed_primes']==[11,13,17,19,23,29,31],'seven small primes')
    low=[p for p in range(11,211) if all(p%d for d in range(2,isqrt(p)+1))]
    need(small['closed_primes']+large['new_full_mark_primes']==low,'thirty-p finite-prime split')
    need(Fraction(small['all_new_regions_common_upper'])<Fraction(233,10000)<Fraction(2,85),'seven-prime cap')
    need(Fraction(large['common_upper'])<Fraction(2,85),'remaining thirty-p cap')
    rows=list(csv.DictReader((root/'thirty_p_large/results/finite_inventory.csv').open()))
    finite=list(csv.DictReader((root/'thirty_p_large/results/finite_bounds.csv').open()))
    expected=[t for t in rows if t['route']=='witness']
    need([int(t['M']) for t in finite]==[int(t['M']) for t in expected],'thirty-p exact witness residual')
    for t,e in zip(finite,expected):
        n={k:int(v) for k,v in t.items()};M=n['M']
        need(n['p']==int(e['p']) and n['states']==M and n['representatives']==M//2+1,'thirty-p state weights')
        need(n['mean_numerator']==2*n['sum_log_upper_ticks']-M*n['log_M_lower_ticks']<0,'thirty-p integer log assembly')
    # Candidate is a positive lower bound only. Its sector result is not
    # promoted to a claim about the global optimizer.
    candidate=parse((root/'radical210/results/finite_certificate.json').read_bytes())
    conclusion=parse((root/'radical210/results/conclusion.json').read_bytes())
    need(conclusion['unique_maximizing_pair']==[210,1] and not conclusion['global_Cstar_value_determined'],'candidate sector scope')
    # Direct definition oracle at this one small modulus, not a re-use of
    # the retained witness upper heights in the other 120 moduli.
    heights=[]
    for s in range(210):
        best=max(1,min(s,210-s))
        for den in range(1,best+1):
            if gcd(den,210)==1:
                residue=s*den%210
                best=min(best,max(den,min(residue,210-residue)))
        heights.append(best)
    product=prod(heights)
    need(product==int(candidate['candidate_P']),'exact candidate product')
    need(product==prod(int(p)**int(e) for p,e in candidate['candidate_factorization'].items()),'candidate factorization')
    need(product**170>210**18690,'candidate exceeds 2/85 by integer comparison')
    for case,key in [('three_prime_large','common_delta_upper'),('three_prime_353','common_delta_upper'),('three_prime_89','common_bound')]:
        name='family_closure.json' if case=='three_prime_large' else 'closure.json'
        need(Fraction(parse((root/case/'results'/name).read_bytes())[key])<Fraction(2,85),'three-prime rational cap')
    return dict(status='RETAINED_EARLY_CLOSURE_INPUTS_PASS',cases=8,case_input_files=inputs,
                result_checks=checks,numeric_import_edges=len(spec['imports']),
                three_prime_partition=[353,89],three_prime_finite_moduli=11911,
                full_geometric_reconstruction_in_this_run=False,
                previous_full_reconstruction_date='2026-10-08',
                current_scope='retained numerical identities, domain joins, rational caps and finite ledger; full routes supplied separately')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    a=p.parse_args()
    try:print(json.dumps(verify(a.root),sort_keys=True))
    except (ValueError,KeyError,TypeError,OSError) as e:p.exit(1,'FAIL '+str(e)+'\n')

