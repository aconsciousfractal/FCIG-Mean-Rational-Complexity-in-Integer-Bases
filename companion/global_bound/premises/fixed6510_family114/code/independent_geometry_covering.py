#!/usr/bin/env python3
"""Portable adapter for the unchanged independent geometry/covering review code."""
from pathlib import Path
import argparse
import shutil
import tempfile
from common import need, read, sha, canonical, module

ROOT=Path(__file__).resolve().parents[1]


def certificate(root, work):
    root=Path(root).resolve();work=Path(work).resolve()
    need(not work.is_relative_to(root), 'independent review outputs outside package')
    work.mkdir(parents=True,exist_ok=True)
    source=root/'reviews/geometry_covering_check.py'
    review=read(root/'reviews/geometry_covering_theory.json')
    need(sha(source)==review['independent_code_sha256'], 'unaltered independent review program')
    before=sha(source)
    imported=module(source,'p72_independent_geometry_covering_portable')
    # The original review was executed in the research workspace. Its byte-exact
    # snapshot stays unchanged; these two runtime bindings make the replay portable.
    imported.ROOT=root
    imported.W=work
    imported.run()
    current=read(work/'geometry_covering_independent.json')
    expected=read(root/'reviews/geometry_covering_independent.json')
    ignored={'runtime_seconds','independent_code_path'}
    normalized=lambda z:{k:v for k,v in z.items() if k not in ignored}
    need(canonical(normalized(current))==canonical(normalized(expected)),
         'complete independent scientific certificate identity')
    need(sha(source)==before and current['all_source_bytes_unchanged']
         and not current['producer_modules_imported'] and not current['new_geometry_raw_replayed'],
         'actual independent arithmetic replay and unchanged source')
    return dict(status='PASS', full_scientific_certificate_identity=True,
                ignored_nonscientific_fields=sorted(ignored),
                current_independent_code_sha256=before,
                wrapper_sha256=sha(Path(__file__)),
                geometry_rows=current['geometry_rows_checked'],
                integer_gap_pairs=current['integer_gap_pairs'],
                integer_endpoint_checks=current['integer_endpoint_checks'],
                covering_rows=current['covering_rows_reconstructed'],
                all_interval_starts=current['covering_interval_start_counts_reconstructed'],
                new_covering_rows=current['new_covering_rows_reconstructed'],
                producer_modules_imported=False, current_geometric_raw_replayed=False,
                all_reviewed_source_bytes_unchanged=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--work',type=Path)
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();report=a.report.resolve()
    need(not report.is_relative_to(root), 'external report')
    own_work=a.work is None
    work=a.work or Path(tempfile.mkdtemp(prefix='p72-independent-geometry-'))
    success=False
    try:
        result=certificate(root,work);report.parent.mkdir(parents=True,exist_ok=True)
        report.write_bytes(canonical(result));print(canonical(result).decode(),end='');success=True
    finally:
        if own_work and success:shutil.rmtree(work)


if __name__=='__main__':main()
