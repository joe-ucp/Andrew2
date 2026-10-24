"""One-command data verification, analysis, estimator tests and independent rebuild."""
import argparse, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--offline',action='store_true',help='Verify supplied raw files and the external author-code cache without downloading')
    ap.add_argument('--source-cache',help='External author-code cache directory; defaults to BETA_CELL_SOURCE_CACHE or ~/.cache/beta-cell-screen-comparison/sources')
    ap.add_argument('--build-manuscript',action='store_true',help='Also regenerate matching Word and TeX sources and compile PDF with pdflatex')
    args=ap.parse_args()
    cache_args=['--source-cache',args.source_cache] if args.source_cache else []
    source_args=cache_args+(['--offline'] if args.offline else [])
    steps=[['code/fetch_inputs.py']+source_args,
           ['-m','unittest','discover','-s','tests','-p','test_source_cache.py'],
           ['tests/run_estimator_validation.py'],['code/analyze.py'],['code/independent_verify.py','--joint'],['tests/run_statistical_crosscheck.py'],
           ['code/audit_source_semantics.py','--defer-revised-checks']+source_args,['code/revise_analysis.py'],['code/audit_source_semantics.py']+source_args,
           ['tests/validate_source_semantics.py']+cache_args,['tests/validate_revision_catalog.py'],
           ['code/audit_identifier_mapping.py'],['tests/validate_identifier_sensitivity.py'],
           ['code/exploratory_diagnostics.py'],['tests/validate_exploratory_diagnostics.py']]
    if args.build_manuscript: steps.append(['code/build_manuscript.py','--compile'])
    for step in steps:
        print('RUN '+' '.join(step),flush=True)
        subprocess.run([sys.executable,*step],cwd=ROOT,check=True)
    print('REPRODUCTION PASSED: revised published-catalog computations verified; S8 significance interpretation rejected, selected S7 support bounded; biological validation not performed.',flush=True)
