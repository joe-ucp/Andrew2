"""Run synthetic exact-estimator validation and save auditable JSON evidence."""
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import platform
import sys
import unittest

import numpy
import scipy


def main():
    root = Path(__file__).resolve().parents[1]
    suite = unittest.defaultTestLoader.discover(str(root / 'tests'), pattern='test_exact_stats.py')
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    evidence = {
        'validated_at_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'passed' if result.wasSuccessful() and not result.skipped else 'failed',
        'synthetic_data_only': True,
        'test_groups_run': result.testsRun,
        'failure_count': len(result.failures),
        'error_count': len(result.errors),
        'skipped_count': len(result.skipped),
        'python_version': platform.python_version(),
        'numpy_version': numpy.__version__,
        'scipy_version': scipy.__version__,
        'independent_comparators': [
            'Full joint enumeration with integer math.comb on synthetic strata',
            'scipy.stats.contingency.odds_ratio(kind=conditional), single-stratum MLE and exact CI',
            'scipy.stats.fisher_exact(alternative=greater), single-stratum inclusive enrichment p',
        ],
        'validation_properties': [
            '20 enumerated stratum-set/odds-ratio distributions',
            'Multi-stratum mean equation and exact equal-tail inversion',
            'Zero and infinite estimates and intervals',
            'Empty and wholly noninformative strata explicitly non-estimable',
            'Outcome-complement reciprocal transformation and stratum permutation',
            'Invalid margins rejected; output serializes with allow_nan=False',
            'Five independent SciPy single-stratum estimates, intervals and p-values',
        ],
        'file_sha256': {
            str(path.relative_to(root)).replace('\\', '/'): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [root/'code'/'exact_stats.py', root/'tests'/'test_exact_stats.py', Path(__file__)]
        },
        'command': 'python tests/run_estimator_validation.py',
        'test_output': stream.getvalue(),
    }
    destination = root / 'evidence' / 'estimator_validation.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(evidence, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(stream.getvalue())
    print(f'Validation evidence: {destination}')
    return 0 if evidence['status'] == 'passed' else 1


if __name__ == '__main__':
    sys.exit(main())
