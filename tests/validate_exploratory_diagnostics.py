"""Validate new exploratory outputs against distinct algorithms and identities."""
from pathlib import Path
import csv
import hashlib
import json
import math
import statistics
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from exact_stats import exact_conditional, conditional_distribution
from exploratory_diagnostics import ConditionalModel, read_inputs, continuous_diagnostics

strata, genes = read_inputs()
info = json.loads((ROOT / 'results/conditional_information.json').read_text())
continuous = json.loads((ROOT / 'results/continuous_summary.json').read_text())
with (ROOT / 'results/stratum_influence.csv').open(newline='') as stream:
    influence = list(csv.DictReader(stream))
comparisons = []
for omitted in [None] + [r['omitted_stratum'] for r in influence]:
    selected = [s for s in strata if s['stratum'] != omitted]
    expected = exact_conditional(selected)
    observed = info['primary_independent_inference'] if omitted is None else next(r for r in influence if r['omitted_stratum'] == omitted)
    diffs = {key: abs(float(observed[key]) - expected[key]) for key in ('conditional_mle', 'ci_lower', 'ci_upper', 'p_enrichment')}
    assert max(diffs.values()) < 1e-8, (omitted, diffs)
    comparisons.append({'omitted': omitted, 'max_abs_difference': max(diffs.values())})

distribution_checks = []
for theta in (1., 1.5, 2.):
    x, expected = conditional_distribution(strata, odds_ratio=theta)
    observed = ConditionalModel(strata).pmf(theta)
    assert np.array_equal(x, np.arange(sum(s['lo'] for s in strata), sum(s['hi'] for s in strata)+1))
    difference = float(np.max(np.abs(expected-observed)))
    assert difference < 1e-9
    distribution_checks.append({'theta': theta, 'max_abs_pmf_difference': difference})

threshold_checks = []
threshold = info['minimum_overlap_for_criterion']
for t in (threshold - 1, threshold):
    # Any feasible allocation with this total has the same conditional likelihood.
    remaining = t - sum(s['lo'] for s in strata)
    synthetic = []
    for s in strata:
        increment = min(remaining, s['hi'] - s['lo'])
        synthetic.append({**s, 'a': s['lo'] + increment})
        remaining -= increment
    assert remaining == 0
    result = exact_conditional(synthetic)
    assert (result['ci_lower'] > 1) == (t == threshold)
    assert (result['p_enrichment'] < .025) == (t == threshold)
    threshold_checks.append({'hypothetical_total': t, 'lower_95': result['ci_lower'], 'null_upper_tail': result['p_enrichment']})

continuous_checks = []
for j, row in enumerate(continuous['summary']):
    # A1 weighted mean must equal its unweighted gene mean under exposed-target weights.
    a1 = statistics.mean(g['values'][j] for g in genes if g['A'])
    a0_terms = []
    variance_terms = []
    for s in strata:
        if not s['m']:
            continue
        controls = [g['values'][j] for g in genes if g['stratum'] == s['stratum'] and not g['A']]
        w = s['m'] / sum(v['m'] for v in strata)
        # Direct inverse-control-cell-size weighting, not a group mean operation.
        a0_terms.extend(w * y / len(controls) for y in controls)
        variance_terms.append(w * statistics.variance(controls))
    a0 = math.fsum(a0_terms)
    delta = a1 - a0
    scaled = delta / math.sqrt(math.fsum(variance_terms))
    assert abs(a1 - row['weighted_mean_lfc_A1']) < 1e-12
    assert abs(delta - row['exposure_standardized_mean_difference']) < 1e-12
    assert abs(scaled - row['sd_scaled_mean_difference']) < 1e-12
    continuous_checks.append({'source_label': row['replicate_label'], 'independent_delta': delta, 'independent_scaled_delta': scaled})

rerun, _, _ = continuous_diagnostics(strata, genes, continuous['bootstrap']['draws'], continuous['bootstrap']['seed'])
assert rerun['summary'] == continuous['summary'], 'Seeded bootstrap not exactly reproducible'
singleton = next(g for g in genes if g['stratum'] == '1|3-5|4+' and g['A'])
singleton_rows = [r for r in continuous['stratum_details'] if r['stratum'] == singleton['stratum']]
assert len(singleton_rows) == 2 and all(r['n_A1'] == 1 and r['singleton_exposure_no_bootstrap_variation'] for r in singleton_rows)

output = {
    'status': 'PASS',
    'scope': 'Computational checks only; no label semantic or biological validation.',
    'estimator_crosschecks': comparisons,
    'distribution_crosschecks': distribution_checks,
    'criterion_boundary_checks': threshold_checks,
    'continuous_algebra_checks': continuous_checks,
    'seeded_bootstrap_exact_reproduction': True,
    'singleton_exposure_flag_verified': True,
    'code_sha256': hashlib.sha256((ROOT / 'code/exploratory_diagnostics.py').read_bytes()).hexdigest(),
    'validation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(ROOT / 'evidence/exploratory_diagnostics_validation.json').write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8')
print(json.dumps(output, indent=2))
