"""Independent V2 catalog regression checks from the archived V1 source ledger.

Does not run or import analyze.py/revise_analysis.py/exact_stats.py. The ledger
was separately reconstructed from raw sources. This test rebuilds all seven
catalog selections and all gene/stratum labels without guide-FDR gating.
Conditional inference uses the independent SciPy-based V1 audit helper.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import bisect
import csv
import hashlib
import json
import math

from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from run_statistical_crosscheck import scipy_conditional

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
NAMES = ['primary', 'no_expression', 'expression_ge1', 'expression_ge100',
         'survival_fdr005', 'survival_harming', 'positive_vs_negative']


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    return float(row[key]) if row[key] else math.nan


def truth(row, key):
    assert row[key] in ('True', 'False'), (key, row[key])
    return row[key] == 'True'


def quantile(ordered, q):
    index = (len(ordered) - 1) * q
    a, b = math.floor(index), math.ceil(index)
    return ordered[a] + (ordered[b] - ordered[a]) * (index - a)


def common_eligible(row):
    return truth(row, 'survival_eligible') and truth(row, 'insulin_covered')


def eligible(row, threshold):
    return common_eligible(row) and (threshold is None or (
        not truth(row, 'expression_ambiguous') and math.isfinite(number(row, 'baseMean'))
        and number(row, 'baseMean') >= threshold))


def favoring(row, fdr=.1):
    return truth(row, 'survival_eligible') and number(row, 'survival_beta') > 0 and number(row, 'survival_fdr') < fdr


def rebuild(ledger, name):
    threshold = {'no_expression': None, 'expression_ge1': 1, 'expression_ge100': 100}.get(name, 10)
    fdr = .05 if name == 'survival_fdr005' else .1
    universe = [r for r in ledger if eligible(r, threshold)]
    logs = sorted(math.log2(number(r, 'baseMean') + 1) for r in universe) if threshold is not None else []
    bounds = sorted(set(quantile(logs, i / 5) for i in range(6))) if logs else []
    levels = [str(i) for i in range(1, len(bounds))] if bounds else ['0']
    labels = ['|'.join(x) for x in product(levels, ['3-5', '6+'], ['3', '4+'])]
    strata = {s: dict(n=0, m=0, k=0, a=0) for s in labels}
    genes = {}
    for row in universe:
        positive = favoring(row, fdr)
        negative = number(row, 'survival_beta') < 0 and number(row, 'survival_fdr') < .1
        if name == 'positive_vs_negative' and not (positive or negative):
            continue
        exposure = negative if name == 'survival_harming' else positive
        # Outcome is source Table S1 only. No stored guide FDR or support count is read.
        outcome = row['published_content_label'] == 'Decrease'
        expr = str(bisect.bisect_left(bounds[1:-1], math.log2(number(row, 'baseMean') + 1)) + 1) if bounds else '0'
        label = '|'.join((expr, '6+' if number(row, 'survival_guides') >= 6 else '3-5',
                          '4+' if number(row, 'covered_insulin_guides') >= 4 else '3'))
        genes[row['gene']] = (exposure, outcome, label)
        counts = strata[label]
        for key, amount in (('n', 1), ('m', exposure), ('k', outcome), ('a', exposure and outcome)):
            counts[key] += int(amount)
    return genes, strata, bounds


def evidence_status(row):
    if row['published_content_label'] == 'Decrease':
        return 'Published lower-content evidence'
    if row['published_content_label'] == 'Increase':
        return 'Published higher-content evidence'
    return 'No published directional call' if row['insulin_guide_source_symbol'] else 'Screen measurement unavailable'


def close(left, right, context, tol=1e-8):
    assert abs(left - right) < tol, (context, left, right)
    return abs(left - right)


def main():
    ledger_path = OUT / 'legacy_protocol/gene_eligibility_ledger.csv'
    ledger = read(ledger_path)
    lookup = {r['gene']: r for r in ledger}
    assert len(lookup) == len(ledger)
    actual = json.loads((OUT / 'results.json').read_text())
    assert actual['version'] == 2 and actual['semantic_status'] == 'GUIDE_SIGNIFICANCE_UNRESOLVED'
    assert [r['name'] for r in actual['analyses']] == NAMES
    definitions = json.loads((OUT / 'stratum_definitions.json').read_text())
    audit, expected_by_name = [], {}
    expected_margins = {
        'primary': (12184, 285, 479, 14), 'no_expression': (17276, 390, 514, 15),
        'expression_ge1': (13718, 311, 480, 14), 'expression_ge100': (10143, 238, 451, 14),
        'survival_fdr005': (12184, 99, 479, 4), 'survival_harming': (12184, 286, 479, 11),
        'positive_vs_negative': (571, 285, 25, 14),
    }
    for result in actual['analyses']:
        name = result['name']
        genes, strata, bounds = rebuild(ledger, name)
        expected_by_name[name] = genes
        observed_rows = read(OUT / f'{name}_genes.csv')
        observed_genes = {r['gene']: (truth(r, 'A'), truth(r, 'L'), r['stratum']) for r in observed_rows}
        assert len(observed_genes) == len(observed_rows), name
        assert genes == observed_genes, f'{name}: membership, A/L, or stratum assignment changed'
        observed_strata = {r['stratum']: {key: int(r[key]) for key in ('n', 'm', 'k', 'a')} for r in read(OUT / f'{name}_strata.csv')}
        assert strata == observed_strata, f'{name}: stratum cells changed'
        assert list(strata) == definitions['all_strata'][name]
        assert len(bounds) == len(definitions['log2_baseMean_plus1_boundaries'][name])
        for x, y in zip(bounds, definitions['log2_baseMean_plus1_boundaries'][name]):
            close(x, y, name + ' boundaries', 1e-12)
        margins = tuple(sum(s[k] for s in strata.values()) for k in ('n', 'm', 'k', 'a'))
        assert margins == expected_margins[name]
        assert margins == tuple(result[k] for k in ('n', 'A', 'L', 'overlap'))
        n, m, k, a = margins
        table = [[a, m-a], [k-a, n-m-k+a]]
        assert table == result['table_2x2']
        inference = scipy_conditional(strata)
        differences = {key: close(inference[key], result[key], (name, key))
                       for key in ('conditional_mle', 'ci_lower', 'ci_upper', 'p_enrichment')}
        for key in ('support_lower', 'support_upper', 'observed_total', 'n_informative_strata'):
            assert inference[key] == result[key], (name, key)
        co = odds_ratio(table, kind='conditional')
        ci = co.confidence_interval(.95)
        for key, value in {'conditional_mle': co.statistic, 'ci_lower': ci.low, 'ci_upper': ci.high,
                           'p_enrichment': fisher_exact(table, alternative='greater').pvalue}.items():
            close(float(value), result['unstratified'][key], (name, 'unstratified', key))
        assert result['conclusion'] == 'unresolved' and inference['ci_lower'] < 1 < inference['ci_upper']
        assert 'guide FDR significance not reconstructed' in result['interpretation']
        audit.append({'name': name, 'n': n, 'A': m, 'L': k, 'overlap': a, **inference,
                      'all_gene_labels_and_stratum_cells_match': True, 'max_abs_inference_difference': max(differences.values())})

    old_primary = {r['gene']: (truth(r, 'A'), truth(r, 'L'), r['stratum'])
                   for r in read(OUT / 'legacy_protocol/primary_genes.csv')}
    assert expected_by_name['primary'] == old_primary
    old_noexpr = {r['gene'] for r in read(OUT / 'legacy_protocol/no_expression_genes.csv')}
    assert set(expected_by_name['no_expression']) - old_noexpr == {'TOMM6'}
    assert not old_noexpr - set(expected_by_name['no_expression'])
    assert not truth(lookup['TOMM6'], 'direction_audit_pass')
    assert expected_by_name['no_expression']['TOMM6'][:2] == (False, True)
    assert 'TOMM6' not in expected_by_name['primary'] and number(lookup['TOMM6'], 'baseMean') == 0
    # Poison historical support diagnostics. The revised definition must not use them.
    poisoned = [{**r, 'direction_audit_pass': 'False', 'stronger_direction_audit_pass': 'False',
                 **{f'{sign}_support_{rep}': '0' for sign in ('positive', 'negative') for rep in ('r0', 'r1')}} for r in ledger]
    assert all(rebuild(poisoned, name)[0] == expected_by_name[name] for name in NAMES)

    expected_candidates = {r['gene']: r for r in ledger if favoring(r)}
    candidate_rows = read(OUT / 'candidate_table_427.csv')
    candidates = {r['gene']: r for r in candidate_rows}
    assert len(candidate_rows) == len(candidates) == 427
    assert set(candidates) == set(expected_candidates)
    assert 'candidate_class' not in candidate_rows[0]
    assert all(f'{sign}_support_{rep}' not in candidate_rows[0] for sign in ('positive', 'negative') for rep in ('r0', 'r1'))
    counters = Counter()
    for gene, expected in expected_candidates.items():
        observed = candidates[gene]
        status = evidence_status(expected)
        primary = eligible(expected, 10)
        assert observed['screen_evidence_status'] == status, gene
        assert truth(observed, 'primary_eligible') == truth(observed, 'in_U') == primary, gene
        assert observed['primary_eligibility_status'] == ('Eligible' if primary else 'Excluded'), gene
        assert 'S8 significance invalidated by Int64 coercion' in observed['guide_significance_status'] and 'selected S7' in observed['guide_significance_status'], gene
        assert truth(observed, 'L') == (expected['published_content_label'] == 'Decrease'), gene
        assert 'failed_published_direction_audit' not in observed['eligibility_reason'], gene
        counters[(status, observed['primary_eligibility_status'])] += 1
    accounting = json.loads((OUT / 'candidate_accounting.json').read_text())
    status_counts = Counter(evidence_status(r) for r in expected_candidates.values())
    assert dict(status_counts) == accounting['screen_evidence_status_counts']
    assert counters == Counter({(r['screen_evidence_status'], r['primary_eligibility_status']): r['n'] for r in accounting['evidence_by_primary_eligibility']})
    assert sum(common_eligible(r) for r in expected_candidates.values()) == accounting['common_coverage'] == 390
    assert sum(eligible(r, 10) for r in expected_candidates.values()) == accounting['primary_eligible'] == 285
    assert accounting['coverage_losses'] == 37 and accounting['additional_primary_losses'] == 105
    assert {g for g, r in expected_candidates.items() if eligible(r, 10) and r['published_content_label'] == 'Decrease'} == set(accounting['primary_lower_content_candidates'])
    assert accounting['excluded_with_published_lower_content'] == [{'gene': 'NUPL1', 'eligibility_reason': 'no_exact_expression_symbol'}]
    assert evidence_status(lookup['NUPL1']) == candidates['NUPL1']['screen_evidence_status'] == 'Published lower-content evidence'
    assert common_eligible(lookup['NUPL1']) and not eligible(lookup['NUPL1'], 10)
    assert not lookup['NUPL1']['expression_source_symbol']
    assert expected_by_name['no_expression']['NUPL1'][:2] == (True, True)
    assert 'NUPL1' not in expected_by_name['primary']
    expected_order = sorted(expected_candidates, key=lambda gene: (
        {'Published lower-content evidence': 0, 'Published higher-content evidence': 1,
         'No published directional call': 2, 'Screen measurement unavailable': 3}[evidence_status(expected_candidates[gene])],
        number(expected_candidates[gene], 'survival_fdr'), gene))
    assert [r['gene'] for r in candidate_rows] == expected_order

    report = {
        'status': 'PASS', 'scope': 'Independent selection, strata, candidate evidence/eligibility axes, and conditional arithmetic from archived source ledger. Raw-source reconstruction is a separate audit. No guide-FDR semantics are validated.',
        'production_modules_imported': [], 'inference_helper': 'tests/run_statistical_crosscheck.py::scipy_conditional (independent SciPy PMFs/Brent roots)',
        'all_seven_analyses': audit,
        'primary_gene_membership_labels_and_strata_unchanged_from_v1': True,
        'no_expression_delta_vs_v1': {'added_gene': 'TOMM6', 'A': False, 'published_L': True, 'legacy_direction_rule_passed': False, 'removed_genes': []},
        'historical_guide_support_poisoning_leaves_expected_catalogs_unchanged': True,
        'candidate_count': 427, 'candidate_status_counts': dict(status_counts),
        'NUPL1': {'published_evidence': 'Published lower-content evidence', 'common_coverage': True, 'primary_eligible': False, 'reason': 'no_exact_expression_symbol', 'included_as_overlap_without_expression_filter': True},
        'hashes': {str(path.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in [ledger_path, Path(__file__), ROOT / 'tests/run_statistical_crosscheck.py', ROOT / 'code/revise_analysis.py', OUT / 'results.json']},
    }
    destination = ROOT / 'evidence/revision_catalog_validation.json'
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'analyses': len(audit), 'largest_inference_difference': max(a['max_abs_inference_difference'] for a in audit),
                      'no_expression': next(a for a in audit if a['name'] == 'no_expression'), 'candidate_status_counts': dict(status_counts)}, indent=2))


if __name__ == '__main__':
    main()
