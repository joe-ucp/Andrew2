"""Rebuild stable-ID sensitivity strata and independently check exact inference.

Mapping/source reconstruction is audited separately. This test validates the
exported gene table, all strata, margins and conditional inference without
importing the production analysis, mapping code or exact estimator.
"""
from itertools import product
from pathlib import Path
import bisect
import csv
import hashlib
import json
import math

from run_statistical_crosscheck import scipy_conditional

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / 'results/identifier_sensitivity'


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def quantile(ordered, probability):
    index = (len(ordered) - 1) * probability
    a, b = math.floor(index), math.ceil(index)
    return ordered[a] + (ordered[b] - ordered[a]) * (index - a)


def main():
    gene_path = FOLDER / 'stable_id_primary_genes.csv'
    stratum_path = FOLDER / 'stable_id_primary_strata.csv'
    source_result_path = ROOT / 'evidence/identifier_mapping_sensitivity.json'
    genes = read(gene_path)
    deposited = json.loads(source_result_path.read_text(encoding='utf-8'))
    assert len({g['HGNC_id'] for g in genes}) == len(genes), 'Duplicate HGNC IDs'
    assert len({g['survival_source_symbol'] for g in genes}) == len(genes), 'Duplicate survival-source rows'
    assert len({g['insulin_source_symbol'] for g in genes}) == len(genes), 'Duplicate insulin-source rows'
    logs = sorted(math.log2(float(g['baseMean']) + 1) for g in genes)
    bounds = sorted(set(quantile(logs, i / 5) for i in range(6)))
    reported_bounds = deposited['log2_baseMean_plus1_quintile_boundaries']
    assert len(bounds) == len(reported_bounds)
    assert all(abs(a-b) < 1e-12 for a, b in zip(bounds, reported_bounds))
    strata = {'|'.join(items): dict(n=0, m=0, k=0, a=0)
              for items in product([str(i) for i in range(1, len(bounds))], ['3-5', '6+'], ['3', '4+'])}
    for g in genes:
        assert math.isfinite(float(g['baseMean'])) and float(g['baseMean']) >= 10
        assert float(g['survival_guides']) >= 3 and int(g['covered_insulin_guides']) >= 3
        assert g['A'] in ('True', 'False') and g['L'] in ('True', 'False')
        a, l = g['A'] == 'True', g['L'] == 'True'
        assert l == (g['published_content_label'] == 'Decrease')
        expression_bin = bisect.bisect_left(bounds[1:-1], math.log2(float(g['baseMean']) + 1)) + 1
        label = '|'.join((str(expression_bin), '6+' if float(g['survival_guides']) >= 6 else '3-5',
                          '4+' if int(g['covered_insulin_guides']) >= 4 else '3'))
        assert label == g['stratum'], g['HGNC_id']
        for key, amount in (('n', 1), ('m', a), ('k', l), ('a', a and l)):
            strata[label][key] += int(amount)
    exported = {r['stratum']: {k: int(r[k]) for k in ('n', 'm', 'k', 'a')} for r in read(stratum_path)}
    assert strata == exported
    margins = {key: sum(s[key] for s in strata.values()) for key in ('n', 'm', 'k', 'a')}
    assert margins == {'n': 12826, 'm': 306, 'k': 510, 'a': 14}
    assert tuple(margins.values()) == tuple(deposited[key] for key in ('n', 'A', 'L', 'overlap'))
    n, m, k, a = margins.values()
    assert deposited['table_2x2'] == [[a, m-a], [k-a, n-m-k+a]]
    expected_overlap = {g['HGNC_id'] for g in genes if g['A'] == 'True' and g['L'] == 'True'}
    assert expected_overlap == {g['HGNC_id'] for g in deposited['overlap_rows']}
    result = scipy_conditional(strata)
    differences = {key: abs(result[key] - deposited['conditional'][key])
                   for key in ('conditional_mle', 'ci_lower', 'ci_upper', 'p_enrichment')}
    assert max(differences.values()) < 1e-8
    for key in ('support_lower', 'support_upper', 'observed_total', 'n_informative_strata'):
        assert result[key] == deposited['conditional'][key]
    assert result['ci_lower'] < 1 < result['ci_upper']
    paths = [gene_path, stratum_path, source_result_path, Path(__file__), ROOT / 'tests/run_statistical_crosscheck.py']
    report = {
        'status': 'PASS',
        'scope': 'Independent reconstruction of strata from the exported stable-ID gene table and independent conditional arithmetic. Identifier mapping and raw-source reconstruction are separately audited. This post-review sensitivity does not alter the exact-symbol primary.',
        'unique_gene_and_source_identifiers': True,
        'all_gene_strata_and_cells_match': True,
        'manual_quintile_boundaries_match': True,
        'published_only_outcome_labels_match': True,
        'margins': margins,
        'independent_scipy': result,
        'max_abs_difference': max(differences.values()),
        'strata_sha256': hashlib.sha256(stratum_path.read_bytes()).hexdigest(),
        'hashes': {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    }
    output = ROOT / 'evidence/identifier_sensitivity_statistical_crosscheck.json'
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'margins': margins, 'independent_scipy': result,
                      'max_abs_difference': max(differences.values())}, indent=2))


if __name__ == '__main__':
    main()
