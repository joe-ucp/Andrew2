"""Post-result diagnostics; does not change the frozen primary analysis.

Read only primary_genes.csv and primary_strata.csv. Conditional calculations
use SciPy Fisher noncentral hypergeometric PMFs and Brent root inversion,
independently of exact_stats.py. Continuous contrasts do not use L labels.
Run from any directory: python code/exploratory_diagnostics.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.stats import hypergeom, nchypergeom_fisher

ROOT = Path(__file__).resolve().parents[1]
CAUTION = (
    'Post-result exploratory diagnostic. The frozen binary labels are retained '
    'for arithmetic only; guide-FDR audit semantics remain unvalidated. '
    'These calculations do not establish biological replication or a survival-function tradeoff.'
)


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_inputs():
    with (ROOT / 'results/primary_strata.csv').open(newline='', encoding='utf-8-sig') as stream:
        strata = list(csv.DictReader(stream))
    for s in strata:
        for key in ('n', 'm', 'k', 'a'):
            s[key] = int(s[key])
        s['lo'], s['hi'] = max(0, s['m'] + s['k'] - s['n']), min(s['m'], s['k'])
        assert s['lo'] <= s['a'] <= s['hi']
    with (ROOT / 'results/primary_genes.csv').open(newline='', encoding='utf-8-sig') as stream:
        genes = list(csv.DictReader(stream))
    assert len({g['gene'] for g in genes}) == len(genes), 'Duplicate primary genes'
    for g in genes:
        assert g['A'] in ('True', 'False') and g['L'] in ('True', 'False')
        g['A'], g['L'] = g['A'] == 'True', g['L'] == 'True'
        g['values'] = [float(g['median_lfc_r0']), float(g['median_lfc_r1'])]
        assert all(math.isfinite(x) for x in g['values']), 'Nonfinite primary LFC'
    assert set(g['stratum'] for g in genes) == set(s['stratum'] for s in strata)
    for s in strata:
        rows = [g for g in genes if g['stratum'] == s['stratum']]
        counts = (len(rows), sum(g['A'] for g in rows), sum(g['L'] for g in rows),
                  sum(g['A'] and g['L'] for g in rows))
        assert counts == tuple(s[k] for k in ('n', 'm', 'k', 'a')), (s['stratum'], counts)
    return strata, genes


class ConditionalModel:
    """Independent product of fixed-margin Fisher noncentral hypergeometrics."""

    def __init__(self, strata):
        self.strata = strata
        self.lo, self.hi = sum(s['lo'] for s in strata), sum(s['hi'] for s in strata)
        self.support = np.arange(self.lo, self.hi + 1)
        self.observed = sum(s['a'] for s in strata)

    def pmf(self, theta):
        result = np.ones(1)
        for s in self.strata:
            x = np.arange(s['lo'], s['hi'] + 1)
            if len(x) == 1:
                p = np.ones(1)
            elif theta == 1:
                p = hypergeom.pmf(x, s['n'], s['k'], s['m'])
            else:
                p = nchypergeom_fisher.pmf(x, s['n'], s['k'], s['m'], theta)
            assert np.all(np.isfinite(p)) and np.all(p >= 0)
            p /= p.sum()
            result = np.convolve(result, p)
        result /= result.sum()
        return result

    def inference(self):
        t, lo, hi = self.observed, self.lo, self.hi
        if lo == hi:
            return {'conditional_mle': None, 'ci_lower': 0.0, 'ci_upper': 'Infinity', 'p_enrichment': 1.0}
        def solve(function):
            return float(math.exp(brentq(lambda z: function(self.pmf(math.exp(z))), -25, 25, xtol=1e-12)))
        mle = 0.0 if t == lo else ('Infinity' if t == hi else solve(lambda p: np.dot(self.support, p) - t))
        lower = 0.0 if t == lo else solve(lambda p: p[self.support >= t].sum() - .025)
        upper = 'Infinity' if t == hi else solve(lambda p: p[self.support <= t].sum() - .025)
        return {'conditional_mle': mle, 'ci_lower': lower, 'ci_upper': upper,
                'p_enrichment': float(self.pmf(1)[self.support >= t].sum())}


def conditional_diagnostics(strata):
    model = ConditionalModel(strata)
    p0 = model.pmf(1)
    tails = np.cumsum(p0[::-1])[::-1]
    eligible = np.flatnonzero(tails < .025)
    threshold = int(model.support[eligible[0]]) if len(eligible) else None
    expectation = sum(s['m'] * s['k'] / s['n'] for s in strata)
    variance = sum(s['m'] * (s['k'] / s['n']) * (1 - s['k'] / s['n']) *
                   (s['n'] - s['m']) / (s['n'] - 1) for s in strata if s['n'] > 1)
    assert abs(np.dot(model.support, p0) - expectation) < 1e-9
    assert abs(np.dot((model.support - expectation) ** 2, p0) - variance) < 1e-8
    probabilities = []
    for theta in (1.0, 1.5, 2.0):
        p = model.pmf(theta)
        probabilities.append({
            'common_odds_ratio': theta,
            'expected_total_overlap': float(np.dot(model.support, p)),
            'criterion_probability': float(p[model.support >= threshold].sum()) if threshold is not None else 0.0,
        })
    primary = model.inference()
    assert primary['ci_lower'] < 1 < primary['ci_upper']
    influence = []
    for s in strata:
        if s['lo'] == s['hi']:
            continue
        reduced = ConditionalModel([v for v in strata if v['stratum'] != s['stratum']])
        inf = reduced.inference()
        influence.append({
            'omitted_stratum': s['stratum'], 'omitted_n': s['n'], 'omitted_A': s['m'],
            'omitted_L': s['k'], 'omitted_overlap': s['a'],
            'omitted_null_expected_overlap': s['m'] * s['k'] / s['n'],
            'remaining_observed_overlap': reduced.observed, 'remaining_n': sum(v['n'] for v in reduced.strata),
            **inf,
            'log_OR_change_vs_primary': math.log(inf['conditional_mle'] / primary['conditional_mle']),
            'exploratory_only': True,
        })
    info = {
        'status': 'post_result_exploratory', 'caution': CAUTION,
        'definition': 'T is total overlap across independent fixed-margin strata; a common theta tilts each stratum Fisher noncentral hypergeometric distribution.',
        'criterion': 'Equal-tail exact 95% lower limit > 1, equivalently inclusive P_theta=1(T >= observed) < 0.025. The inequality is strict.',
        'probability_interpretation': 'Model criterion probability for these already observed fixed margins. This is not prospective experimental, donor-level, or clinical power and does not validate the labels or the common-OR model.',
        'n_strata': len(strata), 'n_informative_strata': len(influence),
        'n_genes': sum(s['n'] for s in strata), 'n_exposed': sum(s['m'] for s in strata),
        'n_labeled_outcome': sum(s['k'] for s in strata), 'observed_overlap': model.observed,
        'null_expected_overlap': expectation, 'null_variance_overlap': variance,
        'null_sd_overlap': math.sqrt(variance),
        'minimum_overlap_for_criterion': threshold,
        'additional_overlaps_relative_to_observed': threshold - model.observed if threshold is not None else None,
        'null_upper_tail_at_threshold': float(tails[threshold - model.lo]) if threshold is not None else None,
        'null_upper_tail_one_below_threshold': float(tails[threshold - model.lo - 1]) if threshold is not None and threshold > model.lo else None,
        'criterion_probabilities': probabilities, 'primary_independent_inference': primary,
        'null_distribution': [{'total_overlap': int(x), 'probability': float(p), 'inclusive_upper_tail': float(q)}
                              for x, p, q in zip(model.support, p0, tails)],
        'influence_interpretation': 'All informative strata are omitted one at a time without redefining remaining strata or margins. Diagnostic selection or a preferred omitted stratum is not justified.',
    }
    return info, influence, model


def continuous_diagnostics(strata, genes, n_boot, seed):
    total_A = sum(s['m'] for s in strata)
    active = [s for s in strata if 0 < s['m'] < s['n']]
    assert sum(s['m'] for s in active) == total_A, 'An exposed stratum lacks a comparator'
    weights = np.array([s['m'] / total_A for s in active])
    assert abs(weights.sum() - 1) < 1e-12
    rng = np.random.default_rng(seed)
    point_A, point_B, point_var = np.zeros(2), np.zeros(2), np.zeros(2)
    boot_A, boot_B, boot_var = np.zeros((n_boot, 2)), np.zeros((n_boot, 2)), np.zeros((n_boot, 2))
    details = []
    for s, w in zip(active, weights):
        groups = [np.array([g['values'] for g in genes if g['stratum'] == s['stratum'] and g['A'] == a]) for a in (True, False)]
        exposed, control = groups
        mean_A, mean_B, var_B = exposed.mean(axis=0), control.mean(axis=0), control.var(axis=0, ddof=1)
        point_A += w * mean_A
        point_B += w * mean_B
        point_var += w * var_B
        for group, accumulator, is_control in ((exposed, boot_A, False), (control, boot_B, True)):
            # Gene vectors retain paired source r0/r1 values; strata and group sizes stay fixed.
            for start in range(0, n_boot, 100):
                stop = min(start + 100, n_boot)
                indices = rng.integers(0, len(group), size=(stop - start, len(group)))
                sampled = group[indices]
                accumulator[start:stop] += w * sampled.mean(axis=1)
                if is_control:
                    boot_var[start:stop] += w * sampled.var(axis=1, ddof=1)
        for j, rep in enumerate(('r0', 'r1')):
            details.append({
                'stratum': s['stratum'], 'replicate_label': rep, 'n_A1': len(exposed), 'n_A0': len(control),
                'target_weight': float(w), 'mean_lfc_A1': float(mean_A[j]), 'mean_lfc_A0': float(mean_B[j]),
                'mean_difference': float(mean_A[j] - mean_B[j]),
                'weighted_contribution': float(w * (mean_A[j] - mean_B[j])),
                'sample_sd_A0': float(math.sqrt(var_B[j])),
                'singleton_exposure_no_bootstrap_variation': len(exposed) == 1,
            })
    contrast, boot_contrast = point_A - point_B, boot_A - boot_B
    denominator, boot_denominator = np.sqrt(point_var), np.sqrt(boot_var)
    assert np.all(denominator > 0) and np.all(boot_denominator > 0)
    standardized, boot_standardized = contrast / denominator, boot_contrast / boot_denominator
    rows = []
    for j, rep in enumerate(('r0', 'r1')):
        ci = np.quantile(boot_contrast[:, j], [.025, .975], method='linear')
        sci = np.quantile(boot_standardized[:, j], [.025, .975], method='linear')
        rows.append({
            'replicate_label': rep, 'n_genes_total': len(genes), 'n_exposed_target': total_A,
            'n_comparator_genes_in_weighted_strata': sum(s['n'] - s['m'] for s in active),
            'n_weighted_strata': len(active), 'weighted_mean_lfc_A1': float(point_A[j]),
            'exposure_standardized_mean_lfc_A0': float(point_B[j]),
            'exposure_standardized_mean_difference': float(contrast[j]),
            'gene_bootstrap_percentile_025': float(ci[0]), 'gene_bootstrap_percentile_975': float(ci[1]),
            'weighted_within_A0_sd': float(denominator[j]),
            'sd_scaled_mean_difference': float(standardized[j]),
            'sd_scaled_gene_bootstrap_percentile_025': float(sci[0]),
            'sd_scaled_gene_bootstrap_percentile_975': float(sci[1]),
            'bootstrap_draws': n_boot, 'bootstrap_seed': seed,
        })
    report = {
        'status': 'post_result_exploratory', 'caution': CAUTION,
        'outcome': 'Per-gene median source insulin LFC, separately for source labels r0 and r1. All finite values are used; no L label or guide FDR cutoff is used in the contrast.',
        'estimand': 'Delta_r = sum_s w_s [mean(Y_r | A=1,s) - mean(Y_r | A=0,s)], with fixed w_s=m_s/sum_s m_s. This standardizes A0 to the observed exposed-gene stratum composition and is descriptive, not causal.',
        'target': 'The 285 frozen survival-positive genes represented by the primary exposure labels, compared with A0 genes in the same frozen strata. The finite-catalog contrast is computed exactly; its resampling interval invokes a hypothetical exchangeable-gene reference.',
        'sign': 'Positive delta means larger low-versus-high insulin LFC among A1 genes after this standardization; in the source orientation this corresponds to lower intracellular insulin content. It does not measure secretion or show a biological tradeoff.',
        'sd_scaling': 'D_r = Delta_r / sqrt(sum_s w_s sample_variance(Y_r | A=0,s)). This descriptive scaling uses within-stratum A0 spread; it is not a meta-analytic or bias-corrected standardized mean difference.',
        'bootstrap': {
            'draws': n_boot, 'seed': seed, 'rng': 'NumPy default_rng (PCG64)',
            'method': 'Resample gene rows with replacement separately within each A-by-stratum cell, retaining each selected gene\'s r0/r1 vector together, and retaining each cell size. Fixed original weights. Recalculate control variance for SD-scaled contrasts in every draw. Report linear-interpolated 2.5th/97.5th percentiles.',
            'limitations': 'Intervals represent an exchangeable-gene resampling reference, not donor or biological-replicate variability, guide measurement error, or source-selection/FDR uncertainty. Pathway/gene dependence is not modeled. A singleton cell has zero empirical resampling variability; sparse exposed cells make these intervals fragile. No hypothesis-test p-values or significance claims are made.',
        },
        'unweighted_strata': [{'stratum': s['stratum'], 'n_A1': s['m'], 'n_A0': s['n'] - s['m'], 'target_weight': 0.0} for s in strata if s not in active],
        'sparse_exposed_cells_n_le_6': [{'stratum': s['stratum'], 'n_A1': s['m'], 'n_A0': s['n'] - s['m'], 'target_weight': s['m'] / total_A} for s in active if s['m'] <= 6],
        'summary': rows, 'stratum_details': details,
    }
    return report, rows, details


def figures(info, model, continuous, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(8.6, 4.7))
    cutoff = info['minimum_overlap_for_criterion']
    x = model.support
    colors = ['#313d4b', '#287e92', '#b15b24']
    for row, color in zip(info['criterion_probabilities'], colors):
        theta = row['common_odds_ratio']
        p = model.pmf(theta)
        ax.plot(x, p, marker='o', ms=3, color=color,
                label=f"OR {theta:g}: criterion probability {row['criterion_probability']:.1%}")
    ax.axvspan(cutoff - .5, 45, color='#d9e8d4', alpha=.55, label=f'Criterion: total overlap ≥ {cutoff}')
    ax.axvline(info['observed_overlap'], color='#8c3b4a', ls='--', lw=1.3, label=f"Observed: {info['observed_overlap']}")
    ax.set(xlim=(-.5, 40.5), ylim=(0, None), xlabel='Total binary overlap T across frozen strata', ylabel='Conditional probability',
           title='Fixed-margin conditional information (post-result diagnostic)')
    ax.legend(frameon=False, fontsize=8, loc='upper right')
    fig.text(.5, .015, 'Model probabilities only; guide-FDR label semantics remain unvalidated.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(out / 'conditional_information.png', dpi=200)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0))
    for j, row in enumerate(continuous['summary']):
        for ax, key, low, high in (
            (axes[0], 'exposure_standardized_mean_difference', 'gene_bootstrap_percentile_025', 'gene_bootstrap_percentile_975'),
            (axes[1], 'sd_scaled_mean_difference', 'sd_scaled_gene_bootstrap_percentile_025', 'sd_scaled_gene_bootstrap_percentile_975')):
            value = row[key]
            ax.errorbar(value, 1-j, xerr=[[value-row[low]], [row[high]-value]], fmt='o', color=colors[j+1], capsize=4)
    for ax in axes:
        ax.axvline(0, color='#7a7a7a', lw=1, ls='--')
        ax.set(yticks=[1, 0], yticklabels=['Source r0', 'Source r1'], ylim=(-.65, 1.65))
        ax.grid(axis='x', alpha=.15)
    axes[0].set(xlabel='Weighted mean difference (source LFC units)', title='Exposure-stratum-standardized contrast')
    axes[1].set(xlabel='Contrast / weighted within-A0 SD', title='Descriptive SD-scaled contrast')
    fig.suptitle('Continuous gene-level insulin LFC comparison', y=.98, fontsize=13)
    fig.text(.5, .04, f"Bars: {continuous['bootstrap']['draws']:,} gene-bootstrap percentile intervals; not biological-replicate uncertainty.\nPost-result exploratory; fixed exposure weights; no guide-FDR outcome labels used.", ha='center', fontsize=8.5)
    fig.tight_layout(rect=(0, .13, 1, .94))
    fig.savefig(out / 'continuous_summary.png', dpi=200)
    plt.close(fig)


def documentation(info, influence, continuous):
    probabilities = {r['common_odds_ratio']: r['criterion_probability'] for r in info['criterion_probabilities']}
    lines = [
        '# Exploratory diagnostics for the frozen primary gene universe', '',
        CAUTION, '',
        'These diagnostics were requested after the primary result was known. They are not part of the frozen confirmatory analysis and do not alter any of its nine analyses, universe definitions, thresholds, or candidate labels.', '',
        '## Conditional information', '',
        f"The frozen primary margins contain {info['n_genes']:,} genes, {info['n_exposed']} exposure labels, {info['n_labeled_outcome']} outcome labels, and {info['observed_overlap']} overlaps in {info['n_informative_strata']} informative strata. Under theta=1, the expected total overlap is **{info['null_expected_overlap']:.6f}**, with SD {info['null_sd_overlap']:.6f}. This is the sum of m_s k_s / n_s, checked against the independently convolved null distribution.", '',
        f"For the prespecified criterion (exact equal-tail 95% lower OR bound >1), the minimum total overlap at these margins is **{info['minimum_overlap_for_criterion']}**. Its inclusive null upper tail is {info['null_upper_tail_at_threshold']:.8f}; one overlap below gives {info['null_upper_tail_one_below_threshold']:.8f}. The strict cutoff is p<0.025. The observed total is {info['additional_overlaps_relative_to_observed']} below this threshold. These hypothetical counts preserve all margins; they are not a claim that adding a few genes would reproduce this scenario.", '',
        f"At a common OR of 1.5, the criterion probability is **{probabilities[1.5]:.2%}**; at OR 2 it is **{probabilities[2.0]:.2%}**. The realized null criterion probability is {probabilities[1.0]:.2%}. These are probabilities under a common-OR fixed-margin conditional model of the already observed design, not prospective experimental, donor-level, or clinical power. Label validity, gene dependence, common-OR heterogeneity, and sampling of the original margins are not represented.", '',
        'Inference is independently implemented with SciPy Fisher noncentral hypergeometric PMFs, direct convolution, and Brent root inversion. It does not import the production exact estimator. Inclusive tails match the frozen interval definition; no mid-p or normal approximation is used.', '',
        '## Leave-one-informative-stratum-out influence', '',
        f"Every one of the {len(influence)} informative strata is omitted in turn, leaving all remaining margins and definitions fixed. The resulting conditional ORs range from {min(r['conditional_mle'] for r in influence):.4f} to {max(r['conditional_mle'] for r in influence):.4f}. All omissions and their exact intervals are provided in results/stratum_influence.csv. They are diagnostics, not alternative selected analyses. A preferred omission, multiplicity-adjusted discovery, or stability guarantee is not claimed.", '',
        '## Continuous estimand and resampling', '',
        continuous['outcome'], '', continuous['estimand'], '', continuous['target'], '', continuous['sign'], '',
        'Here Y is each gene\'s median LFC over its available insulin guides in the relevant source replicate label. The mean is then taken over genes; guides are not treated as independent observations. The per-replicate effects are kept separate. No threshold or weight was optimized for a favorable effect.', '',
        'All 285 exposed genes have same-stratum A0 comparators in 12 strata. Eight A0-only strata have zero target weight. This is an explicit exposed-gene target, not the mean effect over all 12,184 genes. The control mean is standardized to this same target. Each stratum contribution, weight, and group size is disclosed in results/continuous_strata.csv and the JSON.', '',
        continuous['sd_scaling'], '',
        '| Source label | Exposure-standardized difference (LFC units) | Gene-bootstrap percentile interval | SD-scaled difference |',
        '|---|---:|---:|---:|',
    ]
    for row in continuous['summary']:
        lines.append(f"| {row['replicate_label']} | {row['exposure_standardized_mean_difference']:.6f} | [{row['gene_bootstrap_percentile_025']:.6f}, {row['gene_bootstrap_percentile_975']:.6f}] | {row['sd_scaled_mean_difference']:.6f} |")
    lines += ['', f"Bootstrap specification: {continuous['bootstrap']['draws']:,} draws, seed {continuous['bootstrap']['seed']}, {continuous['bootstrap']['rng']}. {continuous['bootstrap']['method']}", '',
              continuous['bootstrap']['limitations'], '',
              'The 1|3-5|4+ exposure cell contains one gene. Its empirical bootstrap mean is fixed in every resample; this does not establish that its true variability is zero. Other small exposure cells (n<=6), with their weights, are listed explicitly in the JSON.', '',
              'The continuous comparison avoids the disputed binary guide-FDR outcome semantics, but retains the original survival exposure definition, source screens, mapping, guide medians, and gene universe. It does not validate the binary endpoint, turn an unresolved primary result into support, or establish a causal survival-versus-function tradeoff.', '',
              '## Reproduction and checks', '',
              'Run `python code/exploratory_diagnostics.py` from any directory with NumPy, SciPy, and matplotlib installed. Defaults are 2,000 bootstrap draws and seed 20261023; command-line overrides are recorded in output. No network access is used. The script reads only the frozen primary gene and stratum files and writes its named diagnostic outputs. SHA256 input and code hashes and package versions are recorded in both JSON outputs.', '',
              'Before calculation, the script verifies unique gene identifiers, finite paired LFCs, legal binary labels, and exact reconstruction of every stratum\'s n,m,k,a. It checks the convolved null mean and variance against analytic hypergeometric formulas, verifies normalized target weights, and retains paired r0/r1 vectors during resampling. Statistical independence here means a separate computational implementation; it is not independent data or independent biological evidence.', '',
              'Run `python tests/validate_exploratory_diagnostics.py` after generation to check all 13 conditional fits against the original recurrence-based estimator, three full conditional distributions, the exact lower-bound criterion at totals 18 and 19, continuous contrasts by independent weighted-sum identities, exact seeded bootstrap reproducibility, and singleton-cell disclosure. Results are saved to evidence/exploratory_diagnostics_validation.json.', '']
    path = ROOT / 'docs/exploratory_diagnostics.md'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bootstrap-draws', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=20261023)
    args = parser.parse_args()
    if args.bootstrap_draws < 100:
        parser.error('At least 100 bootstrap draws are required')
    strata, genes = read_inputs()
    info, influence, model = conditional_diagnostics(strata)
    continuous, rows, detail = continuous_diagnostics(strata, genes, args.bootstrap_draws, args.seed)
    provenance = {
        'inputs': {name: sha(ROOT / 'results' / name) for name in ('primary_genes.csv', 'primary_strata.csv')},
        'code_sha256': sha(Path(__file__)), 'numpy': np.__version__, 'scipy': scipy.__version__,
        'checks': ['Unique genes and finite paired LFCs', 'Every n,m,k,a rebuilt from primary genes',
                   'Convolved null mean and variance match analytic formulas', 'All exposed genes covered by positive fixed weights'],
    }
    info['provenance'] = continuous['provenance'] = provenance
    write_json(ROOT / 'results/conditional_information.json', info)
    write_json(ROOT / 'results/continuous_summary.json', continuous)
    write_csv(ROOT / 'results/stratum_influence.csv', influence)
    write_csv(ROOT / 'results/continuous_summary.csv', rows)
    write_csv(ROOT / 'results/continuous_strata.csv', detail)
    figures(info, model, continuous, ROOT / 'manuscript/figures')
    documentation(info, influence, continuous)
    print(json.dumps({'null_expected_overlap': info['null_expected_overlap'], 'minimum_overlap': info['minimum_overlap_for_criterion'],
                      'criterion_probabilities': info['criterion_probabilities'], 'continuous': rows}, indent=2))


if __name__ == '__main__':
    main()
