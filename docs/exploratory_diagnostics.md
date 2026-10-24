# Exploratory diagnostics for the frozen primary gene universe

Post-result exploratory diagnostic. The frozen binary labels are retained for arithmetic only; guide-FDR audit semantics remain unvalidated. These calculations do not establish biological replication or a survival-function tradeoff.

These diagnostics were requested after the primary result was known. They are not part of the frozen confirmatory analysis and do not alter any of its nine analyses, universe definitions, thresholds, or candidate labels.

## Conditional information

The frozen primary margins contain 12,184 genes, 285 exposure labels, 479 outcome labels, and 14 overlaps in 12 informative strata. Under theta=1, the expected total overlap is **11.518255**, with SD 3.276855. This is the sum of m_s k_s / n_s, checked against the independently convolved null distribution.

For the prespecified criterion (exact equal-tail 95% lower OR bound >1), the minimum total overlap at these margins is **19**. Its inclusive null upper tail is 0.02223366; one overlap below gives 0.04052227. The strict cutoff is p<0.025. The observed total is 5 below this threshold. These hypothetical counts preserve all margins; they are not a claim that adding a few genes would reproduce this scenario.

At a common OR of 1.5, the criterion probability is **31.33%**; at OR 2 it is **75.71%**. The realized null criterion probability is 2.22%. These are probabilities under a common-OR fixed-margin conditional model of the already observed design, not prospective experimental, donor-level, or clinical power. Label validity, gene dependence, common-OR heterogeneity, and sampling of the original margins are not represented.

Inference is independently implemented with SciPy Fisher noncentral hypergeometric PMFs, direct convolution, and Brent root inversion. It does not import the production exact estimator. Inclusive tails match the frozen interval definition; no mid-p or normal approximation is used.

## Leave-one-informative-stratum-out influence

Every one of the 12 informative strata is omitted in turn, leaving all remaining margins and definitions fixed. The resulting conditional ORs range from 0.9485 to 1.4958. All omissions and their exact intervals are provided in results/stratum_influence.csv. They are diagnostics, not alternative selected analyses. A preferred omission, multiplicity-adjusted discovery, or stability guarantee is not claimed.

## Continuous estimand and resampling

Per-gene median source insulin LFC, separately for source labels r0 and r1. All finite values are used; no L label or guide FDR cutoff is used in the contrast.

Delta_r = sum_s w_s [mean(Y_r | A=1,s) - mean(Y_r | A=0,s)], with fixed w_s=m_s/sum_s m_s. This standardizes A0 to the observed exposed-gene stratum composition and is descriptive, not causal.

The 285 frozen survival-positive genes represented by the primary exposure labels, compared with A0 genes in the same frozen strata. The finite-catalog contrast is computed exactly; its resampling interval invokes a hypothetical exchangeable-gene reference.

Positive delta means larger low-versus-high insulin LFC among A1 genes after this standardization; in the source orientation this corresponds to lower intracellular insulin content. It does not measure secretion or show a biological tradeoff.

Here Y is each gene's median LFC over its available insulin guides in the relevant source replicate label. The mean is then taken over genes; guides are not treated as independent observations. The per-replicate effects are kept separate. No threshold or weight was optimized for a favorable effect.

All 285 exposed genes have same-stratum A0 comparators in 12 strata. Eight A0-only strata have zero target weight. This is an explicit exposed-gene target, not the mean effect over all 12,184 genes. The control mean is standardized to this same target. Each stratum contribution, weight, and group size is disclosed in results/continuous_strata.csv and the JSON.

D_r = Delta_r / sqrt(sum_s w_s sample_variance(Y_r | A=0,s)). This descriptive scaling uses within-stratum A0 spread; it is not a meta-analytic or bias-corrected standardized mean difference.

| Source label | Exposure-standardized difference (LFC units) | Gene-bootstrap percentile interval | SD-scaled difference |
|---|---:|---:|---:|
| r0 | -0.039261 | [-0.107276, 0.031040] | -0.060170 |
| r1 | 0.034120 | [-0.103827, 0.163203] | 0.030953 |

Bootstrap specification: 2,000 draws, seed 20261023, NumPy default_rng (PCG64). Resample gene rows with replacement separately within each A-by-stratum cell, retaining each selected gene's r0/r1 vector together, and retaining each cell size. Fixed original weights. Recalculate control variance for SD-scaled contrasts in every draw. Report linear-interpolated 2.5th/97.5th percentiles.

Intervals represent an exchangeable-gene resampling reference, not donor or biological-replicate variability, guide measurement error, or source-selection/FDR uncertainty. Pathway/gene dependence is not modeled. A singleton cell has zero empirical resampling variability; sparse exposed cells make these intervals fragile. No hypothesis-test p-values or significance claims are made.

The 1|3-5|4+ exposure cell contains one gene. Its empirical bootstrap mean is fixed in every resample; this does not establish that its true variability is zero. Other small exposure cells (n<=6), with their weights, are listed explicitly in the JSON.

The continuous comparison avoids the disputed binary guide-FDR outcome semantics, but retains the original survival exposure definition, source screens, mapping, guide medians, and gene universe. It does not validate the binary endpoint, turn an unresolved primary result into support, or establish a causal survival-versus-function tradeoff.

## Reproduction and checks

Run `python code/exploratory_diagnostics.py` from any directory with NumPy, SciPy, and matplotlib installed. Defaults are 2,000 bootstrap draws and seed 20261023; command-line overrides are recorded in output. No network access is used. The script reads only the frozen primary gene and stratum files and writes its named diagnostic outputs. SHA256 input and code hashes and package versions are recorded in both JSON outputs.

Before calculation, the script verifies unique gene identifiers, finite paired LFCs, legal binary labels, and exact reconstruction of every stratum's n,m,k,a. It checks the convolved null mean and variance against analytic hypergeometric formulas, verifies normalized target weights, and retains paired r0/r1 vectors during resampling. Statistical independence here means a separate computational implementation; it is not independent data or independent biological evidence.

Run `python tests/validate_exploratory_diagnostics.py` after generation to check all 13 conditional fits against the original recurrence-based estimator, three full conditional distributions, the exact lower-bound criterion at totals 18 and 19, continuous contrasts by independent weighted-sum identities, exact seeded bootstrap reproducibility, and singleton-cell disclosure. Results are saved to evidence/exploratory_diagnostics_validation.json.
