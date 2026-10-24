# Data dictionary — version 2

## Version and interpretation

Top-level `results/` contains the revised published-catalogue comparison. `results/legacy_protocol/` preserves the original arithmetic and old field names as historical diagnostics within this consolidated project. A historical field named FDR or support does not establish valid guide significance. The complete Table S8 p/FDR fields underwent integer conversion, and the original significance interpretation is rejected. Continuous Table S7 values support only a selected subset. See `source_semantics.md`.

The primary comparison unit is the exact source gene symbol, trimmed only at surrounding whitespace. Primary matching does not change case, resolve aliases or use fuzzy matching. All source-row fields are **one-based physical CSV/worksheet rows including headers**. Missing numeric values mean unavailable evidence, not zero effects. False eligibility flags mean failure of the stated rule, not biological safety.

## Units, identity and signs

| Quantity | Meaning |
|---|---|
| `survival_beta` | Benaglio MAGeCK model beta: positive means relative knockout enrichment in cytokine-exposed versus untreated culture. It is not an absolute survival probability. |
| `survival_fdr` | Original Benaglio gene-level FDR, retained without recalculating it in the shared universe. |
| Insulin `LFC` | Guide log2 enrichment in low versus high intracellular insulin-content bins; positive means lower content following knockout. |
| Table S8 `FDR`, `stored_fdr_r0`, `stored_fdr_r1` | Stored integer values after the saved Power Query `Int64.Type` conversion. Zero does **not** establish FDR <0.10. Historical file column names are preserved where necessary for provenance. |
| Table S7 `S7_FDR` | Continuous FDR as stored in the selected published guide table. The subset does not provide genome-wide continuous probabilities or values for absent guides. |
| `control_count`, `treatment_count` | Source normalized count quantities. Their sum is used as a coverage rule here; it is not an absolute molecular or cell count. |
| `baseMean` | Table S5 proxy from the CALCOCO2-silencing/control RNA experiment. It is not an untreated reference expression or secretion measurement. |
| `published_mean_lfc` | Table S1 source mean log fold change, distinct from the per-replicate median of all valid guide LFCs. |

Survival beta and insulin LFC have different meanings and scales. They must not be subtracted or combined into an efficacy score. Distinct guides are not extra biological replicates.

## Inputs

| File or source | Relevant records |
|---|---|
| D1 `Supplementary_Data_09.CRISPR-KO-screen_results.csv` | 20,708 gene rows; header at row 1. Gene, sgRNA count and cytokine-survival beta/z/p-value/FDR. |
| D2 `41588_2022_1261_MOESM4_ESM.xlsx` | `Table S1-Screen hits ` (trailing space), `Table S5 - RNASeq`, `Table S7-sgRNA_hits`, `Table S8-Screen_all sgRNA`; headers at Excel row 1. |
| D3 `mmc6.xlsx` | Redundant survival check, sheet `KO screening results`, header row 2, data from row 3. Kill-curve replicates do not establish the number of pooled screen replicates. |
| C1 `16_CRISPR_figure.ipynb` | Author notebook pinned to commit `7d10581be640a2a8f4586c53a677b7c6cbbd7b1c`; implementation and denominator provenance. |
| HGNC complete set | Additional pinned post-review source used only for the separate identifier sensitivity. See `identifier_mapping.md`. |
| Author code and metadata | Attributed public snapshots in `data/source_documents`, inventoried in `evidence/source_document_manifest.json`. |

## Current gene ledger and candidate table

`gene_eligibility_ledger.csv` retains source matches, eligibility failures and current evidence axes. `candidate_table_427.csv` retains all 427 source survival-favoring candidates. The candidate table separates screen evidence from primary eligibility rather than combining them into one class.

| Field or group | Current meaning |
|---|---|
| `gene` | Trimmed exact-symbol join key. |
| `survival_source_symbol`, `insulin_guide_source_symbol`, `insulin_hit_source_symbol`, `expression_source_symbol` | Preserved source identifiers. A missing source symbol is not imputed. |
| `survival_source_row`, `hit_source_row`, `expression_source_row` | Physical source row in D1, S1 and S5 respectively. |
| `source_rows_r0`, `source_rows_r1` | Semicolon-separated physical S8 rows supporting the per-gene replicate records. |
| `survival_guides` | Source sgRNA count; eligibility requires at least three and excludes `hsa-` labels. |
| `survival_unique`, `survival_finite`, `survival_eligible` | Source-symbol uniqueness, finite-value and full survival-eligibility checks. |
| `covered_insulin_guides` | Distinct guide keys present once in both replicates, with finite LFC/stored FDR/counts and control+treatment sum >=20 in each replicate. At least three are required; EGFP, LacZ and luciferase are excluded controls. |
| `insulin_covered` | Passes that coverage rule, independent of guide significance. |
| `valid_guides_r0`, `valid_guides_r1` | Valid guide rows per replicate before the separate paired coverage threshold. |
| `median_lfc_r0`, `median_lfc_r1` | Median of all valid guide LFCs per gene and replicate. Includes guides without recoverable significance; descriptive only. |
| `stored_zero_positive_count_r0/r1`, `stored_zero_negative_count_r0/r1` | Counts with positive or negative LFC and stored S8 FDR=0. These are encoded-value diagnostics, not significant-guide counts. |
| `published_content_label` | Table S1 `Decrease`/`Increase`, otherwise `Not a published hit`. Never inferred from median LFC or stored-zero counts. |
| `legacy_stored_zero_rule_2_each`, `legacy_stored_zero_rule_3_each` | Historical published-direction checks using at least two/three stored-zero same-direction guides in each replicate. They do not require the same identities across replicates and do not gate current catalogue membership. |
| `same_positive_sign_guides`, `opposite_signs_guides` | Number of identical paired guide keys with both LFCs positive or with opposite signs. No significance meaning. |
| `same_positive_stored_zero_guides`, `same_negative_stored_zero_guides` | Same-key guide pairs with consistent positive/negative signs and both stored S8 FDR values zero. No calibrated FDR meaning. |
| `GeneID`, `baseMean`, `expression_ambiguous` | S5 identifier, numeric expression proxy and duplicate-name flag. All 75 duplicated source names remain ambiguous in exact-symbol primary matching. |
| `in_U0` | Common coverage before expression filtering: 17,276 genes. No S8 significance rule is applied. |
| `in_U` | Current primary eligibility: common coverage, unique finite expression and baseMean >=10; 12,184 genes. No stored-zero direction gate. |
| `primary_eligible`, `primary_eligibility_status` | Explicit primary-membership Boolean and `Eligible`/`Excluded` label. |
| `eligibility_reason` | `eligible` or semicolon-separated source-match, coverage, expression-ambiguity or expression-threshold exclusions. |
| `screen_evidence_status` | Published lower-content, published higher-content, no published directional call, or screen measurement unavailable. This axis is independent of expression eligibility. |
| `guide_significance_status` | Explicit warning that complete S8 guide significance is unresolved because of binary encoding. Selected continuous S7 evidence is recorded separately. |
| `A` | Primary survival-favoring label: beta>0 and original source FDR<0.10. Named exposure sensitivities redefine A as documented below. Primary A=0 includes nonsignificant and survival-harming genes. |
| `L` | Current published Table S1 `Decrease` label. L=0 means no published lower-content call, not no biological effect or safety. |
| `survival_harming` | Beta<0 and original FDR<0.10. |
| `external_functional_evidence_status` | Refers to bounded candidate annotations and makes explicit that matched survival/content/secretion validation was not performed. |

Across all 427 candidates there are 15 published lower-content, two published higher-content, 384 measured without a published directional call, and 26 without an exact screen gene measurement. Common coverage includes 390 and primary eligibility includes 285. Eleven measured but inadequately covered candidates remain in the no-call category. NUPL1 retains its lower-content label despite exact-symbol expression exclusion. Candidate ordering is evidence category, survival FDR, then source gene; it is not a therapeutic score.

## Guide and selected-S7 evidence

`guide_evidence.csv.gz` retains the legacy source schema: `guide_source_row`, `source_guide_id`, `insulin_guide_source_symbol`, `gene`, `guide_key`, `replicate`, `LFC`, `FDR`, counts, `valid`, `covered`. Its `FDR` is the stored S8 integer. `guide_key` removes only the terminal `_r0`/`_r1` and surrounding whitespace; coordinates and strand remain. `replicate` identifies the two insulin screen replicates.

`candidate_paired_guide_effects.csv` pairs the same guide key across r0/r1, preserving `lfc_r0/r1` and `stored_fdr_r0/r1`. The `same_positive_sign`, `opposite_signs`, `same_positive_stored_zero`, and `same_negative_stored_zero` flags describe signs and encoded values only. All LFCs remain visible, including opposing signs.

`evidence/candidate_S7_selected_guide_evidence.csv` records 14 primary candidates, NUPL1 and TOMM6. `source_s7_positive_fdr_lt010_r0/r1` uses continuous selected S7 FDR<0.10 and positive LFC; `le010` records <=0.10 for comparison. `source_s7_negative_fdr_lt010_r0/r1` is the negative counterpart. `source_s7_positive_fdr_lt010_sameguide` and its keys count identical positive guides across replicates. `S7_rows_r0/r1`, `S8_rows_r0/r1`, and `S7_covers_all_S8_rows_r0/r1` expose the selected subset's incompleteness. `source_s7_any_FDR_exactly010` checks ties. Three-guide flags apply only to this selected evidence. All 14 plus NUPL1 have two matching positive guides but only two of four source guides per replicate represented; no third-guide or genome-wide recovery follows.

`evidence/S7_S8_semantic_linkage.csv` anchors each of the 2,384 S7 rows to its unique S8 row by gene, coordinate, counts and LFC. S7/S8 FDR and other p-value columns remain separately named. `semantic_gene_diagnostics.csv` reports opposing signs and illustrative stored-zero rules without relabeling published decisions.

## Current statistical and exploratory outputs

| Field/file | Meaning |
|---|---|
| `expression_stratum` | Right-closed quintile of log2(baseMean+1), with linear boundaries, duplicate boundaries dropped and lowest endpoint included. Zero in no-expression analysis. |
| `survival_guide_stratum`, `insulin_guide_stratum` | Survival 3–5/6+ guides; insulin three/4+ covered guides. |
| `stratum` | Three component labels joined with `\|`; four possible no-expression strata, otherwise at most twenty. |
| `n`, `m`, `k`, `a` in strata CSVs | Stratum size, A=1 count, L=1 count, and joint count. Empty and noninformative strata are retained. |
| `A1_L0`, `A0_L1`, `A0_L0` | Other cells: m−a, k−a, n−m−k+a. |
| `informative` | Conditional support has more than one feasible joint count. |
| `conditional_mle`, `ci_lower`, `ci_upper` | Common-odds conditional estimate and inclusive equal-tail exact 95% interval. |
| `p_enrichment` | Inclusive upper-tail probability at odds ratio one. |
| `overlap` | Current result JSON count with A=1 and published L=1; old files use `conflicts`, whose biological interpretation is withdrawn. |
| `table_2x2` | Rows A=1/A=0; columns L=1/L=0. Primary: `[[14,271],[465,11434]]`. |
| `fraction_A1`, `fraction_A0` | Published lower-content fractions with catalogue denominators, not patient-risk probabilities. |
| `unstratified` | Cross-product odds ratio and separate conditional Fisher estimate/interval/tail probability. |
| `conclusion` | The interval-based catalogue result; current comparisons remain unresolved. It is not evidence of an individual's safety or clinical effect. |
| `continuous_summary.csv` | Separate replicate-specific exposure-count-weighted mean LFC contrasts, gene-bootstrap percentile bounds, comparator SD-scaled contrasts, draw count and seed. These are exploratory and are not source gene FDRs. |
| `continuous_strata.csv` | Per-stratum component counts and continuous summaries used by the weighted contrast. |
| `stratum_influence.csv` | Leave-one-stratum-out catalogue estimates and uncertainty, with omitted counts and change in log OR. Exploratory only. |
| `conditional_information.json` | Fixed-margin expected overlap, criterion probabilities and information diagnostics; not wet-lab power or independent data. |
| `threshold_free_stratum_diagnostics.csv` | Legacy descriptive r0/r1 within-stratum quartiles/medians; missing summaries mean an empty group. These original-primary distributions do not replace the amended catalogue endpoint. |

Intervals do not propagate all source-call uncertainty, donor variation or dependence between genes. Infinite endpoints are explicit; a non-estimable model is not a null result. Bootstrap intervals concern gene resampling under the stated design, not biological replicate or donor uncertainty.

## Comparison names and historical separation

The seven current comparisons are `primary`, `no_expression` (common coverage; retains TOMM6; coverage strata only), `expression_ge1`, `expression_ge100` (each rebuilds expression quintiles), `survival_fdr005` (primary universe, tighter exposure threshold), `survival_harming` (negative survival exposure), and `positive_vs_negative` (favoring versus harming with primary stratum boundaries).

The two three-guide analyses `three_directional_guides` and `three_guides_fixed_universe` are historical stored-zero diagnostics only. The canonical archived copies are in `results/legacy_protocol/`; any corresponding compatibility output at top level is not one of the seven current comparisons. In historical files `positive_support_*`, `negative_support_*`, `direction_audit_pass`, `stronger_direction_audit_pass`, `stronger_evidence_flag`, and `candidate_class` encode the withdrawn interpretation. They must not be read as current calibrated guide significance or biological validation. Legacy candidate counts 14/2/269/142 combined evidence and eligibility; the current independent axes replace that classification.

`results/identifier_sensitivity/` is a separate post-review HGNC/Ensembl analysis, with explicit source symbols, approved symbols, stable IDs, exclusions and collision handling. Its `mapping_ledger.csv`, `stable_id_primary_genes.csv`, and `stable_id_primary_strata.csv` do not modify the exact-symbol primary. See `identifier_mapping.md` for the rules and membership changes.

## Machine-readable interpretation gate

`evidence/source_semantics_audit.json` records raw-source hashes and `claim_permissions`: the V1 S8 significance claim is `REJECTED`; complete genome-wide continuous q-value reconstruction is false; the published Table S1 catalogue is `PERMITTED_WITH_QUALIFICATION`; selected S7 values are `SOURCE_SUPPORTED_SELECTED_SUBSET_ONLY`. `tests/validate_source_semantics.py` checks these permissions, binary distributions, actual embedded M bytes, linkage counterexample, subset completeness and current catalogue reconstruction. A passing test permits only the stated qualified interpretation; it is not biological validation.
