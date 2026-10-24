# Revised source interpretation and bounded prior-art audit

Review date: 23 October 2026. This source interpretation and literature review was performed after the original results were known. Computational reproduction does not validate continuous guide significance; source interpretation is distinct from the numerical calculation.

## What the two source studies establish

[Benaglio et al., Cell Genomics 2022](https://doi.org/10.1016/j.xgen.2022.100214) supplies a signed cytokine-representation contrast. Normal-function “pro-death” genes correspond to knockout enrichment; “pro-survival” genes correspond to knockout depletion. These terms are not the knockout directions. Relative representation can also reflect growth/selection. The n=3 cell-count time-course is separate from the pooled genome-wide contrast. Article/deposit denominator discrepancies remain visible. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903835/), [deposited data v2](https://zenodo.org/records/7106728).

[Rottner et al., Nature Genetics 2023](https://doi.org/10.1038/s41588-022-01261-2), online December 2022, supplies two independent intracellular-insulin screens and an author-reported directional hit catalog. Positive INSlow/INShigh LFC corresponds to lower content after knockout. Content is not secretion. CALCOCO2 follow-up illustrates that normalization matters; it cannot be imputed to unrelated genes. Supplementary S5 baseMean is from CALCOCO2 silencing/control RNA-seq, so it is an expression proxy rather than an untreated reference. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9839450/).

Detailed library, timing, condition, replication and endpoint differences appear in `source_design_comparison.md` and its table JSON. A common cell-line name does not make the experiments paired or exchangeable.

## Consequential correction: source precision and author calls

Workbook-XML inspection finds only numeric zero and one in all 142,180 S8 guide FDR entries and the associated p-value fields. This is stored data, not display formatting. S7 retains continuous FDRs for its selected rows. All 2,384 S7 rows uniquely match S8 by gene, guide coordinate, counts and LFC; every linked S8 FDR equals nearest-integer rounding of S7 FDR. TOMM6 contains FDR 0.42555 in S7 but zero in S8. Embedded Power Query M converts S8 probability fields to Int64.Type, providing a documented import-coercion explanation; author intent is not known. S8 zero cannot be interpreted as FDR <0.1.

The [official author-code archive, v1.0.0](https://doi.org/10.5281/zenodo.7226348) uses a paired MAGeCK comparison and a hit-selection script involving count filters, FDR and the same guide identity across replicates. `source_semantics.md` and `evidence/source_semantics_audit.json` document the script-level details and exclusions. Full original summary/lookup outputs were absent. V2 therefore separates published Table S1 calls, selected continuous S7 evidence, and historical V1 stored-zero calculations.

The selected S7 subset supplies exactly two positive FDR<0.1 guides in each replicate, with the same identities paired, for each of the fourteen primary candidates and NUPL1. This supports those selected guide statements as source evidence; it is neither a new experiment nor reconstruction of genome-wide q-values. S7 contains two of four guides per replicate for these candidates, so omission does not establish nonsignificance of the others. It does not establish three significant guides per replicate. The selected evidence does not repair the full 142,180-row S8 universe or validate a stored-zero sensitivity. Source-byte checks and numerical reproducibility cannot recover lost precision.

## Closest identified overlap: Okino

[Okino et al., 2026 preprint](https://doi.org/10.64898/2026.03.30.715154) already integrates stress profiling/co-expression with Rottner insulin-content regulators and adds a thapsigargin survival screen in EndoC-betaH1. It expressly discusses secretion-related functions that may reduce ER-stress survival. This precludes claims to the first stress/function connection or first integration of these biological themes. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13060215/), [PubMed record](https://pubmed.ncbi.nlm.nih.gov/41959148/).

The reviewed text does not provide this run's exact coverage/expression-controlled comparison of Benaglio's cytokine-survival shortlist with the published directional content catalog and its exhaustive 427-candidate accounting. That is a distinction in estimand and deliverable, subject to the damaged-guide-field limitation, not a priority certification. The record reviewed is labeled a preprint; this bounded search did not establish a later publication. `source_design_comparison.md` states the concrete assay and integration differences and the missing joint-condition experiment.

## Other adjacent primary work

[Howe et al., Nature Communications 2026, volume 17, article 513](https://www.nature.com/articles/s41467-025-67199-1), first published online 13 December 2025, examines transport between cellular perturbation and human phenotypes using loss-of-function information and explicit assumptions. It supports keeping cellular content distinct from validated clinical function; it is not an analysis of these two catalogs. No GPAT computation or human transport validation is performed here.

[van Tienhoven et al., Cell 2025](https://doi.org/10.1016/j.cell.2025.02.018) associates stress-dependent insulin-mRNA decay with islet vitality/function. Only its [primary abstract](https://pubmed.ncbi.nlm.nih.gov/40112799/) was inspected for this context; it is not an input or candidate validation. A lower insulin-related measurement cannot be interpreted automatically as worse endocrine performance.

## Candidate evidence adds utility, not validation by analogy

`candidate_annotations.json`, `.csv` and `candidate_biology.md` distinguish independent candidate-targeted beta-cell studies, beta-cell observational findings, other-system mechanisms and predicted annotations. Functional beta-cell/islet studies were located for SSR2, NCOA6, CREB3, SLC8B1 and UBE2C; their species, stressors, perturbations and normalizations limit transfer. RAD52 evidence is observational, PEX13 beta-cell work is separately flagged thesis evidence, and other candidates mostly have process-level context. “Not located” is not an absence claim. NUPL1/NUP58 is a fifteenth, excluded illustration and does not silently alter the exact-symbol primary analysis. Speculative alternatives and untested endpoints are explicitly marked.

## Bounded search and retrieval record

Initial query families included Benaglio/Rottner with `trade-off`, `overlap` and `enrichment`; cytokine/insulin-content CRISPR comparisons; and the exact titles of the source and Okino studies. Revision queries combined each candidate with `beta cell`, `islet`, `insulin`, `knockout` or `silencing`, plus process/alias terms including NCLX, GLUT11 and NUP58. Original journal, PubMed/PMC, authoritative gene and institutional repository records were followed. Reviews and commercial summaries served as leads only. No equivalent exact analysis was found, but search-engine coverage and terminology remain incomplete.

Some claims rely on indexed original full text or primary abstracts where original XML was unavailable. Institutional gray literature is labeled. This is not a systematic review, has no completeness guarantee and has not received external human specialist assessment. A later equivalent analysis could change the novelty assessment; corrected full-precision data could change the guide-level interpretation.
