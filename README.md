# Beta-cell survival and insulin content

[Paper (PDF)](manuscript/paper.pdf) · [427-candidate table](results/candidate_table_427.csv) · [Source correction](docs/source_semantics.md) · [Reproduce](#reproduce) · [Citation and license](#citation-and-license)

**Can beta-cell vulnerability be reduced while preserving insulin production and glucose-responsive release?** This study compares published cytokine-survival and intracellular insulin-content knockout screens in human EndoC-betaH1 cells. It identifies candidate perturbations and quantifies their overlap while keeping the two experimental endpoints distinct.

## Three contributions

1. **A documented loss of probability precision.** All 142,180 Table S8 guide FDR entries are stored as 0 or 1. Each of 2,384 uniquely matched Table S7 rows has an S8 FDR equal to nearest-integer rounding of its continuous S7 value; the workbook's saved Power Query explicitly converts the probability columns to `Int64.Type`. Stored zero therefore cannot establish FDR <0.10. The [source analysis](docs/source_semantics.md) traces the conversion and withdraws the original guide-significance interpretation.
2. **A usable, fully accounted candidate resource.** The [427-candidate table](results/candidate_table_427.csv) separates published directional evidence, measured coverage and primary eligibility. Excluded candidates remain visible, and [paired guide data](results/candidate_paired_guide_effects.csv) expose opposing effects. [Biological annotations](docs/candidate_annotations.csv) add bounded, source-linked context for the 14 primary lower-content overlaps and NUPL1.
3. **An unresolved enrichment estimate with explicit uncertainty.** The published-catalog comparison gives a conditional odds ratio of **1.23 (exact 95% CI 0.66–2.13)**. Removing the invalid stored-FDR filter leaves primary membership and estimates unchanged. The result neither establishes enrichment nor demonstrates absence of a survival–function trade-off.

| Primary comparison | Published lower-content call | No published lower-content call |
|---|---:|---:|
| Survival-favoring knockouts | 14 | 271 |
| Remaining eligible genes | 465 | 11,434 |

The corresponding fractions are 14/285 and 465/11,899; the one-sided enrichment p-value is 0.2634. Inference conditions on the declared expression/guide-coverage strata and a common-odds, independent-stratum working model.

## Explore the results

- **Report:** [PDF](manuscript/paper.pdf), [Word](manuscript/paper.docx) and [TeX](manuscript/paper.tex), with source-design comparisons, candidate-specific hypotheses and proposed causal/recovery experiments.
- **Candidate accounting:** [complete table](results/candidate_table_427.csv), [field definitions](docs/data_dictionary.md), [paired guide effects](results/candidate_paired_guide_effects.csv) and [candidate biology](docs/candidate_biology.md). Across all 427 candidates, 15 have published lower-content calls, 2 higher-content calls, 384 measurements without a published directional call, and 26 no exact insulin-screen measurement. Common coverage includes 390 genes; 285 meet primary eligibility.
- **Source evidence:** [probability-field analysis](docs/source_semantics.md), [saved query transformations](evidence/workbook_query_transformations.json), [S7–S8 linkage](evidence/S7_S8_semantic_linkage.csv) and [selected continuous-guide evidence](evidence/candidate_S7_selected_guide_evidence.csv). S7 supports two paired positive FDR <0.10 guides per replicate for the 14 primary candidates and NUPL1; only two of four guides per replicate are represented, so omitted probabilities remain unknown.
- **Interpretation and follow-up:** [assay comparison and proposed experiments](docs/source_design_comparison.md), [bounded literature review](docs/literature_audit.md), [identifier sensitivity](docs/identifier_mapping.md) and [precision/influence diagnostics](docs/exploratory_diagnostics.md).

## Reproduce

Clone [joe-ucp/andrew2](https://github.com/joe-ucp/andrew2) and work from the repository root:

```sh
git clone https://github.com/joe-ucp/andrew2.git
cd andrew2
```

Use **Python 3.12** with the exact packages in [requirements.txt](requirements.txt); [environment.yml](environment.yml) supplies a conda/pip specification. Reuse a suitable environment or install the requirements in an environment outside the repository.

```sh
python -m pip install -r requirements.txt
python reproduce.py
```

The first run requires network access to download the pinned Rottner author-code archive from its original Zenodo source. The archive is **not bundled**. Downloads are checked against the expected SHA-256; a mismatched source is rejected. The default external cache is `~/.cache/beta-cell-screen-comparison/sources`. Set `BETA_CELL_SOURCE_CACHE` or pass `--source-cache PATH` to choose another external directory; the command-line option takes precedence. Cache directories inside the repository are rejected.

To fetch and verify inputs without running the analysis, use `python code/fetch_inputs.py`. Once the required inputs are present, check the cache and reproduce without network access:

```sh
python code/fetch_inputs.py --offline
python reproduce.py --offline
```

Use the same cache setting on subsequent runs. Missing inputs cause the offline check to fail. Source URLs, hashes and terms are documented in [input provenance](evidence/input_manifest.json) and [source notices](docs/data_sources_and_licenses.md).

The driver verifies inputs, reconstructs original and amended calculations, runs independent statistical and source-semantic checks, and generates identifier and continuous/influence diagnostics. Seven amended catalog comparisons and the separate identifier sensitivity reuse source data; they are not independent experimental replications. The two original three-guide calculations remain explicitly labeled historical stored-zero diagnostics in [the historical specification](evidence/protocol.md).

To rebuild synchronized Word, TeX and PDF, also provide `pdflatex` with the packages declared in [paper.tex](manuscript/paper.tex):

```sh
python reproduce.py --offline --build-manuscript
```

[narrative.json](manuscript/narrative.json) is the editable prose source. The manuscript builder combines it with generated results and records dependency checks in [manuscript_consistency.json](evidence/manuscript_consistency.json). [release_manifest.json](evidence/release_manifest.json) indexes the current scientific files and their checksums. Keep environments, installers and temporary renders outside the repository.

## Scope and limitations

Relative representation after cytokine treatment can reflect survival, growth or selection. Intracellular insulin content does not establish insulin synthesis, glucose-responsive secretion or durable functional preservation. A gene appearing in both catalogs is therefore a candidate for matched experiments, not a validated protective or therapeutically useful target. Higher content and absence of a published content call likewise do not establish preserved function.

The screens examine separate perturbation populations and conditions. Conditional intervals do not capture all gene dependence, source-call uncertainty or donor variation. The corrected interpretation and additional diagnostics were developed after the original results were known; they are not a new preregistration. Proposed mediator interventions and injury/recovery studies have not been performed. This analysis establishes no prevention, remission, insulin independence or cure.

## Citation and license

Use [CITATION.cff](CITATION.cff) to cite this report, and cite the original studies and datasets listed in [source provenance and licenses](docs/data_sources_and_licenses.md).

Original analysis code is available under the [MIT License](LICENSE), © 2026 UCP Technology LLC. Third-party data and reference code retain their own terms; the project license does not relicense them.
