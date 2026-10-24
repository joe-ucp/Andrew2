# Source design, comparability and a discriminating follow-up

Reviewed 23 October 2026. Facts below come from the original articles and, for the Rottner call pipeline, the released author code. These designs share a cell-line name and gene identifiers; they do not measure the same edited cells under the same condition.

| Design element | Benaglio cytokine-representation screen | Rottner insulin-content screen |
|---|---|---|
| Model | Human EndoC-betaH1 | Human EndoC-betaH1 |
| Library | GeCKO-V2 | TKOv3; four guides per targeted gene |
| Infection | MOI 0.3 | MOI 0.3; independent transductions |
| Selection and time | 48 h after infection, 3 d puromycin selection; another 14 d culture before 72 h treatment | 7 d puromycin incubation before collection/staining/sorting |
| Condition | IL-1beta 0.5 ng/ml, TNF-alpha 1 ng/ml, IFN-gamma 10 ng/ml; BSA control | Standard culture; intracellular insulin immunostaining, without the cytokine challenge |
| Endpoint | Guide representation after cytokines relative to control | Guide representation in low- versus high-insulin sorted fractions |
| Sampling | 60 million control and 30 million cytokine cells harvested | Live-single-cell gating before fixation/permeabilization; low gate set using isotype control, high gate highest 10% |
| Genome-wide replication | Methods describe the pooled contrast; do not infer multiple independent screens from the n=3 time-course | Two independent screens, with separate transduction, sorting and sequencing |
| Positive effect | Knockout relatively enriched after cytokines | Knockout enriched in low-insulin fraction, consistent with reduced content |

Sources: [Benaglio original article, STAR Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903835/); [Rottner original article, screen methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC9839450/). The Benaglio n=3 refers to a residual-cell-count time-course. Rottner follow-up assays with three or more experiments are likewise separate from its two genome-wide replicates.

The endpoint mismatch is substantive. Relative representation can combine death, growth, cell-cycle effects and prechallenge selection; it does not by itself identify absolute viability. Insulin fluorescence is influenced by biosynthesis, processing, storage, cell state and release. The different libraries, maturation intervals and selection histories also permit perturbation-strength differences. A signed overlap therefore means two separate reported phenotypes, not demonstration that a protected cell has lost glucose-responsive secretion.

The [Rottner author-code release](https://doi.org/10.5281/zenodo.7226348) further distinguishes a published call from V1's reconstructed rule. The released pipeline uses MAGeCK v0.5.9.2 with a paired LOW1,LOW2 versus HIGH1,HIGH2 comparison. Its R hit-selection script applies bin-count and FDR filters and requires the same guide identity across replicates, with at least two such guides per gene. The archived script uses `<= 0.1`, while the article states `< 0.1`; `evidence/source_semantics_audit.json` identifies the script and its hash. All 2,384 uniquely linked S7/S8 FDRs show nearest-integer precision loss; TOMM6 includes 0.42555 rounded to zero. Embedded Power Query M explicitly converts S8 probability columns to Int64.Type, explaining whole-number storage without establishing author intent. Thus S8 zero cannot establish FDR below 0.1. Selected S7 rows support two paired positive guides per replicate for the fourteen candidates and NUPL1, but do not restore full-universe q-values or validate a third-guide claim. Published Table S1 labels remain author-reported calls. This comparison makes no new claim about the authors' statistical error rate.

## Concrete relation to Okino

[Okino et al., 2026 preprint](https://doi.org/10.64898/2026.03.30.715154) already links stress-response networks with Rottner insulin-content regulators. It also performs an EndoC-betaH1 Brunello screen after 72 h thapsigargin exposure, reporting 167 normal-function pro-survival and 47 pro-death genes among 19,114 tested. Its discussion expressly considers secretion-related activities that may worsen ER-stress survival. Thus neither stress/function tension nor integrating insulin-content evidence with stress genomics is new here.

V2 instead interrogates the specified cytokine-survival shortlist against published directional insulin-content calls, retaining coverage, expression and all-candidate accounting. Okino uses a chemical ER-stress contrast and network integration; this run uses the cytokine contrast and a declared finite-catalog enrichment estimand. Neither design directly measures survival, content and secretion after the same perturbation under matched inflammatory conditions. This is a narrow distinction of question and accounting, subject to incomplete recovery of full-precision guide probabilities, not a claim of priority or superior biological validity. See [full preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13060215/).

## Matched-condition discrimination design — proposed, not performed

Use a factorial design within the same experimental batches: targeted perturbation versus non-targeting control, each under basal conditions and the matched inflammatory challenge. Start with independently targeted guides and confirm the intended gene/protein perturbation; include a guide-resistant rescue where feasible. Use the same perturbation timing, starting viable-cell number and sampling schedule across arms. An initial EndoC-betaH1 study can resolve cell-line mechanisms, followed by multiple human donor preparations or independently differentiated stem-cell beta-cell batches for transportability. Wells and cells are technical subsamples; donor or independent preparation is the biological unit.

Measure absolute live-cell recovery and death markers before and after challenge, with parallel proliferation/cell-cycle readouts. In matched aliquots at the same time point, quantify insulin and proinsulin per viable beta cell and per total culture. Measure basal and glucose-stimulated insulin release, retaining absolute amounts, amounts per viable cell, fractions of content released and stimulation ratios. A time-resolved glucose-response assay is especially useful for SLC8B1. Basal and stimulated values must remain visible because a similar fold change can coexist with lower absolute secretion, and content normalization can conceal depleted stores. Measure beta-cell identity and total protein/cell size where they are plausible alternatives.

Predeclare the endpoint hierarchy, donor/batch model and multiplicity treatment. For the initial phenotyping stage, estimate perturbation-by-inflammation interactions and show the four group means with intervals. A preservation claim needs a justified equivalence/noninferiority margin and enough precision to evaluate it; a nonsignificant difference is insufficient. Determine sample size from pilot variability for the relevant biological unit, rather than borrowing the large gene count of this catalog analysis.

| Competing explanation | Pattern that would distinguish it |
|---|---|
| Selective inflammatory protection | Improved absolute live-cell recovery/death markers under inflammation, beyond any basal growth effect |
| Reduced biosynthetic demand | Lower per-cell content/proinsulin production with altered stress response; absolute secretory capacity may also fall |
| Altered release rather than production | Content and released insulin move differently; secretion time course and fraction released discriminate |
| Broad fitness or identity change | Basal growth, cell cycle, cell size or beta-cell identity shifts accompany the representation/content signal |
| Library/guide artifact | Phenotype fails independent targeting or gene rescue |

After joint phenotyping, independently vary a measured candidate mediator while retaining the original perturbation, with mediator-only and direct-injury controls. Attenuation of protection supports contribution; preventing the mediator change and removing protection addresses necessity, while mediator suppression alone tests sufficiency. Functional restoration with meaningful protection retained supports separability in the tested setting; unsuccessful restoration does not establish an unavoidable trade-off. The main report gives illustrative candidate-to-decision examples, including higher-content comparison candidates RAB14 and LRRC61.

A recovery stage must distinguish withdrawal of the injury from withdrawal of the intervention. Permanent knockout followed by cytokine removal cannot establish reversibility of the perturbation; demonstrate that both intervention and affected program can recover before testing temporary suppression. Track survivors versus replacement and retain basal, stimulated, temporal, per-cell and total-population output against prespecified preservation margins.

Candidate-specific readouts in `candidate_annotations.json` refine this common design. They are hypotheses for discrimination, not evidence that any listed mechanism explains the current overlap. No laboratory work or clinical testing has been performed in this run.
