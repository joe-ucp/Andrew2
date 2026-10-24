# Source semantics audit and revised interpretation

The published Table S1 catalogue remains usable as the authors' reported directional decisions. Table S8 cannot reconstruct the authors' guide-level FDR threshold: its four p/FDR columns contain only stored integers 0 and 1. The independent V1 verification established computational agreement, not the scientific validity of treating those stored zeros as FDR <0.10. This semantic audit supersedes that broader interpretation. V2 proceeds as a qualified comparison of published catalogues, with guide significance reconstructed only where continuous selected-guide values are available in Table S7.

Run `python code/audit_source_semantics.py` from the project root. It reads the pinned raw inputs directly and does not import the main analysis construction. The audit checks the workbook SHA-256, inspects actual XML values and types, links S7 to S8, reconstructs the revised universes and all 427 candidate axes, and records script/input/output hashes in `evidence/source_semantics_audit.json`.

## The saved import pipeline identifies whole-number conversion

The workbook contains an embedded Power Query package in `customXml/item1.xml`. Decoding its DataMashup container recovers `Formulas/Section1.m`. The saved `RRApaired_screen sgrna_summary (1)` query explicitly applies `Table.TransformColumnTypes` with `Int64.Type` to **p.low, p.high, p.twosided, and FDR**. LFC, counts, score and other continuous columns retain `type number`. The `sgRNA_hits` query applies `Int64.Type` to p.low only; p.high, p.twosided and FDR retain `type number`. This exactly matches the differing S7/S8 column patterns. `evidence/workbook_query_transformations.json` contains the workbook hash, raw embedded M hash, column declarations and relevant query text with original source paths redacted. The queries are inspected, never executed.

Microsoft documents that column-type conversion normally uses the corresponding `.From` function and that `Int64.From` rounds fractional values to whole numbers, defaulting to ties-to-even. Thus the saved query supplies a concrete conversion stage consistent with the data, beyond a generic rounding hypothesis. It does not identify why the types were chosen, whether a default was accepted, or anyone's intent. The original complete continuous input is still unavailable. [Column conversion](https://learn.microsoft.com/en-us/powerquery-m/table-transformcolumntypes), [integer conversion](https://learn.microsoft.com/en-us/powerquery-m/int64-from).

Producer metadata identifies Microsoft Excel 16.0300. `calcPr` contains `calcId="0"`; it has no explicit `fullPrecision` attribute. Probability cells are stored numeric XML values with default cell style, while high_in_treatment is a separately typed Boolean. This is neither merely Excel display formatting nor floating-point underflow. The saved integer import transformation is the evidenced mechanism; original continuous values cannot be recovered by changing cell formatting.

| Table S8 column | Stored 0 | Stored 1 | Total |
|---|---:|---:|---:|
| p.low | 72,618 | 69,562 | 142,180 |
| p.high | 69,562 | 72,618 | 142,180 |
| p.twosided | 110,998 | 31,182 | 142,180 |
| FDR | 99,564 | 42,616 | 142,180 |

S7 has 2,384 rows. Its p.low is also binary (245 zeros, 2,139 ones), while p.high, p.twosided and FDR contain 2,340, 2,376 and 2,361 distinct stored values respectively. Full distributions and XML examples are in the JSON audit.

## Continuous selected-guide evidence and its limits

All 2,384 S7 rows uniquely match S8 by source gene, genomic coordinate, control count, treatment count and LFC. In every pair the stored S8 FDR equals nearest-integer rounding of the continuous S7 FDR. The decisive counterexample is TOMM6: S7 row 2385 matches S8 row 91896, r0, coordinate chr6:41755542-41755561. The S7 FDR is 0.42555 while S8 stores 0; p.high is 0.86265 versus 1, and p.twosided is 0.2747 versus 0. A stored zero therefore does not establish FDR <0.10. Exact linked rows are in `evidence/S7_S8_semantic_linkage.csv`.

All 14 exact-symbol primary candidates and excluded NUPL1 have **two positive-LFC guides with continuous S7 FDR <0.10 in each replicate, and the same two guide identities across replicates**. Each has only two of its four total guides per replicate represented in S7. No candidate has three such source-supported guides per replicate in this selected table; the six retained by the old three-stored-zero rule do not acquire three-guide support from S7. There are no S7 FDR values exactly equal to 0.10, so `<` versus `<=` does not change these counts. TOMM6 has one positive guide meeting the threshold in r0, two in r1, and one matching identity across replicates.

These are continuous values as stored in a **selected published hit subset**, not complete original MAGeCK probabilities, an independent replication, or a new genome-wide hit analysis. Missing guides are not demonstrated nonsignificant. S7 cannot repair the complete S8 denominator, promote nonpublished genes, or certify absence of opposing effects. Per-gene positive/negative counts, same-guide identities, source-row completeness and strictness checks are in `evidence/candidate_S7_selected_guide_evidence.csv`.

## Encoded-value diagnostics must remain illustrative

In the 12,184-gene original primary universe, 2,421 genes absent from the published directional catalogue satisfy the illustrative rule “positive LFC and stored FDR=0 for at least two guides in each replicate.” That rule therefore cannot be substituted for the published 518 lower-content labels. The complete gene list is in the audit JSON and `evidence/semantic_gene_diagnostics.csv`.

| Group | Genes | Any negative-LFC stored-zero guide | At least one in each replicate | At least two in each replicate |
|---|---:|---:|---:|---:|
| Primary universe | 12,184 | 11,541 | 8,230 | 2,743 |
| Published primary decreases | 479 | 385 | 176 | 15 |
| Fourteen primary candidates | 14 | 10 | 7 | 0 |

These counts describe encoded values and opposing LFC signs. They are not counts of significant opposing guides. Paired identities are also reported: all 14 have at least two identical positive stored-zero guides, but pairing alone does not validate an encoded probability. The selected continuous S7 evidence provides the narrower source-supported assessment above.

## Original code and bounded source recovery

The article links the author code release [Zenodo 7226348](https://doi.org/10.5281/zenodo.7226348), version 1.0.0. Its archived R script uses MAGeCK 0.5.9.2 output, retains **each** low/high count >=10, applies FDR <=0.10, requires same-coordinate support across replicates and at least two such guides per gene in the same direction, removes NDUFS6, and applies a separate expression lookup. Its MAGeCK command compares LOW1,LOW2 against HIGH1,HIGH2 with `--paired`. These details differ from the illustrative legacy count-sum and per-replicate stored-zero rule. The expression lookup is not the S5 baseMean filter used for this cross-screen catalogue.

The 10,474-byte author archive includes scripts but not the original `RRApaired_screen.sgrna_summary (1).txt`, gene summary, count input, separate expression lookup, or output hit table. No rounding or spreadsheet export command was found in that archived analysis code; the workbook itself contains the conversion evidence. The article identifies raw screen reads at ENA PRJEB44712. A complete FASTQ reanalysis and recovery of absent processed files are outside this bounded audit; neither has been claimed. [Rottner et al.](https://doi.org/10.1038/s41588-022-01261-2).

## Independent revised-catalogue checks

Direct raw reconstruction matches V2 membership and A/L assignments exactly for the primary catalogue (N=12,184; A=285; L=479; overlap=14) and coverage-only catalogue (N=17,276; A=390; L=514; overlap=15). TOMM6's published label is retained in coverage-only membership without the failed stored-value direction filter. NUPL1 remains a published lower-content candidate excluded by exact-symbol expression matching.

All 427 candidates match the revised evidence and eligibility axes: 15 published lower-content, two published higher-content, 384 with no published directional call, and 26 with unavailable screen measurements. There are 390 with common coverage and 285 primary eligible. Eleven measured candidates lacking coverage remain “no published directional call,” not “measurement unavailable.” This comparison validates the revised table construction; it does not validate original statistical calibration or a causal survival–content trade-off.

No source value, published decision, or original version was silently corrected. The review disposition is qualified catalogue analysis with guide-significance limits disclosed. See `docs/identifier_mapping.md` for the separate post-review stable-ID sensitivity.
