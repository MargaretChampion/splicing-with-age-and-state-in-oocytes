# State-aware alternative splicing analysis in single-cell mouse oocytes

This project analyzes alternative splicing in mouse oocytes using MARVEL. The dataset originates from work by Kelsey et al, who analyze transcriptome variation of oocytes from reproductively old vs. reproductively young mice. These authors found that, rather than age alone, NSN/SN strongly shapes transcriptomic variation. Thus, the aim of this analysis is to determine if alternative splicing patterns differ based on maternal age and/or chromatin configuration.

Rather than treating NSN/SN structure as background noise, this project uses it as an explicit interpretive constraint for downstream splicing analysis. Chromatin-state assignments were carefully reconstructed and validated before moving into splicing analysis. The main technical contribution of this project was stabilizing the MARVEL skipped-exon (SE) workflow across several incompatible intermediate formats. The major issues were structural rather than biological: malformed `tran_id` strings from `Preprocess_rMATS()`, splice-junction coordinate mismatches between STAR-derived junction tables and MARVEL expectations, and downstream use of pre-quantification rather than quantified PSI objects. After patching SE event strings, rebuilding splice-junction coordinates into MARVEL-compatible form, and correcting downstream object usage, the SE branch became stable. The resulting checkpoint contains 7,552 quantified SE events across 88 samples, together with PSI values, validated splice-feature annotations, and the included/excluded junction counts used for PSI calculation.

With SE quantification working, downstream comparisons were run both across all samples and within an NSN-restricted subset. These analyses identified shared events as well as events detected only in the NSN-restricted comparison, supporting the interpretation that cell-state heterogeneity can partially mask age-associated splicing signal in broader comparisons. Event-level visualization was then completed after fixing a dplyr scoping bug in the plotting helper, allowing candidate PSI plots to be generated cleanly across the SE candidate set.

This project currently focuses on event-level splicing quantification and differential comparison rather than a full end-to-end MARVEL-style integrative analysis. In other words, it establishes a stable lower layer of the workflow: event definition, PSI quantification, differential comparison, and event-level visualization. Higher-level MARVEL analyses, including fuller integration with gene expression, annotation layers, and broader modeling, are intentionally deferred until the underlying event classes are robust.

## Key findings

- This analysis reproduces the original observation that chromatin configuration (NSN/SN), not age alone, is a defining axis of transcriptomic structure in mouse oocytes. Expression-side analyses confirmed that chromatin state dominates global variation, while age remains a weaker but still detectable signal.
- Expression variability analyses showed that age-associated expression variability was most evident in NSN-like oocytes rather than broadly distributed across all cells. Further, age-associated expression variability did not overlap with state-associated differential variability. This supports the interpretation that NSN aging heterogeneity is not simply a restatement of chromatin-state differences.
- At the level of skipped-exon (SE) events, the current analysis identifies a small set of moderate-confidence candidate events, but does not support a strong global shift in splicing patterns with either age or chromatin state. This may reflect partial masking by cell-state heterogeneity in aggregate comparisons, or that skipped-exon variation is not a dominant axis of transcriptomic variation in these oocytes under the conditions analyzed.
- That said, skipped exons represent only one class of alternative splicing events. Analysis of additional event types (A5/A3 and intron retention) is ongoing, and may reveal patterns not captured in the SE analysis.

- ###### additional analysis ###
- I used a published NSN/SN transcriptomic signature defined in an independent dataset and projected it onto the oocyte transcriptomes as a state score prior to downstream classification analyses. The resulting signature scores cleanly recapitulated the previously reported SN/NSN classifications within the dataset, validating the reconstructed chromatin-state assignments. Consistent with both the original study and downstream expression analyses, chromatin configuration emerged as a dominant axis of transcriptomic variation, while age-associated structure remained detectable but substantially weaker.

## Current workflow status

- **SE (skipped exon):** stable and analyzed
- **A5/A3:** appear functional at the preprocessing level, but are not yet deeply characterized
- **RI (intron retention):** still under active debugging due to event-format and modeling challenges

That boundary is deliberate. The immediate goal is to stabilize RI and other event classes before extending the project into a more complete MARVEL-style analysis without carrying upstream instability into downstream interpretation.

## Repository structure

```text
scripts/
├── 01_expression_context/
├── 02_metadata_reconciliation/
├── 03_marvel_preprocessing/
├── 04_se_quantification/
├── 05_differential_splicing/
└── 06_visualization/
