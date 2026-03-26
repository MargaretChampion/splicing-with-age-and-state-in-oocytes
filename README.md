Mouse oocyte single-cell QC + splicing project

Accession PRJNA645959

STAR alignment → GRCm38
featureCounts → GRCm38.102

---
Meta Notes:
Data is single ended
Alignment is against 

Pipeline:
1. Build QC table from STAR outputs
2. Filter low-quality cells
3. PCA + QC visualization
4. Condition-level QC comparison
5. Differential Expression
6. Variability analysis
7. Prepare MARVEL input

Raw data:
- star/ (aligned files)

Outputs:
- data/derived/

----
Pipeline notes:
1.1 SRR12212403 was identified as empty during samtools quickcheck and dropped from analysis
3.1 QC inspection showed that old oocytes had fewer detected genes than young oocytes. This difference was significant by Wilcoxon rank-sum test and remained significant after adjustment for log10(total_counts) in a linear model, suggesting that the age-group shift is not explained solely by library size.
   this is noted in script 05
4.1 old samples show somewhat greater dispersion / heterogeneity, but there is no single unequivocal outlier warranting exclusion based on distance ranking, PCA, and QC metrics.
4.1 sample-level PCA and pairwise distance analysis did not reveal a single additional unequivocal outlier. Samples with the highest mean distances corresponded to expected extremes in library size / complexity rather than obvious technical failures.
   while the source paper noted removing one sample, it is possible that their outlier == the one that I identified as empty and dropped before PCA
4.2 Old and young samples differ in raw junction-support metrics, but these differences are largely explained by overall junction-support depth; after depth adjustment, the residual age-group effect is modest.

5. DESeq2 identified 627 genes differing between young and aged oocytes at FDR < 0.05, of which 103 showed an absolute fold change of at least 1.5. This is broadly consistent with the published analysis, which reported 560 DEGs and 156 genes with fold change ≥ 1.5. In our analysis, large-effect changes were skewed toward higher abundance in young oocytes.
5.1 however DE was performed raw/without accounting for SN/NSN. We will return and account for this
6.1 before moving to variability analysis, we incorporated SN/NSN classifications


Adding SN/NSN chromatin state annotation

A key source of variation in this dataset is oocyte chromatin configuration, classified as SN (surrounded nucleolus) or NSN (non-surrounded nucleolus). This state reflects developmental competence and was reported in the original study to strongly associate with the primary axis of variation (PC1).

To account for this, we integrated the authors’ supplemental chromatin state annotations into our sample metadata.

Because our pipeline operates on SRA run IDs (SRR), while the paper reports chromatin state using author-defined cell IDs (e.g., yGV_34), we constructed a reconciliation crosswalk:

Supplemental table: Cell → Predicted_configuration (SN / NSN / NA)
GEO/SRA metadata: GSM → SRR
Manual reconciliation file: GSM → author cell ID (oGV_XX / yGV_XX)

These were merged to map:

SRR → GSM → Cell → Predicted_configuration

This allows chromatin state to be joined directly onto our analysis-ready metadata.

Final counts:

SN: 66 samples
NSN: 20 samples
Unclassified (NA): 1 sample (oGV_48A, also NA in original paper)

We retain the original labels without imputation. The unclassified sample is kept in the dataset and excluded only where chromatin state is explicitly modeled.

This annotation is critical for downstream analyses (e.g., variability, splicing), as chromatin configuration represents a major biological axis that could otherwise confound age-related effects.

PCA with SN/NSN chromatin-state annotation

To better interpret transcriptomic structure in the oocyte dataset, we integrated the authors’ supplemental predicted chromatin configuration labels (SN, NSN, or NA) into the analysis metadata. This was necessary because the original study identified chromatin state as a major source of heterogeneity among GV oocytes, independent of age.

Because the published supplemental table uses author-defined cell IDs (for example, yGV_34), while our processing pipeline uses SRA run IDs (SRR), we built a reconciliation crosswalk linking:

SRR -> GSM -> author cell ID -> predicted_configuration

This allowed the paper-derived Predicted_configuration field to be added directly to the analysis-ready metadata.

Final chromatin-state counts:

SN: 66
NSN: 20
Unclassified (NA): 1

The unclassified sample was retained as missing rather than imputed, since it was also unclassified in the original paper.

We then recomputed a variance-stabilized expression matrix and performed PCA with chromatin state included in the metadata. In this reprocessed dataset:

PC1 is dominated by chromatin configuration (SN vs NSN-like state)
PC2 captures a weaker but independent age-associated signal

Linear modeling of PCA coordinates showed that predicted_configuration explains most of the variation in PC1, whereas age_group contributes little to PC1 once chromatin state is included. In contrast, age_group is associated with PC2 even after accounting for chromatin state.

This indicates that the major transcriptomic axis in the full dataset reflects developmental/chromatin-state heterogeneity, while aging acts as a secondary axis layered on top of that structure.

This matters for downstream analysis because age-associated expression or variability differences may otherwise reflect a mixture of:

true aging effects
shifts in SN/NSN composition
both

For that reason, chromatin configuration should be considered in downstream modeling, especially for:

differential expression
expression variability analysis
alternative splicing analysis


