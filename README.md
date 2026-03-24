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
5. Prepare MARVEL input

Raw data:
- star/ (aligned files)

Outputs:
- data/derived/

----
Pipeline notes:
1.1 SRR12212403 was identified as empty during samtools quickcheck and dropped from analysis
