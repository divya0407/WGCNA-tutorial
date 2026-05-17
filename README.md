This project performs an integrated Weighted Gene Co-expression Network Analysis (WGCNA) using multiple colorectal cancer (CRC) microarray datasets obtained from the Gene Expression Omnibus (GEO) database. Unlike traditional differential expression analysis, WGCNA enables the identification of groups of highly co-expressed genes that may participate in common biological pathways and regulatory mechanisms.
 The workflow includes preprocessing of raw microarray data, probe annotation, dataset integration, batch correction, co-expression network construction, module–trait relationship analysis, and hub gene identification. 

The following GEO datasets were included in the analysis:

1.GSE4107

2. GSE33113

3. GSE44076

4. GSE113513

Objective:  The objective of this study is to identify biologically relevant co-expression modules and potential hub genes associated with colorectal cancer progression using integrated transcriptomic analysis.

Data Preprocessing: 

•	Raw CEL files were downloaded from GEO and independently preprocessed for each dataset.

Preprocessing Steps:
	Background correction and normalization using the Robust Multi-array Average (RMA) method.

	Probe annotation using platform-specific annotation packages.

	Removal of probes without valid gene symbols.

	Averaging of duplicate probes corresponding to the same gene.

	Export of cleaned gene expression matrices as WGCNA-ready CSV files.

WGCNA Workflow:
   
1. Loading WGCNA-Ready Expression Matrices
Preprocessed CSV files were imported into R.

3. Identification of Common Genes
Only genes shared across all datasets were retained for downstream analysis.

5. Dataset Integration
Expression matrices were merged into a single combined matrix.

7. Batch Effect Correction
Batch effects between GEO datasets were corrected using the ComBat algorithm from the sva package.

9. Variance Filtering
Low-variance genes were removed prior to network construction.
10. Soft Threshold Selection
The optimal soft-thresholding power was determined using the pickSoftThreshold() function of the WGCNA package.

12. Co-expression Network Construction
A weighted gene co-expression network was constructed using the blockwiseModules() function.


14. Module–Trait Relationship Analysis
Module eigengenes were correlated with tumor and normal phenotypes.

16. Hub Gene Identification
Hub genes were identified based on:
•	Gene Significance (GS) 
•	Module Membership (MM) 
•	Statistical significance thresholds

Required Packages:

	WGCNA

	tidyverse

	limma

	sva

	oligo

	GEOquery

	hgu133plus2.db

	hugene20sttranscriptcluster.db

	pd.primeview



