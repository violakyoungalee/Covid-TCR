# Machine learning identifies T cell receptor repertoire signatures associated with COVID-19 severity

Jonathan J. Park, Kyoung A V. Lee, Stanley Z. Lam, Kat Moon, Zhenhao Fang, Sidi Chen

Code and analysis for the manuscript "Machine learning identifies T cell receptor repertoire signatures associated with COVID-19 severity"

## Project Structure
```
.
├── Standardized Data
├── Data preprocessing scripts
├── Immune repertoire statistics scripts
├── Kmer analysis scripts
├── Motif analysis scripts: GLIPH2 and OLGA for specificity analysis
├── Single-cell transcriptome analysis scripts
└── Machine Learning Models

```
## Data acquisition

TCR repertoire data was obtained from datasets published by Adaptive Biotechnologies 49,  ISB-Swedish COVID-19 Biobanking Unit 25, Fifth Medical Center of PLA General Hospital 20, and Wuhan Hankou Hospital China 19. For COVID-19 patients sequenced with Adaptive Biotechnologies immunoSEQ assays, TCR-seq data were obtained from the ImmuneCODE database at https://doi.org/10.21417/ADPT2020COVID; for healthy donor patients, TCR-seq data was obtained at https://doi.org/10.21417/ADPT2020V4CD. Single-cell TCR-seq and gene expression (GEX) data for CD4+ and CD8+ T cell repertoires from COVID-19 patients and healthy donors from the ISB-Swedish COVID-19 Biobanking Unit 25 was obtained from the ArrayExpress database 50 (http://www.ebi.ac.uk/arrayexpress) using the accession number E-MTAB-9357. Single-cell TCR-seq data from COVID-19 patients and healthy donors were also obtained from the Fifth Medical Center of PLA General Hospital, accessed through the supplementary tables of the associated publication 20; and Wuhan Hankou Hospital China, metadata accessed through the supplementary tables of the associated publication 19 and TCR-seq data obtained from the iReceptor platform 51 (http://ireceptor.irmacs.sfu.ca). The Adaptive Biotechnologies dataset comprised bulk TCR-seq data, while the ISB-S, PLA, and WHH datasets comprised single-cell TCR-seq data, with variations in sequencing modalities, patient populations, and sample sizes. Due to the differences in wet lab protocols and the potential presence of batch effects in each of these four datasets, all downstream analyses were performed separately on each individual dataset, and the result from each dataset as well as the consensus findings are reported.

## Data preprocessing

All four TCR-seq datasets were individually but identically pre-processed for standardized analysis with Immunarch v0.6.6 without pooling 52. Data obtained from the Adaptive Biotechnologies ImmuneCODE database were used directly as inputs for Immunarch processing, with 1,475 COVID-19 patient samples and 88 healthy donor patient samples (1,563 samples total) successfully loaded and used for further analysis. For the ISB-Swedish cohort, patients were first filtered by those who had sequencing data available as performed by 10X Genomics. Sequence filtering and processing were performed as follows: for cells with multiple TRA and TRB CDR3 sequences, the first instances, respectively, were selected; only cells with paired TRA and TRB sequences were kept (column chain_pairing = Single pair, Extra alpha, Extra beta or Two chains); sequence files were converted to VDJtools format for input into Immunarch. COVID severity scores were translated from the WHO Ordinal Scale (0-7) to four tiers: healthy donor (0), mild (1-2), moderate (3-4), and severe (5-7). After pre-processing, the CD4 and CD8 datasets were composed of 136,429 and 69,687 clones, represented in a total of 16 healthy donors, 61 mild, 42 moderate, and 24 severe patients, 143 individuals total (16 healthy donors, 108 mild, 93 moderate, and 49 severe repertoires when accounting for patients with samples from two time points, 266 samples total). For the PLA General Hospital and Wuhan Hankou Hospital China cohort, cells with more than one TRA or TRB sequence had the chain with the highest number of reads kept for further analysis, and sequence files were converted to VDJTools format for input into Immunarch. The PLA General Hospital aggregated patient dataset contained 31951 clones across 3 healthy donors (two healthy donors from the original study were excluded for lack of TCR CDR3 amino acid data), 7 moderate, 4 severe, and 6 convalescent patients (of which 4 were the second time point collections of moderate patients – P01, P02, P03, and P04). The Wuhan Hankou Hospital China aggregated patient dataset contained 42001 clones across 5 healthy donors, 5 moderate, and 5 severe patients. Metadata was manually reformatted from supplementary tables.

## Immune repertoire statistics

Clonotype statistics and diversity metrics were calculated using Immunarch v0.6.6 52. For the total number of unique clonotypes, the repExplore function was used with parameter .method = “volume”; for distribution of CDR3 sequence lengths, repExplore function with .method = “len” and .col = “aa”; for the Chao1 estimator, repDiversity function with .method = “chao1”; for Gini-Simpson index, repDiversity function with .method = “gini.simp”; for Inverse Simpson index, repDiversity function with .method = “inv.simp”. Clonal proportion estimates were calculated with the repClonality function with .method = “top”. CDR3, V gene, and J gene usage proportions were calculated and aggregated directly from sample TCR data. Statistical significance testing comparing groups was performed using the two-sided Wilcoxon rank-sum test by the wilcox.test in R.

## K-mer analyses

For K-mer abundance calculations, each VDJtools formatted sample was converted to a vector of CDR3 sequences. The vector was converted to k-mer statistics using the getKmers function from Immunarch, then merged with k-mer statistics of other samples using the R function merge with parameter all = TRUE for full outer join. Empty cells were converted from NAs to 0 counts. The 50,000 top variance unique k-mers were selected for downstream analyses (PCA and machine learning pipelines) with the exception of 3-mers which had 6916 unique k-mers. The selection of 50,000 top variance k-mers was done to keep the data dimensions consistent, increase the efficiency of the analysis pipelines, and use data features that are the most likely to be the most biologically meaningful. Low-variance kmers, either due to lack of representation or due to conservation between healthy and disease samples are unlikely to be significant in the T cell response to COVID given that the signature is not shared or enriched by disease status. K-mer counts were normalized to sum to 1 for each sample prior to downstream analyses. PCA was performed using the prcomp function in R with parameter center = TRUE. 

## Motif analyses

TCR clustering and specificity group analysis was performed using GLIPH2 29. Software executable for analysis was obtained from http://50.255.35.37:8080/ and run with the human v2.0 reference on clonal data for each disease condition and T cell type. Parameters include global_convergence_cutoff=1, local_min_OVE=10, kmer_min_depth=3, simulation_depth=1000, p_depth=1000, ignored_end_length=3, cdr3_length_cutoff=8, motif_distance_cutoff=3, all_aa_interchangeable=1, kmer_sizes=2,3,4, and local_min_pvalue=0.001000. 
Generation probability calculations were performed using OLGA 30. Software installation and setup were performed as described in https://github.com/statbiophys/OLGA and run on clonal data for each disease condition and T cell type. Representative calculations with parameters are as follows: olga-compute_pgen -i input.tsv --humanTRB -o out_pgens.tsv --v_in 1 --j_in 2.

## Single-cell transcriptome analysis

Single-cell transcriptome data from the ISB-S dataset were processed using Seurat v4.0.4. The pipeline included log normalization with a scale factor of 1,000,000, scaling and centering, PCA, nearest-neighbor graph construction, clustering with the Louvain algorithm, UMAP, differential gene expression, and generation of various visualizations. Parameters included: for the FindNeighbors function, dims = 1:10; for FindClusters, resolution = 0.6; for RunUMAP, dims = 1:10; for FindAllMarkers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25. Differential gene expression between clonally expanded clusters and all other cells was performed using a downsampled cell subset (5,000 cells per group) of the data and the FindMarkers function with parameters logfc.threshold = 0.01 and min.pct = 0.1. P-value adjustment was performed using Bonferroni correction. Upregulated or downregulated genes with significance q-value < 1e-4 were then used for functional annotation with DAVID analysis. In addition to default Seurat outputs, custom R scripts were used to generate visualizations including UMAPs associated with CDR3 motifs and disease severity.

## Machine Learning Models

Five ML-based approaches were trained on the k-mer frequency matrix generated from amino acids in the CDR3 region in the T cell repertoires of healthy donors and COVID-19 patients from the ISB-S datasets, using Python v3.8.6 and scikit-learn v0.23.1. These algorithms were: Random Forests (RF), Support Vector Machines (SVM), Naïve Bayes (NB), Gradient Boosting Classifiers (GBC), and K-Nearest Neighbors (KNN). The k-mer frequency matrix dataset was partitioned into subsets to perform binary classification between the healthy donor and the specified disease phenotype, such that models were trained for classification tasks of healthy donor vs moderate disease and healthy donor vs severe disease. The dataset was first partitioned into test and train sets with an 80:20 ratio. Following this test-train partition, to address imbalanced data, healthy donor samples were randomly resampled to be equal to the number of COVID-19 samples represented in the dataset, prior to training. Hyperparameter selection was informed by GridSearchCV, where optimal parameters found by the grid search were adopted when empirical performance on the test set was improved from the default parameters.


## Dependencies

```
Python  3.8.6
numpy 1.12.1
pandas 1.2.4
scikit-learn 0.23.1
Plotly 5.1.0
```

## Contact

Please see manuscript for further details. For any questions, contact jonathan [dot] park [at] yale [dot] edu.
This study was conducted by the Sidi Chen lab at Yale University (https://sidichenlab.org/).
