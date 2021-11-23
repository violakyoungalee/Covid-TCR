# T cell receptor repertoire signatures associated with COVID-19 severity

Code and analysis for the manuscript "T cell receptor repertoire signatures associated with COVID-19 severity"

## Project Structure
```
.
├── code
│   ├── ML_models.py
│   ├── Kmer_matrix.R
│   └── significant_features.py
├── data
│   ├── main
│   ├── intermediate_files
│   └── output
├── figures
└── tables

```
## Data acquisition

CD4+ and CD8+ T cell repertoires for COVID patients and healthy donors (HD) were accessed through Array Express using the accession number E-MTAB-9357 (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9357/). A metadata file, containing the sample ID, blood draw time point, COVID disease severity, sex and age, etc., was compiled using Table S1 (https://data.mendeley.com/datasets/tzydswhhb5/5). COVID patients for which sequences at two time points were taken were labelled according to the time point (e.g. 100_1 or 100_2). Also, COVID severity scores were translated from the WHO Ordinal Scale (0-7) to four tiers: mild (WOS: 1-3), moderate (4-5), severe (6-7), and healthy donor (0) (Marshall 2020, https://www.who.int/docs/default-source/documents/emergencies/minimalcoreoutcomemeasure.pdf). Only patients which were sequenced by 10X Genomics (column tenX with a value of yes) were kept for the analysis. Then, corresponding metadata information was joined to individual repertoire files. After patient filtering, sequences were filtered: for barcodes with TRA and TRB data for more than one sequence, only the first TRA and TRB were taken. After the above filtering, the aggregated CD4+ and CD8+ datasets were composed of 107,397 and 55,994 sequences, represented in a total of 16 healthy, 50 mild, 52 moderate, and 27 severe patients (16 HD, 109 mild, 95 moderate, and 49 severe repertoires when accounting for patients with samples from two time points).

## Dependencies

```
Python  3.8.6
numpy 1.12.1
pandas 1.2.4
scikit-learn 0.23.1
Plotly 5.1.0
```
## Pipeline

#### Kmer Matrix Generation
A list of all the overlapping substrings of k (k = 3, 4, 5, 6) amino acids in the CDR3 region in the CD4 and CD8+ T cell repertoire of healthy donor and covid-19 patients was generated. From each individual patient in our dataset, the frequency of each k-mer in the list was determined. The frequency counts of each kmer of the patients was organized as a data matrix where rows corresponded to patient IDs and columns corresponded to the frequency counts of the kmers found in each patient. A data column representing the patient’s disease phenotype was included.

#### Classical Machine Learning 

Training and Evaluating Machine Learning Models:
ML_models.py contains code that trains classical machine learning models and visualizes their ROC curves using the plotly library. 
Five machine learning algorithms were selected for training on kmer frequency matrix for 32 different permutations: for CD4+ and CD8+ TCR datasets, classification of HD vs. Mild, HD vs. Moderate, HD vs. Severe, HD vs. All COVID were trained for 3mers, 4mers, 5mers and 6mers. These algorithms were Random Forest, Support Vector Machine, Gaussian Naïve Bayes (GNB), Gradient Boosting (GB) and K-Nearest Neighbors (KNN). Testing and training set was partitioned with a 80:20 ratio using sklearn’s test_train_split function. Training features included all kmers represented in the CD8+ TCR repertoire’s alpha chain CDR3 of all patients in the dataset. Five-fold cross-validation was performed within the train set with 100 repetitions (5 x 100) using sklearn’s RepeatedStratifiedKFold function. Performance on the validation set was used to generate ROC plots and performance on test and validation sets were exported. 

#### Interpreting Machine Learning Models: 
To evaluate significant features detected by the machine learning models, a list of top 10,000 kmers conferring the most significance to the support vector machine classifier and the random forest classifier was obtained for each iteration of 100 repetitions of 5-fold cross validation. The top significant features were determined as the kmers commonly present among the top 10,000 kmers across the 500 iterations the models were trained. A separate list of such features was obtained for the support vector machine classifier and the random forest classifier, and the kmers common to the lists from both algorithms were determined to be the kmers with shared importance in distinguishing the immunological profile of healthy donors and covid-19 infected individuals. 

## Contact

Please see manuscript for further details. For any questions, contact jonathan [dot] park [at] yale [dot] edu.
This study was conducted by the Sidi Chen lab at Yale University (https://sidichenlab.org/).
