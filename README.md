# T cell receptor repertoire signatures associated with COVID-19 severity

Code and analysis for the manuscript "T cell receptor repertoire signatures associated with COVID-19 severity"

## Project Structure
```
.
├── code
│   ├── ML_models.py
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

## Contact

Please see manuscript for further details. For any questions, contact jonathan [dot] park [at] yale [dot] edu.
This study was conducted by the Sidi Chen lab at Yale University (https://sidichenlab.org/).
