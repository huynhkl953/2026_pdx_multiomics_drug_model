# PDX Multi-Omics Drug Response Prediction Model

This repository contains the multi-omics drug-response prediction pipeline for Patient-Derived Xenograft (PDX) models[cite: 1, 2]. The workflow evaluates and compares the predictive performance of transcriptomic (RNA-seq) versus proteomic (Mass Spectrometry) data across 1 µM and 10 µM drug concentrations using a ComBat-harmonized ridge-regression framework[cite: 1, 2, 3].

## Repository Structure

```text
2026_pdx_multiomics_drug_model/
├── data/                                      # Input dataset directory
│   ├── 25.11.20_HUMAN_U54_AllSamples_...csv   # RNA-seq gene expression TPM
│   ├── 20251203_DetectedProteins_...csv       # Mass spectrometry proteomics
│   ├── 20251113_All_U54_1uM_inhibition.xlsx   # 1 µM drug inhibition matrix
│   └── 20251113_All_U54_10uM_inhibition.xlsx  # 10 µM drug inhibition matrix
├── function_oncoPredict.R                     # Modified oncoPredict modeling engine
├── predict_drug_score.R                       # Main 9-fold cross-validation pipeline
├── final_revision_Nina.Rmd                    # R Markdown evaluation & figure generator
└── README.md                                  # Repository documentation
