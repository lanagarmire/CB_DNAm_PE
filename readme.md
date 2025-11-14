## Project Summary

This repository contains the code supporting the manuscript:

> **“Cell type proportions rather than DNA methylation in the cord blood show significant associations with severe preeclampsia”**

### Authors
Xiaotong Yang✝, Wenting Liu✝, Zhixin Mao, Yuheng Du, Cameron Lassiter, Fadhl M. AlAkwaa, Paula A. Benny, Lana X. Garmire*  
✝ Equal contribution  
*Corresponding author: Lana X. Garmire, PhD, University of Michigan

### Repository Structure
```
├── scripts/                         
│   ├── in-house/                    # Analysis of Hawaii Biorepository (HiBR) data
│   │   ├── 1_in-house_raw_data_preprocessing.Rmd
│   │   ├── 2_in-house_data_analysis.R
│   │   ├── 2.1_in-house_DMR.R
│   │   ├── 2.2_in-house_DM_on_gene_level.R
│   │   └── 2.3_in-house_Pathifier.R
│
│   ├── public_data/                 # Individual analyses of published datasets
│   │   ├── 3.1_Ching_data_preprocess.R
│   │   ├── 3.2_SOV_Ching_data.R
│   │   ├── 3.3_Ching_data_analysis.R
│   │   ├── 4.1_Herzog_data_analysis.R
│   │   ├── 4.2_SOV_Herzog.R
│   │   └── 5_Kashima_data_analysis_PBMC.R
│
│   ├── pooled_analysis/             # Cross-cohort pooled and harmonization analysis
│   │   └── 6_pooled_data_analysis.R
│
│   ├── preterm_analysis/            # Decoupling preeclampsia and gestational age effects
│   │   └── 7_preterm_control_Fernando_analysis.R
├── LICENSE
├── README.md
```