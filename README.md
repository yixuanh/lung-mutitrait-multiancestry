# Analysis pipeline 


This repository provides a generalized framework for how to derive PRSxtra for your own phenotype of interest and phenotype. Please refer to our manuscript for details "Multi-trait and multi-ancestry genetic analysis of comorbid lung diseases and traits improves genetic discovery and polygenic risk prediction" He Y, Lu W, Jee YH, et al. medRxiv, DOI: https://doi.org/10.1101/2024.08.25.24312558

`regularization.R` provides the script for fitting a regularization regression model on a set of scores for a specific phenotype and population. This script will also output prediction performance for the given phenotype. This script will require the phenotype file (with the headers `ID`, `pheno`), a list of IDs for training and testing (`ID_train.csv` and `ID_test.csv`, and a compiled PRS file (with the headers `ID`, `PRS1`, `PRS2`,...).  

## 1. Meta-analysis of ancestry-specific summary statistics 
If needed, please conduct ancestry-specific meta-analysis. `1_metal.txt` and `1_runmetal.sh` provide the generalized script for running METAL. You can find more information about the software from their [wiki](https://genome.sph.umich.edu/wiki/METAL_Documentation).

For our analysis approach, we used the classical approach using effect size estimates and standard errors. In our project, we meta-analyzed two East Asian biobanks (KCPS2 and TWB) for FEV1, FVC, and FEV1/FVC. 

## 2. Multi-trait analysis within ancestry groups 
`2_mtag.sh` provides the generalized script for running MTAG. You can find more information about the software from their [github page](https://github.com/JonJala/mtag). 

After running MTAG on the set of N traits within each ancestry group, you can generate polygenic risk scores using `PRS-cs.sh`. You can find more information about PRScs from their [github page](https://github.com/getian107/PRScs). Then you can run `regularization.R` on the N scores to get PRS-xt. Repeat for each ancestry group.

## 3. Multi-ancestry polygenic risk score within traits
To generate PRS-xa, run `PRS-csx.sh` on the set of M ancestry-specific summary statistics available for each of your traits. You can find more information about PRScsx from their [github page](https://github.com/getian107/PRScsx). You should get M+1 (M ancestry specific and 1 meta-analyzed) PRS. Then you can run `regularization.R` on the M+1 scores to get PRS-xa. Repeat for each trait.

## 4. Multi-trait and multi-ancestry polygenic risk score (PRSxtra)
First, run `PRS-csx.sh` on the set of summary statistics generated using MTAG from step 2. Then run `regularization.R` on all the output scores to get PRSxtra.

## Additional information
We have included an example of evaluating the performance of PRSxtra to predict COPD in `prediction_COPD_PRSxtra.R`. This framework can be applied to the other scores and other phenotypes. 

In our manuscript, we conducted a pleitropy analysis between smoking behaviors and COPD and lung cancer. The corresponding code can be found at `lung_multitrait_linemodel.R` and `lung_multitrait_linemodel_utils.R`.

We used the LD reference panels generated from the 1000 Genome Project phase 3 samples that were generated by the PRScsx team.

