# mutitrait-multiancestry-lung

Code for "Multi-trait and multi-ancestry genetic analysis of comorbid lung diseases and traits improves genetic discovery and polygenic risk prediction" He Y, Lu W, Jee YH, et al. medRxiv, DOI: https://doi.org/10.1101/2024.08.25.24312558

This repository provides a framework for how to derive PRSxtra for your own phenotype of interest and phenotype.
`regularization.R` provides the script for fitting a regularization regression model on a set of scores. This script will also output prediction performance for the given phenotype.

## 1. Metanalysis of ancestry-specific summary statistics 
`1_metal.txt` and `1_runmetal.sh` provide the generalized script for running METAL. You can find more information about the software from their [wiki](https://genome.sph.umich.edu/wiki/METAL_Documentation).

For our analysis approach, we used the classical approach using effect size estimates and standard errors

## 2. Multi-trait analysis within ancestry groups
`2_mtag.sh` provides the generalized script for running MTAG. You can find more information about the software from their [github page](https://github.com/JonJala/mtag). 

After running MTAG, you can generate polygenic risk scores using `PRS-cs.sh`. Then you can run `regularization.R` to get PRS-xt.

