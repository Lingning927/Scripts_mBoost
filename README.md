# mBoost analysis scripts

This repository contains the core analysis code for the manuscript **"Identifying low quality features and suboptimal models to enhance phenotype prediction"**.

mBoost is a statistical framework for diagnosing and improving machine-learning models for phenotype prediction. It evaluates a model from three training-stage perspectives and one application-stage perspective:

- **Structure improvement (`ps`)**: tests whether the current model structure is still improvable.
- **Feature redundancy (`pr`)**: tests whether the feature set contains too many redundant/noisy variables for the current model.
- **Feature validity (`pv`)**: tests whether selected features are stable beyond what is expected under phenotype permutation.
- **Coverage Ratio (`CR`)**: measures whether a trained model is applicable to a new dataset by comparing the train and test feature distributions.

The manuscript validates mBoost with simulations, a drug-sensitivity prediction example, and 120 gut microbiome case-control cohorts covering 20 diseases. In the microbiome analysis, mBoost identifies suboptimal models and low-quality feature sets, then gives practical guidance such as using a more flexible model, applying feature selection, changing the feature-selection method, or avoiding model application to poorly covered external cohorts.

## Repository scope

This repository provides the core methodological implementation and compact scripts used for the main analyses. It is **not** a complete paper-build pipeline: some intermediate files, high-throughput batch scripts, and the full plotting workflow used to assemble every manuscript figure are not included. The included code and data are intended to make the main mBoost diagnostics reproducible and reusable.

## Directory layout

```text
methods/
  MAC.cpp                  Rcpp implementation of the MAC statistic
  NewMac.R                 R wrapper for MAC
  simulation_methods.R     bootstrap tests for structure improvement

scripts/simulations/
  example1.R               Lasso feature recovery with increasing noise
  example2.R               polynomial regression and structure testing
  example3.R               crescent-shaped classification example
  example4.R               linear-logistic feature redundancy example
  example5.R               nonlinear feature-selection simulation
  example6.R               negative-control simulation

scripts/real_data/
  structure_improvement.R  ps diagnostic for SIAMCAT microbiome models
  feature_redundancy.R     pr diagnostic and cross-cohort AUC checks
  feature_validaty.R       pv diagnostic for feature-selection methods
  coverage_ratio.R         CR calculation between microbiome cohorts

data/
  feat_meta/               processed microbiome feature/meta RData files
  cohort_info.RData        cohort-level metadata

fig/
  figures/                 final main-text figures
  supplemental_figures/    final supplemental figures
```

## Data

The processed microbiome cohort data used by the real-data scripts are included under `data/feat_meta/`. These data are derived from publicly available gut microbiome cohorts and follow the manuscript's real-data analysis setup.

The drug-sensitivity analysis in the manuscript follows the Precily framework. The corresponding external data/code archive is available from Zenodo: [10.5281/zenodo.7024834](https://doi.org/10.5281/zenodo.7024834).

## Requirements

The code was written in R and uses several CRAN/Bioconductor packages, depending on the script:

```r
install.packages(c(
  "Rcpp", "MASS", "glmnet", "ggplot2", "caret", "e1071",
  "randomForest", "ggsci", "patchwork", "Rtsne", "dplyr", "pROC",
  "kSamples", "xgboost", "keras3", "remotes"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("SIAMCAT", "mixOmics"))

remotes::install_github("Lingning927/mBoost")
```

The `parallel` package ships with R. Some scripts use `mclapply()` with multiple cores and are best run on Linux, macOS, or a Unix-like/HPC environment. On Windows, set `mc.cores = 1` or adapt the parallel loops before running large jobs.

Run scripts from the repository root so that relative paths such as `methods/NewMac.R` and `data/feat_meta/` resolve correctly.

## Example usage

```bash
Rscript scripts/simulations/example2.R
Rscript scripts/simulations/example3.R
Rscript scripts/simulations/example5.R 1 RF
```

For real microbiome data, choose the disease, data type, and method near the bottom of each script, then run for example:

```bash
Rscript scripts/real_data/structure_improvement.R
Rscript scripts/real_data/feature_redundancy.R
Rscript scripts/real_data/feature_validaty.R
Rscript scripts/real_data/coverage_ratio.R
```

The manuscript reports repeated simulations and cohort-wide sweeps. The scripts here are compact entry points; extending them to loop over all seeds, diseases, data types, and methods will reproduce the corresponding larger analysis tables.

## Citation

If you use this code, please cite the manuscript:

Shen B, Li M, Hu X, Chen M, Chen W-H, Jiang H. **Identifying low quality features and suboptimal models to enhance phenotype prediction**.
