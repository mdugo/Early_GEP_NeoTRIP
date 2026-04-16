# Bulk RNA-seq analysis of pre- and on-treatment TNBC samples from the NeoTRIP (NCT02620280) trial

## Overview

This repository contains script and data necessary to reproduce the analyses published in Dugo et al., 2026.

## Data Availability

Raw sequencing data are not included due to size and privacy constraints.

Processed gene level data can be accessed at NCBI GEO repository with accession number GSE319641.
A script to retrieve and assemble data from NCBI GEO is provided (script 01-data_download_and_singscore.R).
Clinical data are not publicly available due to privacy restrictions and can be obtained upon reasonable request by contacting the authors at the following email address: translational@fondazionemichelangelo.org

## Reproducibility

To reproduce the analysis, first open the project by clicking on the .Rproj file. This ensures that the correct project root is set and that all file paths defined using the [here package](https://doi.org/10.32614/CRAN.package.here) will work correctly.

The R environment required to reproduce all analyses is provided through the [renv package](https://doi.org/10.32614/CRAN.package.renv), ensuring full computational reproducibility of the workflow across different systems. Restore the R environment by running the following command:

```r
if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv")

renv::restore()
```
Run the 00b-project_setup.R to create all the necessary folders.
Once the environment has been restored and the project structure has been initialized, the remaining scripts can be executed sequentially to fully reproduce the analysis workflow and generate all results and figures.

## Notes
The script 07-Figure6.R that implements XGBoost model training with repeated cross-validation is computationally intensive and may require a significant amount of time to run, depending on the available hardware.

## Author
Matteo Dugo

## License

MIT