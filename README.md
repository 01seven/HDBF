# HDBF

This repository contains the R code accompanying the manuscript:

**Power-Enhanced Cauchy Combination Tests for the Two-Sample High-Dimensional Behrens-Fisher Problem**

The scripts reproduce the Monte Carlo simulation study and the real-data gene ontology (GO) analysis reported in the manuscript.

## Repository contents

- `HDBF_simulation_study.R`: Monte Carlo simulation code for type-I-error and power studies.
- `HDBF_real_data_GO_analysis.R`: real-data GO-set analysis for the ALL leukemia data set and Table 2 summaries.
- `CITATION.cff`: citation metadata for this software repository.

## Tested environment

- R 4.4.3
- Windows

## Required R packages

CRAN packages:

```r
install.packages(c(
  "future", "future.apply", "foreach", "doParallel",
  "highDmean", "highmean", "matrixStats", "mvnfast", "knitr"
))
```

Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "ALL", "genefilter", "hgu95av2.db", "GO.db", "AnnotationDbi"
))
```

## Running the simulation study

The full simulation uses 2000 Monte Carlo replications per covariance structure and can be time-consuming.

```r
source("HDBF_simulation_study.R")
results <- run_all()
```

For a quick smoke check:

```r
source("HDBF_simulation_study.R")
cfg <- hdbf_config()
cfg$n_simulations <- 10
cfg$workers <- 1
cfg$use_progress <- FALSE
results <- run_all(cfg, write_csv = FALSE)
```

When `write_csv = TRUE`, the script writes files named `HDBF_type1_error_results_*.csv` and `HDBF_power_results_*.csv`.

## Running the real-data GO analysis

```r
source("HDBF_real_data_GO_analysis.R")
results <- run_real_data_analysis()
```

To reduce CPU use:

```r
source("HDBF_real_data_GO_analysis.R")
cfg <- hdbf_real_data_config()
cfg$workers <- 1
results <- run_real_data_analysis(cfg, write_report = TRUE)
```

The real-data script writes `WCCT_Advantage_Report.txt` by default.

## Data availability

The real-data analysis uses publicly available acute lymphoblastic leukemia data distributed through the Bioconductor `ALL` package. The manuscript also cites the original data source and related comparison material at:

- https://doi.org/10.1182/blood-2003-09-3243
- https://arxiv.org/abs/2405.02551

## Reproducibility notes

- Both scripts stop early with a clear message if required packages are missing.
- Parallel computation is enabled by default. Set `workers <- 1` in the configuration object for single-core execution.
- Generated CSV and text report outputs are ignored by Git by default.
