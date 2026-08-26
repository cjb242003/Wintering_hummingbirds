# Winter hummingbird analysis code

This repository contains the R code used for the primary analyses, sensitivity analyses, figures, and supplementary tables for the winter hummingbird study.

## Files

- `01_functions.R` — reusable data-preparation, model-fitting, and robustness-analysis functions.
- `02_primary_analysis.R` — reads the BBL data, builds the six species datasets, fits the primary three-species GLMM, and generates the standardized Rufous latitude predictions reported in the manuscript.
- `03_sensitivity_analysis.R` — runs missing-data bounds, January age-classification checks, temporal-confounding analyses, leave-one-high-volume-location-out analyses, primary-model diagnostics, and nonlinear spatial GAMM sensitivity analyses.
- `04_figures.R` — generates Figures 1-4 and saves each figure as both PNG and vector PDF. Also generates Table 1.
- `05_supplement.R` — generates Supplementary Tables S1-S14 and saves them as `Supplementary_Tables.docx`.

## Input data

Place the Bird Banding Laboratory data extract

`NABBP_2025_grp_13.csv`

in the working directory containing the R scripts. The file can be downloaded at [ScienceBase](https://www.sciencebase.gov/catalog/item/68837a85d4be027deba86316).

## Run order

`02_primary_analysis.R` automatically sources `01_functions.R`. From a clean R session, run:

```r
source("02_primary_analysis.R")
source("03_sensitivity_analysis.R")
source("04_figures.R")
source("05_supplement.R")
```

The primary model is fit in `02_primary_analysis.R`. The figure and supplement scripts use objects created by the primary and sensitivity scripts and do not alter the primary model.

## GAMM results cache

The nonlinear spatial GAMMs and `mgcv::k.check()` calculations in `03_sensitivity_analysis.R` are substantially more computationally expensive than the other analyses. Their downstream results are stored in:

`gamm_results_PRIMARY_NOV15_DEC31_EPSG5070.rds`

When this file exists and its stored data signature and settings match the current primary analysis, `03_sensitivity_analysis.R` loads the cached results instead of refitting the GAMMs or rerunning the basis-dimension checks. If the cache is absent, the script fits the GAMMs, runs the checks, and creates the cache. Set `force_recompute_gamms <- TRUE` in `03_sensitivity_analysis.R` only when a complete GAMM rerun is intended.

`05_supplement.R` reads the same cache to construct Supplementary Table S10 and verifies that it matches the current analysis data and GAMM settings.

## Main outputs

Running `04_figures.R` creates:

- `Figure1.png` and `Figure1.pdf`
- `Figure2.png` and `Figure2.pdf`
- `Figure3.png` and `Figure3.pdf`
- `Figure4.png` and `Figure4.pdf`
- `Table1.docx`

The PDF files are vector graphics suitable for publication and remain sharp when enlarged.

Running `05_supplement.R` creates:

- `Supplementary_Tables.docx`

## R packages

The analysis uses the following contributed R packages:

- `data.table`
- `dplyr`
- `tidyr`
- `purrr`
- `lme4`
- `emmeans`
- `DHARMa`
- `mgcv`
- `sf`
- `ggplot2`
- `scales`
- `rnaturalearth`
- `rnaturalearthdata`
- `broom.mixed`
- `flextable`
- `officer`

Package versions used for the final analysis should be reported with the archived code or repository release if exact environment reconstruction is required.
