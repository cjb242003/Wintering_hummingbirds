# Winter hummingbird analysis code

This repository contains the R code used for the primary analyses, sensitivity analyses, figures, and supplementary tables for the winter hummingbird study.

## Citation

If you use this code, please cite the archived release:

> Clarkson, A. E., B. A. Schwab, and C. J. Butler (2026). Code from: Age-specific geographic variation in sex ratios differs among hummingbird species wintering in eastern North America. Zenodo. https://doi.org/10.5281/zenodo.22284800

and the underlying data:

> Nakash, E., M. Malorodova, L.-A. Howes, and A. Celis-Murillo (2025). North American Bird Banding Program dataset 1960–2025 retrieved 2025-07-11. U.S. Geological Survey data release. https://doi.org/10.5066/P1KPZGAR

## License

The code in this repository is released under the MIT License. See [`LICENSE`](LICENSE) for the full text.

## Files

- `01_functions.R` — reusable data-preparation, model-fitting, and robustness-analysis functions.
- `02_primary_analysis.R` — reads the BBL data, builds the six species datasets, fits the primary three-species GLMM, summarizes the structure and record-type composition of the primary dataset, and calculates response-scale effect sizes reported in the manuscript.
- `03_sensitivity_analysis.R` — runs missing-data bounds, January age-classification checks, day-of-window and western-boundary sensitivity analyses, temporal-confounding analyses, leave-one-high-volume-location-out analyses, primary-model diagnostics including residual spatial autocorrelation, geographic variance-retention calculations, and nonlinear spatial GAMM sensitivity analyses.
- `04_figures.R` — generates Figures 1–4 and saves each figure as both PNG and vector PDF. Also generates Table 1.
- `05_supplement.R` — generates Supplementary Tables S1–S10 and saves them as `Supplementary_Tables.docx`.
- `gamm_results_PRIMARY_NOV15_DEC31_EPSG5070.rds` — cached results from the final nonlinear spatial GAMM analyses and basis-dimension checks.

## Input data

The analyses use a single input file, the hummingbird species group (group 13) extract from the North American Bird Banding Program dataset:

`NABBP_2025_grp_13.csv`

Download it from the U.S. Geological Survey data release at https://doi.org/10.5066/P1KPZGAR (the file is also reachable directly through the [ScienceBase catalog item](https://www.sciencebase.gov/catalog/item/68837a85d4be027deba86316)), and place it in the repository root alongside the R scripts.

**Version note.** The analyses reported in the manuscript used the extract retrieved on **2025-07-11**, which contains 1,186,275 rows. USGS data releases are periodically revised, and a later version of this file may not reproduce the reported record counts exactly. If your download differs in row count from the value above, you have a different version of the release.

## Software requirements

The final analyses were run under R 4.6.1. Full version information for every attached and loaded package is given in the [R session information](#r-session-information) section below.

Install the contributed packages used by the analysis with:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "purrr",
  "lme4", "emmeans", "DHARMa", "ape", "mgcv",
  "sf", "ggplot2", "scales",
  "rnaturalearth", "rnaturalearthdata",
  "broom.mixed", "flextable", "officer"
))
```

`sf` requires the system libraries GDAL, GEOS, and PROJ. These are bundled with the macOS and Windows CRAN binaries but must be installed separately on most Linux distributions (e.g. `libgdal-dev`, `libgeos-dev`, `libproj-dev` on Debian/Ubuntu). This is the most common installation failure point.

## Run order

The scripts assume that the working directory is the repository root and use relative paths throughout. No script calls `setwd()`.

`02_primary_analysis.R` automatically sources `01_functions.R`. From a clean R session, run:

```r
source("02_primary_analysis.R")
source("03_sensitivity_analysis.R")
source("04_figures.R")
source("05_supplement.R")
```

**Scripts `02`–`05` must be run sequentially within the same R session.** Scripts `03`, `04`, and `05` depend on objects created by the scripts before them and will fail if sourced on their own in a fresh session. The primary model is fit in `02_primary_analysis.R`; the sensitivity, figure, and supplement scripts read that model but do not alter it.

### Approximate runtimes

On the machine described in the session information below:

- `02_primary_analysis.R` — approximately 5 minutes, including the roughly 1 minute needed to source `01_functions.R`.
- `03_sensitivity_analysis.R` — approximately 10 minutes when the GAMM cache is present, and approximately 90 minutes when the GAMMs are refit from scratch.
- `04_figures.R` and `05_supplement.R` — approximately 1 minute each.

## GAMM results cache

The nonlinear spatial GAMMs and `mgcv::k.check()` calculations in `03_sensitivity_analysis.R` are substantially more computationally expensive than the other analyses. Their downstream results are stored in:

`gamm_results_PRIMARY_NOV15_DEC31_EPSG5070.rds`

The archived repository includes the cache used for the final analyses. When this file exists and its stored data signature and settings match the current primary analysis, `03_sensitivity_analysis.R` loads the cached results instead of refitting the GAMMs or rerunning the basis-dimension checks. If the cache is absent, the script fits the GAMMs, runs the checks, and creates the cache. Set `force_recompute_gamms <- TRUE` in `03_sensitivity_analysis.R` only when a complete GAMM rerun is intended.

`05_supplement.R` reads the same cache to construct Supplementary Table S7 and verifies that it matches the current analysis data and GAMM settings.

## Main outputs

Running `04_figures.R` creates:

- `Figure1.png` and `Figure1.pdf`
- `Figure2.png` and `Figure2.pdf`
- `Figure3.png` and `Figure3.pdf`
- `Figure4.png` and `Figure4.pdf`
- `Table1.docx`

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
- `ape`
- `mgcv`
- `sf`
- `ggplot2`
- `scales`
- `rnaturalearth`
- `rnaturalearthdata`
- `broom.mixed`
- `flextable`
- `officer`

## R session information

The final analyses were run under the following R environment:

```text
R version 4.6.1 (2026-06-24)
Platform: aarch64-apple-darwin23
Running under: macOS Tahoe 26.6.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib; LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] DHARMa_0.5.0            emmeans_2.0.4           broom.mixed_0.2.9.7
 [4] officer_0.7.6           flextable_0.10.0        rnaturalearthdata_1.0.0
 [7] rnaturalearth_1.2.0     scales_1.4.0            ggplot2_4.0.3
[10] purrr_1.2.2             tidyr_1.3.2             sf_1.1-2
[13] mgcv_1.9-4              nlme_3.1-169            lme4_2.0-6
[16] Matrix_1.7-5            dplyr_1.2.1             data.table_1.18.6.1

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1        viridisLite_0.4.3       farver_2.1.2
 [4] S7_0.2.2                fastmap_1.2.0           fontquiver_0.2.1
 [7] promises_1.5.0          digest_0.6.39           mime_0.13
[10] estimability_2.0.0      lifecycle_1.0.5         magrittr_2.0.5
[13] compiler_4.6.1          rlang_1.3.0             tools_4.6.1
[16] utf8_1.2.6              knitr_1.51              labeling_0.4.3
[19] askpass_1.2.1           classInt_0.4-11         plyr_1.8.9
[22] xml2_1.6.0              RColorBrewer_1.1-3      gap.datasets_0.0.6
[25] KernSmooth_2.23-26      withr_3.0.3             grid_4.6.1
[28] gdtools_0.5.1           xtable_1.8-8            e1071_1.7-17
[31] future_1.75.0           iterators_1.0.14        globals_0.19.1
[34] MASS_7.3-65             cli_3.6.6               mvtnorm_1.4-2
[37] rmarkdown_2.31          ragg_1.5.2              reformulas_0.4.4
[40] generics_0.1.4          otel_0.2.0              minqa_1.2.8
[43] DBI_1.3.0               ape_5.8-1               proxy_0.4-29
[46] maps_3.4.3              splines_4.6.1           parallel_4.6.1
[49] vctrs_0.7.3             boot_1.3-32             fontBitstreamVera_0.1.1
[52] qgam_2.0.0              listenv_1.0.0           systemfonts_1.3.2
[55] foreach_1.5.2           units_1.0-1             gap_1.15.2
[58] glue_1.8.1              parallelly_1.48.0       nloptr_2.2.1
[61] codetools_0.2-20        gtable_0.3.6            later_1.4.8
[64] tibble_3.3.1            furrr_0.4.0             pillar_1.11.1
[67] htmltools_0.5.9         openssl_2.4.2           R6_2.6.1
[70] textshaping_1.0.5       Rdpack_2.6.6            doParallel_1.0.17
[73] shiny_1.14.0            evaluate_1.0.5          lattice_0.22-9
[76] rbibutils_2.4.1         backports_1.5.1         broom_1.0.13
[79] httpuv_1.6.17           fontLiberation_0.1.0    class_7.3-23
[82] Rcpp_1.1.2              zip_3.0.2               uuid_1.2-2
[85] coda_0.19-4.1           xfun_0.60               forcats_1.0.1
[88] pkgconfig_2.0.3
```
