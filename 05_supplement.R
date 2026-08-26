# ============================================================
# 05_supplement.R
# Generate Supplementary Tables S1-S13 from objects created by the
# primary and sensitivity analyses. Run 02_primary_analysis.R and
# 03_sensitivity_analysis.R before this script. Cached GAMM results are
# formatted here but are not refit.
#
# Supplementary tables:
#   S1  Demographic composition
#   S2  Primary comparative GLMM fixed effects
#   S3  Primary comparative GLMM random effects
#   S4  Primary age-specific geographic slopes
#   S5  Primary pairwise interspecific contrasts
#   S6  Temporal-confounding omnibus tests
#   S7  Temporal-confounding slopes and contrasts
#   S8  Leave-one-high-volume-location-out sensitivity
#   S9  GAMM nonlinear-spatial sensitivity
#   S10 Extreme-allocation demographic bounds
#   S11 Missing demographic information
#   S12 Primary comparative model diagnostics
#   S13 Geographic sampling support by latitude band
# ============================================================


library(dplyr)
library(tidyr)
library(purrr)
library(broom.mixed)
library(flextable)
library(officer)


# ============================================================
# 1. VALIDATE REQUIRED INPUTS AND LOAD CACHED GAMM RESULTS
# ============================================================

required_objects <- c(
  "results",
  "comparative_primary",
  "comparative_df",
  "m_comparative_full",
  "comparative_temporal",
  "station_influence_top_sites",
  "station_influence_lrt_table",
  "station_influence_age_differences",
  "station_influence_species_contrasts",
  "sensitivity_results",
  "missingness_summary_all",
  "comparative_primary_diagnostics"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Missing required object(s): ",
    paste(missing_objects, collapse = ", "),
    ". Run 01_functions.R, 02_primary_analysis.R, and ",
    "03_sensitivity_analysis.R first."
  )
}


retained_species <- c(
  "Rufous",
  "Black-chinned",
  "Ruby-throated",
  "Calliope",
  "Broad-tailed",
  "Allen's"
)

focal_species_order <- c(
  "Rufous",
  "Black-chinned",
  "Ruby-throated"
)

descriptive_species <- c(
  "Calliope",
  "Broad-tailed",
  "Allen's"
)

lat_band_width <- 2

gamm_results_cache_file <-
  "gamm_results_PRIMARY_NOV15_DEC31_EPSG5070.rds"


# ------------------------------------------------------------
# Load cached GAMM summaries created by 03_sensitivity_analysis.R.
# ------------------------------------------------------------

if (!file.exists(gamm_results_cache_file)) {
  stop(
    "GAMM results cache not found: ",
    gamm_results_cache_file,
    ". Run 03_sensitivity_analysis.R to create it."
  )
}

gamm_results <- readRDS(gamm_results_cache_file)

required_gamm_results <- c(
  "smooth_summary",
  "k_check_summary",
  "model_metadata",
  "data_signature",
  "settings"
)

if (
  !is.list(gamm_results) ||
  !all(required_gamm_results %in% names(gamm_results))
) {
  stop(
    "The GAMM results cache is incomplete or incompatible. ",
    "Recreate it with 03_sensitivity_analysis.R."
  )
}

gamm_sensitivity_summary <- as.data.frame(
  gamm_results$smooth_summary
)

gamm_kcheck_summary <- as.data.frame(
  gamm_results$k_check_summary
)

gamm_model_metadata <- as.data.frame(
  gamm_results$model_metadata
)


# Confirm that the cached GAMM results match the current primary dataset
# and GAMM specification.
expected_gamm_settings <- list(
  species = focal_species_order,
  basis_k = 20L,
  projection = "EPSG:5070",
  spatial_units = "100 km",
  spatial_scale_m = 100000,
  family = "binomial(logit)",
  method = "fREML",
  discrete = TRUE,
  kcheck_seed = 20260823L
)

if (!isTRUE(
  all.equal(
    gamm_results$settings,
    expected_gamm_settings,
    check.attributes = FALSE
  )
)) {
  stop(
    "The GAMM cache was created with settings that do not match the ",
    "published GAMM specification."
  )
}

build_supplement_gamm_signature <- function(df) {
  required <- c(
    "Species",
    "band",
    "winter_year",
    "lat_dd",
    "lon_dd",
    "Female",
    "Immature_f"
  )
  
  missing_cols <- setdiff(required, names(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "`comparative_df` is missing required GAMM signature column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  x <- data.frame(
    Species = as.character(df$Species),
    band = as.character(df$band),
    winter_year = as.character(df$winter_year),
    lat_dd = as.numeric(df$lat_dd),
    lon_dd = as.numeric(df$lon_dd),
    Female = as.integer(df$Female),
    Immature_f = as.character(df$Immature_f),
    stringsAsFactors = FALSE
  )
  
  x <- x[
    do.call(
      order,
      x
    ),
    ,
    drop = FALSE
  ]
  
  rownames(x) <- NULL
  
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  
  saveRDS(
    x,
    tmp,
    version = 2
  )
  
  unname(tools::md5sum(tmp))
}

current_gamm_signature <-
  build_supplement_gamm_signature(
    comparative_df
  )

if (!identical(
  current_gamm_signature,
  gamm_results$data_signature$md5
)) {
  stop(
    "The GAMM results cache does not match the current primary ",
    "comparative dataset. Do not generate Table S9 from this cache."
  )
}

# ============================================================
# 2. FORMATTING HELPERS
# ============================================================

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(
      x < 0.001,
      "<0.001",
      sprintf("%.3f", x)
    )
  )
}


fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    sprintf(paste0("%.", digits, "f"), x)
  )
}


nice_ft <- function(x, font_size = 9) {
  flextable(x) %>%
    theme_booktabs() %>%
    fontsize(size = font_size, part = "all") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    valign(valign = "center", part = "all") %>%
    autofit()
}


# Common formatting for supplementary-table captions.
caption_text_fp <- fp_text(
  font.family = "Times New Roman",
  font.size = 10,
  bold = TRUE,
  italic = TRUE
)

caption_fp <- fp_par(
  text.align = "center",
  padding.top = 0,
  padding.bottom = 0,
  line_spacing = 1,
  keep_with_next = TRUE
)


get_ci_cols <- function(x) {
  lower <- intersect(
    c("asymp.LCL", "lower.CL"),
    names(x)
  )[1]
  
  upper <- intersect(
    c("asymp.UCL", "upper.CL"),
    names(x)
  )[1]
  
  if (
    length(lower) == 0L ||
    length(upper) == 0L ||
    is.na(lower) ||
    is.na(upper)
  ) {
    stop("Could not identify confidence-interval columns.")
  }
  
  c(lower, upper)
}


standardize_emm <- function(x, axis, result_type) {
  x <- as.data.frame(x)
  ci <- get_ci_cols(x)
  
  if ("estimate" %in% names(x)) {
    est <- x$estimate
  } else {
    trend_col <- grep("\\.trend$", names(x), value = TRUE)[1]
    
    if (length(trend_col) == 0L || is.na(trend_col)) {
      stop("Could not identify trend column.")
    }
    
    est <- x[[trend_col]]
  }
  
  tibble(
    Axis = axis,
    Result = result_type,
    Species =
      if ("Species" %in% names(x)) as.character(x$Species) else "",
    Age =
      if ("Immature_f" %in% names(x)) as.character(x$Immature_f) else "",
    Contrast =
      if ("contrast" %in% names(x)) as.character(x$contrast) else "",
    Estimate = est,
    SE = x$SE,
    CI_low = x[[ci[1]]],
    CI_high = x[[ci[2]]],
    z =
      if ("z.ratio" %in% names(x)) x$z.ratio else NA_real_,
    p = x$p.value
  )
}


extract_lrt <- function(x, test_name, axis) {
  a <- as.data.frame(x)
  i <- nrow(a)
  
  tibble(
    Test = test_name,
    Axis = axis,
    `Chi-square` = as.numeric(a$Chisq[i]),
    df = as.integer(a$Df[i]),
    p = as.numeric(a[["Pr(>Chisq)"]][i])
  )
}


clean_fixed_term <- function(x) {
  vapply(
    x,
    function(term) {
      if (term == "(Intercept)") {
        return("Intercept")
      }
      
      parts <- strsplit(term, ":", fixed = TRUE)[[1]]
      
      parts <- vapply(
        parts,
        function(part) {
          dplyr::recode(
            part,
            "SpeciesBlack-chinned" = "Black-chinned",
            "SpeciesRuby-throated" = "Ruby-throated",
            "Immature_fImmature" = "Immature",
            "lat5" = "Latitude",
            "lon5" = "Longitude",
            .default = part
          )
        },
        character(1)
      )
      
      paste(parts, collapse = " × ")
    },
    character(1)
  )
}


clean_species_contrast <- function(x) {
  case_when(
    grepl("Rufous.*Black-chinned", x) ~
      "Rufous – Black-chinned",
    grepl("Rufous.*Ruby-throated", x) ~
      "Rufous – Ruby-throated",
    grepl("Black-chinned.*Ruby-throated", x) ~
      "Black-chinned – Ruby-throated",
    TRUE ~ x
  )
}


format_location <- function(latitude, longitude) {
  paste0(
    sprintf("%.3f", abs(latitude)),
    ifelse(latitude >= 0, "°N", "°S"),
    ", ",
    sprintf("%.3f", abs(longitude)),
    ifelse(longitude >= 0, "°E", "°W")
  )
}


# ============================================================
# TABLE S1. Demographic composition
# ============================================================

table_s1 <- map_dfr(
  retained_species,
  function(sp) {
    d <- results[[sp]]$data_points
    
    n_total <- nrow(d)
    n_immature <- sum(d$Immature == 1, na.rm = TRUE)
    n_adult <- sum(d$Immature == 0, na.rm = TRUE)
    n_female <- sum(d$Female == 1, na.rm = TRUE)
    n_adult_female <- sum(
      d$Immature == 0 & d$Female == 1,
      na.rm = TRUE
    )
    n_immature_female <- sum(
      d$Immature == 1 & d$Female == 1,
      na.rm = TRUE
    )
    
    tibble(
      Species = sp,
      N = n_total,
      `Immature, n (%)` =
        sprintf(
          "%d (%.1f)",
          n_immature,
          100 * n_immature / n_total
        ),
      `Female overall, n (%)` =
        sprintf(
          "%d (%.1f)",
          n_female,
          100 * n_female / n_total
        ),
      `Adult female, n (%)` =
        sprintf(
          "%d (%.1f)",
          n_adult_female,
          100 * n_adult_female / n_adult
        ),
      `Immature female, n (%)` =
        sprintf(
          "%d (%.1f)",
          n_immature_female,
          100 * n_immature_female / n_immature
        )
    )
  }
)

caption_s1 <- paste0(
  "Table S1. Demographic composition of hummingbird records ",
  "within the 15 November–31 December primary sampling window ",
  "for the six species retained in the study. Only records with ",
  "classifiable age and known sex are included. Immature and overall ",
  "female percentages are calculated among all included records; ",
  "adult and immature female percentages are calculated within ",
  "the corresponding age class."
)

ft_s1 <- nice_ft(table_s1)


# ============================================================
# TABLE S2. Primary comparative GLMM fixed effects
# ============================================================

table_s2 <- broom.mixed::tidy(
  m_comparative_full,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
) %>%
  transmute(
    Term = clean_fixed_term(term),
    Estimate = round(estimate, 3),
    SE = round(std.error, 3),
    `95% CI` = sprintf("%.3f to %.3f", conf.low, conf.high),
    z = round(statistic, 3),
    p = fmt_p(p.value)
  )

caption_s2 <- paste0(
  "Table S2. Fixed-effect coefficients from the primary comparative ",
  "generalized linear mixed model of hummingbird sex ratio. Rufous ",
  "Hummingbird and adults are the reference levels. Latitude and ",
  "longitude are expressed per 5° geographic change and centered ",
  "at 31°N and 90°W."
)

ft_s2 <- nice_ft(table_s2)


# ============================================================
# TABLE S3. Primary comparative GLMM random effects
# ============================================================

random_effect_levels <- c(
  species_site_id =
    dplyr::n_distinct(comparative_df$species_site_id),
  species_winter_year =
    dplyr::n_distinct(comparative_df$species_winter_year),
  species_band_id =
    dplyr::n_distinct(comparative_df$species_band_id)
)

table_s3 <- as.data.frame(
  lme4::VarCorr(m_comparative_full)
) %>%
  filter(is.na(var2)) %>%
  transmute(
    `Random effect` = dplyr::recode(
      grp,
      "species_site_id" = "Species × reported location",
      "species_winter_year" = "Species × year",
      "species_band_id" = "Species × band identity",
      .default = grp
    ),
    Levels = as.integer(random_effect_levels[grp]),
    Variance = round(vcov, 3),
    SD = round(sdcor, 3)
  )

caption_s3 <- paste0(
  "Table S3. Random-effect variance components from the primary ",
  "comparative generalized linear mixed model. Random intercepts ",
  "were included for species × reported location, species × year, ",
  "and species × band identity. Variance components are reported ",
  "on the logit scale."
)

ft_s3 <- nice_ft(table_s3)


# ============================================================
# TABLE S4. Primary age-specific geographic slopes
# ============================================================

table_s4 <- bind_rows(
  standardize_emm(
    comparative_primary$latitude_slopes,
    "Latitude",
    "Age-specific slope"
  ),
  standardize_emm(
    comparative_primary$longitude_slopes,
    "Longitude",
    "Age-specific slope"
  )
) %>%
  transmute(
    Axis,
    Species,
    Age,
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    `95% CI` = sprintf("%.3f to %.3f", CI_low, CI_high),
    z = round(z, 3),
    p = fmt_p(p)
  )

caption_s4 <- paste0(
  "Table S4. Age-specific latitude and longitude slopes estimated ",
  "from the primary comparative model. Estimates represent change ",
  "in log-odds of a record being female per 5° geographic change."
)

ft_s4 <- nice_ft(table_s4)


# ============================================================
# TABLE S5. Primary pairwise interspecific contrasts
# ============================================================

table_s5 <- bind_rows(
  standardize_emm(
    comparative_primary$latitude_species_contrasts,
    "Latitude",
    "Species contrast"
  ),
  standardize_emm(
    comparative_primary$longitude_species_contrasts,
    "Longitude",
    "Species contrast"
  )
) %>%
  transmute(
    Axis,
    Contrast = clean_species_contrast(Contrast),
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    `95% CI` = sprintf("%.3f to %.3f", CI_low, CI_high),
    z = round(z, 3),
    p = fmt_p(p)
  )

caption_s5 <- paste0(
  "Table S5. Tukey-adjusted pairwise interspecific comparisons ",
  "of age-dependent latitude and longitude slopes from the primary ",
  "comparative model. Estimates represent the immature–adult slope ",
  "difference in the first species minus that in the second species."
)

ft_s5 <- nice_ft(table_s5)


# ============================================================
# TABLE S6. Temporal-confounding omnibus tests
# ============================================================

table_s6 <- bind_rows(
  extract_lrt(
    comparative_temporal$lrt$overall,
    "Species × Age × Geography",
    "Joint"
  ),
  extract_lrt(
    comparative_temporal$lrt$latitude,
    "Species × Age × Latitude",
    "Latitude"
  ),
  extract_lrt(
    comparative_temporal$lrt$longitude,
    "Species × Age × Longitude",
    "Longitude"
  )
) %>%
  transmute(
    Test,
    Axis,
    `Chi-square` = round(`Chi-square`, 3),
    df,
    p = fmt_p(p)
  )

caption_s6 <- paste0(
  "Table S6. Likelihood-ratio tests from the temporal-confounding ",
  "robustness analysis. Geographic coordinates were centered within ",
  "species and year, and species-, age-, and year-dependent temporal ",
  "terms were included as nuisance effects so geographic gradients ",
  "were estimated primarily from within-year spatial variation."
)

ft_s6 <- nice_ft(table_s6)


# ============================================================
# TABLE S7. Temporal-confounding slopes and contrasts
# ============================================================

table_s7 <- bind_rows(
  standardize_emm(
    comparative_temporal$latitude_slopes,
    "Latitude",
    "Age-specific slope"
  ),
  standardize_emm(
    comparative_temporal$longitude_slopes,
    "Longitude",
    "Age-specific slope"
  ),
  standardize_emm(
    comparative_temporal$latitude_age_differences,
    "Latitude",
    "Immature – Adult"
  ),
  standardize_emm(
    comparative_temporal$longitude_age_differences,
    "Longitude",
    "Immature – Adult"
  ),
  standardize_emm(
    comparative_temporal$latitude_species_contrasts,
    "Latitude",
    "Species contrast"
  ),
  standardize_emm(
    comparative_temporal$longitude_species_contrasts,
    "Longitude",
    "Species contrast"
  )
) %>%
  mutate(
    Contrast = ifelse(
      Result == "Species contrast",
      clean_species_contrast(Contrast),
      Contrast
    ),
    Contrast = gsub(
      " - ",
      " – ",
      Contrast,
      fixed = TRUE
    )
  ) %>%
  transmute(
    Result,
    Axis,
    Species,
    Age,
    Contrast,
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    `95% CI` = sprintf("%.3f to %.3f", CI_low, CI_high),
    z = round(z, 3),
    p = fmt_p(p)
  )

caption_s7 <- paste0(
  "Table S7. Age-specific geographic slopes, within-species ",
  "immature–adult slope differences, and pairwise interspecific ",
  "contrasts from the temporal-confounding robustness analysis. ",
  "Pairwise interspecific p-values are Tukey-adjusted; age-specific ",
  "slopes and prespecified within-species contrasts are not ",
  "multiplicity-adjusted."
)

ft_s7 <- nice_ft(
  table_s7,
  font_size = 8
)


# ============================================================
# TABLE S8. Leave-one-high-volume-location-out sensitivity
# ============================================================

location_lookup <- as.data.frame(
  station_influence_top_sites
) %>%
  transmute(
    Omitted_station = as.character(Station),
    `Omitted location` =
      format_location(
        Mean_latitude,
        Mean_longitude
      ),
    `N omitted` = as.integer(N),
    `% primary omitted` = Percent_of_primary
  )


loo_lrt <- as.data.frame(
  station_influence_lrt_table
) %>%
  mutate(
    Omitted_station = as.character(Omitted_station),
    Axis = case_when(
      grepl("geography", Test, ignore.case = TRUE) ~ "Joint",
      grepl("latitude", Test, ignore.case = TRUE) ~ "Latitude",
      grepl("longitude", Test, ignore.case = TRUE) ~ "Longitude",
      TRUE ~ ""
    )
  ) %>%
  transmute(
    Omitted_station,
    Result = "Omnibus LRT",
    Axis,
    Detail = gsub(
      " x ",
      " × ",
      Test,
      fixed = TRUE
    ),
    Estimate = NA_real_,
    SE = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    Statistic = Chi_square,
    df = as.numeric(df),
    p
  )


loo_age_raw <- as.data.frame(
  station_influence_age_differences
)

loo_age_ci <- get_ci_cols(
  loo_age_raw
)

loo_age <- loo_age_raw %>%
  transmute(
    Omitted_station = as.character(Omitted_station),
    Result = "Within-species age difference",
    Axis,
    Detail = paste0(Species, ": Immature – Adult"),
    Estimate = estimate,
    SE,
    CI_low = .data[[loo_age_ci[1]]],
    CI_high = .data[[loo_age_ci[2]]],
    Statistic = z.ratio,
    df = NA_real_,
    p = p.value
  )


loo_species_raw <- as.data.frame(
  station_influence_species_contrasts
)

loo_species_ci <- get_ci_cols(
  loo_species_raw
)

loo_species <- loo_species_raw %>%
  transmute(
    Omitted_station = as.character(Omitted_station),
    Result = "Species contrast",
    Axis,
    Detail = clean_species_contrast(contrast),
    Estimate = estimate,
    SE,
    CI_low = .data[[loo_species_ci[1]]],
    CI_high = .data[[loo_species_ci[2]]],
    Statistic = z.ratio,
    df = NA_real_,
    p = p.value
  )


table_s8 <- bind_rows(
  loo_lrt,
  loo_age,
  loo_species
) %>%
  left_join(
    location_lookup,
    by = "Omitted_station"
  ) %>%
  select(
    `Omitted location`,
    `N omitted`,
    `% primary omitted`,
    Result,
    Axis,
    Detail,
    Estimate,
    SE,
    CI_low,
    CI_high,
    Statistic,
    df,
    p
  ) %>%
  mutate(
    `% primary omitted` =
      round(`% primary omitted`, 2),
    Estimate = fmt_num(Estimate, 3),
    SE = fmt_num(SE, 3),
    `95% CI` = ifelse(
      is.na(CI_low) | is.na(CI_high),
      "",
      sprintf("%.3f to %.3f", CI_low, CI_high)
    ),
    Statistic = fmt_num(Statistic, 3),
    df = ifelse(
      is.na(df),
      "",
      as.character(as.integer(df))
    ),
    p = fmt_p(p)
  ) %>%
  select(-CI_low, -CI_high)

caption_s8 <- paste0(
  "Table S8. Leave-one-high-volume-reported-location-out sensitivity ",
  "analysis. The four reported coordinate locations contributing the ",
  "largest numbers of complete-case records to the primary comparative ",
  "dataset were omitted one at a time, and the primary comparative ",
  "model was refit without modification. Omnibus rows report likelihood-",
  "ratio tests; within-species rows report immature–adult slope differences; ",
  "and species-contrast rows report Tukey-adjusted pairwise comparisons ",
  "of age-dependent slopes."
)

ft_s8 <- nice_ft(
  table_s8,
  font_size = 7.5
)


# ============================================================
# TABLE S9. GAMM nonlinear-spatial sensitivity
# ============================================================

table_s9_raw <- gamm_sensitivity_summary %>%
  left_join(
    gamm_kcheck_summary,
    by = c("Species", "Term")
  ) %>%
  left_join(
    gamm_model_metadata %>%
      select(
        Species,
        N,
        Model_rank
      ),
    by = "Species"
  )
# Confirm complete GAMM, basis-dimension, and model-metadata output.
stopifnot(
  nrow(table_s9_raw) == 6L,
  !anyNA(
    table_s9_raw[
      c(
        "Species",
        "Term",
        "Smooth",
        "edf",
        "Ref_df",
        "Chi_square",
        "p",
        "k_index",
        "k_p",
        "N",
        "Model_rank"
      )
    ]
  )
)

table_s9 <- table_s9_raw %>%
  transmute(
    Species,
    N,
    Smooth = dplyr::recode(
      Smooth,
      "Immature - Adult spatial difference" =
        "Immature – Adult spatial difference",
      .default = Smooth
    ),
    edf = round(edf, 3),
    `Reference df` = round(Ref_df, 3),
    `Chi-square` = round(Chi_square, 3),
    p = fmt_p(p),
    `k-index` = ifelse(
      is.na(k_index),
      "",
      sprintf("%.2f", k_index)
    ),
    `k-check p` = fmt_p(k_p),
    `Model rank` = Model_rank
  )

# Confirm the expected GAMM table structure.
stopifnot(
  nrow(table_s9) == 6L,
  !anyNA(table_s9$Species)
)

caption_s9 <- paste0(
  "Table S9. Species-specific generalized additive mixed-model ",
  "(GAMM) sensitivity analysis allowing nonlinear two-dimensional ",
  "geographic sex-ratio surfaces. Spatial smooths were fitted using ",
  "EPSG:5070 projected coordinates scaled to 100-km units. The adult ",
  "spatial surface is the baseline smooth, whereas the Immature – Adult ",
  "spatial difference tests whether the geographic surface differs ",
  "between age classes. Models included random-effect smooths for ",
  "reported location, year, and band identity. k-index values and ",
  "associated p-values are from basis-dimension checks calculated once ",
  "in the GAMM sensitivity analysis and stored for reuse."
)

ft_s9 <- nice_ft(
  table_s9,
  font_size = 8
)


# ============================================================
# TABLE S10. Extreme-allocation demographic bounds
# ============================================================

age_bounds_table <- map_dfr(
  descriptive_species,
  function(sp) {
    as.data.frame(
      sensitivity_results[[sp]]$age_bounds
    ) %>%
      filter(Group == "All") %>%
      mutate(Species = sp)
  }
)

sex_bounds_table <- map_dfr(
  descriptive_species,
  function(sp) {
    as.data.frame(
      sensitivity_results[[sp]]$sex_ratio_bounds
    ) %>%
      mutate(Species = sp)
  }
)

table_s10 <- sex_bounds_table %>%
  select(
    Species,
    Age,
    Observed_female_prop,
    Minimum_female_prop,
    Maximum_female_prop
  ) %>%
  left_join(
    age_bounds_table %>%
      select(
        Species,
        Observed_immature_prop,
        Minimum_immature_prop,
        Maximum_immature_prop
      ),
    by = "Species"
  ) %>%
  transmute(
    Species,
    Age,
    `Observed female (%)` =
      round(100 * Observed_female_prop, 1),
    `Minimum female (%)` =
      round(100 * Minimum_female_prop, 1),
    `Maximum female (%)` =
      round(100 * Maximum_female_prop, 1),
    `Observed immature (%)` =
      round(100 * Observed_immature_prop, 1),
    `Minimum immature (%)` =
      round(100 * Minimum_immature_prop, 1),
    `Maximum immature (%)` =
      round(100 * Maximum_immature_prop, 1)
  )

caption_s10 <- paste0(
  "Table S10. Extreme-allocation bounds for age and sex composition ",
  "of Calliope, Broad-tailed, and Allen’s Hummingbirds. Age-composition ",
  "bounds assign all records with unclassifiable age to adults or ",
  "immatures in turn. Age-specific sex-ratio bounds jointly account ",
  "for unclassifiable age and unknown or missing sex by assigning ",
  "uncertain records to the combinations producing the minimum and ",
  "maximum possible female proportions."
)

ft_s10 <- nice_ft(table_s10)


# ============================================================
# TABLE S11. Missing demographic information
# ============================================================

table_s11 <- as.data.frame(
  missingness_summary_all
) %>%
  filter(
    Species %in% retained_species
  ) %>%
  mutate(
    Species = factor(
      Species,
      levels = retained_species
    )
  ) %>%
  arrange(Species) %>%
  transmute(
    Species = as.character(Species),
    `Selected bird-years` = Total_selected_bird_years,
    `Unknown sex N` = Unknown_or_missing_sex_N,
    `Unknown sex (%)` =
      round(Unknown_or_missing_sex_pct, 1),
    `Unclassifiable age N` = Unclassifiable_age_N,
    `Unclassifiable age (%)` =
      round(Unclassifiable_age_pct, 1)
  )

caption_s11 <- paste0(
  "Table S11. Missing sex and unclassifiable-age information among ",
  "the six hummingbird species retained in the study within the ",
  "15 November–31 December primary sampling window. Counts are based ",
  "on the one-record-per-individual-per-year selection used for the ",
  "descriptive uncertainty analysis."
)

ft_s11 <- nice_ft(table_s11)


# ============================================================
# TABLE S12. Primary comparative model diagnostics
# ============================================================

diag_global <- as.data.frame(
  comparative_primary_diagnostics$global
)

table_s12 <- tibble(
  Diagnostic = c(
    "Convergence warning",
    "Optimizer code",
    "Maximum absolute gradient",
    "Singular model",
    "DHARMa dispersion",
    "DHARMa residual uniformity KS",
    "DHARMa outlier test"
  ),
  Statistic = c(
    ifelse(
      is.na(diag_global$Convergence_message[1]) |
        diag_global$Convergence_message[1] == "",
      "None",
      as.character(diag_global$Convergence_message[1])
    ),
    as.character(diag_global$Optimizer_code[1]),
    sprintf("%.6f", diag_global$Max_abs_gradient[1]),
    as.character(diag_global$Singular[1]),
    sprintf("%.3f", diag_global$DHARMa_dispersion[1]),
    sprintf("%.3f", diag_global$DHARMa_uniformity_KS[1]),
    ""
  ),
  p = c(
    "",
    "",
    "",
    "",
    fmt_p(diag_global$DHARMa_dispersion_p[1]),
    fmt_p(diag_global$DHARMa_uniformity_p[1]),
    fmt_p(diag_global$DHARMa_outlier_p[1])
  )
)

caption_s12 <- paste0(
  "Table S12. Diagnostic evaluation of the primary comparative ",
  "generalized linear mixed model. Rows summarize convergence, ",
  "singularity, and simulation-based DHARMa diagnostics for ",
  "dispersion, residual uniformity, and outliers."
)

ft_s12 <- nice_ft(table_s12)


# ============================================================
# TABLE S13. Geographic sampling support by latitude band
# ============================================================

required_support_cols <- c(
  "Species",
  "Immature_f",
  "lat_dd",
  "species_site_id"
)

stopifnot(
  all(required_support_cols %in% names(comparative_df)),
  !anyNA(comparative_df$Species),
  !anyNA(comparative_df$Immature_f),
  !anyNA(comparative_df$lat_dd),
  !anyNA(comparative_df$species_site_id)
)

latitude_support_data <- comparative_df %>%
  mutate(
    Species = as.character(Species),
    Age = as.character(Immature_f),
    Latitude_band_low =
      floor(lat_dd / lat_band_width) * lat_band_width,
    Latitude_band = paste0(
      Latitude_band_low,
      "–<",
      Latitude_band_low + lat_band_width,
      "°N"
    )
  )

all_band_lows <- seq(
  min(latitude_support_data$Latitude_band_low),
  max(latitude_support_data$Latitude_band_low),
  by = lat_band_width
)

latitude_support_long <- latitude_support_data %>%
  group_by(
    Species,
    Age,
    Latitude_band_low,
    Latitude_band
  ) %>%
  summarise(
    Records = n(),
    Reported_locations =
      n_distinct(species_site_id),
    .groups = "drop"
  ) %>%
  complete(
    Species = focal_species_order,
    Age = c("Adult", "Immature"),
    Latitude_band_low = all_band_lows,
    fill = list(
      Records = 0L,
      Reported_locations = 0L
    )
  ) %>%
  mutate(
    Latitude_band = paste0(
      Latitude_band_low,
      "–<",
      Latitude_band_low + lat_band_width,
      "°N"
    ),
    Species = factor(
      Species,
      levels = focal_species_order
    ),
    Age = factor(
      Age,
      levels = c("Adult", "Immature")
    )
  ) %>%
  arrange(
    Species,
    Latitude_band_low,
    Age
  )

table_s13 <- latitude_support_long %>%
  select(
    Species,
    Latitude_band_low,
    Latitude_band,
    Age,
    Records,
    Reported_locations
  ) %>%
  pivot_wider(
    names_from = Age,
    values_from = c(
      Records,
      Reported_locations
    ),
    names_glue = "{Age}_{.value}"
  ) %>%
  transmute(
    Species = as.character(Species),
    Latitude_band_low,
    `Latitude band` = Latitude_band,
    `Adult records` = Adult_Records,
    `Adult reported locations` =
      Adult_Reported_locations,
    `Immature records` = Immature_Records,
    `Immature reported locations` =
      Immature_Reported_locations
  ) %>%
  arrange(
    factor(
      Species,
      levels = focal_species_order
    ),
    Latitude_band_low
  )

# Confirm that latitude-band counts reproduce the primary-model sample size.
n_from_bands <- latitude_support_long %>%
  group_by(Species, Age) %>%
  summarise(
    N = sum(Records),
    .groups = "drop"
  ) %>%
  mutate(
    Species = as.character(Species),
    Age = as.character(Age)
  ) %>%
  arrange(Species, Age)

n_direct <- comparative_df %>%
  transmute(
    Species = as.character(Species),
    Age = as.character(Immature_f)
  ) %>%
  count(
    Species,
    Age,
    name = "N"
  ) %>%
  arrange(Species, Age)

stopifnot(
  isTRUE(
    all.equal(
      n_from_bands,
      n_direct,
      check.attributes = FALSE
    )
  ),
  sum(latitude_support_long$Records) ==
    nrow(comparative_df)
)

table_s13 <- table_s13 %>%
  select(-Latitude_band_low)

caption_s13 <- paste0(
  "Table S13. Geographic sampling support for adult and immature ",
  "Rufous, Black-chinned, and Ruby-throated Hummingbirds in the ",
  "primary comparative analysis. Records are grouped into 2° latitude ",
  "bands. Reported locations are unique species-specific combinations ",
  "of identical reported latitude and longitude. Because coordinate ",
  "precision varies among banding records, reported locations do not ",
  "necessarily correspond to individual physical banding sites."
)

ft_s13 <- nice_ft(
  table_s13,
  font_size = 8
)


# ============================================================
# 3. WRITE SUPPLEMENTARY TABLES TO WORD
#
# The complete supplement is written in landscape orientation.
# Each table begins on a new page, and every caption is inserted as a
# standard Word paragraph immediately before its table. This provides
# consistent caption formatting and stable pagination across tables.
# ============================================================

supplement_tables <- list(
  ft_s1,
  ft_s2,
  ft_s3,
  ft_s4,
  ft_s5,
  ft_s6,
  ft_s7,
  ft_s8,
  ft_s9,
  ft_s10,
  ft_s11,
  ft_s12,
  ft_s13
)

doc <- read_docx()

doc_dims <- docx_dim(
  doc
)

page_short <- min(
  doc_dims$page
)

page_long <- max(
  doc_dims$page
)

section_margins <- page_mar(
  top = unname(doc_dims$margins["top"]),
  bottom = unname(doc_dims$margins["bottom"]),
  left = unname(doc_dims$margins["left"]),
  right = unname(doc_dims$margins["right"]),
  header = unname(doc_dims$margins["header"]),
  footer = unname(doc_dims$margins["footer"])
)

# Define a single landscape section for the complete supplement.
landscape_properties <- prop_section(
  page_size = page_size(
    width = page_short,
    height = page_long,
    orient = "landscape"
  ),
  page_margins = section_margins,
  type = "continuous"
)

landscape_available_width <-
  page_long -
  unname(doc_dims$margins["left"]) -
  unname(doc_dims$margins["right"]) -
  0.10

# Allow long tables to continue across pages and repeat column headers.
supplement_tables <- lapply(
  supplement_tables,
  function(ft) {
    set_table_properties(
      ft,
      opts_word = list(
        split = TRUE,
        keep_with_next = FALSE,
        repeat_headers = TRUE
      )
    )
  }
)

# Scale the widest tables to the usable landscape-page width.
wide_table_numbers <- c(
  7L,   # temporal-confounding slopes and contrasts
  8L,   # leave-one-location-out sensitivity
  9L,   # GAMM sensitivity
  13L   # geographic sampling support
)

supplement_tables[wide_table_numbers] <- lapply(
  supplement_tables[wide_table_numbers],
  fit_to_width,
  max_width = landscape_available_width
)

caption_numbers <- seq_along(
  supplement_tables
)

captions <- list(
  `1` = caption_s1,
  `2` = caption_s2,
  `3` = caption_s3,
  `4` = caption_s4,
  `5` = caption_s5,
  `6` = caption_s6,
  `7` = caption_s7,
  `8` = caption_s8,
  `9` = caption_s9,
  `10` = caption_s10,
  `11` = caption_s11,
  `12` = caption_s12,
  `13` = caption_s13
)

caption_paragraph_properties <- caption_fp

for (i in seq_along(supplement_tables)) {
  
  # Table S1 begins on page 1; each subsequent table starts on a new page.
  if (i > 1L) {
    doc <- body_add_break(
      doc
    )
  }
  
  # Add a consistently formatted caption immediately before the table.
  doc <- body_add_fpar(
    doc,
    fpar(
      ftext(
        captions[[as.character(i)]],
        prop = caption_text_fp
      ),
      fp_p = caption_paragraph_properties
    )
  )
  
  doc <- body_add_flextable(
    doc,
    supplement_tables[[i]]
  )
}

# Apply the landscape section to the document.
doc <- body_set_default_section(
  doc,
  landscape_properties
)

print(
  doc,
  target = "Supplementary_Tables.docx"
)
