# ============================================================
# 03_sensitivity_analysis.R
# Run the sensitivity analyses, model diagnostics, and nonlinear spatial GAMMs
# reported in the manuscript and Supplementary Material.
# Run 02_primary_analysis.R first.
#
# The GAMM section uses a saved results cache when a matching cache is present.
# If the cache is absent, the GAMMs and basis-dimension checks are fit and the
# results are saved for future runs.
# ============================================================


library(data.table)
library(dplyr)
library(lme4)
library(emmeans)
library(DHARMa)
library(mgcv)
library(sf)


# ============================================================
# 0. CHECK REQUIRED INPUTS
# ============================================================

required_primary_objects <- c(
  "results",
  "comparative_primary",
  "comparative_df",
  "m_comparative_full",
  "focal_species",
  "species_ids",
  "csv_file",
  "months_keep",
  "primary_start_month",
  "primary_start_day",
  "primary_end_month",
  "primary_end_day",
  "countries_keep",
  "lon_min",
  "lon_max",
  "lat_min",
  "lat_max"
)

missing_primary_objects <- required_primary_objects[
  !vapply(
    required_primary_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_primary_objects) > 0L) {
  stop(
    "Run 01_functions.R and 02_primary_analysis.R first. Missing object(s): ",
    paste(missing_primary_objects, collapse = ", ")
  )
}

required_functions <- c(
  "in_month_day_window",
  "infer_season_year_fun",
  "classify_immature",
  "prepare_comparative_data",
  "fit_comparative_age_sex_model",
  "fit_temporal_comparative_robustness",
  "run_leave_one_station_out"
)

missing_functions <- required_functions[
  !vapply(
    required_functions,
    exists,
    logical(1),
    mode = "function"
  )
]

if (length(missing_functions) > 0L) {
  stop(
    "The required functions from 01_functions.R are not loaded. Missing: ",
    paste(missing_functions, collapse = ", ")
  )
}

required_comparative_cols <- c(
  "common_station_id",
  "species_site_id",
  "species_winter_year",
  "species_band_id",
  "band",
  "lat_dd",
  "lon_dd",
  "lat5",
  "lon5",
  "Female",
  "Immature",
  "Immature_f",
  "winter_year"
)

missing_comparative_cols <- setdiff(
  required_comparative_cols,
  names(comparative_df)
)

if (length(missing_comparative_cols) > 0L) {
  stop(
    "`comparative_df` is missing required column(s): ",
    paste(missing_comparative_cols, collapse = ", ")
  )
}

if (!requireNamespace("ape", quietly = TRUE)) {
  stop(
    "Package 'ape' is required for species-specific Moran's I diagnostics."
  )
}


# ============================================================
# 1. SENSITIVITY SETTINGS AND HELPERS
# ============================================================

station_leave_out_n <- 4L
restricted_lon_min <- -97
projection_epsg <- 5070L

dharma_simulations <- 1000L
diagnostics_seed <- 123L


# Return the first available column name from a set of candidates.
pick_col <- function(df, candidates, required = TRUE, label = NULL) {
  hit <- intersect(candidates, names(df))
  
  if (length(hit) == 0L) {
    if (required) {
      stop(
        "Could not find column ",
        if (is.null(label)) paste(candidates, collapse = "/") else label,
        "."
      )
    }
    return(NULL)
  }
  
  hit[1L]
}


# Standardize emmeans output used by sensitivity analyses.
tidy_emm <- function(x) {
  d <- as.data.frame(x)
  
  est_col <- pick_col(
    d,
    c(
      "estimate",
      "lat5.trend",
      "lon5.trend",
      "lat5_within.trend",
      "lon5_within.trend"
    ),
    label = "estimate"
  )
  
  lcl_col <- pick_col(
    d,
    c("asymp.LCL", "lower.CL", "LCL"),
    label = "lower confidence limit"
  )
  
  ucl_col <- pick_col(
    d,
    c("asymp.UCL", "upper.CL", "UCL"),
    label = "upper confidence limit"
  )
  
  stat_col <- pick_col(
    d,
    c("z.ratio", "t.ratio"),
    required = FALSE
  )
  
  out <- data.frame(
    Estimate = as.numeric(d[[est_col]]),
    SE = as.numeric(d[["SE"]]),
    CI_low = as.numeric(d[[lcl_col]]),
    CI_high = as.numeric(d[[ucl_col]]),
    Statistic =
      if (is.null(stat_col)) NA_real_ else as.numeric(d[[stat_col]]),
    p = as.numeric(d[["p.value"]]),
    stringsAsFactors = FALSE
  )
  
  for (nm in intersect(c("contrast", "Species", "Immature_f"), names(d))) {
    out[[nm]] <- as.character(d[[nm]])
  }
  
  out
}


# Extract immature-adult geographic slope differences from a fitted model.
age_slope_differences <- function(
    model,
    lat_var = "lat5",
    lon_var = "lon5"
) {
  one_axis <- function(var_name, axis_label) {
    em <- emmeans::emtrends(
      model,
      specs = ~ Immature_f | Species,
      var = var_name
    )
    
    d <- tidy_emm(
      summary(
        emmeans::contrast(
          em,
          method = "revpairwise"
        ),
        infer = c(TRUE, TRUE),
        level = 0.95
      )
    )
    
    d$Axis <- axis_label
    d
  }
  
  dplyr::bind_rows(
    one_axis(lat_var, "Latitude"),
    one_axis(lon_var, "Longitude")
  ) %>%
    dplyr::select(
      dplyr::any_of(
        c(
          "Axis",
          "Species",
          "contrast"
        )
      ),
      Estimate,
      SE,
      CI_low,
      CI_high,
      Statistic,
      p
    )
}


# Extract chi-square, degrees of freedom, and raw probability value from an LRT.
extract_lrt <- function(x, test_name) {
  a <- as.data.frame(x)
  last <- a[nrow(a), , drop = FALSE]
  
  chisq_col <- intersect(
    c("Chisq", "ChiSq", "Chi_square"),
    names(last)
  )[1L]
  
  df_col <- intersect(
    c("Df", "Chi Df", "df"),
    names(last)
  )[1L]
  
  p_col <- intersect(
    c("Pr(>Chisq)", "p", "p.value"),
    names(last)
  )[1L]
  
  if (any(is.na(c(chisq_col, df_col, p_col)))) {
    stop("Could not identify all required likelihood-ratio test columns.")
  }
  
  data.frame(
    Test = test_name,
    Chi_square = as.numeric(last[[chisq_col]]),
    df = as.numeric(last[[df_col]]),
    p = as.numeric(last[[p_col]]),
    stringsAsFactors = FALSE
  )
}


# Retain the earliest record with known sex and classifiable age for each
# species-individual-winter combination. If no complete record exists, retain
# the earliest record so demographic uncertainty can be summarized.
select_first_record_per_bird_winter <- function(d, complete_only = FALSE) {
  d <- data.table::copy(d)
  
  data.table::setorder(
    d,
    Species,
    band,
    winter_year,
    dt
  )
  
  selected <- d[, {
    known_both <- which(
      !is.na(Female) &
        !is.na(Immature)
    )
    
    if (length(known_both) > 0L) {
      .SD[known_both[1L]]
    } else {
      .SD[1L]
    }
  },
  by = .(
    Species,
    band,
    winter_year
  )]
  
  if (isTRUE(complete_only)) {
    selected <- selected[
      !is.na(Female) &
        !is.na(Immature)
    ]
  }
  
  selected
}


# Convert selected focal-species records into the comparative-data structure
# used by the primary GLMM.
records_to_comparative <- function(selected_records) {
  selected_records <- as.data.frame(selected_records)
  
  sensitivity_results <- stats::setNames(
    lapply(
      focal_species,
      function(sp) {
        d <- selected_records[
          selected_records$Species == sp,
          ,
          drop = FALSE
        ]
        
        if (nrow(d) == 0L) {
          stop(
            "No qualifying records for ",
            sp,
            " in the sensitivity dataset."
          )
        }
        
        list(data_points = d)
      }
    ),
    focal_species
  )
  
  prepare_comparative_data(
    results = sensitivity_results,
    comparative_species = focal_species,
    reference_species = "Rufous"
  )
}


# Fit the primary comparative structure to an independently selected set of
# focal-species records and return the quantities used in sensitivity tables.
fit_rebuilt_sensitivity <- function(selected_records, label) {
  d <- records_to_comparative(selected_records)
  
  fit <- fit_comparative_age_sex_model(
    comparative_df = d,
    optimizer = "bobyqa",
    maxfun = 2e5,
    run_lrt = TRUE,
    pairwise_adjust = "tukey"
  )
  
  lrt <- dplyr::bind_rows(
    extract_lrt(
      fit$lrt$overall,
      "Species x age x geography"
    ),
    extract_lrt(
      fit$lrt$latitude,
      "Species x age x latitude"
    ),
    extract_lrt(
      fit$lrt$longitude,
      "Species x age x longitude"
    )
  )
  
  lrt$Analysis <- label
  lrt$N <- nrow(d)
  
  differences <- age_slope_differences(
    fit$model
  )
  differences$Analysis <- label
  differences$N <- nrow(d)
  
  list(
    data = d,
    model = fit$model,
    lrt = lrt,
    differences = differences,
    convergence = fit$convergence
  )
}


# ============================================================
# 2. PRIMARY-WINDOW RECONSTRUCTION
# ============================================================

# Reconstruct the primary-window selection while retaining incomplete
# demographic records for the uncertainty analysis and exact event dates for
# the within-window sensitivity analysis.
raw_primary <- data.table::fread(
  csv_file,
  select = c(
    "band",
    "species_id",
    "event_year",
    "event_month",
    "event_day",
    "iso_country",
    "lat_dd",
    "lon_dd",
    "sex_code",
    "age_code"
  ),
  showProgress = FALSE
)

raw_primary[, event_year := as.integer(event_year)]
raw_primary[, event_month := as.integer(event_month)]
raw_primary[, event_day := as.integer(event_day)]

primary_date_keep <- in_month_day_window(
  month = raw_primary$event_month,
  day = raw_primary$event_day,
  start_month = primary_start_month,
  start_day = primary_start_day,
  end_month = primary_end_month,
  end_day = primary_end_day
)

raw_primary <- raw_primary[
  species_id %in% unname(species_ids) &
    iso_country %in% countries_keep &
    event_month %in% months_keep &
    primary_date_keep &
    is.finite(lat_dd) &
    is.finite(lon_dd) &
    lon_dd >= lon_min &
    lon_dd <= lon_max &
    lat_dd >= lat_min &
    lat_dd <= lat_max &
    !is.na(band)
]

raw_primary[, Species := names(species_ids)[
  match(
    species_id,
    unname(species_ids)
  )
]]

# Preserve missing event days before assigning a value used only for ordering.
raw_primary[, event_day_missing := is.na(event_day)]
raw_primary[, event_day_filled := event_day]
raw_primary[is.na(event_day_filled), event_day_filled := 15L]

raw_primary[, dt := data.table::as.IDate(
  sprintf(
    "%04d-%02d-%02d",
    event_year,
    event_month,
    event_day_filled
  )
)]

raw_primary <- raw_primary[
  !is.na(dt)
]

primary_year_fun <- infer_season_year_fun(
  months_keep
)

raw_primary[, winter_year :=
              primary_year_fun(
                event_year,
                event_month
              )]

raw_primary[, Female := data.table::fifelse(
  sex_code %in% c(5L, 7L),
  1L,
  data.table::fifelse(
    sex_code %in% c(4L, 6L),
    0L,
    NA_integer_
  )
)]

raw_primary[, Immature :=
              classify_immature(
                event_month,
                age_code
              )]

selected_primary_all <- select_first_record_per_bird_winter(
  raw_primary,
  complete_only = FALSE
)

# Verify that the independently reconstructed complete focal-species records
# exactly match the datasets used in the primary model.
primary_reconstruction_matches <- vapply(
  focal_species,
  function(sp) {
    recon <- selected_primary_all[
      Species == sp &
        !is.na(Female) &
        !is.na(Immature)
    ][
      order(
        band,
        winter_year
      )
    ]
    
    primary <- data.table::as.data.table(
      results[[sp]]$data_points
    )[
      order(
        band,
        winter_year
      )
    ]
    
    if (nrow(recon) != nrow(primary)) {
      return(FALSE)
    }
    
    isTRUE(
      all.equal(
        recon[, .(
          band = as.character(band),
          winter_year = as.integer(winter_year),
          lat_dd = as.numeric(lat_dd),
          lon_dd = as.numeric(lon_dd),
          Female = as.integer(Female),
          Immature = as.integer(Immature)
        )],
        primary[, .(
          band = as.character(band),
          winter_year = as.integer(winter_year),
          lat_dd = as.numeric(lat_dd),
          lon_dd = as.numeric(lon_dd),
          Female = as.integer(Female),
          Immature = as.integer(Immature)
        )],
        tolerance = 1e-8,
        check.attributes = FALSE
      )
    )
  },
  logical(1)
)

if (!all(primary_reconstruction_matches)) {
  stop(
    "Primary-window reconstruction does not exactly match the primary ",
    "selected records for at least one focal species."
  )
}


# ============================================================
# 3. DEMOGRAPHIC UNCERTAINTY
# ============================================================

# Calculate extreme-allocation bounds for age composition and age-specific
# sex ratios among records with missing or unclassifiable demographic data.
safe_prop <- function(num, den) {
  if (den > 0L) num / den else NA_real_
}


calc_species_uncertainty <- function(sp) {
  x <- selected_primary_all[
    Species == sp
  ]
  
  n_selected <- nrow(x)
  n_unknown_sex <- sum(is.na(x$Female))
  n_unclassifiable_age <- sum(is.na(x$Immature))
  
  missingness_summary <- data.table::data.table(
    Species = sp,
    Total_selected_bird_years = n_selected,
    Unknown_or_missing_sex_N = n_unknown_sex,
    Unknown_or_missing_sex_pct =
      if (n_selected > 0L) {
        100 * n_unknown_sex / n_selected
      } else {
        NA_real_
      },
    Unclassifiable_age_N = n_unclassifiable_age,
    Unclassifiable_age_pct =
      if (n_selected > 0L) {
        100 * n_unclassifiable_age / n_selected
      } else {
        NA_real_
      }
  )
  
  n_immature <- sum(x$Immature == 1L, na.rm = TRUE)
  n_adult <- sum(x$Immature == 0L, na.rm = TRUE)
  n_age_unknown <- sum(is.na(x$Immature))
  n_total <- n_immature + n_adult + n_age_unknown
  
  age_bounds <- data.table::data.table(
    Group = "All",
    Immature_known = n_immature,
    Adult_known = n_adult,
    Unclassifiable_age = n_age_unknown,
    Total_records = n_total,
    Observed_immature_prop =
      safe_prop(
        n_immature,
        n_immature + n_adult
      ),
    Minimum_immature_prop =
      safe_prop(
        n_immature,
        n_total
      ),
    Maximum_immature_prop =
      safe_prop(
        n_immature + n_age_unknown,
        n_total
      )
  )
  
  # Joint age-sex uncertainty cells.
  AF <- nrow(x[Immature == 0L & Female == 1L])
  AM <- nrow(x[Immature == 0L & Female == 0L])
  IF <- nrow(x[Immature == 1L & Female == 1L])
  IM <- nrow(x[Immature == 1L & Female == 0L])
  AU <- nrow(x[Immature == 0L & is.na(Female)])
  IU <- nrow(x[Immature == 1L & is.na(Female)])
  UF <- nrow(x[is.na(Immature) & Female == 1L])
  UM <- nrow(x[is.na(Immature) & Female == 0L])
  UU <- nrow(x[is.na(Immature) & is.na(Female)])
  
  sex_ratio_bounds <- data.table::data.table(
    Age = c(
      "Adult",
      "Immature"
    ),
    Observed_female_prop = c(
      safe_prop(AF, AF + AM),
      safe_prop(IF, IF + IM)
    ),
    Minimum_female_prop = c(
      safe_prop(AF, AF + AM + AU + UM + UU),
      safe_prop(IF, IF + IM + IU + UM + UU)
    ),
    Maximum_female_prop = c(
      safe_prop(AF + AU + UF + UU, AF + AM + AU + UF + UU),
      safe_prop(IF + IU + UF + UU, IF + IM + IU + UF + UU)
    )
  )
  
  list(
    missingness_summary = missingness_summary,
    age_bounds = age_bounds,
    sex_ratio_bounds = sex_ratio_bounds
  )
}


sensitivity_results <- stats::setNames(
  lapply(
    names(species_ids),
    calc_species_uncertainty
  ),
  names(species_ids)
)

missingness_summary_all <- data.table::rbindlist(
  lapply(
    sensitivity_results,
    function(x) x$missingness_summary
  ),
  use.names = TRUE,
  fill = TRUE
)


# ============================================================
# 4. JANUARY AGE-CLASSIFICATION LOSS
# ============================================================

# Apply the focal-species filters and one-record-per-individual selection to
# January records and summarize the resulting age-classification loss.
raw_january <- data.table::fread(
  csv_file,
  select = c(
    "band",
    "species_id",
    "event_year",
    "event_month",
    "event_day",
    "iso_country",
    "lat_dd",
    "lon_dd",
    "sex_code",
    "age_code"
  ),
  showProgress = FALSE
)

raw_january[, event_year := as.integer(event_year)]
raw_january[, event_month := as.integer(event_month)]
raw_january[, event_day := as.integer(event_day)]

raw_january <- raw_january[
  species_id %in% unname(species_ids[focal_species]) &
    iso_country %in% countries_keep &
    event_month == 1L &
    is.finite(lat_dd) &
    is.finite(lon_dd) &
    lon_dd >= lon_min &
    lon_dd <= lon_max &
    lat_dd >= lat_min &
    lat_dd <= lat_max &
    !is.na(band)
]

raw_january[, Species := names(species_ids)[
  match(
    species_id,
    unname(species_ids)
  )
]]

raw_january[, Female := data.table::fifelse(
  sex_code %in% c(5L, 7L),
  1L,
  data.table::fifelse(
    sex_code %in% c(4L, 6L),
    0L,
    NA_integer_
  )
)]

raw_january[, Immature :=
              classify_immature(
                event_month,
                age_code
              )]

# Missing event days are used only for ordering within January.
raw_january[, event_day_filled := event_day]
raw_january[is.na(event_day_filled), event_day_filled := 15L]

raw_january[, dt := data.table::as.IDate(
  sprintf(
    "%04d-%02d-%02d",
    event_year,
    event_month,
    event_day_filled
  )
)]

raw_january <- raw_january[
  !is.na(dt)
]

data.table::setorder(
  raw_january,
  Species,
  band,
  event_year,
  dt
)

selected_january <- raw_january[, {
  known_both <- which(
    !is.na(Female) &
      !is.na(Immature)
  )
  
  if (length(known_both) > 0L) {
    .SD[known_both[1L]]
  } else {
    .SD[1L]
  }
},
by = .(
  Species,
  band,
  event_year
)]

january_classification_loss_overall <- selected_january[
  ,
  .(
    N = .N,
    Unclassifiable_age_N = sum(
      is.na(Immature)
    ),
    Unclassifiable_age_pct =
      100 * mean(
        is.na(Immature)
      )
  )
]

january_unclassifiable_age_codes <- selected_january[
  is.na(Immature),
  .N,
  by = .(
    Species,
    age_code
  )
][
  order(
    Species,
    -N
  )
]


# ============================================================
# 5. TEMPORAL-CONFOUNDING ROBUSTNESS
# ============================================================

# Center geography within species and winter year and include standardized
# year terms to assess temporal confounding of geographic sex-ratio patterns.
comparative_temporal <- fit_temporal_comparative_robustness(
  comparative_df = comparative_df,
  optimizer = "bobyqa",
  maxfun = 2e5,
  run_lrt = TRUE,
  pairwise_adjust = "tukey"
)

# Quantify the geographic variance retained after within-year centering.
within_year_variance <- comparative_df %>%
  dplyr::group_by(
    Species,
    winter_year
  ) %>%
  dplyr::mutate(
    lat_centered = lat_dd - mean(lat_dd, na.rm = TRUE),
    lon_centered = lon_dd - mean(lon_dd, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Latitude_variance_retained_pct =
      100 * stats::var(lat_centered, na.rm = TRUE) /
      stats::var(lat_dd, na.rm = TRUE),
    Longitude_variance_retained_pct =
      100 * stats::var(lon_centered, na.rm = TRUE) /
      stats::var(lon_dd, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 6. LEAVE-ONE-HIGH-VOLUME-LOCATION-OUT SENSITIVITY
# ============================================================

# Identify the four highest-volume reported locations and refit the primary
# model after omitting each location separately.
station_influence <- run_leave_one_station_out(
  comparative_df = comparative_df,
  top_n = station_leave_out_n,
  station_col = "common_station_id",
  optimizer = "bobyqa",
  maxfun = 2e5,
  run_lrt = TRUE,
  pairwise_adjust = "tukey"
)

station_influence_top_sites <-
  station_influence$top_station_summary

station_influence_lrt_table <-
  station_influence$lrt_summary

station_influence_age_differences <-
  station_influence$age_difference_summary



# ============================================================
# 7. WITHIN-WINDOW DATE SENSITIVITY
# ============================================================

# Reconstruct the primary focal-species dataset with exact event dates retained.
primary_focal_records <- selected_primary_all[
  Species %in% focal_species &
    !is.na(Female) &
    !is.na(Immature)
]

primary_focal_with_dates <- records_to_comparative(
  primary_focal_records
)

# Confirm that adding event-date fields preserves the primary dataset.
primary_date_keys <- primary_focal_with_dates %>%
  dplyr::transmute(
    Species = as.character(Species),
    band = as.character(band),
    winter_year = as.integer(winter_year),
    lat_dd = as.numeric(lat_dd),
    lon_dd = as.numeric(lon_dd),
    Female = as.integer(Female),
    Immature = as.integer(Immature)
  ) %>%
  dplyr::arrange(
    Species,
    band,
    winter_year
  )

primary_model_keys <- comparative_df %>%
  dplyr::transmute(
    Species = as.character(Species),
    band = as.character(band),
    winter_year = as.integer(winter_year),
    lat_dd = as.numeric(lat_dd),
    lon_dd = as.numeric(lon_dd),
    Female = as.integer(Female),
    Immature = as.integer(Immature)
  ) %>%
  dplyr::arrange(
    Species,
    band,
    winter_year
  )

if (!isTRUE(
  all.equal(
    primary_date_keys,
    primary_model_keys,
    tolerance = 1e-8,
    check.attributes = FALSE
  )
)) {
  stop(
    "The date-augmented primary dataset does not match the primary comparative dataset."
  )
}

# Records without an exact event day are excluded only from this analysis.
date_covariate_df <- primary_focal_with_dates %>%
  dplyr::filter(
    !event_day_missing
  ) %>%
  droplevels()

window_start <- as.Date(
  sprintf(
    "%04d-%02d-%02d",
    2000L,
    primary_start_month,
    primary_start_day
  )
)

date_covariate_df$day_in_window <- as.integer(
  as.Date(
    sprintf(
      "%04d-%02d-%02d",
      2000L,
      date_covariate_df$event_month,
      date_covariate_df$event_day
    )
  ) - window_start
)

date_covariate_df$day_win_sc <- as.numeric(
  scale(
    date_covariate_df$day_in_window
  )
)

# Allow the within-window relationship with sex ratio to differ among species
# and age classes while retaining the primary geographic fixed effects and
# random-effect structure.
date_covariate_model <- lme4::glmer(
  Female ~
    Species * Immature_f * lat5 +
    Species * Immature_f * lon5 +
    Species * Immature_f * day_win_sc +
    (1 | species_site_id) +
    (1 | species_winter_year) +
    (1 | species_band_id),
  data = date_covariate_df,
  family = binomial,
  nAGQ = 1,
  control = lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 2e5
    )
  )
)

date_covariate_lrt <- dplyr::bind_rows(
  extract_lrt(
    anova(
      update(
        date_covariate_model,
        . ~ .
        - Species:Immature_f:lat5
        - Species:Immature_f:lon5
      ),
      date_covariate_model,
      test = "Chisq"
    ),
    "Species x age x geography"
  ),
  extract_lrt(
    anova(
      update(
        date_covariate_model,
        . ~ . - Species:Immature_f:lat5
      ),
      date_covariate_model,
      test = "Chisq"
    ),
    "Species x age x latitude"
  ),
  extract_lrt(
    anova(
      update(
        date_covariate_model,
        . ~ . - Species:Immature_f:lon5
      ),
      date_covariate_model,
      test = "Chisq"
    ),
    "Species x age x longitude"
  )
)

date_covariate_lrt$Analysis <- "Day-of-window covariate"
date_covariate_lrt$N <- nrow(date_covariate_df)

date_covariate_differences <- age_slope_differences(
  date_covariate_model
)
date_covariate_differences$Analysis <- "Day-of-window covariate"
date_covariate_differences$N <- nrow(date_covariate_df)


# ============================================================
# 8. WESTERN-BOUNDARY SENSITIVITY
# ============================================================

# Apply the 97-degree W boundary before selecting the first qualifying record
# for each bird-winter so later records within the restricted region remain
# eligible for inclusion.
longitude_boundary_records <- select_first_record_per_bird_winter(
  raw_primary[
    Species %in% focal_species &
      lon_dd >= restricted_lon_min
  ],
  complete_only = TRUE
)

longitude_boundary <- fit_rebuilt_sensitivity(
  longitude_boundary_records,
  "East of 97 W"
)

longitude_boundary_lrt <-
  longitude_boundary$lrt

longitude_boundary_differences <-
  longitude_boundary$differences


# ============================================================
# 9. PRIMARY MODEL DIAGNOSTICS
# ============================================================

# Generate simulation-based residuals from the primary comparative GLMM.
sim_primary <- DHARMa::simulateResiduals(
  fittedModel = m_comparative_full,
  n = dharma_simulations,
  seed = diagnostics_seed,
  plot = FALSE
)

uniformity_test <- DHARMa::testUniformity(
  sim_primary,
  plot = FALSE
)

outlier_test <- DHARMa::testOutliers(
  sim_primary,
  plot = FALSE
)

dharma_observation_level <- data.frame(
  Uniformity_KS = as.numeric(uniformity_test$statistic),
  Uniformity_p = uniformity_test$p.value,
  Outlier_p = outlier_test$p.value,
  stringsAsFactors = FALSE
)

# Project coordinates for residual spatial-autocorrelation diagnostics.
comparative_xy <- comparative_df %>%
  sf::st_as_sf(
    coords = c(
      "lon_dd",
      "lat_dd"
    ),
    crs = 4326,
    remove = FALSE
  ) %>%
  sf::st_transform(
    projection_epsg
  ) %>%
  sf::st_coordinates()

spatial_df <- data.frame(
  x_km = comparative_xy[, 1] / 1000,
  y_km = comparative_xy[, 2] / 1000
)

# Aggregate residuals across species to one value per reported location.
location_group <- factor(
  as.character(
    comparative_df$common_station_id
  )
)

sim_location <- DHARMa::recalculateResiduals(
  sim_primary,
  group = location_group
)

location_coords <- stats::aggregate(
  spatial_df,
  by = list(
    group = location_group
  ),
  FUN = mean
)

if (nrow(location_coords) != length(sim_location$scaledResiduals)) {
  stop(
    "Aggregated residuals and coordinates have different lengths."
  )
}

if (any(duplicated(location_coords[, c("x_km", "y_km")]))) {
  stop(
    "Duplicate coordinates remain after aggregation to reported location."
  )
}

spatial_test_pooled <- DHARMa::testSpatialAutocorrelation(
  sim_location,
  x = location_coords$x_km,
  y = location_coords$y_km,
  plot = FALSE
)

dharma_location_level <- data.frame(
  N_locations = nrow(location_coords),
  Morans_I = as.numeric(
    spatial_test_pooled$statistic["observed"]
  ),
  Morans_I_expected = as.numeric(
    spatial_test_pooled$statistic["expected"]
  ),
  Morans_I_sd = as.numeric(
    spatial_test_pooled$statistic["sd"]
  ),
  Morans_I_p = spatial_test_pooled$p.value,
  stringsAsFactors = FALSE
)

# Aggregate residuals separately for each species and reported location.
species_site_group <- factor(
  as.character(
    comparative_df$species_site_id
  )
)

sim_species_site <- DHARMa::recalculateResiduals(
  sim_primary,
  group = species_site_group
)

species_site_coords <- stats::aggregate(
  spatial_df,
  by = list(
    group = species_site_group
  ),
  FUN = mean
)

species_site_coords$Species <- sub(
  "__.*$",
  "",
  as.character(
    species_site_coords$group
  )
)

if (nrow(species_site_coords) != length(sim_species_site$scaledResiduals)) {
  stop(
    "Species-by-location residuals and coordinates are misaligned."
  )
}

species_site_coords$residual <-
  sim_species_site$scaledResiduals

spatial_test_by_species <- dplyr::bind_rows(
  lapply(
    focal_species,
    function(sp) {
      d <- species_site_coords[
        species_site_coords$Species == sp,
      ]
      
      distances <- as.matrix(
        stats::dist(
          d[, c("x_km", "y_km")]
        )
      )
      
      weights <- 1 / distances
      diag(weights) <- 0
      weights[!is.finite(weights)] <- 0
      
      mi <- ape::Moran.I(
        d$residual,
        weight = weights,
        scaled = FALSE
      )
      
      data.frame(
        Species = sp,
        N_locations = nrow(d),
        Morans_I = mi$observed,
        Morans_I_expected = mi$expected,
        Morans_I_sd = mi$sd,
        Morans_I_p = mi$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

comparative_primary_diagnostics <- list(
  convergence = comparative_primary$convergence,
  observation_level = dharma_observation_level,
  location_level = dharma_location_level,
  spatial_by_species = spatial_test_by_species
)


# ============================================================
# 10. NONLINEAR SPATIAL GAMM SENSITIVITY
# ============================================================

# GAMM settings and cache file.

gamm_species <- c("Rufous", "Black-chinned", "Ruby-throated")
gamm_basis_k <- 20L
gamm_projection_epsg <- 5070L
gamm_spatial_scale_m <- 100000
gamm_kcheck_seed <- 20260823L

gamm_results_cache_file <-
  "gamm_results_PRIMARY_NOV15_DEC31_EPSG5070.rds"

# Set to TRUE to ignore the saved cache and refit the GAMMs and k-checks.
force_recompute_gamms <- FALSE

expected_gamm_settings <- list(
  species = gamm_species,
  basis_k = gamm_basis_k,
  projection = "EPSG:5070",
  spatial_units = "100 km",
  spatial_scale_m = gamm_spatial_scale_m,
  family = "binomial(logit)",
  method = "fREML",
  discrete = TRUE,
  kcheck_seed = gamm_kcheck_seed
)


# ------------------------------------------------------------
# Cache validation helpers
# ------------------------------------------------------------

build_gamm_data_signature <- function(df) {
  
  required <- c(
    "Species", "band", "winter_year",
    "lat_dd", "lon_dd", "Female", "Immature_f"
  )
  
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      "`comparative_df` is missing required GAMM column(s): ",
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
  
  x <- x[do.call(order, x), , drop = FALSE]
  rownames(x) <- NULL
  
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 2)
  
  list(
    md5 = unname(tools::md5sum(tmp)),
    N = nrow(x),
    species_N = as.data.frame(table(x$Species), stringsAsFactors = FALSE)
  )
}


load_gamm_cache <- function(path) {
  
  x <- readRDS(path)
  
  required <- c(
    "smooth_summary",
    "k_check_summary",
    "model_metadata",
    "data_signature",
    "settings"
  )
  
  if (!is.list(x) || !all(required %in% names(x))) {
    stop(
      "The GAMM results cache is incomplete or incompatible. ",
      "Delete it and rerun 03_sensitivity_analysis.R, or restore the matching cache."
    )
  }
  
  if (!isTRUE(all.equal(
    x$settings,
    expected_gamm_settings,
    check.attributes = FALSE
  ))) {
    stop(
      "The GAMM cache was created with different analysis settings. ",
      "Set `force_recompute_gamms <- TRUE` for an intentional rerun."
    )
  }
  
  if (
    exists(
      "comparative_df",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  ) {
    current_signature <- build_gamm_data_signature(
      get(
        "comparative_df",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    )
    
    if (!identical(current_signature$md5, x$data_signature$md5)) {
      stop(
        "The GAMM cache does not match the currently loaded `comparative_df`. ",
        "If the primary data changed intentionally, set ",
        "`force_recompute_gamms <- TRUE` and rerun."
      )
    }
  }
  
  x
}


# ------------------------------------------------------------
# Load matching cached results when available
# ------------------------------------------------------------

if (!force_recompute_gamms && file.exists(gamm_results_cache_file)) {
  
  gamm_results <- load_gamm_cache(gamm_results_cache_file)
  
  gamm_sensitivity_summary <- gamm_results$smooth_summary
  gamm_kcheck_summary <- gamm_results$k_check_summary
  gamm_model_metadata <- gamm_results$model_metadata
  
  message(
    "Loaded cached GAMM results from ", gamm_results_cache_file,
    "; no GAMMs or k-checks were rerun."
  )
  
  
  # ----------------------------------------------------------
  # Otherwise fit GAMMs, run k-checks, and save the results
  # ----------------------------------------------------------
  
} else {
  
  if (!exists("comparative_df", inherits = FALSE)) {
    stop(
      "`comparative_df` was not found. Run 01_functions.R and ",
      "02_primary_analysis.R before fitting the GAMMs."
    )
  }
  
  required <- c(
    "Species", "Female", "Immature_f", "lon_dd", "lat_dd",
    "species_site_id", "species_winter_year", "species_band_id"
  )
  
  missing_cols <- setdiff(required, names(comparative_df))
  if (length(missing_cols) > 0L) {
    stop(
      "`comparative_df` is missing required GAMM column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # Project coordinates and prepare grouping variables
  # ----------------------------------------------------------
  
  comparative_gam <- comparative_df %>%
    mutate(
      Species = factor(Species, levels = gamm_species),
      Immature_f = factor(Immature_f, levels = c("Adult", "Immature")),
      Immature_ord = ordered(
        as.character(Immature_f),
        levels = c("Adult", "Immature")
      ),
      species_site_id = factor(species_site_id),
      species_winter_year = factor(species_winter_year),
      species_band_id = factor(species_band_id)
    )
  
  comparative_gam_sf <- sf::st_as_sf(
    comparative_gam,
    coords = c("lon_dd", "lat_dd"),
    crs = 4326,
    remove = FALSE
  ) %>%
    sf::st_transform(gamm_projection_epsg)
  
  gam_xy <- sf::st_coordinates(comparative_gam_sf)
  comparative_gam$x100 <- gam_xy[, 1] / gamm_spatial_scale_m
  comparative_gam$y100 <- gam_xy[, 2] / gamm_spatial_scale_m
  
  
  # ----------------------------------------------------------
  # Fit one species-specific GAMM
  # ----------------------------------------------------------
  
  fit_species_gamm <- function(species_name) {
    
    dat <- comparative_gam %>%
      filter(Species == species_name) %>%
      droplevels()
    
    mgcv::bam(
      Female ~
        Immature_f +
        s(x100, y100, bs = "tp", k = gamm_basis_k) +
        s(
          x100, y100,
          by = Immature_ord,
          bs = "tp",
          k = gamm_basis_k
        ) +
        s(species_site_id, bs = "re") +
        s(species_winter_year, bs = "re") +
        s(species_band_id, bs = "re"),
      family = binomial(link = "logit"),
      method = "fREML",
      discrete = TRUE,
      data = dat
    )
  }
  
  gamm_models <- setNames(
    lapply(gamm_species, fit_species_gamm),
    gamm_species
  )
  
  
  # ----------------------------------------------------------
  # Extract spatial smooth results
  # ----------------------------------------------------------
  
  extract_gamm_smooth_summary <- function(model, species_name) {
    
    sm <- as.data.frame(summary(model)$s.table)
    sm$Term <- rownames(sm)
    rownames(sm) <- NULL
    
    sm %>%
      filter(
        Term %in% c(
          "s(x100,y100)",
          "s(x100,y100):Immature_ordImmature"
        )
      ) %>%
      transmute(
        Species = species_name,
        Term,
        Smooth = case_when(
          Term == "s(x100,y100)" ~
            "Adult spatial surface",
          Term == "s(x100,y100):Immature_ordImmature" ~
            "Immature - Adult spatial difference",
          TRUE ~ Term
        ),
        edf,
        Ref_df = Ref.df,
        Chi_square = Chi.sq,
        p = `p-value`
      )
  }
  
  gamm_sensitivity_summary <- bind_rows(
    lapply(
      names(gamm_models),
      function(sp) extract_gamm_smooth_summary(gamm_models[[sp]], sp)
    )
  )
  
  
  # ----------------------------------------------------------
  # Extract model metadata for Supplementary Material Table S7
  # ----------------------------------------------------------
  
  gamm_model_metadata <- bind_rows(
    lapply(
      names(gamm_models),
      function(sp) {
        
        model <- gamm_models[[sp]]
        rank <- as.integer(model$rank)
        coefficients <- as.integer(length(coef(model)))
        
        data.frame(
          Species = sp,
          N = as.integer(nobs(model)),
          Rank = rank,
          Coefficients = coefficients,
          Model_rank = paste0(rank, "/", coefficients),
          Formula = paste(deparse(formula(model)), collapse = " "),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  
  # ----------------------------------------------------------
  # Run basis-dimension checks
  # ----------------------------------------------------------
  
  set.seed(gamm_kcheck_seed)
  
  gamm_k_checks <- setNames(
    lapply(gamm_models, mgcv::k.check),
    names(gamm_models)
  )
  
  gamm_kcheck_summary <- bind_rows(
    lapply(
      names(gamm_k_checks),
      function(sp) {
        
        kc <- gamm_k_checks[[sp]]
        
        data.frame(
          Species = sp,
          Term = rownames(kc),
          k_prime = as.numeric(kc[, 1]),
          k_edf = as.numeric(kc[, 2]),
          k_index = as.numeric(kc[, 3]),
          k_p = as.numeric(kc[, 4]),
          row.names = NULL,
          check.names = FALSE
        )
      }
    )
  )
  
  
  # ----------------------------------------------------------
  # Save all downstream GAMM results in one cache
  # ----------------------------------------------------------
  
  gamm_results <- list(
    smooth_summary = gamm_sensitivity_summary,
    k_check_summary = gamm_kcheck_summary,
    model_metadata = gamm_model_metadata,
    data_signature = build_gamm_data_signature(comparative_df),
    settings = expected_gamm_settings
  )
  
  saveRDS(gamm_results, gamm_results_cache_file)
  
  message(
    "Saved GAMM results to ", gamm_results_cache_file,
    ". Future runs will load this cache unless ",
    "`force_recompute_gamms <- TRUE`."
  )
}

