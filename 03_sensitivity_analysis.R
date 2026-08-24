# ============================================================
# 03_sensitivity_analysis.R
# Run the sensitivity analyses, model diagnostics, and nonlinear spatial GAMMs.
# Run 02_primary_analysis.R first.
#
# The GAMM section uses a saved results cache when a matching cache is present.
# If the cache is absent, the GAMMs and basis-dimension checks are fit and the
# results are saved for future runs.
# ============================================================


library(data.table)
library(dplyr)
library(mgcv)
library(sf)


# ============================================================
# 0. CHECK REQUIRED INPUTS
# ============================================================

required_primary_objects <- c(
  "results",
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
  "lat5",
  "lon5"
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


# ============================================================
# 1. SENSITIVITY SETTINGS
# ============================================================

station_leave_out_n <- 4L

dharma_simulations <- 250L
diagnostics_seed <- 123L


# ============================================================
# 2. PRIMARY-WINDOW UNCERTAINTY RECONSTRUCTION
#
# Reconstruct the primary-window selection while retaining records with
# unknown sex or unclassifiable age. Within each species x band x winter, keep
# the earliest complete record when available; otherwise keep the earliest
# record. The complete-case focal subset is checked against the primary data.
# ============================================================

raw_primary <- fread(
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
  )
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
    lat_dd <= lat_max
]

raw_primary[, Species := names(species_ids)[
  match(species_id, unname(species_ids))
]]

# Missing event days are used only for ordering after date-window filtering.
raw_primary[is.na(event_day), event_day := 15L]

raw_primary[, dt := as.IDate(
  sprintf(
    "%04d-%02d-%02d",
    event_year,
    event_month,
    event_day
  )
)]

raw_primary <- raw_primary[
  !is.na(dt) &
    !is.na(band)
]

primary_year_fun <- infer_season_year_fun(months_keep)

raw_primary[, winter_year :=
              primary_year_fun(
                event_year,
                event_month
              )]

raw_primary[, Female := fifelse(
  sex_code %in% c(5L, 7L),
  1L,
  fifelse(
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

setorder(
  raw_primary,
  Species,
  band,
  winter_year,
  dt
)

selected_primary_all <- raw_primary[, {
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

selected_primary_all <- selected_primary_all[, .(
  Species = as.character(Species),
  band = as.character(band),
  winter_year = as.integer(winter_year),
  event_month = as.integer(event_month),
  age_code = as.integer(age_code),
  sex_code = as.integer(sex_code),
  lat_dd = as.numeric(lat_dd),
  lon_dd = as.numeric(lon_dd),
  Female = as.integer(Female),
  Immature = as.integer(Immature)
)]


# Verify that independently reconstructed complete cases match the primary
# focal-species datasets exactly.
primary_reconstruction_matches <- vapply(
  focal_species,
  function(sp) {
    recon <- selected_primary_all[
      Species == sp &
        !is.na(Female) &
        !is.na(Immature)
    ][
      order(band, winter_year)
    ]
    
    primary <- as.data.table(
      results[[sp]]$data_points
    )[
      order(band, winter_year)
    ]
    
    if (nrow(recon) != nrow(primary)) {
      return(FALSE)
    }
    
    isTRUE(
      all.equal(
        recon[, .(
          band,
          winter_year,
          lat_dd,
          lon_dd,
          Female,
          Immature
        )],
        primary[, .(
          band,
          winter_year,
          lat_dd,
          lon_dd,
          Female,
          Immature
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
# 3. MISSINGNESS AND EXTREME-ALLOCATION BOUNDS
#
# Calculate extreme-allocation bounds for age composition and age-specific
# sex ratios among records with missing or unclassifiable demographic data.
# ============================================================

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
  
  missingness_summary <- data.table(
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
  
  # Overall age-composition bounds.
  n_immature <- sum(x$Immature == 1L, na.rm = TRUE)
  n_adult <- sum(x$Immature == 0L, na.rm = TRUE)
  n_age_unknown <- sum(is.na(x$Immature))
  n_total <- n_immature + n_adult + n_age_unknown
  
  age_bounds <- data.table(
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
  
  # Joint age + sex uncertainty cells:
  # AF/AM = adult female/male; IF/IM = immature female/male;
  # AU/IU = known age, unknown sex; UF/UM = unknown age, known female/male;
  # UU = unknown age and unknown sex.
  AF <- nrow(x[Immature == 0L & Female == 1L])
  AM <- nrow(x[Immature == 0L & Female == 0L])
  IF <- nrow(x[Immature == 1L & Female == 1L])
  IM <- nrow(x[Immature == 1L & Female == 0L])
  AU <- nrow(x[Immature == 0L & is.na(Female)])
  IU <- nrow(x[Immature == 1L & is.na(Female)])
  UF <- nrow(x[is.na(Immature) & Female == 1L])
  UM <- nrow(x[is.na(Immature) & Female == 0L])
  UU <- nrow(x[is.na(Immature) & is.na(Female)])
  
  sex_ratio_bounds <- data.table(
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


sensitivity_results <- setNames(
  lapply(
    names(species_ids),
    calc_species_uncertainty
  ),
  names(species_ids)
)

missingness_summary_all <- rbindlist(
  lapply(
    sensitivity_results,
    function(x) x$missingness_summary
  ),
  use.names = TRUE,
  fill = TRUE
)


# ============================================================
# 4. JANUARY AGE-CLASSIFICATION LOSS
#
# Apply the focal-species filters and one-record-per-individual selection to
# January records and summarize the resulting age-classification loss.
# ============================================================

raw_january <- fread(
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
  )
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
  match(species_id, unname(species_ids))
]]

raw_january[, Female := fifelse(
  sex_code %in% c(5L, 7L),
  1L,
  fifelse(
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
raw_january[is.na(event_day), event_day := 15L]

raw_january[, dt := as.IDate(
  sprintf(
    "%04d-%02d-%02d",
    event_year,
    event_month,
    event_day
  )
)]

raw_january <- raw_january[
  !is.na(dt)
]

setorder(
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
    Unclassifiable_age_N =
      sum(
        is.na(Immature)
      ),
    Unclassifiable_age_pct =
      100 *
      mean(
        is.na(Immature)
      )
  )
]

# Summarize the age codes among January records that cannot be assigned to an
# analysis age class.
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
#
# Fit the temporal-confounding robustness model using geography centered
# within species x winter year and standardized year terms.
# ============================================================

comparative_temporal <- fit_temporal_comparative_robustness(
  comparative_df = comparative_df,
  optimizer = "bobyqa",
  maxfun = 2e5,
  run_lrt = TRUE,
  pairwise_adjust = "tukey"
)


# ============================================================
# 6. LEAVE-ONE-HIGH-VOLUME-LOCATION-OUT INFLUENCE ANALYSIS
#
# Identify the four highest-volume reported coordinate locations and refit
# the primary model after omitting each location separately.
# ============================================================

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

station_influence_species_contrasts <-
  station_influence$species_contrast_summary


# ============================================================
# 7. PRIMARY MODEL DIAGNOSTICS
#
# Summarize convergence, singularity, dispersion, residual uniformity, and
# outlier diagnostics for the primary comparative GLMM.
# ============================================================

run_comparative_primary_diagnostics <- function(
    model,
    n_sim = 250L,
    seed = 123L
) {
  conv_message <- model@optinfo$conv$lme4$messages
  optimizer_code <- model@optinfo$conv$opt
  
  max_abs_gradient <-
    if (!is.null(model@optinfo$derivs$gradient)) {
      max(abs(model@optinfo$derivs$gradient))
    } else {
      NA_real_
    }
  
  singular <- lme4::isSingular(
    model,
    tol = 1e-4
  )
  
  set.seed(seed)
  
  sim <- DHARMa::simulateResiduals(
    fittedModel = model,
    n = n_sim,
    plot = FALSE
  )
  
  dispersion_res <- DHARMa::testDispersion(
    sim,
    plot = FALSE
  )
  
  uniformity_res <- DHARMa::testUniformity(
    sim,
    plot = FALSE
  )
  
  outlier_res <- DHARMa::testOutliers(
    sim,
    plot = FALSE
  )
  
  list(
    global = data.table(
      Convergence_message =
        if (is.null(conv_message)) {
          NA_character_
        } else {
          paste(
            conv_message,
            collapse = " | "
          )
        },
      Optimizer_code =
        as.integer(
          optimizer_code
        ),
      Max_abs_gradient =
        as.numeric(
          max_abs_gradient
        ),
      Singular =
        singular,
      DHARMa_dispersion =
        as.numeric(
          dispersion_res$statistic
        ),
      DHARMa_dispersion_p =
        dispersion_res$p.value,
      DHARMa_uniformity_KS =
        as.numeric(
          uniformity_res$statistic
        ),
      DHARMa_uniformity_p =
        uniformity_res$p.value,
      DHARMa_outlier_p =
        outlier_res$p.value
    )
  )
}


comparative_primary_diagnostics <- run_comparative_primary_diagnostics(
  model = m_comparative_full,
  n_sim = dharma_simulations,
  seed = diagnostics_seed
)


# ============================================================
# 8. NONLINEAR SPATIAL GAMM SENSITIVITY
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
  # Extract model metadata for Supplementary Table S10
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