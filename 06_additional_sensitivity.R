# ============================================================
# 06_additional_sensitivity.R
#
# Additional sensitivity analyses, diagnostics, and reporting quantities
# that address issues likely to be raised in review. Nothing in this script
# alters the primary model or any result already reported; every analysis is
# either a robustness check on the primary model or a descriptive quantity
# that should be stated in the manuscript.
#
# Run after 02_primary_analysis.R. Section 10 additionally uses the GAMM
# results from 03_sensitivity_analysis.R when they are available, but the
# rest of the script does not require 03.
#
# Sections
#   1  Settings
#   2  Record composition audit: repeat individuals, sites, years, precision
#   3  Exact likelihood-ratio p-values for the primary model
#   4  Response-scale predicted probabilities (effect sizes)
#   5  Residual spatial autocorrelation and site-level model diagnostics
#   6  Bandings-only sensitivity
#   7  Longitude-boundary sensitivity (east of 97 W)
#   8  Within-window date sensitivity (ageing reliability)
#   9  Geographic variance retained after within-year centering
#  10  GAMM adult-surface nonlinearity
#  11  Write outputs
#
# Outputs are written to `output_dir` as CSV files plus a single RDS holding
# every result object. The CSVs are deliberately tidy so that any of them can
# be dropped into the flextable machinery already used by 05_supplement.R.
# 
# Oh, and I am writing all this without the data, so you will most likely 
# have to adjust it to reflect what is actually in your csv file
# ============================================================


library(data.table)
library(dplyr)
library(lme4)
library(emmeans)
library(DHARMa)
library(sf)


# ============================================================
# 0. CHECK REQUIRED INPUTS
# ============================================================

required_objects <- c(
  "comparative_df",
  "m_comparative_full",
  "comparative_primary",
  "results",
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
    "Run 01_functions.R and 02_primary_analysis.R first. Missing object(s): ",
    paste(missing_objects, collapse = ", ")
  )
}

required_functions <- c(
  "fit_comparative_age_sex_model",
  "classify_immature",
  "in_month_day_window",
  "infer_season_year_fun"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "The required functions from 01_functions.R are not loaded. Missing: ",
    paste(missing_functions, collapse = ", ")
  )
}


# ============================================================
# 1. SETTINGS
# ============================================================

output_dir <- "sensitivity_06"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# DHARMa simulations. The primary analysis used 250; 1,000 is the usual
# default.
dharma_n_sim <- 1000L
dharma_seed <- 123L

# Restricted longitude boundary for the range-edge sensitivity analysis.
# Black-chinned Hummingbirds winter regularly west of roughly this line, so
# the primary 100 W boundary sits inside rather than outside their winter
# range. Refitting east of 97 W tests whether the Black-chinned longitudinal
# result depends on records near the western edge of the study region.
restricted_lon_min <- -97

# Projection used for distance calculations (NAD83 / Conus Albers), matching
# the GAMM sensitivity analysis in 03_sensitivity_analysis.R.
projection_epsg <- 5070L

# Optional columns in the BBL extract. Leave as NULL to auto-detect. Set
# explicitly if auto-detection reports the wrong column or an unexpected set
# of values.
#
#   record_type_col   column distinguishing banding from encounter records
#   banding_values    value(s) in that column denoting an original banding
#   precision_col     column giving reported coordinate precision
#   exact_precision_values  value(s) denoting exact coordinates
record_type_col <- NULL
banding_values <- NULL
precision_col <- NULL
exact_precision_values <- NULL

record_type_candidates <- c(
  "record_source", "event_type", "record_type",
  "band_type", "encounter_type"
)

precision_candidates <- c(
  "coordinates_precision", "coordinate_precision",
  "coord_precision", "location_precision", "precision"
)

message("06_additional_sensitivity.R: writing results to ", output_dir)


# ============================================================
# 1b. SMALL HELPERS
# ============================================================

# emmeans summary column names differ between asymptotic and t-based fits.
# Pull whichever of a set of candidate names is present.
pick_col <- function(df, candidates, required = TRUE, label = NULL) {
  hit <- intersect(candidates, names(df))

  if (length(hit) == 0L) {
    if (required) {
      stop(
        "Could not find column ",
        if (is.null(label)) paste(candidates, collapse = "/") else label,
        " in the emmeans summary."
      )
    }
    return(NULL)
  }

  hit[1L]
}


# Standardize an emmeans summary data frame to a common set of column names.
tidy_emm <- function(x) {
  d <- as.data.frame(x)

  est_col <- pick_col(
    d,
    c(
      "estimate", "lat5.trend", "lon5.trend",
      "lat5_within.trend", "lon5_within.trend"
    ),
    label = "estimate"
  )

  lcl_col <- pick_col(d, c("asymp.LCL", "lower.CL", "LCL"), label = "lower CL")
  ucl_col <- pick_col(d, c("asymp.UCL", "upper.CL", "UCL"), label = "upper CL")
  stat_col <- pick_col(d, c("z.ratio", "t.ratio"), required = FALSE)

  out <- data.frame(
    Estimate = as.numeric(d[[est_col]]),
    SE = as.numeric(d[["SE"]]),
    CI_low = as.numeric(d[[lcl_col]]),
    CI_high = as.numeric(d[[ucl_col]]),
    Statistic = if (is.null(stat_col)) NA_real_ else as.numeric(d[[stat_col]]),
    p = as.numeric(d[["p.value"]]),
    stringsAsFactors = FALSE
  )

  for (nm in intersect(c("contrast", "Species", "Immature_f"), names(d))) {
    out[[nm]] <- as.character(d[[nm]])
  }

  out
}


# Immature - adult geographic slope differences from any model that uses the
# primary fixed-effect structure.
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
        emmeans::contrast(em, method = "revpairwise"),
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
      dplyr::any_of(c("Axis", "Species", "contrast")),
      Estimate, SE, CI_low, CI_high, Statistic, p
    )
}


# Age-specific geographic slopes from any model using the primary structure.
age_specific_slopes <- function(
    model,
    lat_var = "lat5",
    lon_var = "lon5"
) {
  one_axis <- function(var_name, axis_label) {
    d <- tidy_emm(
      summary(
        emmeans::emtrends(
          model,
          specs = ~ Species * Immature_f,
          var = var_name
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
      dplyr::any_of(c("Axis", "Species", "Immature_f")),
      Estimate, SE, CI_low, CI_high, Statistic, p
    )
}


# Pull chi-square, df, and an unrounded p-value out of an anova() comparison.
extract_lrt <- function(x, test_name) {
  if (is.null(x)) {
    return(
      data.frame(
        Test = test_name,
        Chi_square = NA_real_,
        df = NA_real_,
        p = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  a <- as.data.frame(x)
  last <- a[nrow(a), , drop = FALSE]

  chisq_col <- intersect(c("Chisq", "ChiSq", "Chi_square"), names(last))[1]
  df_col <- intersect(c("Df", "Chi Df", "df"), names(last))[1]
  p_col <- intersect(c("Pr(>Chisq)", "p", "p.value"), names(last))[1]

  data.frame(
    Test = test_name,
    Chi_square = if (!is.na(chisq_col)) as.numeric(last[[chisq_col]]) else NA_real_,
    df = if (!is.na(df_col)) as.numeric(last[[df_col]]) else NA_real_,
    p = if (!is.na(p_col)) as.numeric(last[[p_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
}


# Drop unused levels from the grouping factors before refitting a subset.
drop_subset_levels <- function(d) {
  factor_cols <- intersect(
    c(
      "Species", "Immature_f", "site_id", "species_site_id",
      "species_winter_year", "species_band_id", "common_station_id"
    ),
    names(d)
  )

  for (nm in factor_cols) {
    if (is.factor(d[[nm]])) {
      d[[nm]] <- droplevels(d[[nm]])
    }
  }

  d
}


# Refit the primary model on a subset of the primary data and return the
# omnibus tests and the within-species age differences.
refit_subset <- function(d, label) {
  d <- drop_subset_levels(as.data.frame(d))

  fit <- fit_comparative_age_sex_model(
    comparative_df = d,
    optimizer = "bobyqa",
    maxfun = 2e5,
    run_lrt = TRUE,
    pairwise_adjust = "tukey"
  )

  lrt <- dplyr::bind_rows(
    extract_lrt(fit$lrt$overall, "Species x age x geography"),
    extract_lrt(fit$lrt$latitude, "Species x age x latitude"),
    extract_lrt(fit$lrt$longitude, "Species x age x longitude")
  )

  lrt$Analysis <- label
  lrt$N <- nrow(d)

  differences <- age_slope_differences(fit$model)
  differences$Analysis <- label
  differences$N <- nrow(d)

  list(
    fit = fit,
    lrt = lrt,
    differences = differences,
    n_by_species = as.data.frame(dplyr::count(d, Species, name = "N"))
  )
}


# Read the BBL extract once, keeping the columns needed here plus any optional
# columns that are actually present in the file.
read_bbl_with_dates <- function(csv_path, extra_cols = character(0)) {
  header <- names(
    data.table::fread(csv_path, nrows = 0L)
  )

  base_cols <- c(
    "band", "species_id", "event_year", "event_month", "event_day",
    "iso_country", "lat_dd", "lon_dd", "sex_code", "age_code"
  )

  missing_base <- setdiff(base_cols, header)

  if (length(missing_base) > 0L) {
    stop(
      "The BBL extract is missing expected column(s): ",
      paste(missing_base, collapse = ", ")
    )
  }

  keep <- unique(c(base_cols, intersect(extra_cols, header)))

  DT <- data.table::fread(csv_path, select = keep, showProgress = FALSE)

  DT[, event_year := as.integer(event_year)]
  DT[, event_month := as.integer(event_month)]
  DT[, event_day := as.integer(event_day)]

  list(
    data = DT,
    header = header,
    kept_optional = setdiff(keep, base_cols)
  )
}


# Reproduce the primary-window record selection, retaining the event date and
# any optional columns. This mirrors read_winter_first_known_age_sex() in
# 01_functions.R; the reconstruction is verified against the primary data
# below before it is used for anything.
rebuild_primary_records <- function(DT) {
  date_keep <- in_month_day_window(
    month = DT$event_month,
    day = DT$event_day,
    start_month = primary_start_month,
    start_day = primary_start_day,
    end_month = primary_end_month,
    end_day = primary_end_day
  )

  d <- DT[
    species_id %in% unname(species_ids[focal_species]) &
      iso_country %in% countries_keep &
      event_month %in% months_keep &
      date_keep &
      is.finite(lat_dd) &
      is.finite(lon_dd) &
      lon_dd >= lon_min &
      lon_dd <= lon_max &
      lat_dd >= lat_min &
      lat_dd <= lat_max &
      !is.na(band)
  ]

  d[, Species := names(species_ids)[match(species_id, unname(species_ids))]]

  # Missing days are used only for ordering, after date-window filtering.
  d[, event_day_filled := event_day]
  d[is.na(event_day_filled), event_day_filled := 15L]

  d[, dt := as.IDate(
    sprintf("%04d-%02d-%02d", event_year, event_month, event_day_filled)
  )]

  d <- d[!is.na(dt)]

  season_year_fun <- infer_season_year_fun(months_keep)
  d[, winter_year := season_year_fun(event_year, event_month)]

  d[, Female := fifelse(
    sex_code %in% c(5L, 7L),
    1L,
    fifelse(sex_code %in% c(4L, 6L), 0L, NA_integer_)
  )]

  d[, Immature := classify_immature(event_month, age_code)]

  data.table::setorder(d, Species, band, winter_year, dt)

  selected <- d[, {
    known_both <- which(!is.na(Female) & !is.na(Immature))

    if (length(known_both) > 0L) {
      .SD[known_both[1L]]
    } else {
      .SD[1L]
    }
  }, by = .(Species, band, winter_year)]

  selected[!is.na(Female) & !is.na(Immature)]
}


# Choose an optional column by name, or auto-detect from a candidate list.
resolve_optional_column <- function(user_value, candidates, header, label) {
  if (!is.null(user_value)) {
    if (!user_value %in% header) {
      warning(
        "The ", label, " column '", user_value,
        "' is not present in the BBL extract; skipping."
      )
      return(NULL)
    }
    return(user_value)
  }

  hit <- intersect(candidates, header)

  if (length(hit) == 0L) {
    message(
      "No ", label, " column found in the BBL extract. ",
      "Set it explicitly at the top of 06_additional_sensitivity.R if one exists."
    )
    return(NULL)
  }

  message("Using '", hit[1L], "' as the ", label, " column.")
  hit[1L]
}


# ============================================================
# 2. RECORD COMPOSITION AUDIT
#
# Quantities the manuscript should state outright so that reviewers do not
# have to derive them from the supplementary variance components.
# ============================================================

# Work on a plain copy with character keys so that grouping, binding, and
# joining never mix factor and character columns.
audit_df <- data.frame(
  Species = as.character(comparative_df$Species),
  band = as.character(comparative_df$band),
  winter_year = as.integer(comparative_df$winter_year),
  site_id = as.character(comparative_df$site_id),
  station_id = as.character(comparative_df$common_station_id),
  Age = as.character(comparative_df$Immature_f),
  stringsAsFactors = FALSE
)

audit_df$individual_id <- paste(audit_df$Species, audit_df$band, sep = "__")

records_per_individual <- audit_df %>%
  dplyr::count(Species, band, name = "Records") %>%
  dplyr::count(Species, Records, name = "Individuals")

summarise_audit <- function(d, label) {
  data.frame(
    Species = label,
    Records = nrow(d),
    Individuals = length(unique(d$individual_id)),
    Repeat_records = nrow(d) - length(unique(d$individual_id)),
    Repeat_records_pct =
      100 * (nrow(d) - length(unique(d$individual_id))) / nrow(d),
    Individuals_in_multiple_winters = sum(table(d$individual_id) > 1L),
    Individuals_in_both_age_classes = sum(
      tapply(d$Age, d$individual_id, function(a) length(unique(a))) > 1L
    ),
    Reported_locations = length(unique(d$station_id)),
    Winter_years = length(unique(d$winter_year)),
    First_winter_year = min(d$winter_year, na.rm = TRUE),
    Last_winter_year = max(d$winter_year, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

individual_audit <- dplyr::bind_rows(
  lapply(
    focal_species,
    function(sp) {
      summarise_audit(audit_df[audit_df$Species == sp, , drop = FALSE], sp)
    }
  )
)

individual_audit <- dplyr::bind_rows(
  individual_audit,
  summarise_audit(audit_df, "All species")
)

message(
  "Repeat individuals: ",
  sum(individual_audit$Repeat_records[
    individual_audit$Species == "All species"
  ]),
  " of ", nrow(comparative_df), " records."
)


# ============================================================
# 3. EXACT LIKELIHOOD-RATIO P-VALUES FOR THE PRIMARY MODEL
#
# The manuscript reports these rounded to two decimals, which puts the
# omnibus and temporal longitude tests exactly on the 0.05 threshold.
# ============================================================

primary_lrt_exact <- dplyr::bind_rows(
  extract_lrt(comparative_primary$lrt$overall, "Species x age x geography"),
  extract_lrt(comparative_primary$lrt$latitude, "Species x age x latitude"),
  extract_lrt(comparative_primary$lrt$longitude, "Species x age x longitude")
)

primary_lrt_exact$p_rounded_3 <- round(primary_lrt_exact$p, 3)
primary_lrt_exact$Analysis <- "Primary model"

print(primary_lrt_exact)


# ============================================================
# 4. RESPONSE-SCALE PREDICTED PROBABILITIES
#
# Every estimate in the manuscript is on the log-odds scale. These are the
# same effects expressed as predicted probabilities that a record is female,
# averaged over the observed distribution of the other coordinate, with
# random effects excluded.
# ============================================================

predict_along_axis <- function(species_name, axis = c("latitude", "longitude")) {
  axis <- match.arg(axis)

  base <- comparative_df %>%
    dplyr::filter(Species == species_name)

  if (nrow(base) == 0L) {
    return(NULL)
  }

  values <- if (axis == "latitude") {
    stats::quantile(base$lat_dd, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
  } else {
    stats::quantile(base$lon_dd, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
  }

  grid <- expand.grid(
    Age = c("Adult", "Immature"),
    Value = as.numeric(values),
    stringsAsFactors = FALSE
  )

  grid$Quantile <- rep(names(values), each = 2L)

  out <- lapply(
    seq_len(nrow(grid)),
    function(i) {
      newdata <- base

      if (axis == "latitude") {
        newdata$lat_dd <- grid$Value[i]
        newdata$lat5 <- (grid$Value[i] - 31) / 5
      } else {
        newdata$lon_dd <- grid$Value[i]
        newdata$lon5 <- (grid$Value[i] + 90) / 5
      }

      newdata$Immature_f <- factor(
        grid$Age[i],
        levels = levels(comparative_df$Immature_f)
      )

      newdata$Immature <- as.integer(grid$Age[i] == "Immature")

      p <- stats::predict(
        m_comparative_full,
        newdata = newdata,
        type = "response",
        re.form = NA
      )

      data.frame(
        Species = species_name,
        Axis = axis,
        Quantile = grid$Quantile[i],
        Value = grid$Value[i],
        Age = grid$Age[i],
        Predicted_female_probability = mean(p),
        stringsAsFactors = FALSE
      )
    }
  )

  dplyr::bind_rows(out)
}

predicted_probabilities <- dplyr::bind_rows(
  lapply(
    focal_species,
    function(sp) {
      dplyr::bind_rows(
        predict_along_axis(sp, "latitude"),
        predict_along_axis(sp, "longitude")
      )
    }
  )
)

print(predicted_probabilities)


# ============================================================
# 5. RESIDUAL SPATIAL AUTOCORRELATION AND SITE-LEVEL DIAGNOSTICS
#
# For a Bernoulli response, overdispersion is not identifiable and the
# observation-level uniformity test has very little power, so the diagnostics
# currently reported cannot detect the failure that matters here. Residuals
# are therefore aggregated to the reported-location level, where they are
# informative, and tested for spatial autocorrelation.
# ============================================================

comparative_xy <- comparative_df %>%
  sf::st_as_sf(
    coords = c("lon_dd", "lat_dd"),
    crs = 4326,
    remove = FALSE
  ) %>%
  sf::st_transform(projection_epsg) %>%
  sf::st_coordinates()

# Held separately rather than appended to comparative_df so that the object
# created by 02_primary_analysis.R is left exactly as the other scripts expect.
spatial_df <- data.frame(
  x_km = comparative_xy[, 1] / 1000,
  y_km = comparative_xy[, 2] / 1000
)

set.seed(dharma_seed)

sim_primary <- DHARMa::simulateResiduals(
  fittedModel = m_comparative_full,
  n = dharma_n_sim,
  plot = FALSE
)

# ---- Observation-level diagnostics, recomputed at n = 1000 ----

dharma_observation_level <- data.frame(
  Level = "Observation",
  Dispersion = as.numeric(
    DHARMa::testDispersion(sim_primary, plot = FALSE)$statistic
  ),
  Dispersion_p = DHARMa::testDispersion(sim_primary, plot = FALSE)$p.value,
  Uniformity_KS = as.numeric(
    DHARMa::testUniformity(sim_primary, plot = FALSE)$statistic
  ),
  Uniformity_p = DHARMa::testUniformity(sim_primary, plot = FALSE)$p.value,
  Outlier_p = DHARMa::testOutliers(sim_primary, plot = FALSE)$p.value,
  stringsAsFactors = FALSE
)

# ---- Location-level aggregation, pooled across species ----
#
# Aggregating by reported coordinate (rather than by species x coordinate)
# gives one residual per unique location, which is what the spatial test
# requires. DHARMa aggregates in the sorted order of the grouping factor, so
# the site coordinates are aggregated the same way.

location_group <- factor(as.character(comparative_df$common_station_id))

sim_location <- DHARMa::recalculateResiduals(
  sim_primary,
  group = location_group
)

location_coords <- stats::aggregate(
  spatial_df,
  by = list(group = location_group),
  FUN = mean
)

if (nrow(location_coords) != length(sim_location$scaledResiduals)) {
  stop(
    "Aggregated residuals and aggregated coordinates have different lengths; ",
    "the grouping order cannot be assumed to match."
  )
}

if (any(duplicated(location_coords[, c("x_km", "y_km")]))) {
  stop(
    "Duplicate coordinates remain after aggregating to reported location; ",
    "the spatial autocorrelation test requires unique locations."
  )
}

spatial_test_pooled <- DHARMa::testSpatialAutocorrelation(
  sim_location,
  x = location_coords$x_km,
  y = location_coords$y_km,
  plot = FALSE
)

dharma_location_level <- data.frame(
  Level = "Reported location (pooled across species)",
  N_locations = nrow(location_coords),
  Dispersion = as.numeric(
    DHARMa::testDispersion(sim_location, plot = FALSE)$statistic
  ),
  Dispersion_p = DHARMa::testDispersion(sim_location, plot = FALSE)$p.value,
  Uniformity_KS = as.numeric(
    DHARMa::testUniformity(sim_location, plot = FALSE)$statistic
  ),
  Uniformity_p = DHARMa::testUniformity(sim_location, plot = FALSE)$p.value,
  Morans_I = as.numeric(spatial_test_pooled$statistic["observed"]),
  Morans_I_expected = as.numeric(spatial_test_pooled$statistic["expected"]),
  Morans_I_sd = as.numeric(spatial_test_pooled$statistic["sd"]),
  Morans_I_p = spatial_test_pooled$p.value,
  stringsAsFactors = FALSE
)

# ---- Species-specific Moran's I on species x location residuals ----
#
# The pooled test cannot separate species. Aggregating by species x location
# and testing each species separately is more informative, but a DHARMa
# object cannot be subset, so Moran's I is computed directly with ape.

species_site_group <- factor(as.character(comparative_df$species_site_id))

sim_species_site <- DHARMa::recalculateResiduals(
  sim_primary,
  group = species_site_group
)

species_site_coords <- stats::aggregate(
  spatial_df,
  by = list(group = species_site_group),
  FUN = mean
)

species_site_coords$Species <- sub(
  "__.*$",
  "",
  as.character(species_site_coords$group)
)

if (nrow(species_site_coords) != length(sim_species_site$scaledResiduals)) {
  stop("Species x location residuals and coordinates are misaligned.")
}

species_site_coords$residual <- sim_species_site$scaledResiduals

spatial_test_by_species <- NULL

if (requireNamespace("ape", quietly = TRUE)) {

  spatial_test_by_species <- dplyr::bind_rows(
    lapply(
      focal_species,
      function(sp) {
        d <- species_site_coords[species_site_coords$Species == sp, ]

        if (nrow(d) < 10L) {
          return(
            data.frame(
              Species = sp,
              N_locations = nrow(d),
              Morans_I = NA_real_,
              Morans_I_expected = NA_real_,
              Morans_I_sd = NA_real_,
              Morans_I_p = NA_real_,
              stringsAsFactors = FALSE
            )
          )
        }

        distances <- as.matrix(
          stats::dist(d[, c("x_km", "y_km")])
        )

        weights <- 1 / distances
        diag(weights) <- 0
        weights[!is.finite(weights)] <- 0

        # Residuals are uniform on (0, 1); centring is handled by Moran.I.
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

  print(spatial_test_by_species)

} else {
  message(
    "Package 'ape' is not installed, so species-specific Moran's I was not ",
    "computed. install.packages(\"ape\") to enable it. The pooled ",
    "location-level test above was still run."
  )
}


# ============================================================
# 6. BANDINGS-ONLY SENSITIVITY
#
# Encounter records may carry age and sex codes inherited from the original
# banding record rather than a fresh in-hand determination. This refits the
# primary model using original bandings only.
# ============================================================

bbl <- read_bbl_with_dates(
  csv_file,
  extra_cols = unique(c(record_type_candidates, precision_candidates))
)

record_type_col <- resolve_optional_column(
  record_type_col,
  record_type_candidates,
  bbl$header,
  "record type"
)

precision_col <- resolve_optional_column(
  precision_col,
  precision_candidates,
  bbl$header,
  "coordinate precision"
)

primary_records <- rebuild_primary_records(bbl$data)

# Confirm the reconstruction reproduces the primary dataset before it is used.
reconstruction_ok <- vapply(
  focal_species,
  function(sp) {
    recon <- primary_records[
      Species == sp
    ][order(band, winter_year)]

    primary <- as.data.table(
      results[[sp]]$data_points
    )[order(band, winter_year)]

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

if (!all(reconstruction_ok)) {
  stop(
    "The reconstructed primary records do not match the primary datasets; ",
    "the date-aware rebuild in 06_additional_sensitivity.R is out of sync with ",
    "01_functions.R."
  )
}

# Attach the reconstructed event date and optional columns to the primary
# analysis data frame.
join_cols <- c("Species", "band", "winter_year")

attach_cols <- c(
  "event_month", "event_day", "event_day_filled", "dt",
  record_type_col, precision_col
)

attach_cols <- intersect(attach_cols, names(primary_records))

primary_join <- as.data.frame(
  primary_records[, c(join_cols, attach_cols), with = FALSE]
)

primary_join$Species <- as.character(primary_join$Species)
primary_join$band <- as.character(primary_join$band)
primary_join$winter_year <- as.integer(primary_join$winter_year)

comparative_aug <- comparative_df
comparative_aug$Species_chr <- as.character(comparative_aug$Species)
comparative_aug$band <- as.character(comparative_aug$band)
comparative_aug$winter_year <- as.integer(comparative_aug$winter_year)

comparative_aug <- comparative_aug %>%
  dplyr::left_join(
    primary_join,
    by = c("Species_chr" = "Species", "band" = "band",
           "winter_year" = "winter_year")
  )

if (nrow(comparative_aug) != nrow(comparative_df)) {
  stop("Joining event dates duplicated rows; check the join keys.")
}

if (anyNA(comparative_aug$dt)) {
  stop("Some primary records did not receive an event date from the rebuild.")
}

record_type_summary <- NULL
bandings_only <- NULL

if (!is.null(record_type_col)) {

  record_type_summary <- comparative_aug %>%
    dplyr::count(
      Species = Species_chr,
      Record_type = as.character(.data[[record_type_col]]),
      name = "N"
    ) %>%
    dplyr::group_by(Species) %>%
    dplyr::mutate(Percent = 100 * N / sum(N)) %>%
    dplyr::ungroup()

  print(record_type_summary)

  if (is.null(banding_values)) {
    tab <- sort(
      table(as.character(comparative_aug[[record_type_col]])),
      decreasing = TRUE
    )

    banding_values <- names(tab)[1L]

    message(
      "Assuming '", banding_values, "' denotes an original banding record in ",
      "column '", record_type_col, "'. Check record_type_summary and set ",
      "`banding_values` explicitly if this is wrong."
    )
  }

  banding_rows <- as.character(comparative_aug[[record_type_col]]) %in%
    banding_values

  if (sum(banding_rows) < 100L) {
    warning(
      "Fewer than 100 records were classified as bandings; skipping the ",
      "bandings-only refit. Set `banding_values` explicitly."
    )
  } else if (all(banding_rows)) {
    message(
      "All primary records are bandings, so the bandings-only refit is ",
      "identical to the primary model and was skipped."
    )
  } else {
    bandings_only <- refit_subset(
      comparative_aug[banding_rows, , drop = FALSE],
      "Bandings only"
    )

    print(bandings_only$lrt)
    print(bandings_only$differences)
  }
}


# ============================================================
# 7. LONGITUDE-BOUNDARY SENSITIVITY
#
# The 100 W boundary sits inside the regular winter range of Black-chinned
# Hummingbirds, so for that species longitude within the study region is not
# measuring distance from the range edge in the way it is for Rufous. This
# refits the primary model east of 97 W.
# ============================================================

restricted_rows <- comparative_df$lon_dd >= restricted_lon_min

longitude_boundary <- refit_subset(
  comparative_df[restricted_rows, , drop = FALSE],
  paste0("East of ", abs(restricted_lon_min), " W")
)

longitude_boundary_n <- comparative_df %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    N_primary = dplyr::n(),
    N_retained = sum(lon_dd >= restricted_lon_min),
    N_dropped = sum(lon_dd < restricted_lon_min),
    Percent_dropped = 100 * sum(lon_dd < restricted_lon_min) / dplyr::n(),
    .groups = "drop"
  )

print(longitude_boundary_n)
print(longitude_boundary$lrt)
print(longitude_boundary$differences)


# ============================================================
# 8. WITHIN-WINDOW DATE SENSITIVITY
#
# Bill corrugation fades through the winter, so ageing reliability potentially 
# declines across the 15 November to 31 December window. (Although I'm really hoping
# they used more than just corrugation!) Three checks:
#   (a) does the apparent proportion of immatures decline across the window?
#   (b) do the geographic results change when day-of-window is included?
#   (c) do they change when the window is split into late November and
#       December?
# ============================================================

window_start <- as.Date(
  sprintf("%04d-%02d-%02d", 2000L, primary_start_month, primary_start_day)
)

comparative_aug$day_in_window <- as.integer(
  as.Date(
    sprintf(
      "%04d-%02d-%02d",
      2000L,
      comparative_aug$event_month,
      comparative_aug$event_day_filled
    )
  ) - window_start
)

comparative_aug$day_win_sc <- as.numeric(
  scale(comparative_aug$day_in_window)
)

comparative_aug$Half <- ifelse(
  comparative_aug$event_month == 11L,
  "Late November",
  "December"
)

# ---- (a) apparent age composition across the window ----

age_composition_by_date <- comparative_aug %>%
  dplyr::group_by(Species = Species_chr, Half) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Immature_n = sum(Immature_f == "Immature"),
    Immature_pct = 100 * mean(Immature_f == "Immature"),
    .groups = "drop"
  )

age_trend_models <- dplyr::bind_rows(
  lapply(
    focal_species,
    function(sp) {
      d <- comparative_aug[comparative_aug$Species_chr == sp, ]

      m <- stats::glm(
        as.integer(Immature_f == "Immature") ~ day_win_sc,
        data = d,
        family = binomial
      )

      s <- summary(m)$coefficients

      data.frame(
        Species = sp,
        N = nrow(d),
        Slope_per_SD_day = s["day_win_sc", "Estimate"],
        SE = s["day_win_sc", "Std. Error"],
        z = s["day_win_sc", "z value"],
        p = s["day_win_sc", "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  )
)

print(age_composition_by_date)
print(age_trend_models)

# ---- (b) day-of-window as a nuisance covariate in the primary model ----

date_covariate_model <- lme4::glmer(
  Female ~
    Species * Immature_f * lat5 +
    Species * Immature_f * lon5 +
    Species * Immature_f * day_win_sc +
    (1 | species_site_id) +
    (1 | species_winter_year) +
    (1 | species_band_id),
  data = comparative_aug,
  family = binomial,
  nAGQ = 1,
  control = lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

date_covariate_lrt <- dplyr::bind_rows(
  extract_lrt(
    anova(
      update(date_covariate_model, . ~ . - Species:Immature_f:lat5),
      date_covariate_model,
      test = "Chisq"
    ),
    "Species x age x latitude"
  ),
  extract_lrt(
    anova(
      update(date_covariate_model, . ~ . - Species:Immature_f:lon5),
      date_covariate_model,
      test = "Chisq"
    ),
    "Species x age x longitude"
  )
)

date_covariate_lrt$Analysis <- "Day-of-window covariate"

date_covariate_differences <- age_slope_differences(date_covariate_model)
date_covariate_differences$Analysis <- "Day-of-window covariate"

print(date_covariate_lrt)
print(date_covariate_differences)

# ---- (c) split the window ----

window_halves <- lapply(
  c("Late November", "December"),
  function(half) {
    rows <- comparative_aug$Half == half

    if (sum(rows) < 300L) {
      message(
        "Only ", sum(rows), " records in '", half,
        "'; the split refit was skipped for this half."
      )
      return(NULL)
    }

    refit_subset(
      comparative_aug[rows, , drop = FALSE],
      paste0("Window half: ", half)
    )
  }
)

names(window_halves) <- c("Late November", "December")
window_halves <- window_halves[!vapply(window_halves, is.null, logical(1))]

for (nm in names(window_halves)) {
  print(window_halves[[nm]]$lrt)
  print(window_halves[[nm]]$differences)
}


# ============================================================
# 8b. COORDINATE-PRECISION SENSITIVITY (OPTIONAL)
#
# The Methods note that sites do not necessarily correspond to physical
# banding locations because coordinate precision varies. If the extract
# carries a precision field, refit using exact coordinates only.
# ============================================================

precision_summary <- NULL
exact_coordinates_only <- NULL

if (!is.null(precision_col)) {

  precision_summary <- comparative_aug %>%
    dplyr::count(
      Precision = as.character(.data[[precision_col]]),
      name = "N"
    ) %>%
    dplyr::mutate(Percent = 100 * N / sum(N)) %>%
    dplyr::arrange(dplyr::desc(N))

  print(precision_summary)

  if (is.null(exact_precision_values)) {
    message(
      "`exact_precision_values` is not set, so the exact-coordinate refit was ",
      "skipped. Inspect precision_summary and set it, e.g. ",
      "exact_precision_values <- \"exact\"."
    )
  } else {
    exact_rows <- as.character(comparative_aug[[precision_col]]) %in%
      exact_precision_values

    if (sum(exact_rows) < 300L) {
      warning("Too few exact-coordinate records to refit; skipping.")
    } else {
      exact_coordinates_only <- refit_subset(
        comparative_aug[exact_rows, , drop = FALSE],
        "Exact coordinates only"
      )

      print(exact_coordinates_only$lrt)
      print(exact_coordinates_only$differences)
    }
  }
}


# ============================================================
# 9. GEOGRAPHIC VARIANCE RETAINED AFTER WITHIN-YEAR CENTERING
#
# The temporal-confounding analysis centres latitude and longitude within
# species and year, which discards all between-year spatial contrast. The
# loss of significance in that analysis is therefore partly a loss of
# information rather than evidence of confounding. These quantities let the
# manuscript say so with numbers.
# ============================================================

within_year_variance <- comparative_df %>%
  dplyr::group_by(Species, winter_year) %>%
  dplyr::mutate(
    lat_centered = lat_dd - mean(lat_dd, na.rm = TRUE),
    lon_centered = lon_dd - mean(lon_dd, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Latitude_variance_primary = stats::var(lat_dd, na.rm = TRUE),
    Latitude_variance_within_year = stats::var(lat_centered, na.rm = TRUE),
    Latitude_variance_retained_pct =
      100 * stats::var(lat_centered, na.rm = TRUE) /
        stats::var(lat_dd, na.rm = TRUE),
    Longitude_variance_primary = stats::var(lon_dd, na.rm = TRUE),
    Longitude_variance_within_year = stats::var(lon_centered, na.rm = TRUE),
    Longitude_variance_retained_pct =
      100 * stats::var(lon_centered, na.rm = TRUE) /
        stats::var(lon_dd, na.rm = TRUE),
    .groups = "drop"
  )

singleton_species_years <- comparative_df %>%
  dplyr::count(Species, winter_year, name = "N") %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Species_years = dplyr::n(),
    Species_years_with_one_record = sum(N == 1L),
    Records_in_singleton_years = sum(N[N == 1L]),
    Percent_records_uninformative_after_centering =
      100 * sum(N[N == 1L]) / sum(N),
    .groups = "drop"
  )

within_year_variance <- within_year_variance %>%
  dplyr::left_join(singleton_species_years, by = "Species")

print(within_year_variance)


# ============================================================
# 10. GAMM ADULT-SURFACE NONLINEARITY
#
# Table S9 shows a significant, clearly nonlinear adult spatial surface for
# Rufous Hummingbirds (edf near 6). The manuscript currently states that
# adults show little geographic variation, which is true only of linear
# gradients. This section pulls the relevant rows out for reporting.
# ============================================================

gamm_summary_source <- NULL

if (exists("gamm_sensitivity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  gamm_summary_source <- get(
    "gamm_sensitivity_summary",
    envir = .GlobalEnv,
    inherits = FALSE
  )
} else if (
  exists("gamm_results_cache_file", envir = .GlobalEnv, inherits = FALSE) &&
    file.exists(
      get("gamm_results_cache_file", envir = .GlobalEnv, inherits = FALSE)
    )
) {
  cached <- readRDS(
    get("gamm_results_cache_file", envir = .GlobalEnv, inherits = FALSE)
  )
  gamm_summary_source <- cached$smooth_summary
}

gamm_nonlinearity <- NULL

if (!is.null(gamm_summary_source)) {
  gamm_nonlinearity <- as.data.frame(gamm_summary_source) %>%
    dplyr::mutate(
      Effectively_linear = edf < 2.05,
      Note = dplyr::case_when(
        edf < 2.05 & grepl("difference", Smooth, ignore.case = TRUE) ~
          "Difference surface estimated as a plane; supports the linear model",
        edf >= 2.05 & grepl("Adult", Smooth) ~
          "Adult surface is nonlinear; report and interpret in Results",
        TRUE ~ ""
      )
    )

  print(gamm_nonlinearity)
} else {
  message(
    "GAMM results were not found. Run 03_sensitivity_analysis.R first if the ",
    "adult-surface nonlinearity summary is needed."
  )
}


# ============================================================
# 11. WRITE OUTPUTS
# ============================================================

sensitivity_06 <- list(
  individual_audit = individual_audit,
  records_per_individual = records_per_individual,
  primary_lrt_exact = primary_lrt_exact,
  predicted_probabilities = predicted_probabilities,
  dharma_observation_level = dharma_observation_level,
  dharma_location_level = dharma_location_level,
  spatial_test_by_species = spatial_test_by_species,
  record_type_summary = record_type_summary,
  bandings_only_lrt = if (is.null(bandings_only)) NULL else bandings_only$lrt,
  bandings_only_differences =
    if (is.null(bandings_only)) NULL else bandings_only$differences,
  longitude_boundary_n = longitude_boundary_n,
  longitude_boundary_lrt = longitude_boundary$lrt,
  longitude_boundary_differences = longitude_boundary$differences,
  age_composition_by_date = age_composition_by_date,
  age_trend_models = age_trend_models,
  date_covariate_lrt = date_covariate_lrt,
  date_covariate_differences = date_covariate_differences,
  window_half_lrt = dplyr::bind_rows(
    lapply(window_halves, function(x) x$lrt)
  ),
  window_half_differences = dplyr::bind_rows(
    lapply(window_halves, function(x) x$differences)
  ),
  precision_summary = precision_summary,
  exact_coordinates_lrt =
    if (is.null(exact_coordinates_only)) NULL else exact_coordinates_only$lrt,
  exact_coordinates_differences =
    if (is.null(exact_coordinates_only)) {
      NULL
    } else {
      exact_coordinates_only$differences
    },
  within_year_variance = within_year_variance,
  gamm_nonlinearity = gamm_nonlinearity,
  settings = list(
    dharma_n_sim = dharma_n_sim,
    dharma_seed = dharma_seed,
    restricted_lon_min = restricted_lon_min,
    projection_epsg = projection_epsg,
    record_type_col = record_type_col,
    banding_values = banding_values,
    precision_col = precision_col,
    exact_precision_values = exact_precision_values
  ),
  session_info = utils::capture.output(utils::sessionInfo())
)

saveRDS(
  sensitivity_06,
  file.path(output_dir, "sensitivity_06_results.rds")
)

for (nm in names(sensitivity_06)) {
  x <- sensitivity_06[[nm]]

  if (is.data.frame(x) && nrow(x) > 0L) {
    utils::write.csv(
      x,
      file.path(output_dir, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }
}

writeLines(
  sensitivity_06$session_info,
  file.path(output_dir, "sessionInfo.txt")
)

message(
  "06_additional_sensitivity.R complete. ",
  length(list.files(output_dir)),
  " files written to ", output_dir, "."
)


# ============================================================
# 12. WHAT TO DO WITH THESE RESULTS
#
#  individual_audit          -> Methods: state records, individuals, repeat
#                               records, reported locations, and year range.
#  primary_lrt_exact         -> Results: report p to three decimals.
#  predicted_probabilities   -> Abstract and Discussion: give at least one
#                               effect on the probability scale.
#  dharma_*, spatial_test_*  -> Methods and Results: replace the Bernoulli
#                               dispersion test with the location-level
#                               diagnostics and the spatial autocorrelation
#                               test.
#  record_type_summary,
#  bandings_only_*           -> Methods: banding/encounter split; Results or
#                               supplement: bandings-only refit.
#  longitude_boundary_*      -> Discussion and supplement: whether the
#                               Black-chinned longitude result depends on the
#                               western boundary.
#  age_composition_by_date,
#  age_trend_models,
#  date_covariate_*,
#  window_half_*             -> Limitations: ageing reliability across the
#                               sampling window.
#  within_year_variance      -> Results: quantify the information lost to
#                               within-year centering before describing the
#                               temporal analysis as unsupportive.
#  gamm_nonlinearity         -> Results: report the nonlinear adult Rufous
#                               surface and the near-linear difference
#                               surfaces.
# ============================================================
