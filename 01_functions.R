# ============================================================
# 01_functions.R
# Reusable data-preparation, model-fitting, and robustness-analysis functions
# for the winter hummingbird analysis. This file is sourced by
# 02_primary_analysis.R.
# ============================================================

library(data.table)
library(dplyr)
library(lme4)


# ============================================================
# 1. DATA-PREPARATION HELPERS
# ============================================================

# Create a function that assigns calendar records to winter years.
# Windows that cross 31 December are labeled by the ending calendar year.
infer_season_year_fun <- function(months_keep) {
  m <- sort(unique(as.integer(months_keep)))
  m <- m[m >= 1L & m <= 12L]
  
  if (length(m) == 0L) {
    stop("months_keep must contain at least one integer between 1 and 12.")
  }
  
  if (length(m) == 12L) {
    return(function(y, month) y)
  }
  
  L <- length(m)
  consecutive_circular <- FALSE
  start_month <- m[1L]
  
  for (s in seq_along(m)) {
    shifted <- (m - m[s]) %% 12L
    if (identical(sort(shifted), 0:(L - 1L))) {
      consecutive_circular <- TRUE
      start_month <- m[s]
      break
    }
  }
  
  if (!consecutive_circular) {
    warning(
      "months_keep not consecutive on the year-circle; using season year = calendar year."
    )
    return(function(y, month) y)
  }
  
  raw_m <- sort(unique(as.integer(months_keep)))
  crosses_year_boundary <-
    any(raw_m < start_month) && any(raw_m >= start_month)
  
  if (!crosses_year_boundary) {
    return(function(y, month) y)
  }
  
  function(y, month) ifelse(month >= start_month, y + 1L, y)
}


# Convert BBL age codes to the two age classes used in the analysis.
#   Nov-Dec: HY = immature; AHY/SY/ASY/TY/ATY = adult.
#   Jan-Feb: SY = immature; ASY/TY/ATY = adult.
# Age codes that cannot be assigned unambiguously return NA.
classify_immature <- function(month, age_code) {
  month <- as.integer(month)
  age <- as.integer(age_code)
  out <- rep(NA_integer_, length(month))
  
  novdec <- month %in% c(11L, 12L)
  janfeb <- month %in% c(1L, 2L)
  
  out[novdec & age == 2L] <- 1L
  out[janfeb & age == 5L] <- 1L
  
  out[novdec & age %in% c(1L, 5L, 6L, 7L, 8L)] <- 0L
  out[janfeb & age %in% c(6L, 7L, 8L)] <- 0L
  
  out
}


# Identify records within an inclusive month/day window. Missing event days
# are included only when the full month falls inside the requested window.
in_month_day_window <- function(
    month,
    day,
    start_month,
    start_day,
    end_month,
    end_day
) {
  month <- as.integer(month)
  day <- as.integer(day)
  start_month <- as.integer(start_month)
  start_day <- as.integer(start_day)
  end_month <- as.integer(end_month)
  end_day <- as.integer(end_day)
  
  max_day <- c(
    31L, 29L, 31L, 30L, 31L, 30L,
    31L, 31L, 30L, 31L, 30L, 31L
  )
  
  if (
    start_month < 1L || start_month > 12L ||
    end_month < 1L || end_month > 12L ||
    start_day < 1L || start_day > max_day[start_month] ||
    end_day < 1L || end_day > max_day[end_month]
  ) {
    stop("Invalid month/day window.")
  }
  
  keep <- rep(FALSE, length(month))
  valid_month <- !is.na(month) & month >= 1L & month <= 12L
  
  if (start_month == end_month) {
    rows <- valid_month & month == start_month
    
    keep[rows & !is.na(day)] <-
      day[rows & !is.na(day)] >= start_day &
      day[rows & !is.na(day)] <= end_day
    
    if (start_day == 1L && end_day >= max_day[start_month]) {
      keep[rows & is.na(day)] <- TRUE
    }
    
    return(keep)
  }
  
  crosses_year <- start_month > end_month
  
  if (!crosses_year) {
    keep[valid_month & month > start_month & month < end_month] <- TRUE
  } else {
    keep[valid_month & (month > start_month | month < end_month)] <- TRUE
  }
  
  start_rows <- valid_month & month == start_month
  end_rows <- valid_month & month == end_month
  
  keep[start_rows & !is.na(day)] <-
    day[start_rows & !is.na(day)] >= start_day
  
  keep[end_rows & !is.na(day)] <-
    day[end_rows & !is.na(day)] <= end_day
  
  if (start_day == 1L) {
    keep[start_rows & is.na(day)] <- TRUE
  }
  
  if (end_day >= max_day[end_month]) {
    keep[end_rows & is.na(day)] <- TRUE
  }
  
  keep
}


# Assign site IDs from exact reported longitude-latitude coordinates.
assign_site_id_exact_coordinates <- function(df) {
  stopifnot(all(c("lon_dd", "lat_dd") %in% names(df)))
  
  if (any(!is.finite(df$lon_dd)) || any(!is.finite(df$lat_dd))) {
    stop(
      "All records must have finite longitude and latitude before site IDs are assigned."
    )
  }
  
  df$site_id <- interaction(
    df$lon_dd,
    df$lat_dd,
    drop = TRUE,
    lex.order = TRUE
  )
  
  df
}


# Assign a common reported-location ID across species from exact coordinates.
# `common_station_id` refers to a reported coordinate location and is used by
# the leave-one-high-volume-location-out analysis.
assign_common_location_id_exact_coordinates <- function(
    df,
    id_col = "common_station_id"
) {
  required <- c("lon_dd", "lat_dd")
  missing_cols <- setdiff(required, names(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "`df` is missing required coordinate column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (nrow(df) == 0L) {
    stop("Cannot assign reported-location IDs from an empty dataset.")
  }
  
  if (any(!is.finite(df$lon_dd)) || any(!is.finite(df$lat_dd))) {
    stop(
      "All records must have finite longitude and latitude before common location IDs are assigned."
    )
  }
  
  df[[id_col]] <- interaction(
    df$lon_dd,
    df$lat_dd,
    drop = TRUE,
    lex.order = TRUE
  )
  
  df
}


# Read one species from the BBL extract, apply analysis filters, classify
# age and sex, and retain the first qualifying record per individual per winter.
read_winter_first_known_age_sex <- function(
    csv_path,
    species_id_target,
    months_keep = c(12L, 1L, 2L),
    start_month = NULL,
    start_day = NULL,
    end_month = NULL,
    end_day = NULL,
    countries_keep = c("US", "CA"),
    states_keep = NULL,
    male_codes = c(4L, 6L),
    female_codes = c(5L, 7L),
    season_year_fun = NULL,
    lon_min = -Inf,
    lon_max = Inf,
    lat_min = -Inf,
    lat_max = Inf,
    keep_cols = c(
      "band", "species_id", "event_year", "event_month", "event_day",
      "iso_country", "iso_subdivision", "lat_dd", "lon_dd",
      "sex_code", "age_code"
    )
) {
  if (is.null(season_year_fun)) {
    season_year_fun <- infer_season_year_fun(months_keep)
  }
  
  DT <- fread(csv_path, select = keep_cols, showProgress = TRUE)
  
  DT[, event_year := as.integer(event_year)]
  DT[, event_month := as.integer(event_month)]
  DT[, event_day := as.integer(event_day)]
  
  use_date_window <- any(
    !vapply(
      list(start_month, start_day, end_month, end_day),
      is.null,
      logical(1)
    )
  )
  
  if (use_date_window) {
    if (any(vapply(
      list(start_month, start_day, end_month, end_day),
      is.null,
      logical(1)
    ))) {
      stop(
        "When using a partial-month date window, supply start_month, ",
        "start_day, end_month, and end_day."
      )
    }
    
    date_keep <- in_month_day_window(
      month = DT$event_month,
      day = DT$event_day,
      start_month = start_month,
      start_day = start_day,
      end_month = end_month,
      end_day = end_day
    )
  } else {
    date_keep <- rep(TRUE, nrow(DT))
  }
  
  DT <- DT[
    species_id == species_id_target &
      (is.null(countries_keep) | iso_country %in% countries_keep) &
      event_month %in% months_keep &
      date_keep &
      is.finite(lat_dd) &
      is.finite(lon_dd) &
      lon_dd >= lon_min &
      lon_dd <= lon_max &
      lat_dd >= lat_min &
      lat_dd <= lat_max
  ]
  
  if (!is.null(states_keep)) {
    DT <- DT[iso_subdivision %in% states_keep]
  }
  
  # Missing day is used only to order records after date-window filtering.
  DT[is.na(event_day), event_day := 15L]
  
  DT[, dt := as.IDate(sprintf(
    "%04d-%02d-%02d",
    event_year,
    event_month,
    event_day
  ))]
  
  DT <- DT[!is.na(dt) & !is.na(band)]
  
  DT[, winter_year := season_year_fun(event_year, event_month)]
  
  DT[, Female := fifelse(
    sex_code %in% female_codes,
    1L,
    fifelse(sex_code %in% male_codes, 0L, NA_integer_)
  )]
  
  DT[, Immature := classify_immature(event_month, age_code)]
  
  setorder(DT, band, winter_year, dt)
  
  # Prefer the earliest record with both known sex and classifiable age.
  # If no such record exists, keep the earliest record temporarily; unresolved
  # records are removed below from the primary complete-case dataset.
  selected <- DT[, {
    known_row <- which(!is.na(Female) & !is.na(Immature))
    
    if (length(known_row) > 0L) {
      .SD[known_row[1L]]
    } else {
      .SD[1L]
    }
  }, by = .(band, winter_year)]
  
  complete <- selected[
    !is.na(Female) &
      !is.na(Immature)
  ]
  
  out <- complete[, .(
    band = as.character(band),
    lon_dd = as.numeric(lon_dd),
    lat_dd = as.numeric(lat_dd),
    Female = as.integer(Female),
    Immature = as.integer(Immature),
    winter_year = as.integer(winter_year)
  )]
  
  as.data.frame(out)
}


# Run the species-level data-preparation pipeline and return a consistent
# results object.
winter_age_sex_pipeline <- function(
    csv_path,
    species_id,
    months_keep = c(12L, 1L, 2L),
    start_month = NULL,
    start_day = NULL,
    end_month = NULL,
    end_day = NULL,
    countries_keep = c("US", "CA"),
    states_keep = NULL,
    male_codes = c(4L, 6L),
    female_codes = c(5L, 7L),
    season_year_fun = NULL,
    lon_min = -Inf,
    lon_max = Inf,
    lat_min = -Inf,
    lat_max = Inf
) {
  list(
    data_points = read_winter_first_known_age_sex(
      csv_path = csv_path,
      species_id_target = species_id,
      months_keep = months_keep,
      start_month = start_month,
      start_day = start_day,
      end_month = end_month,
      end_day = end_day,
      countries_keep = countries_keep,
      states_keep = states_keep,
      male_codes = male_codes,
      female_codes = female_codes,
      season_year_fun = season_year_fun,
      lon_min = lon_min,
      lon_max = lon_max,
      lat_min = lat_min,
      lat_max = lat_max
    )
  )
}


# ============================================================
# 2. PRIMARY COMPARATIVE GLMM
# ============================================================

# Combine the three focal species and create geographic predictors and
# species-specific grouping factors for the comparative GLMM.
prepare_comparative_data <- function(
    results,
    comparative_species = c(
      "Rufous",
      "Black-chinned",
      "Ruby-throated"
    ),
    reference_species = "Rufous"
) {
  if (!is.list(results)) {
    stop("`results` must be a named list of species results.")
  }
  
  missing_species <- setdiff(comparative_species, names(results))
  
  if (length(missing_species) > 0L) {
    stop(
      "Missing comparative species in `results`: ",
      paste(missing_species, collapse = ", ")
    )
  }
  
  if (!reference_species %in% comparative_species) {
    stop("`reference_species` must be included in `comparative_species`.")
  }
  
  pieces <- lapply(
    comparative_species,
    function(sp) {
      x <- results[[sp]]
      
      if (is.null(x$data_points)) {
        stop(
          "`results[[\"", sp, "\"]]$data_points` is missing."
        )
      }
      
      d <- as.data.frame(x$data_points)
      
      required <- c(
        "Female",
        "Immature",
        "lat_dd",
        "lon_dd",
        "winter_year",
        "band"
      )
      
      missing_cols <- setdiff(required, names(d))
      
      if (length(missing_cols) > 0L) {
        stop(
          "Comparative data for ", sp,
          " are missing required column(s): ",
          paste(missing_cols, collapse = ", ")
        )
      }
      
      d <- assign_site_id_exact_coordinates(d)
      d$Species <- sp
      d
    }
  )
  
  comparative_df <- dplyr::bind_rows(pieces)
  
  # Center at 31 N and 90 W; one unit equals 5 degrees.
  comparative_df$lat5 <- (comparative_df$lat_dd - 31) / 5
  comparative_df$lon5 <- (comparative_df$lon_dd + 90) / 5
  
  comparative_df$Species <- factor(
    comparative_df$Species,
    levels = c(
      reference_species,
      setdiff(comparative_species, reference_species)
    )
  )
  
  comparative_df$Immature_f <- factor(
    comparative_df$Immature,
    levels = c(0, 1),
    labels = c("Adult", "Immature")
  )
  
  comparative_df$species_site_id <- interaction(
    comparative_df$Species,
    as.character(comparative_df$site_id),
    drop = TRUE,
    sep = "__"
  )
  
  comparative_df$species_winter_year <- interaction(
    comparative_df$Species,
    as.character(comparative_df$winter_year),
    drop = TRUE,
    sep = "__"
  )
  
  comparative_df$species_band_id <- interaction(
    comparative_df$Species,
    as.character(comparative_df$band),
    drop = TRUE,
    sep = "__"
  )
  
  comparative_df <- assign_common_location_id_exact_coordinates(
    comparative_df,
    id_col = "common_station_id"
  )
  
  comparative_df
}


# Fit the primary three-species GLMM and calculate the slopes, contrasts,
# likelihood-ratio tests, and convergence summaries used downstream.
fit_comparative_age_sex_model <- function(
    comparative_df,
    optimizer = "bobyqa",
    maxfun = 2e5,
    run_lrt = TRUE,
    pairwise_adjust = "tukey"
) {
  required <- c(
    "Female",
    "Species",
    "Immature_f",
    "lat5",
    "lon5",
    "species_site_id",
    "species_winter_year",
    "species_band_id"
  )
  
  missing_cols <- setdiff(required, names(comparative_df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "`comparative_df` is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  full_model <- lme4::glmer(
    Female ~
      Species * Immature_f * lat5 +
      Species * Immature_f * lon5 +
      (1 | species_site_id) +
      (1 | species_winter_year) +
      (1 | species_band_id),
    data = comparative_df,
    family = binomial,
    nAGQ = 1,
    control = lme4::glmerControl(
      optimizer = optimizer,
      optCtrl = list(maxfun = maxfun)
    )
  )
  
  conv_message <- full_model@optinfo$conv$lme4$messages
  
  convergence <- data.frame(
    Convergence_message =
      if (is.null(conv_message)) {
        NA_character_
      } else {
        paste(conv_message, collapse = " | ")
      },
    Optimizer_code = full_model@optinfo$conv$opt,
    Max_abs_gradient =
      if (!is.null(full_model@optinfo$derivs$gradient)) {
        max(abs(full_model@optinfo$derivs$gradient))
      } else {
        NA_real_
      },
    Singular = lme4::isSingular(full_model, tol = 1e-4),
    row.names = NULL
  )
  
  lrt <- list(
    overall = NULL,
    latitude = NULL,
    longitude = NULL
  )
  
  if (isTRUE(run_lrt)) {
    no_lat3 <- update(
      full_model,
      . ~ . - Species:Immature_f:lat5
    )
    
    no_lon3 <- update(
      full_model,
      . ~ . - Species:Immature_f:lon5
    )
    
    no_threeway <- update(
      full_model,
      . ~ .
      - Species:Immature_f:lat5
      - Species:Immature_f:lon5
    )
    
    lrt$latitude <- anova(no_lat3, full_model, test = "Chisq")
    lrt$longitude <- anova(no_lon3, full_model, test = "Chisq")
    lrt$overall <- anova(no_threeway, full_model, test = "Chisq")
  }
  
  lat_emtrends <- emmeans::emtrends(
    full_model,
    specs = ~ Species * Immature_f,
    var = "lat5"
  )
  
  lon_emtrends <- emmeans::emtrends(
    full_model,
    specs = ~ Species * Immature_f,
    var = "lon5"
  )
  
  latitude_slopes <- as.data.frame(
    summary(lat_emtrends, infer = c(TRUE, TRUE), level = 0.95)
  )
  
  longitude_slopes <- as.data.frame(
    summary(lon_emtrends, infer = c(TRUE, TRUE), level = 0.95)
  )
  
  lat_age_emm <- emmeans::emtrends(
    full_model,
    specs = ~ Immature_f | Species,
    var = "lat5"
  )
  
  lon_age_emm <- emmeans::emtrends(
    full_model,
    specs = ~ Immature_f | Species,
    var = "lon5"
  )
  
  latitude_age_difference_emm <- emmeans::contrast(
    lat_age_emm,
    method = "revpairwise"
  )
  
  longitude_age_difference_emm <- emmeans::contrast(
    lon_age_emm,
    method = "revpairwise"
  )
  
  latitude_age_differences <- as.data.frame(
    summary(
      latitude_age_difference_emm,
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  longitude_age_differences <- as.data.frame(
    summary(
      longitude_age_difference_emm,
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  latitude_species_contrasts <- as.data.frame(
    summary(
      pairs(
        update(latitude_age_difference_emm, by = NULL),
        adjust = pairwise_adjust
      ),
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  longitude_species_contrasts <- as.data.frame(
    summary(
      pairs(
        update(longitude_age_difference_emm, by = NULL),
        adjust = pairwise_adjust
      ),
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  sample_sizes <- as.data.frame(
    dplyr::count(
      comparative_df,
      Species,
      name = "N"
    )
  )
  
  list(
    data = comparative_df,
    model = full_model,
    convergence = convergence,
    sample_sizes = sample_sizes,
    lrt = lrt,
    latitude_slopes = latitude_slopes,
    longitude_slopes = longitude_slopes,
    latitude_age_differences = latitude_age_differences,
    longitude_age_differences = longitude_age_differences,
    latitude_species_contrasts = latitude_species_contrasts,
    longitude_species_contrasts = longitude_species_contrasts
  )
}


run_comparative_age_sex_analysis <- function(
    results,
    comparative_species = c(
      "Rufous",
      "Black-chinned",
      "Ruby-throated"
    ),
    reference_species = "Rufous",
    optimizer = "bobyqa",
    maxfun = 2e5,
    run_lrt = TRUE,
    pairwise_adjust = "tukey"
) {
  comparative_df <- prepare_comparative_data(
    results = results,
    comparative_species = comparative_species,
    reference_species = reference_species
  )
  
  fit_comparative_age_sex_model(
    comparative_df = comparative_df,
    optimizer = optimizer,
    maxfun = maxfun,
    run_lrt = run_lrt,
    pairwise_adjust = pairwise_adjust
  )
}


# ============================================================
# 3. TEMPORAL-CONFOUNDING ROBUSTNESS
# ============================================================

# Center latitude and longitude within species x winter year and standardize
# winter year for the temporal-confounding robustness model.
prepare_within_winter_geography <- function(comparative_df) {
  required <- c(
    "Species",
    "winter_year",
    "lat_dd",
    "lon_dd"
  )
  
  missing_cols <- setdiff(required, names(comparative_df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "`comparative_df` is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- comparative_df %>%
    dplyr::mutate(
      winter_year_num = as.numeric(as.character(winter_year))
    ) %>%
    dplyr::group_by(
      Species,
      winter_year
    ) %>%
    dplyr::mutate(
      lat5_within = (lat_dd - mean(lat_dd, na.rm = TRUE)) / 5,
      lon5_within = (lon_dd - mean(lon_dd, na.rm = TRUE)) / 5
    ) %>%
    dplyr::ungroup()
  
  year_mean <- mean(out$winter_year_num, na.rm = TRUE)
  year_sd <- stats::sd(out$winter_year_num, na.rm = TRUE)
  
  if (!is.finite(year_sd) || year_sd == 0) {
    stop("Winter year has zero or undefined standard deviation.")
  }
  
  out$year_sc <- (out$winter_year_num - year_mean) / year_sd
  
  attr(out, "winter_year_center") <- year_mean
  attr(out, "winter_year_scale") <- year_sd
  
  out
}


fit_temporal_comparative_robustness <- function(
    comparative_df,
    optimizer = "bobyqa",
    maxfun = 2e5,
    run_lrt = TRUE,
    pairwise_adjust = "tukey"
) {
  temporal_df <- prepare_within_winter_geography(comparative_df)
  
  full_model <- lme4::glmer(
    Female ~
      Species * Immature_f * lat5_within +
      Species * Immature_f * lon5_within +
      Species * Immature_f * year_sc +
      (1 | species_site_id) +
      (1 | species_winter_year) +
      (1 | species_band_id),
    data = temporal_df,
    family = binomial,
    nAGQ = 1,
    control = lme4::glmerControl(
      optimizer = optimizer,
      optCtrl = list(maxfun = maxfun)
    )
  )
  
  conv_message <- full_model@optinfo$conv$lme4$messages
  
  convergence <- data.frame(
    Convergence_message =
      if (is.null(conv_message)) {
        NA_character_
      } else {
        paste(conv_message, collapse = " | ")
      },
    Optimizer_code = full_model@optinfo$conv$opt,
    Max_abs_gradient =
      if (!is.null(full_model@optinfo$derivs$gradient)) {
        max(abs(full_model@optinfo$derivs$gradient))
      } else {
        NA_real_
      },
    Singular = lme4::isSingular(full_model, tol = 1e-4),
    row.names = NULL
  )
  
  lrt <- list(
    overall = NULL,
    latitude = NULL,
    longitude = NULL
  )
  
  if (isTRUE(run_lrt)) {
    no_lat3 <- update(
      full_model,
      . ~ . - Species:Immature_f:lat5_within
    )
    
    no_lon3 <- update(
      full_model,
      . ~ . - Species:Immature_f:lon5_within
    )
    
    no_threeway <- update(
      full_model,
      . ~ .
      - Species:Immature_f:lat5_within
      - Species:Immature_f:lon5_within
    )
    
    lrt$latitude <- anova(no_lat3, full_model, test = "Chisq")
    lrt$longitude <- anova(no_lon3, full_model, test = "Chisq")
    lrt$overall <- anova(no_threeway, full_model, test = "Chisq")
  }
  
  lat_emtrends <- emmeans::emtrends(
    full_model,
    specs = ~ Species * Immature_f,
    var = "lat5_within"
  )
  
  lon_emtrends <- emmeans::emtrends(
    full_model,
    specs = ~ Species * Immature_f,
    var = "lon5_within"
  )
  
  latitude_slopes <- as.data.frame(
    summary(lat_emtrends, infer = c(TRUE, TRUE), level = 0.95)
  )
  
  longitude_slopes <- as.data.frame(
    summary(lon_emtrends, infer = c(TRUE, TRUE), level = 0.95)
  )
  
  lat_age_emm <- emmeans::emtrends(
    full_model,
    specs = ~ Immature_f | Species,
    var = "lat5_within"
  )
  
  lon_age_emm <- emmeans::emtrends(
    full_model,
    specs = ~ Immature_f | Species,
    var = "lon5_within"
  )
  
  latitude_age_difference_emm <- emmeans::contrast(
    lat_age_emm,
    method = "revpairwise"
  )
  
  longitude_age_difference_emm <- emmeans::contrast(
    lon_age_emm,
    method = "revpairwise"
  )
  
  latitude_age_differences <- as.data.frame(
    summary(
      latitude_age_difference_emm,
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  longitude_age_differences <- as.data.frame(
    summary(
      longitude_age_difference_emm,
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  latitude_species_contrasts <- as.data.frame(
    summary(
      pairs(
        update(latitude_age_difference_emm, by = NULL),
        adjust = pairwise_adjust
      ),
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  longitude_species_contrasts <- as.data.frame(
    summary(
      pairs(
        update(longitude_age_difference_emm, by = NULL),
        adjust = pairwise_adjust
      ),
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  list(
    data = temporal_df,
    model = full_model,
    convergence = convergence,
    lrt = lrt,
    latitude_slopes = latitude_slopes,
    longitude_slopes = longitude_slopes,
    latitude_age_differences = latitude_age_differences,
    longitude_age_differences = longitude_age_differences,
    latitude_species_contrasts = latitude_species_contrasts,
    longitude_species_contrasts = longitude_species_contrasts
  )
}


# ============================================================
# 4. LEAVE-ONE-HIGH-VOLUME-LOCATION-OUT ROBUSTNESS
# ============================================================

# Identify the highest-volume reported coordinate locations and refit the
# primary model after omitting each location separately.
run_leave_one_station_out <- function(
    comparative_df,
    top_n = 4,
    station_col = "common_station_id",
    optimizer = "bobyqa",
    maxfun = 2e5,
    run_lrt = TRUE,
    pairwise_adjust = "tukey"
) {
  required <- c(
    "Female",
    "Species",
    "Immature_f",
    "lat5",
    "lon5",
    "lat_dd",
    "lon_dd",
    "species_site_id",
    "species_winter_year",
    "species_band_id",
    station_col
  )
  
  missing_cols <- setdiff(required, names(comparative_df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "`comparative_df` is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1) {
    stop("`top_n` must be a single positive integer.")
  }
  
  top_n <- as.integer(top_n)
  
  d_rank <- comparative_df
  d_rank$.station_id_internal <- as.character(d_rank[[station_col]])
  
  station_summary <- d_rank %>%
    dplyr::group_by(.station_id_internal) %>%
    dplyr::summarise(
      N = dplyr::n(),
      Percent_of_primary = 100 * dplyr::n() / nrow(d_rank),
      Mean_latitude = mean(lat_dd, na.rm = TRUE),
      Mean_longitude = mean(lon_dd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(N),
      .station_id_internal
    ) %>%
    dplyr::rename(Station = .station_id_internal)
  
  if (nrow(station_summary) == 0L) {
    stop("No reported locations were available for the influence analysis.")
  }
  
  top_n <- min(top_n, nrow(station_summary))
  top_station_summary <- station_summary[seq_len(top_n), , drop = FALSE]
  top_stations <- as.character(top_station_summary$Station)
  
  extract_lrt_row <- function(x, test_name) {
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
      Chi_square =
        if (!is.na(chisq_col)) as.numeric(last[[chisq_col]]) else NA_real_,
      df =
        if (!is.na(df_col)) as.numeric(last[[df_col]]) else NA_real_,
      p =
        if (!is.na(p_col)) as.numeric(last[[p_col]]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  fits <- setNames(
    lapply(
      top_stations,
      function(location_to_drop) {
        keep <- as.character(comparative_df[[station_col]]) != location_to_drop
        d <- comparative_df[keep, , drop = FALSE]
        
        factor_cols <- intersect(
          c(
            "Species",
            "Immature_f",
            "site_id",
            "species_site_id",
            "species_winter_year",
            "species_band_id",
            station_col
          ),
          names(d)
        )
        
        for (nm in factor_cols) {
          if (is.factor(d[[nm]])) {
            d[[nm]] <- droplevels(d[[nm]])
          }
        }
        
        fit_comparative_age_sex_model(
          comparative_df = d,
          optimizer = optimizer,
          maxfun = maxfun,
          run_lrt = run_lrt,
          pairwise_adjust = pairwise_adjust
        )
      }
    ),
    top_stations
  )
  
  lrt_summary <- dplyr::bind_rows(
    lapply(
      top_stations,
      function(st) {
        fit <- fits[[st]]
        idx <- match(st, top_station_summary$Station)
        
        dplyr::bind_rows(
          extract_lrt_row(
            fit$lrt$overall,
            "Overall species x age x geography"
          ),
          extract_lrt_row(
            fit$lrt$latitude,
            "Species x age x latitude"
          ),
          extract_lrt_row(
            fit$lrt$longitude,
            "Species x age x longitude"
          )
        ) %>%
          dplyr::mutate(
            Omitted_station = st,
            N_omitted = top_station_summary$N[idx],
            Percent_omitted = top_station_summary$Percent_of_primary[idx],
            N_remaining = nrow(fit$data),
            .before = 1
          )
      }
    )
  )
  
  age_difference_summary <- dplyr::bind_rows(
    lapply(
      top_stations,
      function(st) {
        fit <- fits[[st]]
        
        lat <- as.data.frame(fit$latitude_age_differences)
        lat$Axis <- "Latitude"
        
        lon <- as.data.frame(fit$longitude_age_differences)
        lon$Axis <- "Longitude"
        
        dplyr::bind_rows(lat, lon) %>%
          dplyr::mutate(Omitted_station = st) %>%
          dplyr::relocate(Omitted_station, Axis)
      }
    )
  )
  
  species_contrast_summary <- dplyr::bind_rows(
    lapply(
      top_stations,
      function(st) {
        fit <- fits[[st]]
        
        lat <- as.data.frame(fit$latitude_species_contrasts)
        lat$Axis <- "Latitude"
        
        lon <- as.data.frame(fit$longitude_species_contrasts)
        lon$Axis <- "Longitude"
        
        dplyr::bind_rows(lat, lon) %>%
          dplyr::mutate(Omitted_station = st) %>%
          dplyr::relocate(Omitted_station, Axis)
      }
    )
  )
  
  list(
    top_station_summary = top_station_summary,
    lrt_summary = lrt_summary,
    age_difference_summary = age_difference_summary,
    species_contrast_summary = species_contrast_summary
  )
}
