# ============================================================
# 02_primary_analysis.R
# Build the analysis datasets, fit the primary comparative GLMM, summarize
# the structure of the primary dataset, and calculate response-scale effect sizes.
#
# Primary comparative model:
#   Female ~
#     Species * Immature_f * lat5 +
#     Species * Immature_f * lon5 +
#     (1 | species_site_id) +
#     (1 | species_winter_year) +
#     (1 | species_band_id)
#
# Geographic scaling:
#   lat5 = (latitude - 31) / 5
#   lon5 = (longitude + 90) / 5
# ============================================================


# ============================================================
# 0. LOAD FUNCTIONS
# ============================================================

source("01_functions.R")


# ============================================================
# 1. ANALYSIS SETTINGS
# ============================================================

# Raw USGS Bird Banding Laboratory data file.
csv_file <- "NABBP_2025_grp_13.csv"

# Primary analytical window: 15 November through 31 December.
months_keep <- c(11L, 12L)
primary_start_month <- 11L
primary_start_day <- 15L
primary_end_month <- 12L
primary_end_day <- 31L

countries_keep <- c("US", "CA")

# Eastern North America: at or east of 100 degrees W.
lon_min <- -100
lon_max <- Inf
lat_min <- -Inf
lat_max <- Inf


# ============================================================
# 2. SPECIES DEFINITIONS
# ============================================================

species_ids <- c(
  "Rufous"        = 4330,
  "Black-chinned" = 4290,
  "Ruby-throated" = 4280,
  "Calliope"      = 4360,
  "Broad-tailed"  = 4320,
  "Allen's"       = 4340
)

# Focal species included in the comparative GLMM.
focal_species <- c(
  "Rufous",
  "Black-chinned",
  "Ruby-throated"
)

# Additional species included in descriptive demographic summaries.
descriptive_species <- c(
  "Calliope",
  "Broad-tailed",
  "Allen's"
)

species_order <- c(
  focal_species,
  descriptive_species
)


# ============================================================
# 3. BUILD SPECIES DATASETS
# ============================================================

run_species <- function(species_id) {
  winter_age_sex_pipeline(
    csv_path = csv_file,
    species_id = species_id,
    months_keep = months_keep,
    start_month = primary_start_month,
    start_day = primary_start_day,
    end_month = primary_end_month,
    end_day = primary_end_day,
    countries_keep = countries_keep,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max
  )
}

results <- purrr::map(
  species_ids,
  run_species
)

# Verify that all expected species are present and contain at least one
# qualifying record after filtering.
stopifnot(
  identical(names(results), species_order),
  all(
    vapply(
      results,
      function(x) nrow(x$data_points) > 0L,
      logical(1)
    )
  )
)


# ============================================================
# 4. PRIMARY THREE-SPECIES COMPARATIVE GLMM
# ============================================================

comparative_primary <- run_comparative_age_sex_analysis(
  results = results,
  comparative_species = focal_species,
  reference_species = "Rufous",
  optimizer = "bobyqa",
  maxfun = 2e5,
  run_lrt = TRUE,
  pairwise_adjust = "tukey"
)

# Store the fitted model and comparative dataset under short names used by
# the sensitivity, figure, and supplement scripts.
comparative_df <- comparative_primary$data
m_comparative_full <- comparative_primary$model


# ============================================================
# 5. PRIMARY DATASET STRUCTURE
# ============================================================

# Summarize the number of retained records, unique individuals, repeated
# individuals across winters, reported locations, and winter-year coverage.
audit_df <- data.frame(
  Species = as.character(comparative_df$Species),
  band = as.character(comparative_df$band),
  winter_year = as.integer(comparative_df$winter_year),
  reported_location = as.character(comparative_df$common_station_id),
  stringsAsFactors = FALSE
)

audit_df$individual_id <- paste(
  audit_df$Species,
  audit_df$band,
  sep = "__"
)

summarize_primary_dataset <- function(d, label) {
  n_individuals <- length(unique(d$individual_id))
  
  data.frame(
    Species = label,
    Records = nrow(d),
    Individuals = n_individuals,
    Additional_records = nrow(d) - n_individuals,
    Individuals_in_multiple_winters = sum(table(d$individual_id) > 1L),
    Reported_locations = length(unique(d$reported_location)),
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
      summarize_primary_dataset(
        audit_df[audit_df$Species == sp, , drop = FALSE],
        sp
      )
    }
  ),
  summarize_primary_dataset(
    audit_df,
    "All species"
  )
)


# ============================================================
# 6. PRIMARY RECORD-TYPE COMPOSITION
# ============================================================

# Reconstruct the retained focal-species records with the raw record-type field
# so the contribution of original bandings and subsequent encounters can be
# summarized for the primary dataset.
raw_header <- names(
  data.table::fread(
    csv_file,
    nrows = 0L,
    showProgress = FALSE
  )
)

record_type_candidates <- c(
  "event_type",
  "record_source",
  "record_type",
  "band_type",
  "encounter_type"
)

record_type_col <- intersect(
  record_type_candidates,
  raw_header
)[1L]

if (length(record_type_col) == 0L || is.na(record_type_col)) {
  stop(
    "No recognized record-type column was found in the BBL extract."
  )
}

raw_record_data <- data.table::fread(
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
    "age_code",
    record_type_col
  ),
  showProgress = FALSE
)

raw_record_data[, event_year := as.integer(event_year)]
raw_record_data[, event_month := as.integer(event_month)]
raw_record_data[, event_day := as.integer(event_day)]

primary_date_keep <- in_month_day_window(
  month = raw_record_data$event_month,
  day = raw_record_data$event_day,
  start_month = primary_start_month,
  start_day = primary_start_day,
  end_month = primary_end_month,
  end_day = primary_end_day
)

raw_record_data <- raw_record_data[
  species_id %in% unname(species_ids[focal_species]) &
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

raw_record_data[, Species := names(species_ids)[
  match(
    species_id,
    unname(species_ids)
  )
]]

# Missing event days are used only to order records after date-window filtering.
raw_record_data[is.na(event_day), event_day := 15L]

raw_record_data[, dt := data.table::as.IDate(
  sprintf(
    "%04d-%02d-%02d",
    event_year,
    event_month,
    event_day
  )
)]

raw_record_data <- raw_record_data[
  !is.na(dt)
]

primary_year_fun <- infer_season_year_fun(
  months_keep
)

raw_record_data[, winter_year :=
                  primary_year_fun(
                    event_year,
                    event_month
                  )]

raw_record_data[, Female := data.table::fifelse(
  sex_code %in% c(5L, 7L),
  1L,
  data.table::fifelse(
    sex_code %in% c(4L, 6L),
    0L,
    NA_integer_
  )
)]

raw_record_data[, Immature :=
                  classify_immature(
                    event_month,
                    age_code
                  )]

data.table::setorder(
  raw_record_data,
  Species,
  band,
  winter_year,
  dt
)

selected_primary_records <- raw_record_data[, {
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

selected_primary_records <- selected_primary_records[
  !is.na(Female) &
    !is.na(Immature)
]

# Confirm that the independently reconstructed records are identical to the
# records used in the primary comparative model.
reconstructed_primary <- selected_primary_records[, .(
  Species = as.character(Species),
  band = as.character(band),
  winter_year = as.integer(winter_year),
  lat_dd = as.numeric(lat_dd),
  lon_dd = as.numeric(lon_dd),
  Female = as.integer(Female),
  Immature = as.integer(Immature)
)]

primary_model_records <- comparative_df %>%
  dplyr::transmute(
    Species = as.character(Species),
    band = as.character(band),
    winter_year = as.integer(winter_year),
    lat_dd = as.numeric(lat_dd),
    lon_dd = as.numeric(lon_dd),
    Female = as.integer(Female),
    Immature = as.integer(Immature)
  )

reconstructed_primary <- reconstructed_primary[
  do.call(
    order,
    reconstructed_primary
  ),
]

primary_model_records <- primary_model_records[
  do.call(
    order,
    primary_model_records
  ),
]

rownames(primary_model_records) <- NULL

if (!isTRUE(
  all.equal(
    as.data.frame(reconstructed_primary),
    as.data.frame(primary_model_records),
    tolerance = 1e-8,
    check.attributes = FALSE
  )
)) {
  stop(
    "The record-type reconstruction does not match the primary comparative dataset."
  )
}

record_type_values <- toupper(
  trimws(
    as.character(
      selected_primary_records[[record_type_col]]
    )
  )
)

record_type_labels <- dplyr::case_when(
  record_type_values %in% c("B", "BANDING") ~ "Banding",
  record_type_values %in% c("E", "ENCOUNTER") ~ "Encounter",
  TRUE ~ NA_character_
)

if (anyNA(record_type_labels)) {
  stop(
    "Unrecognized record-type value(s) in the retained primary dataset: ",
    paste(
      sort(
        unique(
          record_type_values[
            is.na(record_type_labels)
          ]
        )
      ),
      collapse = ", "
    )
  )
}

record_type_summary <- data.frame(
  Record_type = record_type_labels,
  stringsAsFactors = FALSE
) %>%
  dplyr::count(
    Record_type,
    name = "Records"
  ) %>%
  dplyr::mutate(
    Percent = 100 * Records / sum(Records)
  )


# ============================================================
# 7. RESPONSE-SCALE EFFECT SIZES
# ============================================================

# Calculate predicted female probabilities at the 10th, 50th, and 90th
# percentiles of observed latitude and longitude for each focal species and
# age class. Predictions exclude random effects and average over the observed
# distribution of the coordinate not being varied.
predict_along_axis <- function(
    species_name,
    axis = c(
      "latitude",
      "longitude"
    )
) {
  axis <- match.arg(axis)
  
  base <- comparative_df %>%
    dplyr::filter(
      Species == species_name
    )
  
  if (nrow(base) == 0L) {
    stop(
      "No primary-analysis records found for ",
      species_name,
      "."
    )
  }
  
  values <- if (axis == "latitude") {
    stats::quantile(
      base$lat_dd,
      probs = c(
        0.10,
        0.50,
        0.90
      ),
      na.rm = TRUE
    )
  } else {
    stats::quantile(
      base$lon_dd,
      probs = c(
        0.10,
        0.50,
        0.90
      ),
      na.rm = TRUE
    )
  }
  
  grid <- expand.grid(
    Age = c(
      "Adult",
      "Immature"
    ),
    Value = as.numeric(values),
    stringsAsFactors = FALSE
  )
  
  grid$Quantile <- rep(
    names(values),
    each = 2L
  )
  
  out <- lapply(
    seq_len(
      nrow(grid)
    ),
    function(i) {
      newdata <- base
      
      if (axis == "latitude") {
        newdata$lat_dd <- grid$Value[i]
        newdata$lat5 <-
          (grid$Value[i] - 31) / 5
      } else {
        newdata$lon_dd <- grid$Value[i]
        newdata$lon5 <-
          (grid$Value[i] + 90) / 5
      }
      
      newdata$Immature_f <- factor(
        grid$Age[i],
        levels = levels(
          comparative_df$Immature_f
        )
      )
      
      newdata$Immature <-
        as.integer(
          grid$Age[i] == "Immature"
        )
      
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
        predict_along_axis(
          sp,
          "latitude"
        ),
        predict_along_axis(
          sp,
          "longitude"
        )
      )
    }
  )
)
