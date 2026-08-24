# ============================================================
# 02_primary_analysis.R
# Build the analysis datasets, fit the primary comparative GLMM, and generate
# the standardized Rufous latitude predictions reported in the manuscript.
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

# Raw USGS bird-banding data file. Update this path if the file is stored
# elsewhere.
csv_file <- "NABBP_2025_grp_13.csv"

# Primary analytical window: 15 November through 31 December.
months_keep <- c(11L, 12L)
primary_start_month <- 11L
primary_start_day <- 15L
primary_end_month <- 12L
primary_end_day <- 31L

countries_keep <- c("US", "CA")

# Eastern North America: east of 100 degrees W.
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
# 5. STANDARDIZED RUFOUS LATITUDE PREDICTIONS FOR THE ABSTRACT
#
# Generate population-level predictions at 30 and 36 degrees N, averaging
# response-scale predictions across the observed Rufous longitude distribution.
# Random effects are excluded.
# ============================================================

rufous_prediction_base <- comparative_df %>%
  dplyr::filter(
    Species == "Rufous"
  )

rufous_latitude_predictions <- dplyr::bind_rows(
  lapply(
    c("Adult", "Immature"),
    function(age_class) {
      dplyr::bind_rows(
        lapply(
          c(30, 36),
          function(latitude) {
            newdata <- rufous_prediction_base
            newdata$lat_dd <- latitude
            newdata$lat5 <- (latitude - 31) / 5
            newdata$Immature_f <- factor(
              age_class,
              levels = levels(comparative_df$Immature_f)
            )
            newdata$Immature <- as.integer(
              age_class == "Immature"
            )
            
            predicted_probability <- predict(
              m_comparative_full,
              newdata = newdata,
              type = "response",
              re.form = NA
            )
            
            data.frame(
              Species = "Rufous",
              Age = age_class,
              Latitude = latitude,
              Predicted_female_probability =
                mean(predicted_probability),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)
