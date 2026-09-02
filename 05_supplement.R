# ============================================================
# 05_supplement.R
#
# Generate Supplementary Material Tables S1-S10 for the hummingbird
# demographic and geographic analyses. Run 01_functions.R,
# 02_primary_analysis.R, and 03_sensitivity_analysis.R before this script.
#
# Supplementary tables:
#   S1  Demographic composition
#   S2  Supporting-species demographic uncertainty
#   S3  Primary age-specific geographic slopes
#   S4  Geographic sampling support by latitude band
#   S5  Primary pairwise interspecific contrasts
#   S6  Primary comparative GLMM fixed and random effects
#   S7  GAMM nonlinear-spatial sensitivity
#   S8  Temporal-confounding robustness analysis
#   S9  Leave-one-high-volume-location-out sensitivity
#   S10 Sampling-date and western-boundary robustness analyses
#
# Objects named table_s#, caption_s#, and ft_s# correspond directly to the
# final Supplementary Material table number used in the manuscript.
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
  "sensitivity_results",
  "missingness_summary_all",
  "within_year_variance",
  "date_covariate_lrt",
  "date_covariate_differences",
  "longitude_boundary_lrt",
  "longitude_boundary_differences"
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

lat_band_width <- 4
lat_band_origin <- 22

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
    "GAMM specification used in the analysis."
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
    "comparative dataset. Recreate the GAMM cache before generating the supplement."
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
      ifelse(
        x < 0.01,
        sprintf("%.3f", x),
        sprintf("%.2f", x)
      )
    )
  )
}


# Numerical precision follows the final manuscript throughout the supplement:
#   - model estimates, SEs, confidence limits, z statistics, variances, and SDs: 2 decimals
#   - chi-square statistics and reference df: 1 decimal
#   - EDF and k-index: 2 decimals
#   - tabulated percentages: 1 decimal
#   - P values: 2 decimals unless raw P < 0.01, then 3 decimals; P < 0.001 as the floor
# A true minus sign is used for negative displayed values.
fmt_num <- function(x, digits = 2) {
  # Prevent display artifacts such as −0.00 after rounding.
  threshold <- 0.5 * 10^(-digits)
  x_clean <- ifelse(!is.na(x) & abs(x) < threshold, 0, x)
  out <- ifelse(
    is.na(x_clean),
    "",
    sprintf(paste0("%.", digits, "f"), x_clean)
  )
  sub("^-", "−", out)
}


fmt_ci <- function(lower, upper, digits = 2) {
  paste0(
    fmt_num(lower, digits),
    " to ",
    fmt_num(upper, digits)
  )
}


fmt_count <- function(x) {
  ifelse(
    is.na(x),
    "",
    format(
      as.integer(x),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  )
}


fmt_model_rank <- function(x) {
  vapply(
    strsplit(as.character(x), "/", fixed = TRUE),
    function(parts) {
      if (length(parts) != 2L || !all(grepl("^[0-9]+$", parts))) {
        return(paste(parts, collapse = "/"))
      }
      paste(fmt_count(as.integer(parts)), collapse = "/")
    },
    character(1)
  )
}


nice_ft <- function(x) {
  font_size <- 12
  
  ft <- flextable(x) %>%
    border_remove() %>%
    font(
      fontname = "Times New Roman",
      part = "all"
    ) %>%
    fontsize(
      size = font_size,
      part = "all"
    ) %>%
    bold(
      part = "header"
    ) %>%
    align(
      align = "left",
      part = "all"
    ) %>%
    valign(
      valign = "top",
      part = "all"
    ) %>%
    line_spacing(
      space = 1,
      part = "all"
    ) %>%
    padding(
      padding = 1,
      part = "all"
    ) %>%
    set_table_properties(
      layout = "autofit",
      width = 1,
      align = "left",
      opts_word = list(
        split = FALSE,
        keep_with_next = FALSE,
        repeat_headers = TRUE
      )
    )
  
  # AOS uses italic statistical symbols. Keep column keys unchanged for the
  # data workflow but display the appropriate symbols in the header.
  if ("p" %in% names(x)) {
    ft <- compose(
      ft,
      j = "p",
      part = "header",
      value = as_paragraph(as_i("P"))
    )
  }
  
  if ("z" %in% names(x)) {
    ft <- compose(
      ft,
      j = "z",
      part = "header",
      value = as_paragraph(as_i("z"))
    )
  }
  
  if ("n" %in% names(x)) {
    ft <- compose(
      ft,
      j = "n",
      part = "header",
      value = as_paragraph(as_i("n"))
    )
  }
  
  ft
}


# Set common typography, single spacing, and alignment for all supplementary tables.
flextable::set_flextable_defaults(
  font.family = "Times New Roman",
  hansi.family = "Times New Roman",
  cs.family = "Times New Roman",
  eastasia.family = "Times New Roman",
  font.size = 12,
  text.align = "left",
  padding = 1,
  line_spacing = 1,
  table.layout = "autofit"
)


# Format table captions and panel headings.
caption_text_fp <- fp_text(
  font.family = "Times New Roman",
  font.size = 12,
  bold = TRUE,
  italic = FALSE
)

caption_stat_text_fp <- fp_text(
  font.family = "Times New Roman",
  font.size = 12,
  bold = TRUE,
  italic = TRUE
)

caption_fp <- fp_par(
  text.align = "left",
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
    `χ²` = as.numeric(a$Chisq[i]),
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

table_s1_demographic_composition <- map_dfr(
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
      Species = ifelse(sp == "Allen's", "Allen’s", sp),
      n = fmt_count(n_total),
      `Immature, n (%)` =
        sprintf(
          "%s (%.1f)",
          fmt_count(n_immature),
          100 * n_immature / n_total
        ),
      `Female overall, n (%)` =
        sprintf(
          "%s (%.1f)",
          fmt_count(n_female),
          100 * n_female / n_total
        ),
      `Adult female, n (%)` =
        sprintf(
          "%s (%.1f)",
          fmt_count(n_adult_female),
          100 * n_adult_female / n_adult
        ),
      `Immature female, n (%)` =
        sprintf(
          "%s (%.1f)",
          fmt_count(n_immature_female),
          100 * n_immature_female / n_immature
        )
    )
  }
)

caption_s1_demographic_composition <- paste0(
  "Supplementary Material Table S1. Demographic composition of hummingbird records ",
  "for the six study species. Only records with classifiable age and known sex are included. ",
  "Female percentages are shown overall and within age classes."
)

ft_s1_demographic_composition <- nice_ft(table_s1_demographic_composition)


# ============================================================
# TABLE S2. Supporting-species demographic uncertainty
# ============================================================

supporting_species_age_bounds <- map_dfr(
  descriptive_species,
  function(sp) {
    as.data.frame(
      sensitivity_results[[sp]]$age_bounds
    ) %>%
      filter(Group == "All") %>%
      mutate(Species = sp)
  }
)

supporting_species_sex_bounds <- map_dfr(
  descriptive_species,
  function(sp) {
    as.data.frame(
      sensitivity_results[[sp]]$sex_ratio_bounds
    ) %>%
      mutate(Species = sp)
  }
)

supporting_species_missingness <- as.data.frame(
  missingness_summary_all
) %>%
  filter(
    Species %in% descriptive_species
  ) %>%
  transmute(
    Species = as.character(Species),
    `Selected bird-years` = fmt_count(Total_selected_bird_years),
    `Unknown sex, n (%)` =
      sprintf(
        "%s (%.1f)",
        fmt_count(Unknown_or_missing_sex_N),
        Unknown_or_missing_sex_pct
      ),
    `Unclassifiable age, n (%)` =
      sprintf(
        "%s (%.1f)",
        fmt_count(Unclassifiable_age_N),
        Unclassifiable_age_pct
      )
  )

supporting_species_age_bounds_compact <- supporting_species_age_bounds %>%
  transmute(
    Species = as.character(Species),
    `Immature, observed (range %)` =
      sprintf(
        "%.1f (%.1f–%.1f)",
        100 * Observed_immature_prop,
        100 * Minimum_immature_prop,
        100 * Maximum_immature_prop
      )
  )

supporting_species_sex_bounds_compact <- supporting_species_sex_bounds %>%
  mutate(
    Age = as.character(Age),
    Sex_bound = sprintf(
      "%.1f (%.1f–%.1f)",
      100 * Observed_female_prop,
      100 * Minimum_female_prop,
      100 * Maximum_female_prop
    )
  ) %>%
  select(
    Species,
    Age,
    Sex_bound
  ) %>%
  pivot_wider(
    names_from = Age,
    values_from = Sex_bound,
    names_glue = "{Age} female, observed (range %)"
  )

required_sex_bound_cols <- c(
  "Adult female, observed (range %)",
  "Immature female, observed (range %)"
)

stopifnot(
  all(required_sex_bound_cols %in% names(supporting_species_sex_bounds_compact))
)

table_s2_supporting_species_uncertainty <- supporting_species_missingness %>%
  left_join(
    supporting_species_age_bounds_compact,
    by = "Species"
  ) %>%
  left_join(
    supporting_species_sex_bounds_compact,
    by = "Species"
  ) %>%
  mutate(
    Species = factor(
      Species,
      levels = descriptive_species
    )
  ) %>%
  arrange(Species) %>%
  mutate(
    Species = as.character(Species),
    Species = ifelse(Species == "Allen's", "Allen’s", Species)
  )

stopifnot(nrow(table_s2_supporting_species_uncertainty) == length(descriptive_species))

caption_s2_supporting_species_uncertainty <- paste0(
  "Supplementary Material Table S2. Missing demographic information and extreme-allocation ",
  "bounds for Calliope, Broad-tailed, and Allen’s Hummingbirds. Ranges show the minimum and ",
  "maximum percentages under extreme allocation of records with unclassifiable age or unknown sex."
)

ft_s2_supporting_species_uncertainty <- nice_ft(table_s2_supporting_species_uncertainty)


# ============================================================
# TABLE S3. Primary age-specific geographic slopes
# ============================================================

table_s3_age_specific_slopes <- bind_rows(
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
    Estimate = fmt_num(Estimate, 2),
    SE = fmt_num(SE, 2),
    `95% CI` = fmt_ci(CI_low, CI_high, 2),
    z = fmt_num(z, 2),
    p = fmt_p(p)
  )

caption_s3_age_specific_slopes <- paste0(
  "Supplementary Material Table S3. Age-specific latitude and longitude slopes from the ",
  "primary comparative model. Estimates are changes in log-odds of a record being female ",
  "per 5° geographic change."
)

ft_s3_age_specific_slopes <- nice_ft(table_s3_age_specific_slopes)


# ============================================================
# TABLE S4. Geographic sampling support by 4° latitude band
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

# Use 4° latitude bands with 34°N as a band boundary.
latitude_support_data <- comparative_df %>%
  mutate(
    Species = as.character(Species),
    Age = as.character(Immature_f),
    Latitude_band_low =
      floor((lat_dd - lat_band_origin) / lat_band_width) *
      lat_band_width + lat_band_origin,
    Latitude_band = paste0(
      Latitude_band_low,
      "–<",
      Latitude_band_low + lat_band_width,
      "°N"
    )
  )

latitude_band_lower_bounds <- seq(
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
    Latitude_band_low = latitude_band_lower_bounds,
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

table_s4_geographic_sampling_support <- latitude_support_long %>%
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
    `Adult records` = fmt_count(Adult_Records),
    `Adult reported locations` =
      fmt_count(Adult_Reported_locations),
    `Immature records` = fmt_count(Immature_Records),
    `Immature reported locations` =
      fmt_count(Immature_Reported_locations)
  ) %>%
  arrange(
    factor(
      Species,
      levels = focal_species_order
    ),
    Latitude_band_low
  )

# Verify that latitude-band counts reproduce the primary-model sample size.
sampling_support_counts_from_bands <- latitude_support_long %>%
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

sampling_support_counts_direct <- comparative_df %>%
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
      sampling_support_counts_from_bands,
      sampling_support_counts_direct,
      check.attributes = FALSE
    )
  ),
  sum(latitude_support_long$Records) ==
    nrow(comparative_df)
)

table_s4_geographic_sampling_support <- table_s4_geographic_sampling_support %>%
  select(-Latitude_band_low)

caption_s4_geographic_sampling_support <- paste0(
  "Supplementary Material Table S4. Geographic sampling support for adult and immature ",
  "Rufous, Black-chinned, and Ruby-throated Hummingbirds by 4° latitude band."
)

ft_s4_geographic_sampling_support <- nice_ft(table_s4_geographic_sampling_support)


# ============================================================
# TABLE S5. Primary pairwise interspecific contrasts
# ============================================================

table_s5_interspecific_contrasts <- bind_rows(
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
    Estimate = fmt_num(Estimate, 2),
    SE = fmt_num(SE, 2),
    `95% CI` = fmt_ci(CI_low, CI_high, 2),
    z = fmt_num(z, 2),
    p = fmt_p(p)
  )

caption_s5_interspecific_contrasts <- paste0(
  "Supplementary Material Table S5. Tukey-adjusted pairwise comparisons of age-dependent ",
  "latitude and longitude slopes among species. Estimates are the immature–adult slope ",
  "difference in the first species minus that in the second species."
)

ft_s5_interspecific_contrasts <- nice_ft(table_s5_interspecific_contrasts)


# ============================================================
# TABLE S6. Primary comparative GLMM fixed and random effects
# ============================================================

table_s6_fixed_effects <- broom.mixed::tidy(
  m_comparative_full,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
) %>%
  transmute(
    Term = clean_fixed_term(term),
    Estimate = fmt_num(estimate, 2),
    SE = fmt_num(std.error, 2),
    `95% CI` = fmt_ci(conf.low, conf.high, 2),
    z = fmt_num(statistic, 2),
    p = fmt_p(p.value)
  )

primary_glmm_random_effect_levels <- c(
  species_site_id =
    dplyr::n_distinct(comparative_df$species_site_id),
  species_winter_year =
    dplyr::n_distinct(comparative_df$species_winter_year),
  species_band_id =
    dplyr::n_distinct(comparative_df$species_band_id)
)

table_s6_random_effects <- as.data.frame(
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
    Levels = fmt_count(primary_glmm_random_effect_levels[grp]),
    Variance = fmt_num(vcov, 2),
    SD = fmt_num(sdcor, 2)
  )

caption_s6_primary_glmm <- paste0(
  "Supplementary Material Table S6. Fixed-effect coefficients and random-effect variance ",
  "components from the primary comparative generalized linear mixed model. Rufous Hummingbird ",
  "and adults are reference levels; latitude and longitude are scaled per 5° and centered at ",
  "31°N and 90°W."
)

ft_s6_fixed_effects <- nice_ft(table_s6_fixed_effects)
ft_s6_random_effects <- nice_ft(table_s6_random_effects)


# ============================================================
# TABLE S7. GAMM nonlinear-spatial sensitivity
# ============================================================

table_s7_gamm_joined <- gamm_sensitivity_summary %>%
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

stopifnot(
  nrow(table_s7_gamm_joined) == 6L,
  !anyNA(
    table_s7_gamm_joined[
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

table_s7_gamm_summary <- table_s7_gamm_joined %>%
  transmute(
    Species,
    n = fmt_count(N),
    Smooth = dplyr::recode(
      Smooth,
      "Immature - Adult spatial difference" =
        "Immature–Adult spatial difference",
      .default = Smooth
    ),
    EDF = fmt_num(edf, 2),
    `Reference df` = fmt_num(Ref_df, 1),
    `χ²` = fmt_num(Chi_square, 1),
    p = fmt_p(p),
    `k-index` = ifelse(
      is.na(k_index),
      "",
      fmt_num(k_index, 2)
    ),
    `k-check p` = fmt_p(k_p),
    `Model rank` = fmt_model_rank(Model_rank)
  )

stopifnot(
  nrow(table_s7_gamm_summary) == 6L,
  !anyNA(table_s7_gamm_summary$Species)
)

caption_s7_gamm_sensitivity <- paste0(
  "Supplementary Material Table S7. Species-specific generalized additive mixed-model ",
  "(GAMM) sensitivity analysis allowing nonlinear two-dimensional geographic sex-ratio surfaces. ",
  "Adult spatial surfaces are the baseline; immature–adult difference smooths test whether ",
  "geographic surfaces differ between age classes. k-index and k-check probability values assess ",
  "basis dimension."
)

caption_s7_gamm_sensitivity_runs <- list(
  ftext(
    paste0(
      "Supplementary Material Table S7. Species-specific generalized additive mixed-model ",
      "(GAMM) sensitivity analysis allowing nonlinear two-dimensional geographic sex-ratio surfaces. ",
      "Adult spatial surfaces are the baseline; immature–adult difference smooths test whether ",
      "geographic surfaces differ between age classes. "
    ),
    prop = caption_text_fp
  ),
  ftext("k", prop = caption_stat_text_fp),
  ftext("-index and ", prop = caption_text_fp),
  ftext("k", prop = caption_stat_text_fp),
  ftext("-check probability values assess basis dimension.", prop = caption_text_fp)
)

ft_s7_gamm_summary <- nice_ft(table_s7_gamm_summary)


# ============================================================
# TABLE S8. Temporal-confounding robustness analysis
# ============================================================

# Panel A: omnibus tests.
table_s8_temporal_lrt <- bind_rows(
  extract_lrt(
    comparative_temporal$lrt$overall,
    "Species × age × geography",
    "Joint"
  ),
  extract_lrt(
    comparative_temporal$lrt$latitude,
    "Species × age × latitude",
    "Latitude"
  ),
  extract_lrt(
    comparative_temporal$lrt$longitude,
    "Species × age × longitude",
    "Longitude"
  )
) %>%
  transmute(
    Test,
    `χ²` = fmt_num(`χ²`, 1),
    df,
    p = fmt_p(p)
  )

# Panel B: within-species age differences and longitude species contrasts.
table_s8_temporal_effects <- bind_rows(
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
  # Latitude species contrasts are omitted because the corresponding
  # species × age × latitude interaction is unsupported.
  standardize_emm(
    comparative_temporal$longitude_species_contrasts,
    "Longitude",
    "Species contrast"
  )
) %>%
  mutate(
    Comparison = case_when(
      Result == "Species contrast" ~ clean_species_contrast(Contrast),
      Result == "Immature – Adult" ~ Species,
      TRUE ~ Contrast
    )
  ) %>%
  transmute(
    Axis,
    Comparison,
    Estimate = fmt_num(Estimate, 2),
    SE = fmt_num(SE, 2),
    `95% CI` = fmt_ci(CI_low, CI_high, 2),
    z = fmt_num(z, 2),
    p = fmt_p(p)
  )

# Panel C: geographic variance retained after within-year centering.
table_s8_temporal_variance <- as.data.frame(
  within_year_variance
) %>%
  mutate(
    Species = factor(
      Species,
      levels = focal_species_order
    )
  ) %>%
  arrange(Species) %>%
  transmute(
    Species = as.character(Species),
    `Latitude variance retained (%)` =
      fmt_num(Latitude_variance_retained_pct, 1),
    `Longitude variance retained (%)` =
      fmt_num(Longitude_variance_retained_pct, 1)
  )

caption_s8_temporal_confounding <- paste0(
  "Supplementary Material Table S8. Temporal-confounding robustness analysis using coordinates ",
  "centered within species and year. Panel A shows omnibus likelihood-ratio tests; Panel B shows ",
  "within-species immature–adult slope differences for single-species rows and Tukey-adjusted ",
  "longitude contrasts for species-pair rows; Panel C shows geographic variance retained after ",
  "centering. For species-pair rows, estimates are the immature–adult slope difference in the first ",
  "species minus that in the second species. Within-species contrasts are not multiplicity-adjusted."
)

ft_s8_temporal_lrt <- nice_ft(table_s8_temporal_lrt)
ft_s8_temporal_effects <- nice_ft(table_s8_temporal_effects)
ft_s8_temporal_variance <- nice_ft(table_s8_temporal_variance)


# ============================================================
# TABLE S9. Leave-one-high-volume-location-out sensitivity
# ============================================================

leave_one_location_lookup <- as.data.frame(
  station_influence_top_sites
) %>%
  transmute(
    Omitted_station = as.character(Station),
    `Omitted location` =
      format_location(
        Mean_latitude,
        Mean_longitude
      ),
    `n omitted` = as.integer(N),
    `% primary omitted` = fmt_num(Percent_of_primary, 1)
  )

# Omnibus tests for each omitted location.
leave_one_location_lrt_compact <- as.data.frame(
  station_influence_lrt_table
) %>%
  transmute(
    Omitted_station = as.character(Omitted_station),
    Test = case_when(
      grepl("geography", Test, ignore.case = TRUE) ~ "Joint P",
      grepl("latitude", Test, ignore.case = TRUE) ~ "Latitude P",
      grepl("longitude", Test, ignore.case = TRUE) ~ "Longitude P",
      TRUE ~ NA_character_
    ),
    p = fmt_p(p)
  ) %>%
  filter(!is.na(Test)) %>%
  distinct(Omitted_station, Test, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = Test,
    values_from = p
  )

# Within-species age differences for each omitted location.
leave_one_location_age_compact <- as.data.frame(
  station_influence_age_differences
) %>%
  transmute(
    Omitted_station = as.character(Omitted_station),
    Species = as.character(Species),
    Axis = as.character(Axis),
    Result = paste0(
      fmt_num(estimate, 2),
      " (",
      fmt_p(p.value),
      ")"
    )
  ) %>%
  filter(
    Species %in% focal_species_order,
    Axis %in% c("Latitude", "Longitude")
  ) %>%
  mutate(
    Column = paste(Species, Axis, sep = "__")
  ) %>%
  select(
    Omitted_station,
    Column,
    Result
  ) %>%
  distinct(Omitted_station, Column, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = Column,
    values_from = Result
  )

required_leave_one_location_columns <- c(
  "Rufous__Latitude",
  "Rufous__Longitude",
  "Black-chinned__Latitude",
  "Black-chinned__Longitude",
  "Ruby-throated__Latitude",
  "Ruby-throated__Longitude"
)

stopifnot(
  all(c("Joint P", "Latitude P", "Longitude P") %in% names(leave_one_location_lrt_compact)),
  all(required_leave_one_location_columns %in% names(leave_one_location_age_compact))
)

table_s9_leave_one_location_out <- leave_one_location_lookup %>%
  left_join(
    leave_one_location_lrt_compact,
    by = "Omitted_station"
  ) %>%
  left_join(
    leave_one_location_age_compact,
    by = "Omitted_station"
  ) %>%
  transmute(
    `Omitted location`,
    `n omitted` = fmt_count(`n omitted`),
    `% omitted` = `% primary omitted`,
    `Joint P`,
    `Latitude P`,
    `Longitude P`,
    `Rufous lat. β (P)` = `Rufous__Latitude`,
    `Rufous lon. β (P)` = `Rufous__Longitude`,
    `Black-chinned lat. β (P)` = `Black-chinned__Latitude`,
    `Black-chinned lon. β (P)` = `Black-chinned__Longitude`,
    `Ruby-throated lat. β (P)` = `Ruby-throated__Latitude`,
    `Ruby-throated lon. β (P)` = `Ruby-throated__Longitude`
  )

stopifnot(nrow(table_s9_leave_one_location_out) == 4L)

caption_s9_leave_one_location_out <- paste0(
  "Supplementary Material Table S9. Leave-one-high-volume-location-out sensitivity analysis. ",
  "The four highest-volume reported locations were omitted individually. Joint, latitude, and ",
  "longitude columns show likelihood-ratio test probability values; species-specific cells show ",
  "immature–adult slope differences with probability values in parentheses."
)

ft_s9_leave_one_location_out <- nice_ft(table_s9_leave_one_location_out)


# ============================================================
# TABLE S10. Sampling-date and 97°W western-boundary robustness analyses
# ============================================================

# Panel A: joint, latitude, and longitude interaction tests under each
# robustness check.
date_covariate_lrt_table <- as.data.frame(
  date_covariate_lrt
) %>%
  filter(
    grepl("geography|latitude|longitude", Test, ignore.case = TRUE)
  ) %>%
  transmute(
    Analysis = "Day-of-window covariate",
    Test = case_when(
      grepl("geography", Test, ignore.case = TRUE) ~
        "Species × age × geography",
      grepl("latitude", Test, ignore.case = TRUE) ~
        "Species × age × latitude",
      grepl("longitude", Test, ignore.case = TRUE) ~
        "Species × age × longitude",
      TRUE ~ Test
    ),
    `χ²` = fmt_num(Chi_square, 1),
    df = as.integer(df),
    p = fmt_p(p)
  )

boundary_97w_lrt_table <- as.data.frame(
  longitude_boundary_lrt
) %>%
  filter(
    grepl("geography|latitude|longitude", Test, ignore.case = TRUE)
  ) %>%
  transmute(
    Analysis = "East of 97°W",
    Test = case_when(
      grepl("geography", Test, ignore.case = TRUE) ~
        "Species × age × geography",
      grepl("latitude", Test, ignore.case = TRUE) ~
        "Species × age × latitude",
      grepl("longitude", Test, ignore.case = TRUE) ~
        "Species × age × longitude",
      TRUE ~ Test
    ),
    `χ²` = fmt_num(Chi_square, 1),
    df = as.integer(df),
    p = fmt_p(p)
  )

stopifnot(
  nrow(date_covariate_lrt_table) == 3L,
  nrow(boundary_97w_lrt_table) == 3L
)

table_s10_interaction_tests <- bind_rows(
  date_covariate_lrt_table,
  boundary_97w_lrt_table
)

stopifnot(nrow(table_s10_interaction_tests) == 6L)

# Panel B: the two prespecified within-species contrasts emphasized in the
# manuscript under each robustness check.
select_focal_robustness_differences <- function(x, analysis_label) {
  as.data.frame(x) %>%
    filter(
      (Species == "Rufous" & Axis == "Latitude") |
        (Species == "Black-chinned" & Axis == "Longitude")
    ) %>%
    transmute(
      Analysis = analysis_label,
      Species = as.character(Species),
      Axis = as.character(Axis),
      Estimate = fmt_num(Estimate, 2),
      SE = fmt_num(SE, 2),
      `95% CI` = fmt_ci(CI_low, CI_high, 2),
      p = fmt_p(p)
    )
}

table_s10_focal_contrasts <- bind_rows(
  select_focal_robustness_differences(
    date_covariate_differences,
    "Day-of-window covariate"
  ),
  select_focal_robustness_differences(
    longitude_boundary_differences,
    "East of 97°W"
  )
)

stopifnot(nrow(table_s10_focal_contrasts) == 4L)

caption_s10_additional_robustness <- paste0(
  "Supplementary Material Table S10. Sampling-date and 97°W western-boundary robustness analyses. ",
  "Panel A shows joint, latitude, and longitude likelihood-ratio tests; Panel B shows the focal ",
  "Rufous latitude and Black-chinned longitude immature–adult contrasts."
)

ft_s10_interaction_tests <- nice_ft(table_s10_interaction_tests)
ft_s10_focal_contrasts <- nice_ft(table_s10_focal_contrasts)


# ============================================================
# 3. FORMAT SUPPLEMENTARY TABLES FOR WORD
# ============================================================

# Store each final supplementary table in a named entry matching its final
# Supplementary Material table number, caption, panel labels, and orientation.
supplement_entries <- list(
  S1 = list(
    caption = caption_s1_demographic_composition,
    tables = list(ft_s1_demographic_composition),
    panel_labels = NULL,
    orientation = "landscape"
  ),
  S2 = list(
    caption = caption_s2_supporting_species_uncertainty,
    tables = list(ft_s2_supporting_species_uncertainty),
    panel_labels = NULL,
    orientation = "landscape"
  ),
  S3 = list(
    caption = caption_s3_age_specific_slopes,
    tables = list(ft_s3_age_specific_slopes),
    panel_labels = NULL,
    orientation = "portrait"
  ),
  S4 = list(
    caption = caption_s4_geographic_sampling_support,
    tables = list(ft_s4_geographic_sampling_support),
    panel_labels = NULL,
    orientation = "portrait"
  ),
  S5 = list(
    caption = caption_s5_interspecific_contrasts,
    tables = list(ft_s5_interspecific_contrasts),
    panel_labels = NULL,
    orientation = "portrait"
  ),
  S6 = list(
    caption = caption_s6_primary_glmm,
    tables = list(ft_s6_fixed_effects, ft_s6_random_effects),
    panel_labels = c("Fixed effects", "Random effects"),
    orientation = "portrait"
  ),
  S7 = list(
    caption = caption_s7_gamm_sensitivity,
    caption_runs = caption_s7_gamm_sensitivity_runs,
    tables = list(ft_s7_gamm_summary),
    panel_labels = NULL,
    orientation = "landscape"
  ),
  S8 = list(
    caption = caption_s8_temporal_confounding,
    tables = list(ft_s8_temporal_lrt, ft_s8_temporal_effects, ft_s8_temporal_variance),
    panel_labels = c(
      "Panel A. Omnibus tests",
      "Panel B. Slope differences and longitude species contrasts",
      "Panel C. Geographic variance retained after within-year centering"
    ),
    orientation = "landscape"
  ),
  S9 = list(
    caption = caption_s9_leave_one_location_out,
    tables = list(ft_s9_leave_one_location_out),
    panel_labels = NULL,
    orientation = "landscape"
  ),
  S10 = list(
    caption = caption_s10_additional_robustness,
    tables = list(ft_s10_interaction_tests, ft_s10_focal_contrasts),
    panel_labels = c(
      "Panel A. Interaction tests",
      "Panel B. Focal within-species contrasts"
    ),
    orientation = "portrait"
  )
)

stopifnot(length(supplement_entries) == 10L)

# Use consistent margins across portrait and landscape sections.
aos_margins <- page_mar(
  top = 1,
  bottom = 1,
  left = 1,
  right = 1,
  header = 0.5,
  footer = 0.5
)

# Add page numbers in the footer.
page_number_footer <- block_list(
  fpar(
    run_word_field(
      field = "PAGE",
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 12
      )
    ),
    fp_p = fp_par(
      text.align = "center",
      line_spacing = 1
    )
  )
)

portrait_properties <- prop_section(
  page_size = page_size(
    width = 8.5,
    height = 11,
    orient = "portrait"
  ),
  page_margins = aos_margins,
  type = "nextPage",
  footer_default = page_number_footer
)

landscape_properties <- prop_section(
  page_size = page_size(
    width = 11,
    height = 8.5,
    orient = "landscape"
  ),
  page_margins = aos_margins,
  type = "nextPage",
  footer_default = page_number_footer
)

# Define the final portrait section without an additional page break.
portrait_final_properties <- prop_section(
  page_size = page_size(
    width = 8.5,
    height = 11,
    orient = "portrait"
  ),
  page_margins = aos_margins,
  type = "continuous",
  footer_default = page_number_footer
)

# Panel labels are separate paragraphs above their tables.
panel_text_fp <- fp_text(
  font.family = "Times New Roman",
  font.size = 12,
  bold = TRUE,
  italic = FALSE
)

panel_fp <- fp_par(
  text.align = "left",
  padding.top = 0,
  padding.bottom = 0,
  line_spacing = 1,
  keep_with_next = TRUE
)

# Italicize statistical symbols inside compound headers.
compose_n_pct_header <- function(ft, column, prefix) {
  compose(
    ft,
    j = column,
    part = "header",
    value = as_paragraph(
      prefix,
      as_i("n"),
      " (%)"
    )
  )
}

compose_p_header <- function(ft, column, prefix, suffix = "") {
  compose(
    ft,
    j = column,
    part = "header",
    value = as_paragraph(
      prefix,
      as_i("P"),
      suffix
    )
  )
}

ft_s1_demographic_composition <- compose_n_pct_header(ft_s1_demographic_composition, "Immature, n (%)", "Immature, ")
ft_s1_demographic_composition <- compose_n_pct_header(ft_s1_demographic_composition, "Female overall, n (%)", "Female overall, ")
ft_s1_demographic_composition <- compose_n_pct_header(ft_s1_demographic_composition, "Adult female, n (%)", "Adult female, ")
ft_s1_demographic_composition <- compose_n_pct_header(ft_s1_demographic_composition, "Immature female, n (%)", "Immature female, ")

ft_s2_supporting_species_uncertainty <- compose_n_pct_header(ft_s2_supporting_species_uncertainty, "Unknown sex, n (%)", "Unknown sex, ")
ft_s2_supporting_species_uncertainty <- compose_n_pct_header(
  ft_s2_supporting_species_uncertainty,
  "Unclassifiable age, n (%)",
  "Unclassifiable age, "
)

# Use merged spanning headers for the leave-one-location-out table.
ft_s9_leave_one_location_out <- set_header_labels(
  ft_s9_leave_one_location_out,
  `Joint P` = "Joint P",
  `Latitude P` = "Lat. P",
  `Longitude P` = "Lon. P",
  `Rufous lat. β (P)` = "Lat. β (P)",
  `Rufous lon. β (P)` = "Lon. β (P)",
  `Black-chinned lat. β (P)` = "Lat. β (P)",
  `Black-chinned lon. β (P)` = "Lon. β (P)",
  `Ruby-throated lat. β (P)` = "Lat. β (P)",
  `Ruby-throated lon. β (P)` = "Lon. β (P)"
)

ft_s9_leave_one_location_out <- compose(
  ft_s9_leave_one_location_out,
  j = "n omitted",
  part = "header",
  value = as_paragraph(as_i("n"), " omitted")
)

ft_s9_leave_one_location_out <- compose_p_header(ft_s9_leave_one_location_out, "Joint P", "Joint ")
ft_s9_leave_one_location_out <- compose_p_header(ft_s9_leave_one_location_out, "Latitude P", "Lat. ")
ft_s9_leave_one_location_out <- compose_p_header(ft_s9_leave_one_location_out, "Longitude P", "Lon. ")
ft_s9_leave_one_location_out <- compose_p_header(ft_s9_leave_one_location_out, "Rufous lat. β (P)", "Lat. β (", ")")
ft_s9_leave_one_location_out <- compose_p_header(ft_s9_leave_one_location_out, "Rufous lon. β (P)", "Lon. β (", ")")
ft_s9_leave_one_location_out <- compose_p_header(
  ft_s9_leave_one_location_out,
  "Black-chinned lat. β (P)",
  "Lat. β (",
  ")"
)
ft_s9_leave_one_location_out <- compose_p_header(
  ft_s9_leave_one_location_out,
  "Black-chinned lon. β (P)",
  "Lon. β (",
  ")"
)
ft_s9_leave_one_location_out <- compose_p_header(
  ft_s9_leave_one_location_out,
  "Ruby-throated lat. β (P)",
  "Lat. β (",
  ")"
)
ft_s9_leave_one_location_out <- compose_p_header(
  ft_s9_leave_one_location_out,
  "Ruby-throated lon. β (P)",
  "Lon. β (",
  ")"
)

ft_s9_leave_one_location_out <- add_header_row(
  ft_s9_leave_one_location_out,
  values = c(
    "", "", "",
    "Omnibus LRT",
    "Rufous",
    "Black-chinned",
    "Ruby-throated"
  ),
  colwidths = c(1, 1, 1, 3, 2, 2, 2),
  top = TRUE
) %>%
  bold(part = "header") %>%
  align(align = "left", part = "header")

# Format k and P as statistical symbols in the GAMM diagnostic headings.
ft_s7_gamm_summary <- set_header_labels(
  ft_s7_gamm_summary,
  `Reference df` = "Ref. df"
)

ft_s7_gamm_summary <- compose(
  ft_s7_gamm_summary,
  j = "k-index",
  part = "header",
  value = as_paragraph(as_i("k"), "-index")
)

ft_s7_gamm_summary <- compose(
  ft_s7_gamm_summary,
  j = "k-check p",
  part = "header",
  value = as_paragraph(as_i("k"), "-check ", as_i("P"))
)

# Update stored table objects after composing statistical headers.
supplement_entries$S1$tables[[1L]] <- ft_s1_demographic_composition
supplement_entries$S2$tables[[1L]] <- ft_s2_supporting_species_uncertainty
supplement_entries$S7$tables[[1L]] <- ft_s7_gamm_summary
supplement_entries$S9$tables[[1L]] <- ft_s9_leave_one_location_out

# Reapply common typography after composing headers.
all_supplement_flextables <- unlist(
  lapply(supplement_entries, function(entry) entry$tables),
  recursive = FALSE
)

all_supplement_flextables <- lapply(
  all_supplement_flextables,
  function(ft) {
    ft %>%
      font(fontname = "Times New Roman", part = "all") %>%
      fontsize(size = 12, part = "all") %>%
      bold(part = "header") %>%
      align(align = "left", part = "all") %>%
      valign(valign = "top", part = "all") %>%
      line_spacing(space = 1, part = "all") %>%
      padding(padding = 1, part = "all") %>%
      border_remove() %>%
      set_table_properties(
        layout = "autofit",
        width = 1,
        align = "left",
        opts_word = list(
          split = FALSE,
          keep_with_next = FALSE,
          repeat_headers = TRUE
        )
      )
  }
)

# Update stored tables after formatting.
flextable_index <- 1L
for (entry_index in seq_along(supplement_entries)) {
  entry_table_count <- length(supplement_entries[[entry_index]]$tables)
  supplement_entries[[entry_index]]$tables <- all_supplement_flextables[
    flextable_index:(flextable_index + entry_table_count - 1L)
  ]
  flextable_index <- flextable_index + entry_table_count
}

# Set stable column widths for portrait and landscape pages.
set_fixed_widths <- function(ft, widths) {
  stopifnot(length(widths) == ncol(ft$body$dataset))
  for (column_index in seq_along(widths)) {
    ft <- width(ft, j = column_index, width = widths[column_index])
  }
  set_table_properties(
    ft,
    layout = "fixed",
    width = 1,
    align = "left",
    opts_word = list(
      split = FALSE,
      keep_with_next = FALSE,
      repeat_headers = TRUE
    )
  )
}

# S1. Demographic composition.
supplement_entries$S1$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S1$tables[[1L]],
  c(1.15, 0.55, 1.10, 1.25, 1.20, 1.25)
)

# S2. Supporting-species demographic uncertainty.
supplement_entries$S2$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S2$tables[[1L]],
  c(1.10, 1.20, 1.15, 1.30, 1.45, 1.40, 1.40)
)

# S3. Primary age-specific geographic slopes.
supplement_entries$S3$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S3$tables[[1L]],
  c(0.70, 1.10, 0.75, 0.75, 0.55, 1.30, 0.55, 0.55)
)

# S4. Geographic sampling support.
supplement_entries$S4$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S4$tables[[1L]],
  c(1.20, 1.00, 0.95, 1.20, 1.00, 1.15)
)

# S5. Primary pairwise interspecific contrasts.
supplement_entries$S5$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S5$tables[[1L]],
  c(0.70, 1.90, 0.80, 0.60, 1.45, 0.55, 0.50)
)

# S6. Primary comparative GLMM fixed and random effects.
supplement_entries$S6$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S6$tables[[1L]],
  c(2.55, 0.75, 0.60, 1.40, 0.60, 0.60)
)
supplement_entries$S6$tables[[2L]] <- set_fixed_widths(
  supplement_entries$S6$tables[[2L]],
  c(3.00, 0.70, 1.35, 1.35)
)
supplement_entries$S6$tables[[2L]] <- set_table_properties(
  supplement_entries$S6$tables[[2L]],
  layout = "fixed",
  width = 1,
  align = "left",
  opts_word = list(
    split = FALSE,
    keep_with_next = TRUE,
    repeat_headers = TRUE
  )
)

# S7. GAMM nonlinear-spatial sensitivity.
supplement_entries$S7$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S7$tables[[1L]],
  c(1.05, 0.55, 1.75, 0.55, 0.70, 0.65, 0.45, 0.65, 0.80, 1.00)
)

# S8. Temporal-confounding robustness.
supplement_entries$S8$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S8$tables[[1L]],
  c(4.30, 1.20, 0.80, 0.80)
)
supplement_entries$S8$tables[[2L]] <- set_fixed_widths(
  supplement_entries$S8$tables[[2L]],
  c(0.90, 3.05, 0.80, 0.60, 1.60, 0.65, 0.65)
)
supplement_entries$S8$tables[[3L]] <- set_fixed_widths(
  supplement_entries$S8$tables[[3L]],
  c(1.50, 3.00, 3.00)
)
supplement_entries$S8$tables[[3L]] <- set_table_properties(
  supplement_entries$S8$tables[[3L]],
  layout = "fixed",
  width = 1,
  align = "left",
  opts_word = list(
    split = FALSE,
    keep_with_next = TRUE,
    repeat_headers = TRUE
  )
)

# S9. Leave-one-high-volume-location-out sensitivity.
supplement_entries$S9$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S9$tables[[1L]],
  c(1.15, 0.50, 0.55, 0.50, 0.50, 0.50, 0.80, 0.80, 0.90, 0.90, 0.90, 0.90)
)

# S10. Sampling-date and western-boundary robustness.
supplement_entries$S10$tables[[1L]] <- set_fixed_widths(
  supplement_entries$S10$tables[[1L]],
  c(1.70, 2.50, 0.90, 0.55, 0.55)
)
supplement_entries$S10$tables[[2L]] <- set_fixed_widths(
  supplement_entries$S10$tables[[2L]],
  c(1.55, 1.00, 0.70, 0.70, 0.55, 1.40, 0.50)
)

# ------------------------------------------------------------
# Emphasize row-label columns; numerical/statistical cells remain roman.
# ------------------------------------------------------------
bold_body_column <- function(ft, column) {
  if (!column %in% names(ft$body$dataset)) {
    stop("Cannot bold missing body column: ", column)
  }
  bold(ft, j = column, part = "body")
}

# S1. Demographic composition.
supplement_entries$S1$tables[[1L]] <- bold_body_column(
  supplement_entries$S1$tables[[1L]], "Species"
)

# S2. Supporting-species demographic uncertainty.
supplement_entries$S2$tables[[1L]] <- bold_body_column(
  supplement_entries$S2$tables[[1L]], "Species"
)

# S3. Primary age-specific geographic slopes.
supplement_entries$S3$tables[[1L]] <- bold_body_column(
  supplement_entries$S3$tables[[1L]], "Species"
)

# S4. Geographic sampling support.
supplement_entries$S4$tables[[1L]] <- bold_body_column(
  supplement_entries$S4$tables[[1L]], "Species"
)

# S5. Primary pairwise interspecific contrasts.
supplement_entries$S5$tables[[1L]] <- bold_body_column(
  supplement_entries$S5$tables[[1L]], "Contrast"
)

# S6. Primary comparative GLMM fixed and random effects.
supplement_entries$S6$tables[[1L]] <- bold_body_column(
  supplement_entries$S6$tables[[1L]], "Term"
)
supplement_entries$S6$tables[[2L]] <- bold_body_column(
  supplement_entries$S6$tables[[2L]], "Random effect"
)

# S7. GAMM nonlinear-spatial sensitivity.
supplement_entries$S7$tables[[1L]] <- bold_body_column(
  supplement_entries$S7$tables[[1L]], "Species"
)

# S8. Temporal-confounding robustness.
supplement_entries$S8$tables[[1L]] <- bold_body_column(
  supplement_entries$S8$tables[[1L]], "Test"
)
supplement_entries$S8$tables[[2L]] <- bold_body_column(
  supplement_entries$S8$tables[[2L]], "Comparison"
)
supplement_entries$S8$tables[[3L]] <- bold_body_column(
  supplement_entries$S8$tables[[3L]], "Species"
)

# S9. Leave-one-high-volume-location-out sensitivity.
supplement_entries$S9$tables[[1L]] <- bold_body_column(
  supplement_entries$S9$tables[[1L]], "Omitted location"
)

# S10. Sampling-date and western-boundary robustness.
supplement_entries$S10$tables[[1L]] <- bold_body_column(
  supplement_entries$S10$tables[[1L]], "Test"
)
supplement_entries$S10$tables[[2L]] <- bold_body_column(
  supplement_entries$S10$tables[[2L]], "Species"
)

# ============================================================
# 4. FINAL TABLE ORDER AND PAGE GROUPING
#
# Tables are stored and written in final numerical order. S1 and S2 share the
# first landscape page; each remaining table begins on a new page.
# ============================================================

supplement_page_groups <- list(
  list(
    orientation = "landscape",
    blocks = list(
      supplement_entries$S1,
      supplement_entries$S2
    )
  ),
  list(
    orientation = "portrait",
    blocks = list(supplement_entries$S3)
  ),
  list(
    orientation = "portrait",
    blocks = list(supplement_entries$S4)
  ),
  list(
    orientation = "portrait",
    blocks = list(supplement_entries$S5)
  ),
  list(
    orientation = "portrait",
    blocks = list(supplement_entries$S6)
  ),
  list(
    orientation = "landscape",
    blocks = list(supplement_entries$S7)
  ),
  list(
    orientation = "landscape",
    blocks = list(supplement_entries$S8)
  ),
  list(
    orientation = "landscape",
    blocks = list(supplement_entries$S9)
  ),
  list(
    orientation = "portrait",
    blocks = list(supplement_entries$S10)
  )
)

stopifnot(
  length(supplement_page_groups) == 9L,
  sum(vapply(
    supplement_page_groups,
    function(x) length(x$blocks),
    integer(1)
  )) == 10L
)

# ============================================================
# 5. WRITE FINAL SUPPLEMENT TO WORD
# ============================================================

supplement_doc <- read_docx()

for (page_group_index in seq_along(supplement_page_groups)) {
  page_group <- supplement_page_groups[[page_group_index]]
  
  for (block_index in seq_along(page_group$blocks)) {
    table_block <- page_group$blocks[[block_index]]
    
    caption_runs <- if (!is.null(table_block$caption_runs)) {
      table_block$caption_runs
    } else {
      list(ftext(table_block$caption, prop = caption_text_fp))
    }
    
    supplement_doc <- body_add_fpar(
      supplement_doc,
      do.call(
        fpar,
        c(caption_runs, list(fp_p = caption_fp))
      )
    )
    
    for (table_index in seq_along(table_block$tables)) {
      if (!is.null(table_block$panel_labels)) {
        supplement_doc <- body_add_fpar(
          supplement_doc,
          fpar(
            ftext(
              table_block$panel_labels[table_index],
              prop = panel_text_fp
            ),
            fp_p = panel_fp
          )
        )
      }
      
      supplement_doc <- body_add_flextable(
        supplement_doc,
        table_block$tables[[table_index]]
      )
    }
  }
  
  # Start each later page group on a new page, using a section break when
  # the page orientation changes.
  if (page_group_index < length(supplement_page_groups)) {
    next_page_orientation <-
      supplement_page_groups[[page_group_index + 1L]]$orientation
    
    if (identical(page_group$orientation, next_page_orientation)) {
      supplement_doc <- body_add_break(supplement_doc, pos = "after")
    } else {
      current_section_properties <-
        if (identical(page_group$orientation, "landscape")) {
          landscape_properties
        } else {
          portrait_properties
        }
      
      supplement_doc <- body_end_block_section(
        supplement_doc,
        block_section(current_section_properties)
      )
    }
  }
}

# End with the final portrait section.
supplement_doc <- body_set_default_section(
  supplement_doc,
  portrait_final_properties
)

print(
  supplement_doc,
  target = "Supplementary_Tables.docx"
)

