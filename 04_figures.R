# ============================================================
# 04_figures.R
# Generate Figures 1-4 from the fitted primary analysis objects.
# Run 02_primary_analysis.R first.
#
# Figure 1: age x sex composition for all six study species
# Figure 2: geographic sampling coverage for the three focal species
# Figure 3: model-predicted P(Female) by species and age class
# Figure 4: age-specific latitude and longitude slopes
# ============================================================


library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)


# ============================================================
# 0. CHECK REQUIRED INPUTS AND DEFINE PLOT SETTINGS
# ============================================================

required_objects <- c(
  "results",
  "comparative_primary",
  "comparative_df",
  "m_comparative_full",
  "focal_species",
  "species_order",
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
    "Missing required object(s): ",
    paste(missing_objects, collapse = ", "),
    ". Run 01_functions.R and 02_primary_analysis.R first."
  )
}


focal_display_order <- c(
  "Rufous",
  "Black-chinned",
  "Ruby-throated"
)

retained_display_order <- c(
  focal_display_order,
  "Calliope",
  "Broad-tailed",
  "Allen's"
)

stopifnot(
  identical(focal_species, focal_display_order),
  identical(species_order, retained_display_order)
)


# Load map boundaries shared by Figures 2 and 3.
us_states <- map_data("state")

world_map <- map_data("world") %>%
  filter(
    region %in% c("USA", "Canada")
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
  
  c(lower = lower, upper = upper)
}


# ============================================================
# FIGURE 1
# Age x sex composition across all six study species
# ============================================================

class_levels <- c(
  "Adult Female",
  "Adult Male",
  "Immature Female",
  "Immature Male"
)

composition_df <- purrr::imap_dfr(
  results[retained_display_order],
  function(res, sp) {
    res$data_points %>%
      transmute(
        Species = sp,
        Class = case_when(
          Immature == 0 & Female == 1 ~ "Adult Female",
          Immature == 0 & Female == 0 ~ "Adult Male",
          Immature == 1 & Female == 1 ~ "Immature Female",
          Immature == 1 & Female == 0 ~ "Immature Male"
        )
      ) %>%
      count(
        Species,
        Class,
        name = "N"
      )
  }
) %>%
  complete(
    Species = retained_display_order,
    Class = class_levels,
    fill = list(N = 0)
  ) %>%
  group_by(Species) %>%
  mutate(
    Total = sum(N),
    Prop = N / Total
  ) %>%
  ungroup() %>%
  mutate(
    Species = factor(
      Species,
      levels = rev(retained_display_order)
    ),
    Class = factor(
      Class,
      levels = class_levels
    )
  )


species_labels <- composition_df %>%
  distinct(
    Species,
    Total
  ) %>%
  mutate(
    label = paste0(
      as.character(Species),
      "  (n = ",
      scales::comma(Total),
      ")"
    )
  )

label_vector <- setNames(
  species_labels$label,
  species_labels$Species
)


class_colors <- c(
  "Immature Male" = "#56B4E9",
  "Immature Female" = "#E69F00",
  "Adult Male" = "#0072B2",
  "Adult Female" = "#D55E00"
)


fig1 <- ggplot(
  composition_df,
  aes(
    x = Prop,
    y = Species,
    fill = Class
  )
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = ifelse(
        Prop >= 0.07,
        scales::percent(
          Prop,
          accuracy = 1
        ),
        ""
      )
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 3.1
  ) +
  scale_fill_manual(
    values = class_colors,
    breaks = c(
      "Immature Male",
      "Immature Female",
      "Adult Male",
      "Adult Female"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    labels = scales::percent_format(
      accuracy = 10
    ),
    breaks = seq(
      0,
      1,
      0.2
    ),
    limits = c(
      0,
      1
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_discrete(
    labels = label_vector
  ) +
  labs(
    x = "Proportion of classifiable records",
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(
      size = 11,
      face = "plain",
      margin = margin(
        b = 10
      )
    ),
    legend.text = element_text(
      size = 9
    ),
    legend.box.margin = margin(
      b = 4,
      unit = "mm"
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.title.x = element_text(
      size = 11
    ),
    panel.grid.major.x = element_line(
      linewidth = 0.2,
      color = "grey90"
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20,
      unit = "pt"
    )
  )


# ============================================================
# FIGURE 2
# Geographic coverage of primary-analysis records
# ============================================================

# Build the common 100-km grid used for sampling counts.
land <- rnaturalearth::ne_countries(
  country = c(
    "United States of America",
    "Canada"
  ),
  returnclass = "sf"
)

land_5070 <- st_transform(
  land,
  5070
)

land_union_5070 <- st_union(
  land_5070
)

grid_geom <- st_make_grid(
  land_union_5070,
  cellsize = c(
    100000,
    100000
  ),
  what = "polygons",
  square = TRUE
)

grid_5070 <- st_sf(
  cell_id = seq_along(
    grid_geom
  ),
  geometry = grid_geom
)

grid_5070 <- grid_5070[
  lengths(
    st_intersects(
      grid_5070,
      land_union_5070
    )
  ) > 0,
  ,
  drop = FALSE
]


make_sampling_cells <- function(
    res,
    species_name
) {
  pts <- st_as_sf(
    res$data_points,
    coords = c(
      "lon_dd",
      "lat_dd"
    ),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(
      5070
    )
  
  counts <- st_join(
    pts,
    grid_5070,
    join = st_intersects,
    left = FALSE
  ) %>%
    st_drop_geometry() %>%
    count(
      cell_id,
      name = "N"
    )
  
  # Every primary-analysis record must be assigned exactly once.
  if (sum(counts$N) != nrow(res$data_points)) {
    stop(
      "Figure 2 grid assignment did not preserve the exact number of ",
      "primary-analysis records for ", species_name, "."
    )
  }
  
  grid_5070 %>%
    inner_join(
      counts,
      by = "cell_id"
    ) %>%
    mutate(
      Species = species_name
    ) %>%
    st_transform(
      4326
    )
}


sampling_sf <- do.call(
  rbind,
  lapply(
    focal_display_order,
    function(sp) {
      make_sampling_cells(
        results[[sp]],
        sp
      )
    }
  )
)

sampling_sf$Species <- factor(
  sampling_sf$Species,
  levels = focal_display_order
)


# Set the map extent to include every occupied 100-km cell.
fig2_bbox <- st_bbox(
  sampling_sf
)

fig2_extent_pad <- 0.5

fig2_xlim <- c(
  min(-102, unname(fig2_bbox["xmin"]) - fig2_extent_pad),
  max(-72.5, unname(fig2_bbox["xmax"]) + fig2_extent_pad)
)

fig2_ylim <- c(
  min(24, unname(fig2_bbox["ymin"]) - fig2_extent_pad),
  max(43.5, unname(fig2_bbox["ymax"]) + fig2_extent_pad)
)


fig2 <- ggplot() +
  geom_path(
    data = world_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    inherit.aes = FALSE,
    colour = "grey50",
    linewidth = 0.35
  ) +
  geom_path(
    data = us_states,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    inherit.aes = FALSE,
    colour = "grey75",
    linewidth = 0.25
  ) +
  geom_sf(
    data = sampling_sf,
    aes(
      fill = N
    ),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.1
  ) +
  facet_wrap(
    ~ Species,
    ncol = 2
  ) +
  coord_sf(
    xlim = fig2_xlim,
    ylim = fig2_ylim,
    expand = FALSE,
    default_crs = st_crs(
      4326
    )
  ) +
  scale_fill_viridis_c(
    trans = "log10",
    breaks = c(
      1,
      3,
      10,
      30,
      100,
      300,
      1000
    ),
    labels = scales::comma,
    name = "Records per 100 x 100 km cell"
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 12,
      margin = margin(
        t = 0,
        b = 11
      )
    ),
    strip.background = element_blank(),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    legend.position = c(
      0.63,
      0.13
    ),
    legend.justification = c(
      0,
      0
    ),
    legend.title = element_text(
      size = 11,
      face = "plain",
      margin = margin(
        b = 10
      )
    ),
    legend.text = element_text(
      size = 9
    ),
    legend.background = element_rect(
      fill = scales::alpha(
        "white",
        0.9
      ),
      colour = NA
    ),
    panel.spacing.x = grid::unit(
      3,
      "mm"
    ),
    panel.spacing.y = grid::unit(
      3,
      "mm"
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20,
      unit = "pt"
    )
  )


# ============================================================
# FIGURE 3
# Model-predicted P(Female)
#
# Predictions are shown only in 100 × 100 km cells containing at least one
# record of the corresponding species × age class in the primary comparative
# analysis. Cell color represents only predicted P(Female).
#
# Sampling-support symbols:
#   no symbol    = 1–4 records
#   open circle  = 5–9 records
#   filled circle = ≥10 records
# ============================================================

fig3_cell_km <- 100


# ------------------------------------------------------------
# Build a reproducible common grid from all qualifying focal-species records.
# ------------------------------------------------------------

fig3_raw <- fread(
  csv_file,
  select = c(
    "band",
    "species_id",
    "event_year",
    "event_month",
    "event_day",
    "iso_country",
    "lat_dd",
    "lon_dd"
  )
)

fig3_lookup <- data.table(
  Species = names(species_ids),
  species_id = as.integer(
    species_ids
  )
)

fig3_raw <- merge(
  fig3_raw,
  fig3_lookup,
  by = "species_id"
)

fig3_raw <- fig3_raw[
  Species %in%
    focal_display_order
]

fig3_raw[, event_year :=
           as.integer(
             event_year
           )]

fig3_raw[, event_month :=
           as.integer(
             event_month
           )]

fig3_raw[, event_day :=
           as.integer(
             event_day
           )]


fig3_date_keep <- in_month_day_window(
  month = fig3_raw$event_month,
  day = fig3_raw$event_day,
  start_month = primary_start_month,
  start_day = primary_start_day,
  end_month = primary_end_month,
  end_day = primary_end_day
)

fig3_raw <- fig3_raw[
  iso_country %in% countries_keep &
    event_month %in% months_keep &
    fig3_date_keep &
    is.finite(lat_dd) &
    is.finite(lon_dd) &
    lon_dd >= lon_min &
    lon_dd <= lon_max &
    lat_dd >= lat_min &
    lat_dd <= lat_max
]

# Missing day is used only for ordering after date-window filtering.
fig3_raw[
  is.na(event_day),
  event_day := 15L
]

fig3_raw[, dt :=
           as.IDate(
             sprintf(
               "%04d-%02d-%02d",
               event_year,
               event_month,
               event_day
             )
           )]

fig3_raw <- fig3_raw[
  !is.na(dt) &
    !is.na(band)
]

fig3_winter_fun <- infer_season_year_fun(
  months_keep
)

fig3_raw[, winter_year :=
           fig3_winter_fun(
             event_year,
             event_month
           )]

setorder(
  fig3_raw,
  Species,
  band,
  winter_year,
  dt
)

fig3_records <- fig3_raw[
  ,
  .SD[1L],
  by = .(
    Species,
    band,
    winter_year
  )
]


fig3_pts_5070 <- st_as_sf(
  fig3_records,
  coords = c(
    "lon_dd",
    "lat_dd"
  ),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(
    5070
  )


cell_m <- fig3_cell_km * 1000
fig3_bbox <- st_bbox(
  fig3_pts_5070
)

x0 <- floor(
  fig3_bbox["xmin"] /
    cell_m
) * cell_m

y0 <- floor(
  fig3_bbox["ymin"] /
    cell_m
) * cell_m

fig3_grid_geom <- st_make_grid(
  fig3_pts_5070,
  cellsize = cell_m,
  what = "polygons",
  square = TRUE,
  offset = c(
    x0,
    y0
  )
)

fig3_grid <- st_sf(
  cell_id = seq_along(
    fig3_grid_geom
  ),
  geometry = fig3_grid_geom
)


# ------------------------------------------------------------
# Count age-specific sampling support from the primary model data.
# ------------------------------------------------------------

fig3_model_pts <- st_as_sf(
  comparative_df,
  coords = c(
    "lon_dd",
    "lat_dd"
  ),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(
    5070
  )

fig3_age_counts <- st_join(
  fig3_model_pts,
  fig3_grid,
  join = st_within,
  left = FALSE
) %>%
  st_drop_geometry() %>%
  mutate(
    Species = as.character(
      Species
    ),
    Age = as.character(
      Immature_f
    )
  ) %>%
  count(
    Species,
    Age,
    cell_id,
    name = "N_age"
  ) %>%
  mutate(
    Support = case_when(
      N_age <= 4 ~ "1–4",
      N_age <= 9 ~ "5–9",
      TRUE ~ "≥10"
    ),
    Species = factor(
      Species,
      levels = focal_display_order
    ),
    Age = factor(
      Age,
      levels = c(
        "Adult",
        "Immature"
      )
    ),
    Support = factor(
      Support,
      levels = c(
        "1–4",
        "5–9",
        "≥10"
      )
    )
  )

fig3_support_sf <- fig3_grid %>%
  inner_join(
    fig3_age_counts,
    by = "cell_id"
  )


# ------------------------------------------------------------
# Calculate population-level fixed-effect predictions at cell centers.
# ------------------------------------------------------------

predict_fig3_species <- function(
    species_name
) {
  
  predict_age <- function(
    age_name
  ) {
    cells <- fig3_support_sf %>%
      filter(
        Species == species_name,
        Age == age_name
      )
    
    if (nrow(cells) == 0L) {
      return(NULL)
    }
    
    centers_ll <- st_centroid(
      cells
    ) %>%
      st_transform(
        4326
      )
    
    xy <- st_coordinates(
      centers_ll
    )
    
    cells$lon_dd <- xy[, 1]
    cells$lat_dd <- xy[, 2]
    cells$lat5 <- (
      cells$lat_dd - 31
    ) / 5
    cells$lon5 <- (
      cells$lon_dd + 90
    ) / 5
    
    newdata <- data.frame(
      Species = factor(
        rep(
          species_name,
          nrow(cells)
        ),
        levels = levels(
          comparative_df$Species
        )
      ),
      Immature_f = factor(
        rep(
          age_name,
          nrow(cells)
        ),
        levels = levels(
          comparative_df$Immature_f
        )
      ),
      lat5 = cells$lat5,
      lon5 = cells$lon5
    )
    
    fixed_terms <- delete.response(
      terms(
        lme4::nobars(
          formula(
            m_comparative_full
          )
        )
      )
    )
    
    X <- model.matrix(
      fixed_terms,
      data = newdata
    )
    
    beta <- lme4::fixef(
      m_comparative_full
    )
    
    X <- X[
      ,
      names(beta),
      drop = FALSE
    ]
    
    cells$pred <- plogis(
      as.numeric(
        X %*% beta
      )
    )
    
    cells
  }
  
  rbind(
    predict_age(
      "Adult"
    ),
    predict_age(
      "Immature"
    )
  )
}


predicted_sf <- do.call(
  rbind,
  lapply(
    focal_display_order,
    predict_fig3_species
  )
)

predicted_sf$Species <- factor(
  predicted_sf$Species,
  levels = focal_display_order
)

predicted_sf$Age <- factor(
  predicted_sf$Age,
  levels = c(
    "Adult",
    "Immature"
  )
)

predicted_sf$Support <- factor(
  predicted_sf$Support,
  levels = c(
    "1–4",
    "5–9",
    "≥10"
  )
)


# ------------------------------------------------------------
# Verify that grid-cell counts reproduce the primary model sample sizes.
# ------------------------------------------------------------

stopifnot(
  all(
    predicted_sf$N_age >= 1
  )
)

fig3_support_record_check <- fig3_age_counts %>%
  group_by(
    Species,
    Age
  ) %>%
  summarise(
    N_from_cells = sum(
      N_age
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Species = as.character(
      Species
    ),
    Age = as.character(
      Age
    )
  ) %>%
  arrange(
    Species,
    Age
  )

fig3_model_record_check <- comparative_df %>%
  transmute(
    Species = as.character(
      Species
    ),
    Age = as.character(
      Immature_f
    )
  ) %>%
  count(
    Species,
    Age,
    name = "N_model"
  ) %>%
  arrange(
    Species,
    Age
  )

stopifnot(
  isTRUE(
    all.equal(
      fig3_support_record_check %>%
        select(
          Species,
          Age,
          N_model = N_from_cells
        ),
      fig3_model_record_check,
      check.attributes = FALSE
    )
  )
)


# ------------------------------------------------------------
# Sampling-support symbols and plotting extent.
# ------------------------------------------------------------

fig3_support_points <- predicted_sf %>%
  filter(
    Support %in% c(
      "5–9",
      "≥10"
    )
  ) %>%
  st_centroid()

fig3_support_points$Support_symbol <- factor(
  as.character(
    fig3_support_points$Support
  ),
  levels = c(
    "5–9",
    "≥10"
  )
)


# Set the plotting extent to include every occupied species x age cell.
fig3_extent_all <- predicted_sf %>%
  st_transform(
    4326
  )

fig3_bbox_all <- st_bbox(
  fig3_extent_all
)

fig3_extent_pad <- 1.0

fig3_xlim <- c(
  min(
    -100,
    unname(
      fig3_bbox_all["xmin"]
    ) -
      fig3_extent_pad
  ),
  max(
    -74,
    unname(
      fig3_bbox_all["xmax"]
    ) +
      fig3_extent_pad
  )
)

fig3_ylim <- c(
  min(
    25,
    unname(
      fig3_bbox_all["ymin"]
    ) -
      fig3_extent_pad
  ),
  max(
    42,
    unname(
      fig3_bbox_all["ymax"]
    ) +
      fig3_extent_pad
  )
)


fig3 <- ggplot() +
  geom_path(
    data = world_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    inherit.aes = FALSE,
    colour = "grey50",
    linewidth = 0.30
  ) +
  geom_path(
    data = us_states,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    inherit.aes = FALSE,
    colour = "grey75",
    linewidth = 0.22
  ) +
  geom_sf(
    data = predicted_sf,
    aes(
      fill = pred
    ),
    inherit.aes = FALSE,
    color = "darkgray",
    linewidth = 0.20
  ) +
  geom_sf(
    data = fig3_support_points,
    aes(
      shape = Support_symbol
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 0.75,
    stroke = 0.28
  ) +
  facet_grid(
    Species ~ Age
  ) +
  coord_sf(
    xlim = fig3_xlim,
    ylim = fig3_ylim,
    expand = FALSE,
    default_crs = st_crs(
      4326
    )
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "gray98",
    high = "#B2182B",
    midpoint = 0.50,
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    labels = scales::percent,
    name = "Predicted P(Female)"
  ) +
  scale_shape_manual(
    values = c(
      "5–9" = 1,
      "≥10" = 16
    ),
    breaks = c(
      "5–9",
      "≥10"
    ),
    labels = expression(
      5-9,
      "" >= 10
    ),
    name = "Cell records"
  ) +
  guides(
    fill = guide_colorbar(
      order = 1
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        size = 2.2,
        stroke = 0.45
      )
    )
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    strip.text.x = element_text(
      face = "bold",
      size = 12,
      hjust = 0.33,
      margin = margin(
        b = 11
      )
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 12
    ),
    strip.background = element_blank(),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    legend.position = "right",
    legend.box.spacing = grid::unit(
      26,
      "pt"
    ),
    legend.box = "vertical",
    legend.spacing.y = grid::unit(
      25,
      "pt"
    ),
    legend.title = element_text(
      margin = margin(
        b = 10
      )
    ),
    panel.spacing.x = grid::unit(
      3,
      "mm"
    ),
    panel.spacing.y = grid::unit(
      3,
      "mm"
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20,
      unit = "pt"
    )
  )


# ============================================================
# FIGURE 4
# Age-specific geographic slopes from the primary comparative GLMM
#
# Slopes are changes in log-odds of a record being female per 5° increase
# in latitude or longitude.
# ============================================================

latitude_slopes <- as.data.frame(
  comparative_primary$latitude_slopes
)

longitude_slopes <- as.data.frame(
  comparative_primary$longitude_slopes
)

lat_ci <- get_ci_cols(
  latitude_slopes
)

lon_ci <- get_ci_cols(
  longitude_slopes
)


fig4_df <- bind_rows(
  latitude_slopes %>%
    transmute(
      Species = as.character(
        Species
      ),
      Age = as.character(
        Immature_f
      ),
      Axis = "Latitude",
      estimate = lat5.trend,
      conf.low = .data[[lat_ci["lower"]]],
      conf.high = .data[[lat_ci["upper"]]]
    ),
  longitude_slopes %>%
    transmute(
      Species = as.character(
        Species
      ),
      Age = as.character(
        Immature_f
      ),
      Axis = "Longitude",
      estimate = lon5.trend,
      conf.low = .data[[lon_ci["lower"]]],
      conf.high = .data[[lon_ci["upper"]]]
    )
) %>%
  mutate(
    Species = factor(
      Species,
      levels = rev(
        focal_display_order
      )
    ),
    Age = factor(
      Age,
      levels = c(
        "Adult",
        "Immature"
      )
    ),
    Axis = factor(
      Axis,
      levels = c(
        "Latitude",
        "Longitude"
      ),
      labels = c(
        "Latitude slope",
        "Longitude slope"
      )
    )
  )


fig4_xmax <- max(
  abs(
    c(
      fig4_df$conf.low,
      fig4_df$conf.high
    )
  ),
  na.rm = TRUE
)

fig4_xmax <- ceiling(
  fig4_xmax * 10
) / 10


fig4 <- ggplot(
  fig4_df,
  aes(
    x = estimate,
    y = Species
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.45,
    color = "grey50"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.65
  ) +
  geom_point(
    size = 2.7
  ) +
  facet_grid(
    Axis ~ Age
  ) +
  coord_cartesian(
    xlim = c(
      -fig4_xmax,
      fig4_xmax
    )
  ) +
  labs(
    x = "Change in log-odds of female per 5 degrees geographic change (95% CI)",
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    strip.text.x = element_text(
      face = "bold",
      size = 12,
      margin = margin(
        b = 11
      )
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 12
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.title.x = element_text(
      size = 11
    ),
    panel.spacing.x = grid::unit(
      3,
      "mm"
    ),
    panel.spacing.y = grid::unit(
      3,
      "mm"
    ),
    plot.margin = margin(
      t = 20,
      r = 20,
      b = 20,
      l = 20,
      unit = "pt"
    )
  )


# ============================================================
# SAVE FIGURES AS PNG AND VECTOR PDF
# ============================================================

ggsave(
  "Figure1.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure2.png",
  fig2,
  width = 9,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure3.png",
  fig3,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure4.png",
  fig4,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure1.pdf",
  fig1,
  width = 8,
  height = 5,
  bg = "white"
)

ggsave(
  "Figure2.pdf",
  fig2,
  width = 9,
  height = 8,
  bg = "white"
)

ggsave(
  "Figure3.pdf",
  fig3,
  width = 8,
  height = 6,
  bg = "white"
)

ggsave(
  "Figure4.pdf",
  fig4,
  width = 8,
  height = 6,
  bg = "white"
)

