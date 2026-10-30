# ==============================================================================
# 014_ExposureTrend.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#   To describe the temporal trend (2000-2019) of annual exposure to extreme
#   climatic events (ECEs), this script
#     (1) collapses each exposure's monthly grid-cell series into one annual
#         metric per calendar year,
#     (2) computes the global annual series over all low- and middle-income
#         country (LMIC) grid cells,
#     (3) computes the region-specific annual series for the five World Bank
#         regions,
#     (4) exports the wide annual-trend tables and the global / regional
#         trend figures.
#
# Methodology Overview:
#   1. Annual aggregation: for each ECE and each calendar year, the 12 monthly
#      grid-cell values are summed, then summarised across grid cells either by
#      averaging (continuous metrics: heatwave days, extreme-precipitation
#      volume) or as the proportion of exposed grid-cell-months (binary month
#      indicators: flood / drought / tropical-cyclone / wildfire months).
#   2. Grid: ~0.1-deg grid cells over 128 LMICs; monthly columns span 2000-2020
#      (wildfire 2002-2020), stored after four identifier columns
#      (OBJECTID, COUNTRY, lon, lat).
#   3. Regional split: grid cells are mapped to a region through
#      mortality_country_merge.csv (COUNTRY -> REGION).
#   4. Trend lines: cubic-polynomial GLM smooths are overlaid on the annual
#      point series for visualisation only.
#
#
# DEPENDENCIES
#   R        >= 4.0.3
#   tibble    == 3.1.8
#   ggplot2   >= 3.3.0
#   gridExtra  == 2.3
#
# INPUT
#   frame_hw_90_2_days_LIMCs_month_10.rda  -- heatwave (>= 90th pct, >= 2 days)
#   frame_EP_99_cum_LIMCs_month_10.rda     -- extreme precipitation (>= 99th pct)
#   frame_flood_total_LIMCs_month_10.rda   -- total flood
#   frame_dr_6m_LIMCs_month_10.rda         -- 6-month drought (SPEI <= -1)
#   frame_total_TC_LIMCs_month_10.rda      -- total tropical-cyclone exposure
#   frame_wf_exposure_LIMCs_month_10.rda   -- wildfire exposure
#   mortality_country_merge.csv            -- country -> region mapping
#
# OUTPUT
#   frame_EWE_year_global.csv             -- global annual series (wide, 6 x 21)
#   frame_EWE_year_regional_<code>.csv    -- regional annual series (wide, 6 x 21)
#   EWE_exposure_year_global.jpg          -- global figure (6 ECE panels)
#   EWE_<fig>_year_regional.jpg           -- regional figure (global + 5 regions)
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(tibble);library(ggplot2);library(gridExtra)

# --- paths (edit these to point at your own copies) ---------------------------
data_path <- file.path(
  "Repository/Data",
  "frame_u5mr_ECE_month.rda"
)

out_dir <- file.path(
  "Repository/Output", "Result"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

years <- 2000:2019


# ==============================================================================
# 1. Exposure configuration
# ==============================================================================
# One row per exposure definition. Column semantics:
#   family     : ECE family, used as the wide-table row label and figure title
#   file       : grid-cell frame (rda) in `burden_dir`; the object name equals
#                the file name without the ".rda" suffix
#   agg        : annual aggregation rule
#                 "mean" = mean annual days / mm across grid cells (continuous)
#                 "prop" = proportion of exposed grid-cell-months (%) (binary)
#   start_year : first year with monthly columns in the frame
#   fig        : short tag for the regional figure file name
#   yvar       : long-table column name used by the plotting code
#   ylab       : y-axis label
#   y_lo,y_hi,y_step : regional-figure y-axis limits and tick step
ece_config <- tribble(
  ~family,                 ~file,                                    ~agg,   ~start_year, ~fig,    ~yvar,                  ~ylab,                                    ~y_lo, ~y_hi, ~y_step,
  "Heatwave",              "frame_hw_90_2_days_LIMCs_month_10.rda", "mean", 2000,        "HW",    "Heatwave",             "Heatwave days (day)",                       5.0,  60.0,  10.0,
  "Extreme precipitation", "frame_EP_99_cum_LIMCs_month_10.rda",    "mean", 2000,        "EP",    "Extreme.Precipitation","Cumulative extreme precipitation (mm)",    30.0, 200.0,  20.0,
  "Flood",                 "frame_flood_total_LIMCs_month_10.rda",  "prop", 2000,        "Flood", "Flood",                "Proportion of flood month (%)",            -1.0,  10.0,   1.0,
  "Drought",               "frame_dr_6m_LIMCs_month_10.rda",        "prop", 2000,        "SPEI", "Drought",              "Proportion of drought month (%)",           5.0,  60.0,  10.0,
  "Tropical cyclone",      "frame_total_TC_LIMCs_month_10.rda",     "prop", 2000,        "TC",    "Tropical.Cyclone",     "Proportion of tropical cyclone month (%)", -0.1,   1.0,   0.2,
  "Wildfire",              "frame_wf_exposure_LIMCs_month_10.rda",  "prop", 2002,        "WF",    "Wildfire",             "Proportion of wildfire month (%)",          0.0,  25.0,   5.0
)

n_ece <- nrow(ece_config)

# Plot order of the five regions (global panel is always panel "A").
region_info <- tribble(
  ~region,                                ~code,
  "Sub-Saharan Africa",                    "SSA",
  "Central & South Asia",                  "CSA",
  "Europe, Middle East & North Africa",    "EMN",
  "Latin America & Caribbean",             "LAC",
  "East Asia & Pacific",                   "EAP"
)

n_region <- nrow(region_info)

# Country -> region lookup (named character vector; unmatched -> NA).
region_df <- read.csv(file.path(burden_dir, "mortality_country_merge.csv"))
region_map <- setNames(as.character(region_df$REGION),
                       as.character(region_df$COUNTRY))


# ==============================================================================
# 2. Core functions
# ==============================================================================

#' Sum the 12 calendar months of one year for every grid cell.
#'
#' Monthly columns start at position 5 (after OBJECTID, COUNTRY, lon, lat), and
#' column 5 corresponds to January of `start_year`. Column `j` of a year is
#' therefore at position `5 + (year - start_year) * 12 + (j - 1)`.
#'
#' @param dat        grid-cell frame (data.frame).
#' @param year       calendar year (integer).
#' @param start_year first year of the monthly columns (integer).
#' @return A numeric vector of length nrow(dat), one annual total per cell.
yearly_sum <- function(dat, year, start_year) {
  cols <- 5 + (year - start_year) * 12 + 0:11
  rowSums(dat[, cols, drop = FALSE], na.rm = TRUE)
}

#' Annual exposure series (2000-2019) for a set of grid cells.
#'
#' @param dat        grid-cell frame (data.frame, possibly a regional subset).
#' @param start_year first year with monthly columns in `dat`.
#' @param agg        "mean" (average days / mm per cell) or "prop" (% of
#'                   exposed cell-months).
#' @return A numeric vector of length 20 (NA for years before `start_year`).
annual_series <- function(dat, start_year, agg) {
  out <- rep(NA_real_, length(years))
  for (j in seq_along(years)) {
    y <- years[j]
    if (y < start_year) next
    x <- yearly_sum(dat, y, start_year)
    out[j] <- if (agg == "mean") round(mean(x, na.rm = TRUE), 1)
              else               round(sum(x, na.rm = TRUE) / (nrow(dat) * 12) * 100, 2)
  }
  out
}

#' Convert a wide (ECE x year) matrix into a long plotting table.
#'
#' @param mat numeric matrix with `n_ece` rows and one column per year.
#' @return A data.frame with columns `year` plus one column per ECE (`yvar`).
make_long <- function(mat) {
  d <- as.data.frame(t(mat))
  names(d) <- ece_config$yvar
  d$year <- years
  d[, c("year", ece_config$yvar)]
}

# Panel styling shared by all figures.
panel_theme <- theme_bw() +
  theme(
    plot.title   = element_text(face = "bold", color = "black", size = 22),
    axis.title.x = element_text(face = "bold", color = "black", size = 20),
    axis.title.y = element_text(face = "bold", color = "black", size = 20),
    axis.text.x  = element_text(face = "bold", color = "black", size = 15),
    axis.text.y  = element_text(face = "bold", color = "black", size = 15)
  )

#' Draw one annual-trend panel (points + line + cubic-polynomial smooth).
#'
#' @param dat      long table (must contain `year` and `yvar`).
#' @param yvar     column name to plot on the y axis.
#' @param title    panel title (e.g. "(B) Sub-Saharan Africa").
#' @param ylab     y-axis label.
#' @param x_breaks x-axis tick positions.
#' @param y_limits optional y-axis limits; if NULL the scale is left free.
#' @param y_breaks optional y-axis tick positions (used with `y_limits`).
trend_panel <- function(dat, yvar, title, ylab, x_breaks,
                        y_limits = NULL, y_breaks = NULL) {
  p <- ggplot(dat, aes(x = year, y = .data[[yvar]])) +
    geom_point() +
    geom_line(size = 1) +
    stat_smooth(method = "glm", formula = y ~ poly(x, 3),
                se = TRUE, size = 1, colour = "royalblue3",
                fill = "darkslategray3", lty = 2) +
    labs(title = title, x = "Year", y = ylab) +
    scale_x_continuous(breaks = x_breaks, labels = as.character(x_breaks)) +
    panel_theme
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, breaks = y_breaks)
  }
  p
}

#' Open a JPEG device with the figure settings used for every export.
jpeg_open <- function(path) {
  jpeg(path, res = 300, width = 4000, height = 5000)
}


# ==============================================================================
# 3. Compute the global and regional annual series
# ==============================================================================
global_mat   <- matrix(NA_real_, nrow = n_ece, ncol = length(years))
regional_mat <- replicate(n_region,
                          matrix(NA_real_, nrow = n_ece, ncol = length(years)),
                          simplify = FALSE)

for (i in seq_len(n_ece)) {
  cfg <- ece_config[i, ]
  obj <- sub("\\.rda$", "", cfg$file)

  # Load the frame into a local environment and pull out the data.frame.
  env <- new.env()
  load(file.path(burden_dir, cfg$file), envir = env)
  frame <- as.data.frame(env[[obj]])
  rm(env)

  # Global series: every LMIC grid cell.
  global_mat[i, ] <- annual_series(frame, cfg$start_year, cfg$agg)

  # Regional series: split cells by their region (unmatched -> NA, excluded).
  cell_region <- region_map[as.character(frame$COUNTRY)]
  for (k in seq_len(n_region)) {
    sub <- frame[!is.na(cell_region) & cell_region == region_info$region[k], ]
    regional_mat[[k]][i, ] <- annual_series(sub, cfg$start_year, cfg$agg)
  }

  rm(frame, sub, cell_region)
  gc()
}


# ==============================================================================
# 4. Build long tables for plotting
# ==============================================================================
long_global   <- make_long(global_mat)
long_regional <- lapply(regional_mat, make_long)

# x-axis ticks: the global wildfire panel starts in 2002; all other panels 2000.
x_breaks_2000 <- c(2000, 2005, 2010, 2015, 2019)
x_breaks_2002 <- c(2002, 2005, 2010, 2015, 2019)


# ==============================================================================
# 5. Global figure (six ECE panels, free y scales)
# ==============================================================================
jpeg_open(file.path(out_dir, "FigS3_EWE_exposure_year_Global.jpeg"))

global_panels <- lapply(seq_len(n_ece), function(i) {
  cfg <- ece_config[i, ]
  dat <- long_global[long_global$year >= cfg$start_year, ]
  xb  <- if (cfg$start_year == 2002) x_breaks_2002 else x_breaks_2000
  trend_panel(dat, cfg$yvar, paste0("(", LETTERS[i], ") ", cfg$family),
              cfg$ylab, xb)
})

grid.arrange(grobs = global_panels, nrow = 3)

dev.off()


# ==============================================================================
# 6. Regional figures (global + five regions, shared y scale per ECE)
# ==============================================================================
for (i in seq_len(n_ece)) {
  cfg <- ece_config[i, ]
  y_limits <- c(cfg$y_lo, cfg$y_hi)
  y_breaks <- seq(cfg$y_lo, cfg$y_hi, by = cfg$y_step)

  dat_global <- long_global[long_global$year >= cfg$start_year, ]

  panels <- c(
    list(trend_panel(dat_global, cfg$yvar, "(A) Global", cfg$ylab,
                     x_breaks_2000, y_limits, y_breaks)),
    lapply(seq_len(n_region), function(k) {
      dat <- long_regional[[k]][long_regional[[k]]$year >= cfg$start_year, ]
      trend_panel(dat, cfg$yvar,
                  paste0("(", LETTERS[k + 1], ") ", region_info$region[k]),
                  cfg$ylab, x_breaks_2000, y_limits, y_breaks)
    })
  )

  jpeg_open(file.path(out_dir, paste0("FigS",i+3,"_EWE_", cfg$fig, "_year_Regional.jpeg")))
  grid.arrange(grobs = panels, nrow = 3)
  dev.off()
}


# ==============================================================================
# 7. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
