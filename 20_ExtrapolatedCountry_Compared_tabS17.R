# ==============================================================================
# 019_ExtrapolatedCountry_Compared_tabS17.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
# (exposure summary statistics, DHS-covered vs. extrapolated LMICs)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    Summarise the extreme climatic event (ECE) exposure across the study
#    countries, reporting the mean +/- standard deviation within the total, DHS
#    and extrapolated LMIC country groups over 2000-2019.

# Exposure definition:
#    For each ECE, the gridded monthly time series is summed to an ANNUAL value
#    per grid-year (heatwave days, extreme-precipitation mm, or the number of
#    flood / drought / tropical-cyclone / wildfire months).  The summary is the
#    mean +/- SD of these annual values over the EXPOSED grid-years only
#    (annual value > 0).
#
# Country classification:
#    `mortality_country_merge.csv` carries the `DHS` flag (1 = DHS survey
#    available, 0 = extrapolated for global coverage).
#
# INPUT
#    frame_<key>_LIMCs_month_10.rda -- gridded monthly ECE exposure (703714 grids)
#    mortality_country_merge.csv    -- COUNTRY / DHS flag / REGION
#
# OUTPUT
#    TabS17_ExtrapolatedCountry_Compared.csv -- exposure rows of Table S17

# ==============================================================================

# --- paths (edit these to point at your own copies) ---------------------------
data_path <- file.path(
  "Repository/Data",
  "frame_u5mr_ECE_month.rda"
)

out_dir <- file.path(
  "Repository/Output", "Result"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. Country classification (DHS vs. extrapolated)
# ==============================================================================
country_merge <- read.csv(file.path(src_dir, "mortality_country_merge.csv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
names(country_merge) <- c("Num", "DHS", "COUNTRY", "REGION")
cc <- country_merge$COUNTRY
dhs_flag <- setNames(country_merge$DHS, country_merge$COUNTRY)

cat("Countries in mortality_country_merge.csv:", nrow(country_merge), "\n")
cat("  DHS (=1):", sum(country_merge$DHS == 1),
    " | Extrapolated (=0):", sum(country_merge$DHS == 0), "\n")


# ==============================================================================
# 2. Exposure configuration (six ECE, in Table S17 row order)
# ==============================================================================
exposure_config <- data.frame(
  frame_file = c("frame_hw_90_2_days_LIMCs_month_10", "frame_EP_99_cum_LIMCs_month_10",
                 "frame_flood_total_LIMCs_month_10",  "frame_dr_6m_LIMCs_month_10",
                 "frame_total_TC_LIMCs_month_10",     "frame_wf_exposure_LIMCs_month_10"),
  label      = c("Heatwave (Days)", "Extreme precipitation (mm)", "Flood (Months)",
                 "Drought (Months)", "Tropical cyclone (Months)", "Wildfire (Months)"),
  start_year = c(2000L, 2000L, 2000L, 2000L, 2000L, 2002L),
  stringsAsFactors = FALSE
)

monthly_year <- function(cn) {
  as.integer(substr(sub("^X", "", cn), 1, 4))
}

#' Annual exposure values over exposed grid-years, tagged with the country.
#' Returns a data.frame(country, value) with one row per exposed grid-year
#' (annual value > 0).
aggregate_exposure <- function(frame_file, start_year) {
  e <- new.env()
  load(file.path(src_dir, paste0(frame_file, ".rda")), envir = e)
  frame <- e[[ls(e)[1]]]
  rm(e)
  frame <- frame[order(frame$OBJECTID), ]
  valid <- frame$COUNTRY != "South Korea"

  monthly_cols <- grep("^X?20[0-9]{2}[._-][0-9]{2}", colnames(frame))
  yrs  <- monthly_year(colnames(frame)[monthly_cols])
  keep <- yrs %in% seq.int(start_year, 2019)

  m   <- as.matrix(frame[valid, monthly_cols[keep]])
  yk  <- yrs[keep]
  uy  <- sort(unique(yk))
  country <- frame$COUNTRY[valid]

  # Monthly -> annual sums (heatwave days / precipitation mm / exposed months).
  ann <- sapply(uy, function(y) rowSums(m[, yk == y, drop = FALSE], na.rm = TRUE))

  exposed <- ann > 0
  exposed[is.na(exposed)] <- FALSE
  vals <- as.vector(ann)[as.vector(exposed)]
  cnt  <- rep(country, length(uy))[as.vector(exposed)]

  data.frame(country = cnt, value = vals, stringsAsFactors = FALSE)
}

#' Mean +/- SD of `value` over the observations in each country group.
summarise_values <- function(df) {
  g <- list(
    Total        = df$country %in% cc,
    DHS          = df$country %in% cc[dhs_flag == 1],
    Extrapolated = df$country %in% cc[dhs_flag == 0]
  )
  sapply(g, function(idx) {
    v <- df$value[idx]
    c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE))
  })
}

exposure_stats <- lapply(seq_len(nrow(exposure_config)), function(i) {
  cat("Aggregating:", exposure_config$label[i], "\n")
  summarise_values(aggregate_exposure(exposure_config$frame_file[i],
                                      exposure_config$start_year[i]))
})
names(exposure_stats) <- exposure_config$label


# ==============================================================================
# 3. Format the exposure rows of Table S17 and export
# ==============================================================================
fmt_msd <- function(x) paste0(sprintf("%.1f", x["mean"]), "±", sprintf("%.1f", x["sd"]))
group_cols <- c("Total LMICs", "DHS LMICs", "Extrapolated LMICs")

tabS17 <- data.frame(
  Indicator = c("Extreme climatic events*", paste0("  ", exposure_config$label)),
  stringsAsFactors = FALSE, check.names = FALSE
)

for (g in seq_along(group_cols)) {
  cells <- sapply(exposure_stats, function(st) fmt_msd(st[, g]))
  tabS17[[group_cols[g]]] <- c("", cells)
}

write.csv(tabS17, file.path(out_dir, "TabS17_ExtrapolatedCountry_Compared.csv"),
          row.names = FALSE)

cat("\n=== Table S17 exposure: DHS-covered vs. extrapolated LMICs (mean +/- SD, 2000-2019) ===\n\n")
print(tabS17, row.names = FALSE)
