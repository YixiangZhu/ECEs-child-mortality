# ==============================================================================
# 14_GlobalRegionalBurdenPlot_region_figS31_tabS22.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
# (Region-stratified burden post-processing -- Figure S31, Table S22)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To tabulate and visualise the global / regional burden of under-5 mortality
#    attributable to the six extreme climatic events, using the REGION-STRATIFIED
#    results produced by 13_BurdenCalculation_region.R (Table S22 / Figure S31).
#
# Methodology Overview:
#    This script does NOT recompute anything -- it only READS the per-exposure
#    region-level process files already written to burden_result_new2_region
#    (frame_<key>_LIMCs_year_10_case[_low/_high]_region.csv) and produces:
#      TabS22  -- combined region-level burden summary table (annual attributable
#                 deaths / deaths per 100,000 under-5 person-year / PAF, %), with
#                 thousand separators in the annual-death cells.
#      FigS31  -- region stacked bar of the attributable child mortality rate
#                 (per 100,000 under-5 person-years).
#    The wildfire burden is annualised over 18 years (2002-2019), all others 20.
#
#
# REQUIREMENTS
#    R       >= 4.0.3
#    ggplot2      == 3.4.1
#    dplyr        == 1.1.1
#    RColorBrewer == 1.1.3
#    stringr      == 1.5.0
#
# INPUT
#    frame_<key>_LIMCs_year_10_case[_low/_high]_region.csv
#
# OUTPUT
#    FigS31_BurdenRate_ByRegion_region.png
#    TabS22_GlobalRegionalBurden_region.csv
# ==============================================================================

library(ggplot2);library(dplyr);library(RColorBrewer);library(stringr)

# ------------------------------------------------------------------------------
# 1. Paths and configuration
# ------------------------------------------------------------------------------
src_dir <- "Repository/Output/burden_result_regional"
out_dir <- "Repository/Output/burden_result_regional"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Six exposures: key (file prefix), display name, number of reportable years.
ece <- data.frame(
  key   = c("hw_90_2_days", "tp_99", "flood", "spei_6", "tc_R34", "wf_exposure"),
  ewe   = c("Heatwave", "Extreme precipitation", "Flood", "Drought",
            "Tropical Cyclone", "Wildfire"),
  years = c(20L, 20L, 20L, 20L, 20L, 18L),
  stringsAsFactors = FALSE
)

ewe_table <- c("Heatwave", "Extreme Precipitation", "Flood", "Drought",
               "Tropical Cyclone", "Wildfire")

ordered_regions <- c("Global", "Sub-Saharan Africa", "Central & South Asia",
                     "Europe, Middle East & North Africa", "Latin America & Caribbean",
                     "East Asia & Pacific")

ewe_cols <- brewer.pal(9, "RdBu")[c(2:4, 6:8)]

# ------------------------------------------------------------------------------
# 2. Read region-level results (point / lower / upper) for all six exposures
# ------------------------------------------------------------------------------
read_region <- function(key, ewe, years) {
  est  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_region.csv")),
                   check.names = FALSE)
  low  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_low_region.csv")),
                   check.names = FALSE)
  high <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_high_region.csv")),
                   check.names = FALSE)

  # Columns (by position): REGION, <var>_case, case, pop_sum, <var>_burden, <var>_PAF
  colnames(est)  <- c("region", "case_total", "case", "pop_sum", "burden", "PAF")
  colnames(low)  <- c("region", "case_low",  "case", "pop_sum", "burden_low",  "PAF_low")
  colnames(high) <- c("region", "case_high", "case", "pop_sum", "burden_high", "PAF_high")

  out <- est
  out$EWE            <- ewe
  out$years          <- years
  out$case_mean      <- out$case_total / years        # annual attributable deaths
  out$case_mean_low  <- low$case_low  / years
  out$case_mean_high <- high$case_high / years
  out$case_low       <- low$case_low
  out$case_high      <- high$case_high
  out$burden_low     <- low$burden_low
  out$burden_high    <- high$burden_high
  out$PAF_low        <- low$PAF_low
  out$PAF_high       <- high$PAF_high
  out
}

region_all <- do.call(rbind, lapply(seq_len(nrow(ece)),
                                    function(i) read_region(ece$key[i], ece$ewe[i], ece$years[i])))

region_all <- region_all[!is.na(region_all$region), ]
region_all$region <- factor(region_all$region, levels = ordered_regions)
region_all$EWE    <- factor(region_all$EWE, levels = ece$ewe)
region_all <- region_all[order(region_all$region, region_all$EWE), ]

# ------------------------------------------------------------------------------
# 3. FigS31 -- region stacked bar of attributable mortality rate (per 100k)
# ------------------------------------------------------------------------------
region_rate <- do.call(rbind, lapply(seq_len(nrow(ece)), function(i) {
  d <- read.csv(file.path(src_dir, paste0("frame_", ece$key[i], "_LIMCs_year_10_case_region.csv")),
                check.names = FALSE)
  colnames(d) <- c("region", "case_total", "case", "pop_sum", "burden", "PAF")
  data.frame(region = d$region, EWE = ece$ewe[i], burden = d$burden)
}))
region_rate <- region_rate[!is.na(region_rate$region), ]
region_rate$region <- factor(region_rate$region, levels = ordered_regions)
region_rate$EWE    <- factor(region_rate$EWE, levels = ece$ewe)

p_bar <- ggplot(region_rate, aes(x = region, y = burden, fill = EWE)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8,
           color = "black", linewidth = 0.3) +
  theme_bw() +
  scale_fill_manual(values = ewe_cols) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "Region",
       y = "Attributable child mortality rate\n(per 100,000 under-5 person-years)",
       fill = "Extreme climatic events") +
  theme(axis.title   = element_text(face = "bold", color = "black", size = 18, lineheight = 1),
        axis.text    = element_text(face = "bold", color = "black", size = 13),
        legend.text  = element_text(face = "bold", color = "black", size = 16, lineheight = 1),
        legend.title = element_text(face = "bold", color = "black", size = 18, lineheight = 1),
        legend.position = "bottom")

ggsave(file.path(out_dir, "FigS31_Burden_ByRegion_region.png"),
       plot = p_bar, width = 12.5, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 4. TabS22 -- combined region-level burden summary table
# ------------------------------------------------------------------------------
fmt_annual <- function(est, low, high) {
  f <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  paste0(f(est), " (", f(low), ", ", f(high), ")")
}
fmt_1d <- function(est, low, high) {
  paste0(sprintf("%.1f", est), " (", sprintf("%.1f", low), ", ", sprintf("%.1f", high), ")")
}
fmt_2d <- function(est, low, high) {
  paste0(sprintf("%.2f", est), " (", sprintf("%.2f", low), ", ", sprintf("%.2f", high), ")")
}

build_section <- function(est_col, low_col, high_col, fmt) {
  wide <- matrix("", nrow = length(ordered_regions), ncol = nrow(ece) + 2)
  colnames(wide) <- c("Region", ewe_table, "Total")
  for (r in seq_along(ordered_regions)) {
    wide[r, 1] <- ordered_regions[r]
    sub_all <- region_all[region_all$region == ordered_regions[r], ]
    for (e in seq_len(nrow(ece))) {
      sub <- sub_all[sub_all$EWE == ece$ewe[e], ]
      if (nrow(sub) == 1)
        wide[r, e + 1] <- fmt(sub[[est_col]], sub[[low_col]], sub[[high_col]])
    }
    if (nrow(sub_all) == nrow(ece)) {
      wide[r, nrow(ece) + 2] <- fmt(sum(sub_all[[est_col]],  na.rm = TRUE),
                                    sum(sub_all[[low_col]],  na.rm = TRUE),
                                    sum(sub_all[[high_col]], na.rm = TRUE))
    }
  }
  as.data.frame(wide, stringsAsFactors = FALSE, check.names = FALSE)
}

section_header <- function(label, ncol) {
  row <- as.list(rep("", ncol))
  row[[1]] <- label
  names(row) <- c("Region", ewe_table, "Total")
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
}

tab_annual <- build_section("case_mean", "case_mean_low", "case_mean_high", fmt_annual)
tab_burden <- build_section("burden",    "burden_low",    "burden_high",    fmt_1d)
tab_paf    <- build_section("PAF",       "PAF_low",       "PAF_high",       fmt_2d)

tabS22 <- rbind(
  section_header("Annual attributable death", nrow(ece) + 2), tab_annual,
  section_header("Attributable deaths per 100,000 under-5 person-year", nrow(ece) + 2), tab_burden,
  section_header("PAF (%)", nrow(ece) + 2), tab_paf
)

write.csv(tabS22, file.path(out_dir, "TabS22_GlobalRegionalBurden_region.csv"), row.names = FALSE)

cat("\n=== Table S22: global and regional burden (region-stratified estimates) ===\n\n")
print(tabS22, row.names = FALSE)
cat("\nFiles written to", out_dir, "\n")
