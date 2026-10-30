# ==============================================================================
# 11_GlobalRegionalBurdenPlot_fig4_tabS19_figsS29_30.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
# (Global / regional burden plots and tables -- Figure 4, Table S19, Figures S29-S30)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To visualise and tabulate the global / regional burden of under-5 mortality
#    attributable to the six extreme climatic events, using the region-level
#    results produced by 10_BurdenCalculation.R.
#
# Methodology Overview:
#    1. Read the region-level point / lower / upper burden results
#       (frame_<key>_LIMCs_year_10_case[_low/_high]_region.csv) and the
#       PAF / burden time-trend files (frame_<key>_LIMCs_year_10_{PAF,burden}_trend.csv).
#    2. Drop the NA region row (unmapped grids, 0 attributable deaths) before
#       plotting / summarising.
#    3. Generate:
#         Fig4A  -- region stacked bar of the attributable child mortality rate
#                   (per 100,000 under-5 person-years) by ECE and region.
#         (Fig4B -- global spatial distribution maps of attributable child mortality rates 
#                    across LIMCs, generated via ArcGIS Pro.)
#         Fig4C  -- temporal trend of the population attributable fraction (PAF, %)
#                   for each ECE, 2002-2019 (stacked area).
#         TabS19 -- combined region-level burden summary table (annual attributable
#                   deaths / deaths per 100,000 under-5 person-year / PAF, %),
#                   values reported as "estimate (95% CI)" with fixed decimals.
#         FigS29 -- 7-panel temporal trend of PAF (%) for each ECE + total.
#         FigS30 -- 7-panel temporal trend of the attributable mortality rate for
#                   each ECE + total.
#
# Differences from the original 011 script:
#   (1) Input  path  -> burden_result_new2  (Revised_3 / Reanalysis).
#   (2) Exposure keys -> hw_90_2_days / spei_6 / wf_exposure (were hw / dr_6m / wf).
#   (3) mortality_country_merge.csv is now 4 columns (COUNTRY = 3, REGION = 4).
#   (4) Output path  -> result_raw.
#   (5) TabS19 is exported as one combined table (annual / burden / PAF).
#
# REQUIREMENTS
#    R            >= 4.0.3
#    ggplot2      == 3.4.1
#    dplyr        == 1.1.1
#    RColorBrewer == 1.1.3
#    stringr      == 1.5.0
#    patchwork    == 1.2.0
#
# INPUT
#    frame_<key>_LIMCs_year_10_case[_low/_high]_region.csv
#    frame_<key>_LIMCs_year_10_{PAF,burden}_trend.csv
#
# OUTPUT
#    Fig4A_BurdenRate_ByRegion.png
#    Fig4C_PAF_Trend.png
#    TabS19_GlobalRegionalBurden.csv
#    FigS29_PAF_Trend.png
#    FigS30_Burden_Trend.png
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------------------------
library(ggplot2);library(dplyr);library(RColorBrewer);library(stringr);library(patchwork)

# ------------------------------------------------------------------------------
# 2. Paths and configuration
# ------------------------------------------------------------------------------
src_dir  <- "Repository/Output/burden_result/"
out_dir  <- "Repository/Output/burden_result/"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Six exposures: key (output prefix), internal display name, number of reportable years.
ece <- data.frame(
  key   = c("hw_90_2_days", "tp_99", "flood", "spei_6", "tc_R34", "wf_exposure"),
  ewe   = c("Heatwave", "Extreme precipitation", "Flood", "Drought",
            "Tropical Cyclone", "Wildfire"),
  years = c(20L, 20L, 20L, 20L, 20L, 18L),
  stringsAsFactors = FALSE
)

# Display names used in Table S19 column headers.
ewe_table <- c("Heatwave", "Extreme Precipitation", "Flood", "Drought",
               "Tropical Cyclone", "Wildfire")

ordered_regions <- c("Global", "Sub-Saharan Africa", "Central & South Asia",
                     "Europe, Middle East & North Africa", "Latin America & Caribbean",
                     "East Asia & Pacific")

ewe_cols <- brewer.pal(9, "RdBu")[c(2:4, 6:8)]   # 6 diverging colours

# ------------------------------------------------------------------------------
# 3. Read region-level results (point / lower / upper) for all six exposures
# ------------------------------------------------------------------------------
read_region <- function(key, ewe, years) {
  est  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_region.csv")),
                   check.names = FALSE)
  low  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_low_region.csv")),
                   check.names = FALSE)
  high <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_high_region.csv")),
                   check.names = FALSE)

  # New CSVs: REGION, <var>_case, case, pop_sum, <var>_burden, <var>_PAF
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

# Drop the NA region row (unmapped grids, 0 attributable deaths).
region_all <- region_all[!is.na(region_all$region), ]
region_all$region <- factor(region_all$region, levels = ordered_regions)
region_all$EWE    <- factor(region_all$EWE, levels = ece$ewe)
region_all <- region_all[order(region_all$region, region_all$EWE), ]

# ------------------------------------------------------------------------------
# 4. Fig4A -- region stacked bar of attributable mortality rate (per 100k)
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

p_4A <- ggplot(region_rate, aes(x = region, y = burden, fill = EWE)) +
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

ggsave(file.path(out_dir, "Fig4A_Burden_ByRegion.png"),
       plot = p_4A, width = 12.5, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Fig4C -- temporal trend of the population attributable fraction (PAF, %)
#            for each extreme climatic event, 2002-2019.
# ------------------------------------------------------------------------------
# Each PAF-trend CSV has 6 rows (one per EWE) x 21 cols (EWE + year_2000..2019);
# only the exposure's own row is filled. Take row i from exposure i.
paf_long <- data.frame()
for (i in seq_len(nrow(ece))) {
  f <- file.path(src_dir, paste0("frame_", ece$key[i], "_LIMCs_year_10_PAF_trend.csv"))
  d <- read.csv(f, check.names = FALSE)
  yrs <- 2002:2019
  paf_long <- rbind(paf_long,
                    data.frame(year = yrs,
                               EWE  = ece$ewe[i],
                               PAF  = as.numeric(d[i, paste0("year_", yrs)])))
}
paf_long$EWE <- factor(paf_long$EWE, levels = ece$ewe)

p_4C <- ggplot(paf_long, aes(x = year, y = PAF, fill = EWE)) +
  geom_area() +
  scale_fill_manual(values = ewe_cols) +
  scale_x_continuous(name = "Year",
                     breaks = c(2002, 2005, 2010, 2015, 2019),
                     labels = c("2002", "2005", "2010", "2015", "2019")) +
  labs(x = "Year", y = "Population attributable fraction (%)",
       fill = "Extreme climatic events") +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", color = "black", size = 18, lineheight = 1),
        axis.text    = element_text(face = "bold", color = "black", size = 18),
        legend.text  = element_text(face = "bold", color = "black", size = 18, lineheight = 1),
        legend.title = element_text(face = "bold", color = "black", size = 18, lineheight = 1),
        legend.position = "bottom")

ggsave(file.path(out_dir, "Fig4C_PAF_Trend.png"),
       plot = p_4C, width = 12, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 6. TabS19 -- combined region-level burden summary table
# ------------------------------------------------------------------------------
# Three stacked sections, one per metric:
#   "Annual attributable death"                     -> case_mean (annual)
#   "Attributable deaths per 100,000 under-5 ..."   -> burden (per 100k)
#   "PAF (%)"                                       -> PAF (%)
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

tab_annual <- build_section("case_mean",      "case_mean_low",  "case_mean_high",  fmt_annual)
tab_burden <- build_section("burden",         "burden_low",     "burden_high",     fmt_1d)
tab_paf    <- build_section("PAF",            "PAF_low",        "PAF_high",        fmt_2d)

tabS19 <- rbind(
  section_header("Annual attributable death",                       nrow(ece) + 2), tab_annual,
  section_header("Attributable deaths per 100,000 under-5 person-year", nrow(ece) + 2), tab_burden,
  section_header("PAF (%)",                                         nrow(ece) + 2), tab_paf
)

write.csv(tabS19, file.path(out_dir, "TabS19_GlobalRegionalBurden.csv"), row.names = FALSE)

cat("\n=== Table S19: global and regional burden of under-5 mortality attributable to ECE ===\n\n")
print(tabS19, row.names = FALSE)


# ------------------------------------------------------------------------------
# 7. FigS29 and FigS30 -- temporal trend of PAF (%) and burden for each ECE (7-panel)
# ------------------------------------------------------------------------------
read_trend_metric <- function(metric) {
  # metric: "PAF" or "burden"
  long <- data.frame()
  for (i in seq_len(nrow(ece))) {
    f <- file.path(src_dir, paste0("frame_", ece$key[i], "_LIMCs_year_10_", metric, "_trend.csv"))
    d <- read.csv(f, check.names = FALSE)
    yrs <- 2000:2019
    long <- rbind(long,
                  data.frame(year = yrs,
                             EWE  = ece$ewe[i],
                             value = as.numeric(d[i, paste0("year_", yrs)])))
  }
  # Total = sum of the six exposures (wildfire is NA before 2002).
  tot <- long %>%
    group_by(year) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(EWE = "Total")
  long <- rbind(long, tot)
  long$EWE <- factor(long$EWE, levels = c(ece$ewe, "Total"))
  long
}

# Widen to year x EWE for the per-panel plots.
to_wide <- function(long) {
  w <- long %>% tidyr::pivot_wider(names_from = EWE, values_from = value) %>% as.data.frame()
  w[order(w$year), ]
}

panel_theme <- theme(
  axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
  axis.text.x  = element_text(size = 13), axis.text.y  = element_text(size = 13)
)

panel_plot <- function(wide, col, ylab, start_year) {
  d <- wide[wide$year >= start_year, ]
  breaks <- if (start_year > 2000) c(2002, 2005, 2010, 2015, 2019) else c(2000, 2005, 2010, 2015, 2019)
  ggplot(d) +
    geom_point(aes(x = year, y = .data[[col]])) +
    geom_line(aes(x = year, y = .data[[col]]), linewidth = 1) +
    theme_bw() +
    labs(title = col) +
    stat_smooth(aes(x = year, y = .data[[col]]), method = "glm",
                formula = y ~ poly(x, 3), se = TRUE, linewidth = 1,
                col = "royalblue3", fill = "darkslategray3", lty = 2) +
    xlab("Year") + ylab(ylab) + panel_theme +
    scale_x_continuous(breaks = breaks, labels = as.character(breaks))
}

make_trend_figure <- function(metric, ylab, outfile) {
  long <- read_trend_metric(metric)
  wide <- to_wide(long)

  fig_names <- ece$ewe
  # panels p1..p6 use 2000-2019 (wildfire p6 uses 2002-2019); Total p7 uses 2002-2019.
  ps <- lapply(seq_along(fig_names), function(k) {
    panel_plot(wide, fig_names[k], ylab, start_year = if (k == 6) 2002 else 2000)
  })
  p_total <- panel_plot(wide, "Total", ylab, start_year = 2002)

  bottom_grid <- (ps[[1]] | ps[[2]]) / (ps[[3]] | ps[[4]]) / (ps[[5]] | ps[[6]])
  final <- p_total / bottom_grid + plot_layout(heights = c(0.8, 2))

  jpeg(file.path(out_dir, outfile), res = 300, width = 2600, height = 3400)
  print(final)
  dev.off()
  cat("written:", outfile, "\n")
}

make_trend_figure("PAF",    "PAF (%)",                    "FigS29_PAF_Trend.png")
make_trend_figure("burden", "Deaths per 100,000",         "FigS30_Burden_Trend.png")

cat("\nAll figures and tables written to:", out_dir, "\n")
