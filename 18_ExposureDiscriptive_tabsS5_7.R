# ==============================================================================
# 18_ExposureDiscriptive_tabsS5_7.R
#
#
# In addition to Table S5 (annual exposure summary + within-family divergence),
# this version also generates two further descriptive tables from the same
# individual-level frame (frame_u5mr_ECE_month.rda):
#
#   Table S6 -- Summary of missing exposure data across cases and controls for
#               the 20 extreme-climatic-event definitions.
#   Table S7 -- Summary of exposure concordance / discordance status for
#               maternal (sibling) pairs exposed to each event.
#
# Definitions used here (see also the "DATA NOTICE" in 017_ExposureDiscriptive_1.R:
# the shipped frame is a ~1% microdata sample, so the counts approximate but do
# not exactly reproduce the published tables, which used the full dataset).
#
#   * "missing" exposure (Table S6): a child is counted as missing for an event
#     when ANY of its 13 monthly lag columns (lag_0m_<prefix> ... lag_12m_<prefix>)
#     is NA -- i.e. the gridded exposure could not be determined for at least one
#     of the 13 months of the child's exposure window.  Counts are reported
#     separately for all children (Total), deceased cases (u5mr == 1) and alive
#     controls (u5mr == 0).
#
#   * "exposed" status (Table S7): the child's annual (lag 0-12 month) exposure
#     metric is > 0 (i.e. the child experienced >= 1 event month / day).  The
#     annual metric follows make_annual() below (sum for heatwave / precipitation /
#     flood / cyclone / wildfire-months; SPEI <= -1 for drought; R34 excluding the
#     R64 core for "TC_34 knots").
#
#   * "pair" (Table S7): a mother (caseid) with at least one deceased case and at
#     least one alive control, where exposure is non-missing for both.
#     "Exposure-concordant" = the case and control(s) have the SAME exposure status
#     (both exposed, or both unexposed).  "Exposure-discordant" = they differ.
#
# OUTPUT (written to out_dir)
#   Table_S6_MissingExposure.csv             -- Table S6 (formatted "N (%)")
#   Table_S6_MissingExposure_raw.csv         -- Table S6 (numeric, long form)
#   Table_S7_ConcordanceDiscordance.csv      -- Table S7 (formatted "N (%)")
#   Table_S7_ConcordanceDiscordance_raw.csv  -- Table S7 (numeric, incl. 2x2 breakdown)
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(dplyr);library(tibble)

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
# 1. Load the analysis frame
# ==============================================================================
load(data_path)
data <- as.data.frame(frame_u5mr_ECE_month)


# ==============================================================================
# 2. Exposure configuration
# ==============================================================================
# One row per exposure definition. Column semantics:
#   family     : grouping used as the Table S5 section header
#   label      : event name shown in Table S5
#   unit       : measurement unit appended to the event name
#   prefix     : monthly lag-column prefix (lag_0m_<prefix> ... lag_12m_<prefix>)
#   type       : how the 13 monthly values are collapsed to one annual metric
#                "sum"   = sum over months (days / mm / exposed months)
#                "spei"  = months with SPEI <= -1 (drought), then sum
#                "storm" = tropical-storm months (R34 exposure excluding the
#                          hurricane core, i.e. R64), then sum
#                "mean"  = mean over exposed months (0 -> no event; used for
#                          wildfire speed / duration / spread)
#   drop_zero  : whether the Case/Control columns are restricted to children
#                with >= 1 exposed month (TRUE for count-type "months"
#                exposures). The within-family divergence is always computed
#                over all children in exposed families, irrespective of this
#                flag (matching the published Table S5).
ece_config <- tribble(
  ~family,                 ~label,              ~prefix,       ~type,    ~unit,    ~drop_zero,
  "Heatwave",              "HW_90%_2D",         "hw_90_2",     "sum",    "Days",   FALSE,
  "Heatwave",              "HW_90%_3D",         "hw_90_3",     "sum",    "Days",   FALSE,
  "Heatwave",              "HW_95%_2D",         "hw_95_2",     "sum",    "Days",   FALSE,
  "Heatwave",              "HW_95%_3D",         "hw_95_3",     "sum",    "Days",   FALSE,
  "Extreme precipitation", "EP_95%",            "extr_tp_95",  "sum",    "mm",     FALSE,
  "Extreme precipitation", "EP_99%",            "extr_tp_99",  "sum",    "mm",     FALSE,
  "Flood",                 "Total flood",       "flood_total", "sum",    "Months", TRUE,
  "Flood",                 "Large flood",       "flood_1",     "sum",    "Months", TRUE,
  "Flood",                 "Extreme flood",     "flood_2",     "sum",    "Months", TRUE,
  "Drought",               "1-Month Drought",   "spei_1",      "spei",   "Months", TRUE,
  "Drought",               "3-Month Drought",   "spei_3",      "spei",   "Months", TRUE,
  "Drought",               "6-Month Drought",   "spei_6",      "spei",   "Months", TRUE,
  "Drought",               "9-Month Drought",   "spei_9",      "spei",   "Months", TRUE,
  "Tropical cyclone",      "Total TC",          "total_TC",    "sum",    "Months", TRUE,
  "Tropical cyclone",      "TC_34 knots",       "R34exp",      "storm",  "Months", TRUE,
  "Tropical cyclone",      "TC_64 knots",       "R64exp",      "sum",    "Months", TRUE,
  "Wildfire",              "Wildfire exposure", "wf_exposure", "sum",    "Months", TRUE,
  "Wildfire",              "Wildfire speed",    "wf_speed",    "mean",   "km/Day", FALSE,
  "Wildfire",              "Wildfire duration", "wf_duration", "mean",   "Days",   FALSE,
  "Wildfire",              "Wildfire spread",   "wf_spread",   "mean",   "km2/Day", FALSE
)

# Lag-column suffixes: lag_0m (current month) ... lag_12m (12 months before).
lag_names <- c(paste0("lag_", 0:12, "m"))

# Order of the six families in Table S5.
family_order <- c("Heatwave", "Extreme precipitation", "Flood", "Drought",
                  "Tropical cyclone", "Wildfire")


# ==============================================================================
# 3. Core functions
# ==============================================================================

#' Collapse the 13 monthly lag columns of one exposure into an annual metric.
#'
#' @param prefix lag-column prefix (see `ece_config`).
#' @param type   aggregation rule ("sum" / "spei" / "storm" / "mean").
#' @param dat    analysis frame (must contain lag_0m_<prefix> ... lag_12m_<prefix>).
#' @return A numeric vector of length nrow(dat), one annual value per child.
make_annual <- function(prefix, type, dat) {
  cols <- paste0(lag_names, "_", prefix)
  stopifnot(all(cols %in% names(dat)))
  x <- as.matrix(dat[cols])

  if (type == "spei") {
    return(rowSums(1 * (x <= -1), na.rm = TRUE))        # drought months
  }
  if (type == "storm") {
    r64 <- as.matrix(dat[paste0(lag_names, "_R64exp")]) # hurricane core
    x <- x * (r64 != 1)                                 # exclude hurricane months
    return(rowSums(x, na.rm = TRUE))
  }
  if (type == "mean") {
    x[x == 0] <- NA                                     # no event -> excluded
    return(rowMeans(x, na.rm = TRUE))                   # mean over exposed months
  }
  rowSums(x, na.rm = TRUE)                              # "sum": days / mm / months
}

#' Case/control summary and within-family exposure divergence for one exposure.
#'
#' Case = deceased child (`u5mr == 1`); controls = alive siblings (`u5mr == 0`).
#' Families without any control, or without any exposed child, contribute no
#' usable sibling contrast and are dropped before summarising.
describe_exposure <- function(prefix, type, drop_zero, dat) {
  x <- make_annual(prefix, type, dat)
  d <- data.frame(caseid = dat$caseid, u5mr = dat$u5mr, x = x)
  d <- d[!is.na(d$x), ]                     # drop never-exposed children ("mean")

  # Family-level filters: keep families with >= 1 control and >= 1 exposed child.
  fam <- d %>%
    group_by(caseid) %>%
    summarise(n_control = sum(u5mr == 0, na.rm = TRUE),
              total_x    = sum(x, na.rm = TRUE),
              .groups    = "drop")
  keep <- fam$caseid[fam$n_control >= 1 & fam$total_x > 0]
  d <- d[d$caseid %in% keep, ]

  # Within-family divergence: case minus the mean of its matched controls,
  # computed over ALL children in each exposed family (no zero restriction).
  div <- d %>%
    group_by(caseid) %>%
    summarise(case_x = mean(x[u5mr == 1], na.rm = TRUE),
              ctrl_x = mean(x[u5mr == 0], na.rm = TRUE),
              .groups = "drop") %>%
    mutate(divergence = case_x - ctrl_x)

  # Case/control mean +/- SD. For count-type ("months") exposures, restrict the
  # summary to children with >= 1 exposed month.
  if (drop_zero) d <- d[d$x > 0, ]
  case_x    <- d$x[d$u5mr == 1]
  control_x <- d$x[d$u5mr == 0]

  data.frame(
    case_mean    = mean(case_x, na.rm = TRUE),
    case_sd      = sd(case_x, na.rm = TRUE),
    control_mean = mean(control_x, na.rm = TRUE),
    control_sd   = sd(control_x, na.rm = TRUE),
    div_median   = median(div$divergence, na.rm = TRUE),
    div_iqr      = IQR(div$divergence, na.rm = TRUE)
  )
}


# ==============================================================================
# 4. Compute descriptive statistics for every exposure
# ==============================================================================
stats <- lapply(seq_len(nrow(ece_config)), function(i) {
  r <- ece_config[i, ]
  out <- describe_exposure(r$prefix, r$type, r$drop_zero, data)
  out$family <- r$family
  out$label  <- r$label
  out$unit   <- r$unit
  out
})
result <- do.call(rbind, stats)
result <- result[, c("family", "label", "unit", "case_mean", "case_sd",
                     "control_mean", "control_sd", "div_median", "div_iqr")]

# Order rows by the Table S5 family order (within-family order is preserved).
result$family <- factor(result$family, levels = family_order)
result <- result[order(result$family), ]


# ==============================================================================
# 5. Format Table S5
# ==============================================================================
result$event    <- paste0(result$label, " (", result$unit, ")")
result$case_txt <- sprintf("%.1f±%.1f", result$case_mean,    result$case_sd)
result$ctrl_txt <- sprintf("%.1f±%.1f", result$control_mean, result$control_sd)
result$div_txt  <- sprintf("%.2f, %.2f",    result$div_median,   result$div_iqr)

table_s5 <- result[, c("family", "event", "case_txt", "ctrl_txt", "div_txt")]
names(table_s5) <- c("Extreme climatic event (family)",
                     "Extreme climatic event",
                     "Case",
                     "Control",
                     "Exposure divergence (Median, IQR)")


# ==============================================================================
# 6. Export Table S5
# ==============================================================================
write.csv(table_s5, file.path(out_dir, "TabS5_ExposureSummary.csv"), row.names = FALSE)

cat("\n=== Table S5: annual ECE exposure among cases and controls ===\n\n")
print(table_s5, row.names = FALSE)
cat("\nCase and Control are reported as mean ± SD.\n")
cat("Exposure divergence is the grand median (IQR) of within-family\n")
cat("case-minus-control differences, pooled across all strata.\n")


# ==============================================================================
# 7. Table S6 -- Summary of missing exposure data across cases and controls
# ==============================================================================
# A child is "missing" for an event when any of its 13 monthly lag columns is NA.
missing_flag <- function(prefix, dat) {
  cols <- paste0(lag_names, "_", prefix)
  stopifnot(all(cols %in% names(dat)))
  x <- as.matrix(dat[cols])
  rowSums(is.na(x)) > 0
}

# Table S6 display names (same 20-row order as ece_config).
label_s6 <- c(
  "HW_Month_90%&2D", "HW_Month_90%&3D", "HW_Month_95%&2D", "HW_Month_95%&3D",
  "EP_Month_95%", "EP_Month_99%",
  "Total Flood", "Large Flood", "Extreme Flood",
  "1-Month Drought", "3-Month Drought", "6-Month Drought", "9-Month Drought",
  "Total TC", "TC_34 knots", "TC_64 knots",
  "Wildfire", "Wildfire Duration", "Wildfire Speed", "Wildfire Spread"
)

fmt_n_pct <- function(n, denom) {
  pct <- if (denom == 0) 0 else n / denom * 100
  sprintf("%s (%s%%)", prettyNum(n, big.mark = ",", scientific = FALSE),
          round(pct, 1))
}

fmt_n <- function(n) prettyNum(n, big.mark = ",", scientific = FALSE)

s6_raw <- lapply(seq_len(nrow(ece_config)), function(i) {
  miss <- missing_flag(ece_config$prefix[i], data)
  data.frame(
    label        = label_s6[i],
    total_miss   = sum(miss),
    case_miss    = sum(miss & data$u5mr == 1),
    control_miss = sum(miss & data$u5mr == 0),
    stringsAsFactors = FALSE
  )
})
s6_raw <- do.call(rbind, s6_raw)

# Group denominators (constant across events).
n_total   <- nrow(data)
n_case    <- sum(data$u5mr == 1)
n_control <- sum(data$u5mr == 0)

table_s6 <- data.frame(
  label   = s6_raw$label,
  total   = mapply(fmt_n_pct, s6_raw$total_miss, n_total),
  case    = mapply(fmt_n_pct, s6_raw$case_miss, n_case),
  control = mapply(fmt_n_pct, s6_raw$control_miss, n_control),
  stringsAsFactors = FALSE
)
names(table_s6) <- c("Extreme climatic event",
                     sprintf("Total (n = %s)", prettyNum(n_total, big.mark = ",")),
                     sprintf("Case (n = %s)", prettyNum(n_case, big.mark = ",")),
                     sprintf("Control (n = %s)", prettyNum(n_control, big.mark = ",")))

write.csv(table_s6, file.path(out_dir, "TabS6_MissingExposure.csv"), row.names = FALSE)

cat("\n=== Table S6: missing exposure data (cases vs controls) ===\n\n")
print(table_s6, row.names = FALSE)


# ==============================================================================
# 8. Table S7 -- Exposure concordance / discordance for maternal pairs
# ==============================================================================
# Table S7 display names. The four wildfire rows all report the "Wildfire
# exposure" status (speed / duration / spread share the same exposure window).
label_s7 <- c(
  "HW_Month_90%&2D", "HW_Month_90%&3D", "HW_Month_95%&2D", "HW_Month_95%&3D",
  "EP_Month_95%", "EP_Month_99%",
  "Total Flood", "Large Flood", "Extreme Flood",
  "1-Month Drought", "3-Month Drought", "6-Month Drought", "9-Month Drought",
  "Total TC", "TC_34 knots", "TC_64 knots",
  "Wildfire Exposure", "Wildfire Exposure", "Wildfire Exposure", "Wildfire Exposure"
)

# For the concordance table the four wildfire rows all use wf_exposure (sum type).
s7_prefix <- ece_config$prefix
s7_type   <- ece_config$type
s7_prefix[17:20] <- "wf_exposure"
s7_type[17:20]   <- "sum"

s7_raw <- lapply(seq_len(nrow(ece_config)), function(i) {
  ann <- make_annual(s7_prefix[i], s7_type[i], data)
  miss <- missing_flag(s7_prefix[i], data)     # any of the 13 lag months is NA
  d2  <- data.frame(caseid = data$caseid, u5mr = data$u5mr, x = ann)
  d2  <- d2[!miss, ]                           # drop children with missing exposure
  d2$exposed <- d2$x > 0

  fam <- d2 %>%
    group_by(caseid) %>%
    summarise(
      ncase        = sum(u5mr == 1),
      nctrl        = sum(u5mr == 0),
      case_exp     = any(exposed[u5mr == 1]),
      ctrl_exp     = any(exposed[u5mr == 0]),
      .groups      = "drop"
    )
  # Maternal pairs: >= 1 case and >= 1 control.
  fam <- fam[fam$ncase >= 1 & fam$nctrl >= 1, ]

  data.frame(
    label          = label_s7[i],
    pairs          = nrow(fam),
    both_exposed   = sum(fam$case_exp &  fam$ctrl_exp),
    both_unexposed = sum(!fam$case_exp & !fam$ctrl_exp),
    case_only      = sum(fam$case_exp & !fam$ctrl_exp),
    control_only   = sum(!fam$case_exp & fam$ctrl_exp),
    stringsAsFactors = FALSE
  )
})
s7_raw <- do.call(rbind, s7_raw)

# Concordant = same exposure status; discordant = different status.
s7_raw$concordant <- s7_raw$both_exposed + s7_raw$both_unexposed
s7_raw$discordant <- s7_raw$case_only + s7_raw$control_only

table_s7 <- data.frame(
  label     = s7_raw$label,
  pairs     = fmt_n(s7_raw$pairs),
  concordant = mapply(fmt_n_pct, s7_raw$concordant, s7_raw$pairs),
  discordant = mapply(fmt_n_pct, s7_raw$discordant, s7_raw$pairs),
  stringsAsFactors = FALSE
)
names(table_s7) <- c("Extreme climatic event", "Pairs (mother)",
                     "Exposure-concordant pairs", "Exposure-discordant pairs")

write.csv(table_s7, file.path(out_dir, "TabS7_ConcordanceDiscordance.csv"), row.names = FALSE)

cat("\n=== Table S7: exposure concordance / discordance (maternal pairs) ===\n\n")
print(table_s7, row.names = FALSE)


# ==============================================================================
# 9. Session information (for reproducibility)
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
