# ==============================================================================
# 015_IndividualDescriptive.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#   To summarise the individual characteristics of deceased children (cases) and
#   their matched alive siblings (controls), this script
#     (1) recodes the individual-level covariates into the Table S4 categories,
#     (2) computes the frequency and column percentage of each category among
#         cases and controls,
#     (3) obtains the between-group difference P value from a sibling-matched
#         conditional logistic regression; for multi-category variables an
#         overall Wald test (Type III) is conducted,
#     (4) exports Table S4 (Summary statistics and the between-group differences
#         for deceased children and their matched alive children).
#
# Methodology Overview:
#   1. Sibling-matched case/control split: within each maternal lineage
#      (`caseid`), the deceased child is the case (`u5mr == 1`) and its alive
#      siblings are the controls (`u5mr == 0`).
#   2. Categorical recoding: born order and twin count are top-coded (>=6 and
#      >=2, respectively); mother's age at the child's death is cut into three
#      bands (<18, 18-35, >35). The remaining covariates (region, residence,
#      education, water / sanitation infrastructure, floor / roof / wall
#      material, birth interval) are already coded as the Table S4 categories in
#      the analysis frame.
#   3. Between-group difference: for covariates that can vary within a household
#      (infant's gender, born order, twin, mother's age, birth interval), a
#      conditional logistic regression stratified on `caseid` is fitted and the
#      overall Wald Type III P value is extracted. Covariates invariant within a
#      household (region, residence, mother's education, water, sanitation,
#      floor, roof, wall material) are balanced by design through the
#      within-household matching and are therefore not tested (P value "-").
#
# DATA NOTICE
#   This script ships with a representative ~1% microdata sample for each of the
#   168 included surveys. The summary statistics therefore approximate, but do
#   not exactly reproduce, the published Table S4 (which used the full dataset).
#
# DEPENDENCIES
#   R        >= 4.0.3
#   dplyr    >= 1.0.0
#   tibble    == 3.1.8
#   survival >= 3.3
#   car      >= 3.0
#
# INPUT
#   frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      the outcome `u5mr` (1 = deceased case, 0 = alive control), the matching
#      stratum `caseid`, and the individual / household covariates.
#
# OUTPUT
#   Table_S4_IndividualSummary.csv     -- formatted Table S4 (count (%) per
#      category with the between-group difference P value).
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(dplyr);library(tibble);library(survival);library(car)

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

# Cases / controls denominator (used for the column percentages).
n_case  <- sum(data$u5mr == 1, na.rm = TRUE)
n_ctrl  <- sum(data$u5mr == 0, na.rm = TRUE)


# ==============================================================================
# 2. Recode continuous / semi-continuous covariates into the Table S4 categories
# ==============================================================================

# Born order: top-code at 6+ (raw values range 1 .. 17).
data$bord6 <- ifelse(data$bord >= 6, 6, data$bord)

# Twin: top-code at 2+ (raw values 0 .. 3).
data$twin2 <- ifelse(data$twin >= 2, 2, data$twin)

# Mother's age at the child's death (for cases) or at the matched time point
# (for controls): three bands (<18, 18-35, >35).
data$mother_age_level <- ifelse(data$maternal_age_death < 18, 1,
                         ifelse(data$maternal_age_death <= 35, 2, 3))


# ==============================================================================
# 3. Table S4 configuration
# ==============================================================================
# One entry per Table S4 row-block. Column semantics:
#   col    : analysis-frame column.
#   label  : section title shown in the table.
#   tested : whether the case/control difference is statistically tested with a
#            conditional logistic regression. FALSE for household-invariant
#            covariates that are balanced by design (P value reported as "-").
#   levels : named vector <display-label> = <data value>, in display order.
# ==============================================================================
tab_config <- list(
  region         = list(col = "region", label = "Region", tested = FALSE,
    levels = c("Sub-Saharan Africa"                 = "Sub-Saharan Africa",
               "Central & South Asia"               = "Central & South Asia",
               "Europe, Middle East & North Africa" = "Europe, Middle East & North Africa",
               "Latin America & Caribbean"          = "Latin America & Caribbean",
               "East Asia & Pacific"                = "East Asia & Pacific")),
  residence      = list(col = "URBAN_RURA", label = "Residence", tested = FALSE,
    levels = c("Urban" = "U", "Rural" = "R")),
  sex            = list(col = "sex", label = "Infant's gender", tested = TRUE,
    levels = c("Male" = "1", "Female" = "2")),
  bord           = list(col = "bord6", label = "Born order", tested = TRUE,
    levels = c("1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5", "≥6" = "6")),
  twin           = list(col = "twin2", label = "Twin", tested = TRUE,
    levels = c("0" = "0", "1" = "1", "≥2" = "2")),
  mother_age     = list(col = "mother_age_level", label = "Mother's age", tested = TRUE,
    levels = c("<18" = "1", "18~35" = "2", ">35" = "3")),
  birth_interval = list(col = "birth_interval_level",
    label = "Maternal average birth interval", tested = FALSE,
    levels = c("< 24 months" = "1", "24~48 months" = "2", "> 48 months" = "3")),
  maternal_edu   = list(col = "maternal_edu", label = "Mother's highest education", tested = FALSE,
    levels = c("No education" = "0", "Primary" = "1", "Secondary" = "2", "Higher" = "3")),
  water_sour     = list(col = "water_sour", label = "Water infrastructures", tested = FALSE,
    levels = c("Basic (surface water)" = "1", "Intermediate (well)" = "2",
               "Modern (piped/tap)" = "3", "Others" = "4")),
  toilet_type    = list(col = "toilet_type", label = "Sanitation infrastructures", tested = FALSE,
    levels = c("Basic (no facility)" = "1", "Intermediate (latrine)" = "2",
               "Modern (flush toilet)" = "3", "Others" = "4")),
  floor_mat      = list(col = "floor_mat", label = "Floor material", tested = FALSE,
    levels = c("Natural" = "1", "Rudimentary" = "2", "Finished" = "3", "Others" = "4")),
  roof_mat       = list(col = "roof_mat", label = "Roof material", tested = FALSE,
    levels = c("Natural" = "1", "Rudimentary" = "2", "Finished" = "3", "Others" = "4")),
  wall_mat       = list(col = "wall_mat", label = "Wall material", tested = FALSE,
    levels = c("Natural" = "1", "Rudimentary" = "2", "Finished" = "3", "Others" = "4"))
)


# ==============================================================================
# 4. Core functions
# ==============================================================================

#' Frequency and column percentage of each category among cases and controls.
#'
#' @param col        analysis-frame column name.
#' @param levels_vec named vector <display-label> = <data value>.
#' @param dat        analysis frame.
#' @return data.frame with one row per category: category, case_n, control_n,
#'         case_pct, control_pct (percentages use the full case/control total).
tabulate_var <- function(col, levels_vec, dat) {
  x     <- as.character(dat[[col]])
  labs  <- names(levels_vec)
  vals  <- unname(levels_vec)

  case_n <- vapply(vals, function(v) sum(x == v & dat$u5mr == 1, na.rm = TRUE),
                   integer(1))
  ctrl_n <- vapply(vals, function(v) sum(x == v & dat$u5mr == 0, na.rm = TRUE),
                   integer(1))

  data.frame(
    category    = labs,
    case_n      = case_n,
    control_n   = ctrl_n,
    case_pct    = case_n / n_case * 100,
    control_pct = ctrl_n / n_ctrl * 100,
    stringsAsFactors = FALSE
  )
}

#' Overall Wald (Type III) P value from a sibling-matched conditional logistic
#' regression. Returns NA_real_ if the model cannot be estimated (e.g. the
#' covariate is invariant within strata in a small sample).
#'
#' @param col covariate column (entered as as.factor()).
#' @param dat analysis frame.
wald_pvalue <- function(col, dat) {
  f <- as.formula(sprintf("u5mr ~ as.factor(%s) + strata(caseid)", col))
  m <- tryCatch(
    clogit(f, method = "breslow", data = dat, na.action = na.omit),
    error = function(e) NULL
  )
  if (is.null(m)) return(NA_real_)

  a <- tryCatch(
    Anova(m, type = "III", test.statistic = "Wald"),
    error = function(e) NULL
  )
  if (is.null(a)) return(NA_real_)

  pcol <- grep("Pr\\(>", colnames(a))
  if (length(pcol) == 0) return(NA_real_)
  as.numeric(a[1, pcol[1]])
}

#' Format a P value for display ("<0.001" below the threshold, otherwise 3 dp).
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}


# ==============================================================================
# 5. Compute the case/control summary and P value for every variable
# ==============================================================================
rows <- lapply(tab_config, function(cfg) {
  t <- tabulate_var(cfg$col, cfg$levels, data)

  p <- NA_real_
  if (isTRUE(cfg$tested)) p <- wald_pvalue(cfg$col, data)

  data.frame(
    variable = cfg$label,
    category = t$category,
    case_n   = t$case_n,
    case_pct = t$case_pct,
    control_n   = t$control_n,
    control_pct = t$control_pct,
    p_value     = p,
    stringsAsFactors = FALSE
  )
})

result <- do.call(rbind, rows)
rownames(result) <- NULL

# Display P value: "" for untested (household-invariant) variables; for tested
# variables, show the formatted P value once per variable (first category row).
result$p_show <- "-"
for (nm in names(tab_config)) {
  cfg <- tab_config[[nm]]
  idx <- which(result$variable == cfg$label)
  if (isTRUE(cfg$tested)) {
    result$p_show[idx]    <- ""
    result$p_show[idx[1]] <- fmt_p(result$p_value[idx[1]])
  }
}


# ==============================================================================
# 6. Format Table S4 and export
# ==============================================================================
# Formatted cell: "N (xx.xx%)".
result$case_txt  <- sprintf("%s (%.2f%%)", format(result$case_n, big.mark = ","),
                            result$case_pct)
result$ctrl_txt  <- sprintf("%s (%.2f%%)", format(result$control_n, big.mark = ","),
                            result$control_pct)

table_s4 <- result[, c("variable", "category", "case_txt", "ctrl_txt", "p_show")]
names(table_s4) <- c("Variables", "Category", "Cases", "Controls",
                     "P value for difference")

# Raw (numeric) export for downstream use.
raw_out <- result[, c("variable", "category", "case_n", "case_pct",
                      "control_n", "control_pct", "p_value")]
write.csv(table_s4, file.path(out_dir, "TabS4_IndividualSummary.csv"),
          row.names = FALSE)


# ==============================================================================
# 7. Print
# ==============================================================================
cat("\n=== Table S4: individual characteristics of cases and controls ===\n\n")
cat(sprintf("Cases (N = %s)   Controls (N = %s)\n",
            format(n_case, big.mark = ","), format(n_ctrl, big.mark = ",")))
print(table_s4, row.names = FALSE)
cat("\nFootnotes:\n")
cat(" 1 Mother's age at the child's death (for cases) or at the corresponding matched time point (for controls).\n")
cat(" 2 P values were derived from conditional logistic regression; for multi-category variables an overall Wald test (Type III) was conducted.\n")
cat("   Variables invariant within households (region, residence, education, water/sanitation, floor/roof/wall material)\n")
cat("   were not tested because they are balanced by design through within-household matching.\n")


# ==============================================================================
# 8. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
