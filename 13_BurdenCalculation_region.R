# ==============================================================================
# 13_BurdenCalculation_region.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
# (Region-stratified burden calculation -- Table S22)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To recompute the burden of under-5 mortality attributable to the six extreme
#    climatic events using REGION-SPECIFIC effect estimates (from the regional
#    EWE-related child-mortality risk analysis), so that the results correspond to
#    Table S22 ("Global and regional burden ... based on regional EWE-related
#    child-mortality risk estimation").
#
# Methodology Overview:
#    1. Same grid-level attrdl() framework as 10_BurdenCalculation.R (the global
#       run), with the same deterministic per-grid Monte-Carlo seed and checkpoint
#       resume -- but each grid uses the effect estimate of ITS OWN region:
#           coef_<exposure>_<region>.rda / vcov_<exposure>_<region>.rda /
#           cbTEM_<exposure>_<region>.rda
#    2. Each region also uses its own cumulative-lag window (the region-specific
#       "maximum significant lag" from the regional analysis), passed to attrdl()
#       as sublag = c(0, cum_lag).
#    3. Tropical cyclone has no Europe / MENA model, so TC x MENA burden is set to
#       0 (matches the "0 (0, 0)" cell in Table S22).
#    4. The "Global" row is the SUM of the five regional estimates (each region
#       computed with its own regional effect), NOT a separately pooled model.
#
# REQUIREMENTS
#    R        >= 4.0.3
#    dlnm     == 2.4.7
#    dplyr    == 1.1.1

# INPUT
#    frame_<key>_LIMCs_month_10.rda -- gridded monthly ECE exposure
#    coef_<key>_<region>.rda / vcov_<key>_<region>.rda / cbTEM_<key>_<region>.rda
#    mortality_country_merge.csv    -- COUNTRY -> REGION mapping
#    mortality_2000_2019_info.csv   -- country under-5 baseline mortality rate
#    pop_data_info_10.rda           -- grid under-5 population
#
# OUTPUT
#    list_<key>_LIMCs_month_AN_region.rda -- grid-level monthly AN (checkpoint)
#    frame_<key>_LIMCs_year_10_case[_low/_high]_{grid,country,region}.csv
#    TableS22_region_burden.csv           -- the summary table (Table S22)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------
library(dlnm);library(dplyr)

# ------------------------------------------------------------------------------
# 2. Custom attributable-risk function attrdl (identical to the global run)
# ------------------------------------------------------------------------------
attrdl <- function(x, basis, cases, model = NULL, coef = NULL, vcov = NULL,
                   type = "af", dir = "back", tot = TRUE, cen, range = NULL,
                   sim = FALSE, nsim = 1000, sublag = NULL) {

  if (packageVersion("dlnm") < "2.2.0") stop("Update dlnm package to version >= 2.2.0")

  name <- deparse(substitute(basis))
  type <- match.arg(type, c("an", "af"))
  dir  <- match.arg(dir, c("back", "forw"))

  if (missing(cen) && is.null(cen <- attr(basis, "argvar")$cen)) stop("'cen' must be provided")
  if (!is.numeric(cen) && length(cen) > 1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL

  if (!is.null(range)) x[x < range[1] | x > range[2]] <- cen

  lag <- attr(basis, "lag")
  if (NCOL(x) == 1L) {
    at <- if (dir == "back") tsModel:::Lag(x, seq(lag[1], lag[2])) else
      matrix(rep(x, diff(lag) + 1), length(x))
  } else {
    if (dir == "forw") stop("'x' must be a vector when dir='forw'")
    if (ncol(at <- x) != diff(lag) + 1) stop("Dimension of 'x' not compatible with 'basis'")
  }

  if (NROW(cases) != NROW(at)) stop("'x' and 'cases' not consistent")
  if (NCOL(cases) > 1L) {
    if (dir == "back") stop("'cases' must be a vector if dir='back'")
    if (ncol(cases) != diff(lag) + 1) stop("Dimension of 'cases' not compatible")
    den <- sum(rowMeans(cases, na.rm = TRUE), na.rm = TRUE)
    cases <- rowMeans(cases)
  } else {
    den <- sum(cases, na.rm = TRUE)
    if (dir == "forw") cases <- rowMeans(as.matrix(tsModel:::Lag(cases, -seq(lag[1], lag[2]))))
  }

  if (!is.null(model)) {
    model.class <- class(model)
    all_coef <- dlnm:::getcoef(model, model.class)
    all_vcov <- dlnm:::getvcov(model, model.class)
    cond <- paste0("^", name, "([[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}|[0-9]+)?$")
    ind  <- grep(cond, names(all_coef))
    if (length(ind) == 0L) ind <- grep(name, names(all_coef), fixed = TRUE)
    if (length(ind) == 0L) stop(paste("No coefficients matching the basis name", name, "were found in the model."))
    coef <- all_coef[ind]
    vcov <- all_vcov[ind, ind, drop = FALSE]
    model.link <- tryCatch(dlnm:::getlink(model, model.class), error = function(e) "log")
    if (!is.null(model.link) && !(model.link %in% c("log", "logit", "cox")))
      warning("Model link function is not explicitly 'log'.")
  }

  typebasis <- ifelse(length(coef) != ncol(basis), "one", "cb")
  predvar <- if (typebasis == "one") x else seq(NROW(at))
  predlag <- if (typebasis == "one") 0 else dlnm:::seqlag(lag)

  if (!is.null(sublag)) {
    if (length(sublag) != 2L || sublag[1] > sublag[2])
      stop("'sublag' must be a length-2 vector giving the first and last lag, e.g. c(0, 3)")
    if (sublag[1] < lag[1] || sublag[2] > lag[2])
      stop(paste0("'sublag' (", sublag[1], "-", sublag[2], ") exceeds the lag range"))
  }

  if (typebasis == "cb") {
    Xpred <- dlnm:::mkXpred(typebasis, basis, at, predvar, predlag, cen)
    Xpredall <- 0
    for (i in seq_along(predlag)) {
      current_lag <- predlag[i]
      if (is.null(sublag) || (current_lag >= sublag[1] && current_lag <= sublag[2])) {
        ind_i <- seq_along(predvar) + length(predvar) * (i - 1)
        Xpredall <- Xpredall + Xpred[ind_i, , drop = FALSE]
      }
    }
  } else {
    basis_one <- do.call(onebasis, c(list(x = x), attr(basis, "argvar")))
    Xpredall  <- dlnm:::mkXpred(typebasis, basis_one, x, predvar, predlag, cen)
  }

  if (length(coef) == 0L) stop("Extracted coefficient vector 'coef' is empty; check model and cross-basis.")
  if (length(coef) != ncol(Xpredall)) stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if (any(dim(vcov) != c(length(coef), length(coef)))) stop("arguments 'coef' and 'vcov' do not match")
  if (typebasis == "one" && dir == "back") stop("only dir='forw' allowed for reduced estimates")

  af <- 1 - exp(-drop(as.matrix(Xpredall %*% coef)))
  an <- af * cases

  if (tot) {
    isna <- is.na(an)
    af <- sum(an[!isna]) / sum(cases[!isna])
    an <- af * den
  }

  if (sim) {
    k <- length(coef)
    eigen_decomp <- eigen(vcov)
    X <- matrix(rnorm(length(coef) * nsim), nsim)
    coefsim <- coef + eigen_decomp$vectors %*% diag(sqrt(pmax(eigen_decomp$values, 0)), k) %*% t(X)
    if (tot) {
      afsim <- apply(coefsim, 2, function(coefi) {
        ani <- (1 - exp(-drop(Xpredall %*% coefi))) * cases
        sum(ani[!is.na(ani)]) / sum(cases[!is.na(ani)])
      })
      ansim <- afsim * den
      sim_vec <- if (type == "an") ansim else afsim
      est <- if (type == "an") an else af
      ci  <- quantile(sim_vec, probs = c(0.025, 0.975), na.rm = TRUE)
      res <- c(est = est, low = ci[1], high = ci[2])
      names(res) <- c("est", "2.5%", "97.5%")
    } else {
      afsim <- 1 - exp(-(Xpredall %*% coefsim))
      ansim <- afsim * cases
      sim_mat <- if (type == "an") ansim else afsim
      ci_low  <- apply(sim_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
      ci_high <- apply(sim_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
      est     <- if (type == "an") an else af
      res <- data.frame(est = est, low = ci_low, high = ci_high, row.names = NULL)
      colnames(res) <- c("est", "2.5%", "97.5%")
    }
  } else {
    res <- if (type == "an") an else af
  }
  return(res)
}

# ------------------------------------------------------------------------------
# 3. Configuration
# ------------------------------------------------------------------------------
src_dir <- "Repository/Data/Burden_Source"
out_dir <- "Repository/Output/burden_result_regional"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Burden_Source subfolders
exposure_dir <- file.path(src_dir, "Global_exposure")
pop_dir      <- file.path(src_dir, "Pop_data")
coef_dir     <- file.path(src_dir, "coef_vcov")
cbTEM_dir    <- file.path(src_dir, "cbTEM_regional")

# Exposure config: short prefix (used for the region coef/vcov/cbTEM files) mapped
# to the global burden key / var / frame / number-of-years.
exposure_config <- data.frame(
  exp        = c("hw",    "EP",    "flood",  "spei",  "TC",       "wf"),
  key        = c("hw_90_2_days", "tp_99", "flood", "spei_6", "tc_R34", "wf_exposure"),
  var        = c("hw_90_2_days", "tp_99", "flood", "drought", "tc", "wf_exposure"),
  label      = c("Heatwave", "Extreme precipitation", "Flood", "Drought",
                 "Tropical cyclone", "Wildfire"),
  frame_file = c("frame_hw_90_2_days_LIMCs_month_10", "frame_EP_99_cum_LIMCs_month_10",
                 "frame_flood_total_LIMCs_month_10", "frame_dr_6m_LIMCs_month_10",
                 "frame_total_TC_LIMCs_month_10", "frame_wf_exposure_LIMCs_month_10"),
  num_years  = c(20L, 20L, 20L, 20L, 20L, 18L),
  start_year = c(2000L, 2000L, 2000L, 2000L, 2000L, 2002L),
  stringsAsFactors = FALSE
)

# Region config: short code -> data label (match the REGION column produced by
# merging mortality_country_merge.csv, and the label used by the region effect files).
region_config <- data.frame(
  code  = c("SSA",  "CSA",  "MENA", "LAC",  "EAP"),
  label = c("Sub-Saharan Africa", "Central & South Asia",
            "Europe, Middle East & North Africa", "Latin America & Caribbean",
            "East Asia & Pacific"),
  stringsAsFactors = FALSE
)

# Region-specific cumulative-lag window (months)
cumlag_region <- list(
  hw    = c(SSA = 6,  CSA = 8,  MENA = 3,  LAC = 1,  EAP = 7),
  EP    = c(SSA = 10, CSA = 7,  MENA = 8,  LAC = 3,  EAP = 6),
  flood = c(SSA = 6,  CSA = 6,  MENA = 6,  LAC = 7,  EAP = 6),
  spei  = c(SSA = 3,  CSA = 3,  MENA = 3,  LAC = 3,  EAP = 3),
  TC    = c(SSA = 4,  CSA = 1,  LAC = 3,  EAP = 1),
  wf    = c(SSA = 11, CSA = 11, MENA = 1,  LAC = 2,  EAP = 4)
)

# Run control
nsim_burden <- 1000L
run_keys    <- exposure_config$key          # all six exposures
max_grids   <- NULL                          # debug only (NULL = all grids)
mc_seed_base <- 20240201L

# Command-line debug override:
#   Rscript 010_BurdenCalculation_region.R <key> <max_grids>
# <key> limits the run to one exposure; <max_grids> processes only the first N grids.
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && cli_args[1] %in% exposure_config$key) run_keys <- cli_args[1]
if (length(cli_args) >= 2 && grepl("^[0-9]+$", cli_args[2])) max_grids <- as.integer(cli_args[2])

# ------------------------------------------------------------------------------
# 4. Shared data (population, country->region mapping, baseline mortality)
# ------------------------------------------------------------------------------
region_country_merge <- read.csv(file.path(pop_dir, "mortality_country_merge.csv"),
                                 check.names = FALSE)[, c(3, 4)]   # COUNTRY, REGION

load(file.path(pop_dir, "pop_data_info_10.rda"))                   # -> pop_data_info_10

mortality_2000_2019_info <- read.csv(file.path(pop_dir, "mortality_2000_2019_info.csv"),
                                     check.names = FALSE)

mortality_grid <- merge(
  pop_data_info_10[, 1:4],
  mortality_2000_2019_info[, c(2, 4:23)],
  by = "COUNTRY", all.x = TRUE
)
mortality_grid <- mortality_grid[order(mortality_grid$pointid), ]

mortality_case_grid <- mortality_grid
for (k in 1:20) {
  mortality_case_grid[, 4 + k] <- pop_data_info_10[, 4 + k] * mortality_case_grid[, 4 + k]
}

frame_pop_LIMCs_year_10_s             <- pop_data_info_10[, 5:24]
frame_mortality_case_LIMCs_month_10_s <- mortality_case_grid[, 5:24]

if (!is.null(max_grids)) {
  frame_pop_LIMCs_year_10_s             <- frame_pop_LIMCs_year_10_s[seq_len(max_grids), , drop = FALSE]
  frame_mortality_case_LIMCs_month_10_s <- frame_mortality_case_LIMCs_month_10_s[seq_len(max_grids), , drop = FALSE]
  pop_data_info_10                      <- pop_data_info_10[seq_len(max_grids), , drop = FALSE]
}

# ------------------------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------------------------
monthly_cols_of <- function(frame) grep("^X?20[0-9]{2}[._-][0-9]{2}", colnames(frame))

parse_monthly_dates <- function(cn) {
  s <- sub("^X", "", cn)
  d1 <- as.Date(s, format = "%Y.%m.%d"); d2 <- as.Date(s, format = "%Y-%m-%d")
  d3 <- as.Date(s, format = "%Y.%m");       d4 <- as.Date(s, format = "%Y_%m")
  out <- d1
  out[is.na(out)] <- d2[is.na(out)]
  out[is.na(out)] <- d3[is.na(out)]
  out[is.na(out)] <- d4[is.na(out)]
  out
}

build_cases <- function(yearly_cases, num_months) {
  n_core <- length(yearly_cases) * 12
  if (num_months >= n_core) {
    c(rep(yearly_cases, each = 12), rep(yearly_cases[length(yearly_cases)], num_months - n_core))
  } else {
    rep(yearly_cases, each = 12)[1:num_months]
  }
}

# ------------------------------------------------------------------------------
# 6. Per-exposure region-stratified computation
# ------------------------------------------------------------------------------
label_to_code <- setNames(region_config$code, region_config$label)

for (row_i in which(exposure_config$key %in% run_keys)) {
  cfg <- exposure_config[row_i, ]
  exp       <- cfg$exp
  key       <- cfg$key
  var       <- cfg$var
  num_years <- cfg$num_years
  start_year<- cfg$start_year
  col0 <- start_year - 1999L

  cat("\n================ Computing (region):", key, " ================\n")

  # ---- 6.1 Load exposure frame, merge region, align with pointid ----
  load(file.path(exposure_dir, paste0(cfg$frame_file, ".rda")))
  frame <- get(cfg$frame_file)
  frame$COUNTRY[frame$COUNTRY == "Taiwan"] <- "China"
  frame <- merge(frame, region_country_merge, by = "COUNTRY", all.x = TRUE)
  frame <- frame[order(frame$OBJECTID), ]
  if (!is.null(max_grids)) frame <- frame[seq_len(max_grids), , drop = FALSE]

  monthly_cols <- monthly_cols_of(frame)
  num_months   <- length(monthly_cols)
  date_seq     <- parse_monthly_dates(colnames(frame)[monthly_cols])
  attr_cols    <- seq_len(min(monthly_cols) - 1)
  region_col   <- which(colnames(frame) == "REGION")
  grid_info    <- frame[, c(attr_cols, region_col)]
  grid_region_code <- as.character(label_to_code[as.character(frame$REGION)])

  cat("  grids:", nrow(frame), "| monthly cols:", num_months,
      "| years:", num_years, "\n")

  # ---- 6.2 Load region-specific coef / vcov / cbTEM (once per region) ----
  region_effects <- list()
  for (j in seq_len(nrow(region_config))) {
    code  <- region_config$code[j]
    label <- region_config$label[j]
    cl <- cumlag_region[[exp]][code]
    if (exp == "TC" && code == "MENA") { region_effects[[code]] <- NULL; next }
    if (is.na(cl)) { region_effects[[code]] <- NULL; next }

    e <- new.env(); load(file.path(coef_dir, paste0("coef_", exp, "_", code, ".rda")), envir = e); coef_r <- e[[ls(e)[1]]]
    e <- new.env(); load(file.path(coef_dir, paste0("vcov_", exp, "_", code, ".rda")), envir = e); vcov_r <- e[[ls(e)[1]]]
    e <- new.env(); load(file.path(cbTEM_dir, paste0("cbTEM_", exp, "_", code, ".rda")), envir = e); cbTEM_r <- e[[ls(e)[1]]]
    if (length(coef_r) != ncol(cbTEM_r))
      warning(sprintf("  [%s|%s] coef length %d != cbTEM cols %d", exp, code, length(coef_r), ncol(cbTEM_r)))

    region_effects[[code]] <- list(coef = coef_r, vcov = vcov_r, cbTEM = cbTEM_r,
                                   sublag = c(0, cl))
  }

  # ---- 6.3 Grid-level attributable deaths (region-specific effect) ----
  output_file <- file.path(out_dir, paste0("list_", key, "_LIMCs_month_AN_region.rda"))
  total_grids <- nrow(frame)
  if (file.exists(output_file)) {
    load(output_file)          # -> list_an
    cat("  checkpoint found, resuming ...\n")
  } else {
    # Pre-extract the per-grid numeric data (compact), then loop over all grids
    # sequentially in a single thread.
    exposure_mat  <- as.matrix(frame[, monthly_cols])                       # total_grids x num_months
    cases_year_mat<- as.matrix(frame_mortality_case_LIMCs_month_10_s[, col0:(col0 + num_years - 1)])
    idx           <- seq_len(total_grids)

    compute_grids <- function(cd) {
      out <- vector("list", length(cd$ids))
      for (k in seq_along(cd$ids)) {
        i <- cd$ids[k]
        exposure     <- as.numeric(cd$exposure[k, ])
        yearly_cases <- cd$cases_year[k, ] / 12
        cases_vector <- build_cases(yearly_cases, num_months)
        code_i <- cd$region_code[k]
        eff <- region_effects[[code_i]]
        set.seed(mc_seed_base + row_i * 1000000L + i)
        if (is.null(eff)) {
          # No regional model for this grid (e.g. TC in MENA): attributable deaths = 0.
          out[[k]] <- data.frame(date = date_seq, exposure = exposure, case = cases_vector,
                                 AN_v = 0, AN_l = 0, AN_h = 0)
        } else {
          an <- attrdl(x = exposure, basis = eff$cbTEM, cases = cases_vector,
                       coef = eff$coef, vcov = eff$vcov,
                       type = "an", dir = "forw", tot = FALSE, cen = 0,
                       sim = TRUE, nsim = nsim_burden, sublag = eff$sublag)
          out[[k]] <- data.frame(date = date_seq, exposure = exposure, case = cases_vector,
                                 AN_v = an$est, AN_l = an[, "2.5%"], AN_h = an[, "97.5%"])
        }
      }
      out
    }

    cat("  computing", total_grids, "grids (single-threaded) ...\n")
    list_an <- compute_grids(list(ids = idx,
                                  exposure = exposure_mat,
                                  cases_year = cases_year_mat,
                                  region_code = grid_region_code))
    rm(exposure_mat, cases_year_mat, compute_grids); gc()

    save(list_an, file = output_file)
  }
  cat("  grid-level attributable deaths computed.\n")

  # ---- 6.4 Monthly -> annual aggregation ----
  mat_AN_v_month <- t(sapply(list_an, function(df) df$AN_v[1:num_months]))
  mat_AN_l_month <- t(sapply(list_an, function(df) df$AN_l[1:num_months]))
  mat_AN_h_month <- t(sapply(list_an, function(df) df$AN_h[1:num_months]))

  mat_AN_v_year <- matrix(NA, nrow = total_grids, ncol = num_years)
  mat_AN_l_year <- matrix(NA, nrow = total_grids, ncol = num_years)
  mat_AN_h_year <- matrix(NA, nrow = total_grids, ncol = num_years)
  for (y in 1:num_years) {
    cols_idx <- ((y - 1) * 12 + 1):(y * 12)
    mat_AN_v_year[, y] <- rowSums(mat_AN_v_month[, cols_idx], na.rm = TRUE)
    mat_AN_l_year[, y] <- rowSums(mat_AN_l_month[, cols_idx], na.rm = TRUE)
    mat_AN_h_year[, y] <- rowSums(mat_AN_h_month[, cols_idx], na.rm = TRUE)
  }
  colnames(mat_AN_v_year) <- paste0("year_", start_year:(start_year + num_years - 1))
  colnames(mat_AN_l_year) <- paste0("year_", start_year:(start_year + num_years - 1))
  colnames(mat_AN_h_year) <- paste0("year_", start_year:(start_year + num_years - 1))

  frame_year_case      <- data.frame(grid_info, mat_AN_v_year)
  frame_year_case_low  <- data.frame(grid_info, mat_AN_l_year)
  frame_year_case_high <- data.frame(grid_info, mat_AN_h_year)

  # ---- 6.5 Spatial aggregation (grid / country / region) ----
  calc_spatial_burden <- function(df_year_case) {
    nattr <- ncol(grid_info)
    val_cols <- (nattr + 1):ncol(df_year_case)
    grid_df <- data.frame(df_year_case[, seq_len(nattr)], check.names = FALSE)
    grid_df$ewe_case <- rowSums(df_year_case[, val_cols], na.rm = TRUE)
    grid_df$case     <- rowSums(frame_mortality_case_LIMCs_month_10_s[, col0:(col0 + num_years - 1)], na.rm = TRUE)
    grid_df$pop_sum  <- rowSums(frame_pop_LIMCs_year_10_s[, col0:(col0 + num_years - 1)], na.rm = TRUE)
    grid_df$ewe_case_pop <- (grid_df$ewe_case / grid_df$pop_sum) * 100000
    grid_df$ewe_case_pop[is.na(grid_df$ewe_case_pop)] <- 0

    country_df <- grid_df %>%
      group_by(COUNTRY) %>%
      summarise(ewe_case = sum(ewe_case, na.rm = TRUE),
                case     = sum(case, na.rm = TRUE),
                pop_sum  = sum(pop_sum, na.rm = TRUE)) %>%
      rbind(data.frame(COUNTRY = "Global",
                       ewe_case = sum(.$ewe_case, na.rm = TRUE),
                       case     = sum(.$case, na.rm = TRUE),
                       pop_sum  = sum(.$pop_sum, na.rm = TRUE)))

    region_df <- grid_df %>%
      group_by(REGION) %>%
      summarise(ewe_case = sum(ewe_case, na.rm = TRUE),
                case     = sum(case, na.rm = TRUE),
                pop_sum  = sum(pop_sum, na.rm = TRUE)) %>%
      rbind(data.frame(REGION = "Global",
                       ewe_case = sum(.$ewe_case, na.rm = TRUE),
                       case     = sum(.$case, na.rm = TRUE),
                       pop_sum  = sum(.$pop_sum, na.rm = TRUE)))

    names(grid_df)[names(grid_df) == "ewe_case"]     <- paste0(var, "_case")
    names(grid_df)[names(grid_df) == "ewe_case_pop"] <- paste0(var, "_case_pop")
    names(country_df)[names(country_df) == "ewe_case"] <- paste0(var, "_case")
    names(region_df)[names(region_df) == "ewe_case"]   <- paste0(var, "_case")

    country_df[[paste0(var, "_burden")]] <- (country_df[[paste0(var, "_case")]] / country_df$pop_sum) * 100000
    country_df[[paste0(var, "_PAF")]]    <- (country_df[[paste0(var, "_case")]] / country_df$case) * 100
    region_df[[paste0(var, "_burden")]]  <- (region_df[[paste0(var, "_case")]] / region_df$pop_sum) * 100000
    region_df[[paste0(var, "_PAF")]]     <- (region_df[[paste0(var, "_case")]] / region_df$case) * 100

    list(grid = grid_df, country = country_df, region = region_df)
  }

  res_est  <- calc_spatial_burden(frame_year_case)
  res_low  <- calc_spatial_burden(frame_year_case_low)
  res_high <- calc_spatial_burden(frame_year_case_high)

  filter_sk <- function(df) if ("COUNTRY" %in% colnames(df)) df[df$COUNTRY != "South Korea", ] else df
  res_est$grid  <- filter_sk(res_est$grid)
  res_low$grid  <- filter_sk(res_low$grid)
  res_high$grid <- filter_sk(res_high$grid)

  # ---- 6.6 Export per-exposure files ----
  p <- function(x) file.path(out_dir, x)
  write.csv(res_est$grid,      p(paste0("frame_", key, "_LIMCs_year_10_case_grid.csv")), row.names = FALSE)
  write.csv(res_est$country,   p(paste0("frame_", key, "_LIMCs_year_10_case_country.csv")), row.names = FALSE)
  write.csv(res_est$region,    p(paste0("frame_", key, "_LIMCs_year_10_case_region.csv")), row.names = FALSE)
  write.csv(res_low$grid,      p(paste0("frame_", key, "_LIMCs_year_10_case_low_grid.csv")), row.names = FALSE)
  write.csv(res_low$country,   p(paste0("frame_", key, "_LIMCs_year_10_case_low_country.csv")), row.names = FALSE)
  write.csv(res_low$region,    p(paste0("frame_", key, "_LIMCs_year_10_case_low_region.csv")), row.names = FALSE)
  write.csv(res_high$grid,     p(paste0("frame_", key, "_LIMCs_year_10_case_high_grid.csv")), row.names = FALSE)
  write.csv(res_high$country,  p(paste0("frame_", key, "_LIMCs_year_10_case_high_country.csv")), row.names = FALSE)
  write.csv(res_high$region,   p(paste0("frame_", key, "_LIMCs_year_10_case_high_region.csv")), row.names = FALSE)

  cat("  done", key, "-> files exported\n")
  rm(list_an, frame, region_effects, cbTEM_r, coef_r, vcov_r, e); gc()
}

# ------------------------------------------------------------------------------
