# ==============================================================================
# 10_BurdenCalculation.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
# (Global burden calculation -- grid-level attributable under-5 deaths, the data
#  behind Figure 4 / 6 and Tables S19-S22)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To estimate the number of under-5 deaths attributable to each of six extreme
#    climatic events (ECE) -- heatwave, extreme precipitation, flood, drought,
#    tropical cyclone and wildfire -- across the 126 low- and middle-income
#    countries (LMICs) for 2000-2019, at the grid, country, region and global
#    level, together with 95% confidence intervals from Monte-Carlo simulation.
#
# Methodology Overview:
#    1. Cross-basis coefficients (coef), covariance matrix (vcov) and the
#       cross-basis template (cbTEM) are loaded and passed to attrdl().
#    2. For every grid the forward attributable number (AN) is computed with
#       Monte-Carlo simulation (1000 draws) under a deterministic per-grid seed
#       (mc_seed_base + exposure index * 1e6 + grid index), so the point estimate
#       and 95% CI are reproducible.
#    3. The six exposures share one population / baseline-mortality dataset and
#       are processed in a single loop with independent per-exposure checkpoints
#       (saved every 10,000 grids), allowing a clean resume after interruption.
#    4. Monthly columns and the grid attribute / REGION columns are detected.
#    5. Grid-level monthly AN -> annual -> spatial (grid / country / region) ->
#       time trend (case / PAF / burden, point estimate and 95% CI).
#
# REQUIREMENTS
#    R       >= 4.0.3
#    dlnm    == 2.4.7
#    dplyr   == 1.1.1
#
# INPUT
#    frame_<key>_LIMCs_month_10.rda -- gridded monthly ECE exposure (703714 grids)
#    coef_<key>.rda / vcov_<key>.rda -- pre-extracted cross-basis coefficients / covariance
#    cbTEM_<key>.rda                -- cross-basis template (argvar / arglag / lag)
#    mortality_country_merge.csv    -- COUNTRY -> REGION mapping
#    mortality_2000_2019_info.csv   -- country under-5 baseline mortality rate, 2000-2019
#    pop_data_info_10.rda           -- grid under-5 population, 2000-2019
#
# OUTPUT
#    list_<key>_LIMCs_month_AN.rda  -- grid-level monthly AN (checkpoint file)
#    frame_<key>_LIMCs_year_10_case[_low/_high]_{grid,country,region}.csv
#    frame_<key>_LIMCs_year_10_{case,PAF,burden}[_low/_high]_trend.csv
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------------------------
library(dlnm);library(dplyr)

# ------------------------------------------------------------------------------
# 2. Custom attributable-risk function attrdl
# ------------------------------------------------------------------------------
attrdl <- function(x, basis, cases, model = NULL, coef = NULL, vcov = NULL,
                   type = "af", dir = "back", tot = TRUE, cen, range = NULL,
                   sim = FALSE, nsim = 1000, sublag = NULL) {

  # Check dlnm package version
  if (packageVersion("dlnm") < "2.2.0") {
    stop("Update dlnm package to version >= 2.2.0")
  }

  # Extract name and validate arguments
  name <- deparse(substitute(basis))
  type <- match.arg(type, c("an", "af"))
  dir  <- match.arg(dir, c("back", "forw"))

  # Define the centering value
  if (missing(cen) && is.null(cen <- attr(basis, "argvar")$cen)) {
    stop("'cen' must be provided")
  }
  if (!is.numeric(cen) && length(cen) > 1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL

  # Truncate exposure to the specified range (values outside are set to cen, i.e. zero risk)
  if (!is.null(range)) x[x < range[1] | x > range[2]] <- cen

  # Build the lag matrix
  lag <- attr(basis, "lag")
  if (NCOL(x) == 1L) {
    at <- if (dir == "back") tsModel:::Lag(x, seq(lag[1], lag[2])) else
      matrix(rep(x, diff(lag) + 1), length(x))
  } else {
    if (dir == "forw") stop("'x' must be a vector when dir='forw'")
    if (ncol(at <- x) != diff(lag) + 1)
      stop("Dimension of 'x' not compatible with 'basis'")
  }

  # Validate cases dimensions and compute the denominator / forward moving average
  if (NROW(cases) != NROW(at)) stop("'x' and 'cases' not consistent")
  if (NCOL(cases) > 1L) {
    if (dir == "back") stop("'cases' must be a vector if dir='back'")
    if (ncol(cases) != diff(lag) + 1) stop("Dimension of 'cases' not compatible")
    den <- sum(rowMeans(cases, na.rm = TRUE), na.rm = TRUE)
    cases <- rowMeans(cases)
  } else {
    den <- sum(cases, na.rm = TRUE)
    if (dir == "forw")
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases, -seq(lag[1], lag[2]))))
  }

  # Extract coefficients and covariance from the model when provided
  if (!is.null(model)) {
    model.class <- class(model)
    all_coef <- dlnm:::getcoef(model, model.class)
    all_vcov <- dlnm:::getvcov(model, model.class)

    cond <- paste0("^", name, "([[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}|[0-9]+)?$")
    ind  <- grep(cond, names(all_coef))

    if (length(ind) == 0L) {
      ind <- grep(name, names(all_coef), fixed = TRUE)
    }

    if (length(ind) == 0L) {
      stop(paste("No coefficients matching the basis name", name, "were found in the model."))
    }

    coef <- all_coef[ind]
    vcov <- all_vcov[ind, ind, drop = FALSE]

    model.link <- tryCatch(dlnm:::getlink(model, model.class), error = function(e) "log")
    if (!is.null(model.link) && !(model.link %in% c("log", "logit", "cox"))) {
      warning("Model link function is not explicitly 'log', but clogit models operate on the log-RR/odds scale.")
    }
  }

  typebasis <- ifelse(length(coef) != ncol(basis), "one", "cb")

  # Prepare basis-transform arguments
  predvar <- if (typebasis == "one") x else seq(NROW(at))
  predlag <- if (typebasis == "one") 0 else dlnm:::seqlag(lag)

  # Validate the sublag range
  if (!is.null(sublag)) {
    if (length(sublag) != 2L || sublag[1] > sublag[2]) {
      stop("'sublag' must be a length-2 vector giving the first and last lag, e.g. c(0, 3)")
    }
    if (sublag[1] < lag[1] || sublag[2] > lag[2]) {
      stop(paste0("'sublag' (", sublag[1], "-", sublag[2],
                  ") exceeds the lag range of the cross-basis (", lag[1], "-", lag[2], ")"))
    }
  }

  # Build the centred basis matrix
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

  # Dimension checks
  if (length(coef) == 0L) stop("Extracted coefficient vector 'coef' is empty; check model and cross-basis.")
  if (length(coef) != ncol(Xpredall)) stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if (any(dim(vcov) != c(length(coef), length(coef)))) stop("arguments 'coef' and 'vcov' do not match")
  if (typebasis == "one" && dir == "back") stop("only dir='forw' allowed for reduced estimates")

  # Compute AF and AN
  af <- 1 - exp(-drop(as.matrix(Xpredall %*% coef)))
  an <- af * cases

  # Summary statistics
  if (tot) {
    isna <- is.na(an)
    af <- sum(an[!isna]) / sum(cases[!isna])
    an <- af * den
  }

  # Monte-Carlo simulation for confidence intervals
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

      res <- data.frame(
        est  = est,
        low  = ci_low,
        high = ci_high,
        row.names = NULL
      )
      colnames(res) <- c("est", "2.5%", "97.5%")
    }
  } else {
    res <- if (type == "an") an else af
  }

  return(res)
}

# ------------------------------------------------------------------------------
# 3. Configuration and shared data
# ------------------------------------------------------------------------------

src_dir <- "Repository/Data/Burden_Source"
out_dir <- "Repository/Output/burden_result"

# Burden_Source subfolders
exposure_dir <- file.path(src_dir, "Global_exposure")
pop_dir      <- file.path(src_dir, "Pop_data")
coef_dir     <- file.path(src_dir, "coef_vcov")
cbTEM_dir    <- file.path(src_dir, "cbTEM_global")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Single configuration source for the six exposures.
#   key        : output-file prefix
#   var        : output-column prefix
#   frame_file : exposure time-series frame in burden_source
#   coef_file  : cross-basis coefficients (object `coef`)
#   vcov_file  : cross-basis covariance (object `vcov`)
#   cb_file    : cross-basis template (object `cbTEM`; rgvar/arglag/lag attributes used)
#   sort_key   : column used to align row order with pointid
#   num_years  : number of years aggregated (20 years 2000-2019; wildfire 18 years 2002-2019)
#   start_year : first reportable year (wildfire exposure starts 2002; others 2000)
ece_burden <- data.frame(
  key        = c("hw_90_2_days", "tp_99", "flood", "spei_6", "tc_R34", "wf_exposure"),
  var        = c("hw_90_2_days", "tp_99", "flood", "drought", "tc", "wf_exposure"),
  frame_file = c("frame_hw_90_2_days_LIMCs_month_10", "frame_EP_99_cum_LIMCs_month_10",
                 "frame_flood_total_LIMCs_month_10", "frame_dr_6m_LIMCs_month_10",
                 "frame_total_TC_LIMCs_month_10", "frame_wf_exposure_LIMCs_month_10"),
  coef_file  = c("coef_hw_90_2_days", "coef_EP_99_cum", "coef_flood_total",
                 "coef_spei_6", "coef_TC_total", "coef_wf_exposure"),
  vcov_file  = c("vcov_hw_90_2_days", "vcov_EP_99_cum", "vcov_flood_total",
                 "vcov_spei_6", "vcov_TC_total", "vcov_wf_exposure"),
  cb_file    = c("cbTEM_hw_90_2_days", "cbTEM_EP_99_cum", "cbTEM_flood_total",
                 "cbTEM_spei_6", "cbTEM_total_TC", "cbTEM_wf_exposure"),
  sort_key   = c("OBJECTID", "OBJECTID", "OBJECTID", "OBJECTID", "OBJECTID", "OBJECTID"),
  num_years  = c(20L, 20L, 20L, 20L, 20L, 18L),
  start_year = c(2000L, 2000L, 2000L, 2000L, 2000L, 2002L),
  stringsAsFactors = FALSE
)

# Lag window over which the attributable number is accumulated (0-based).
sublag_map <- list(
  hw_90_2_days = c(0, 8),
  tp_99        = c(0, 6),
  flood        = c(0, 6),
  spei_6       = c(0, 7),
  tc_R34       = c(0, 0),
  wf_exposure  = c(0, 11)
)

# Run control:
#   nsim_burden : number of Monte-Carlo draws for CI
#   run_keys    : subset of exposures to run; full set =
#                 c("hw_90_2_days","tp_99","flood","spei_6","tc_R34","wf_exposure")
#   max_grids   : debug only -- process only the first N grids
nsim_burden <- 1000L
run_keys    <- c("hw_90_2_days", "tp_99", "flood", "spei_6", "tc_R34", "wf_exposure")
max_grids   <- NULL

# Fixed random-seed base: the Monte-Carlo CI (sim = TRUE) is the only source of
# randomness (rnorm). Each grid uses a deterministic seed of the form
# "base + exposure index + grid index", so the point estimate and CI are exactly reproducible.
mc_seed_base <- 20240101L

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) > 0) run_keys <- cli_args

# Six-row trend summary table (fixed row order: heatwave / extreme precipitation /
# flood / drought / tropical cyclone / wildfire)
ewe_names <- c("Heatwave", "Extreme precipitation", "Flood", "Drought",
               "Tropical cyclone", "Wildfire")
create_trend_frame <- function() {
  df <- data.frame(matrix(NA, 6, 21))
  colnames(df) <- c("EWE", paste0("year_", 2000:2019))
  df[, 1] <- ewe_names
  return(df)
}
trend_case  <- create_trend_frame()
trend_paf   <- create_trend_frame()
trend_burden<- create_trend_frame()
trend_case_low   <- create_trend_frame()
trend_paf_low    <- create_trend_frame()
trend_burden_low <- create_trend_frame()
trend_case_high   <- create_trend_frame()
trend_paf_high    <- create_trend_frame()
trend_burden_high <- create_trend_frame()

# ---- Shared data (population, country->region mapping, baseline mortality) ----
region_country_merge <- read.csv(file.path(pop_dir, "mortality_country_merge.csv"),
                                 check.names = FALSE)[, c(3, 4)]   # COUNTRY, REGION

load(file.path(pop_dir, "pop_data_info_10.rda"))                   # -> pop_data_info_10

mortality_2000_2019_info <- read.csv(file.path(pop_dir, "mortality_2000_2019_info.csv"),
                                     check.names = FALSE)

# Grid-level baseline deaths: population x mortality rate
mortality_grid <- merge(
  pop_data_info_10[, 1:4],                        # pointid, COUNTRY, lon, lat
  mortality_2000_2019_info[, c(2, 4:23)],         # COUNTRY, mortality rate 2000-2019
  by = "COUNTRY", all.x = TRUE
)
mortality_grid <- mortality_grid[order(mortality_grid$pointid), ]

mortality_case_grid <- mortality_grid
for (k in 1:20) {
  mortality_case_grid[, 4 + k] <-
    pop_data_info_10[, 4 + k] * mortality_case_grid[, 4 + k]
}

# Shared analysis matrices (all exposures); column 1 of each corresponds to year 2000
frame_pop_LIMCs_year_10_s             <- pop_data_info_10[, 5:24]      # 2000-2019 population
frame_mortality_case_LIMCs_month_10_s <- mortality_case_grid[, 5:24]   # 2000-2019 deaths

# Debug: keep only the first max_grids grids (aligned with the exposure frames)
if (!is.null(max_grids)) {
  frame_pop_LIMCs_year_10_s             <- frame_pop_LIMCs_year_10_s[seq_len(max_grids), , drop = FALSE]
  frame_mortality_case_LIMCs_month_10_s <- frame_mortality_case_LIMCs_month_10_s[seq_len(max_grids), , drop = FALSE]
  pop_data_info_10                      <- pop_data_info_10[seq_len(max_grids), , drop = FALSE]
}

# ---- Helper functions ----
# Detect monthly columns
monthly_cols_of <- function(frame) {
  grep("^X?20[0-9]{2}[._-][0-9]{2}", colnames(frame))
}

# Parse monthly column names to Date
parse_monthly_dates <- function(cn) {
  s <- sub("^X", "", cn)
  d1 <- as.Date(s, format = "%Y.%m.%d")
  d2 <- as.Date(s, format = "%Y-%m-%d")
  d3 <- as.Date(s, format = "%Y.%m")
  d4 <- as.Date(s, format = "%Y_%m")
  out <- d1
  out[is.na(out)] <- d2[is.na(out)]
  out[is.na(out)] <- d3[is.na(out)]
  out[is.na(out)] <- d4[is.na(out)]
  out
}

# Baseline deaths vector (length matching the exposure months).
# yearly_cases holds the annual deaths for the reportable years (wildfire has only
# 18 years, 2002-2019); months beyond the reportable period are filled with the
# last year's value to serve as the forward-lag buffer.
build_cases <- function(yearly_cases, num_months) {
  n_core <- length(yearly_cases) * 12
  if (num_months >= n_core) {
    c(rep(yearly_cases, each = 12),
      rep(yearly_cases[length(yearly_cases)], num_months - n_core))
  } else {
    rep(yearly_cases, each = 12)[1:num_months]
  }
}

# ------------------------------------------------------------------------------
# 4. Per-exposure computation
# ------------------------------------------------------------------------------
for (row_i in which(ece_burden$key %in% run_keys)) {
  cfg    <- ece_burden[row_i, ]
  key    <- cfg$key
  var    <- cfg$var
  sublag <- sublag_map[[key]]
  num_years  <- cfg$num_years
  start_year <- cfg$start_year
  col0 <- start_year - 1999L   # mortality/pop matrices: column 1 = 2000, so start_year sits at start_year - 1999

  cat("\n================ Computing:", key, " ================\n")

  # ---- 4.1 Load exposure frame, merge region, align row order with pointid ----
  load(file.path(exposure_dir, paste0(cfg$frame_file, ".rda")))
  frame <- get(cfg$frame_file)
  frame <- merge(frame, region_country_merge, by = "COUNTRY", all.x = TRUE)
  frame <- frame[order(frame[[cfg$sort_key]]), ]
  if (!is.null(max_grids)) frame <- frame[seq_len(max_grids), , drop = FALSE]

  monthly_cols <- monthly_cols_of(frame)
  num_months   <- length(monthly_cols)
  date_seq     <- parse_monthly_dates(colnames(frame)[monthly_cols])
  attr_cols    <- seq_len(min(monthly_cols) - 1)
  region_col   <- which(colnames(frame) == "REGION")
  grid_info    <- frame[, c(attr_cols, region_col)]

  cat("  grids:", nrow(frame), "| monthly cols:", length(monthly_cols),
      "| years:", num_years, "\n")

  # ---- 4.2 Load pre-extracted coefficients, covariance and cross-basis template ----
  load(file.path(coef_dir, paste0(cfg$coef_file, ".rda")))   # -> coef
  load(file.path(coef_dir, paste0(cfg$vcov_file, ".rda")))   # -> vcov
  load(file.path(cbTEM_dir, paste0(cfg$cb_file, ".rda")))    # -> cbTEM

  cb_coef <- coef
  cb_vcov <- vcov

  # ---- 4.3 Grid-level attributable deaths ----
  output_file <- file.path(out_dir, paste0("list_", key, "_LIMCs_month_AN.rda"))
  total_grids <- nrow(frame)
  grid_start <- 1L
  grid_end   <- total_grids
  if (file.exists(output_file)) {
    load(output_file)          # loads object list_an
    cat("  checkpoint file found, resuming ...\n")
  } else {
    list_an <- vector("list", total_grids)
  }

  save_batch_size <- 10000
  for (i in grid_start:grid_end) {
    if (!is.null(list_an[[i]])) next
    if (i %% 100 == 0 || i == grid_end) {
      cat(sprintf("    %d / %d (%.1f%%)\n", i, total_grids, 100 * i / total_grids))
    }

    exposure <- as.numeric(frame[i, monthly_cols])
    yearly_cases <- as.numeric(frame_mortality_case_LIMCs_month_10_s[i, col0:(col0 + num_years - 1)]) / 12
    cases_vector <- build_cases(yearly_cases, num_months)

    # Deterministic per-grid seed: the Monte-Carlo CI is reproducible, including after resume.
    set.seed(mc_seed_base + row_i * 1000000L + i)

    an <- attrdl(x = exposure, basis = cbTEM, cases = cases_vector,
                 coef = cb_coef, vcov = cb_vcov,
                 type = "an", dir = "forw", tot = FALSE, cen = 0,
                 sim = TRUE, nsim = nsim_burden, sublag = sublag)

    list_an[[i]] <- data.frame(
      date     = date_seq,
      exposure = exposure,
      case     = cases_vector,
      AN_v     = an$est,
      AN_l     = an[, "2.5%"],
      AN_h     = an[, "97.5%"]
    )

    if (i %% save_batch_size == 0 || i == grid_end) {
      save(list_an, file = output_file)
    }
  }
  cat("  grid-level attributable deaths computed.\n")

  # ---- 4.4 Monthly -> annual aggregation (point / lower / upper) ----
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

  # ---- 4.5 Spatial aggregation (grid / country / region) ----
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

    # Rename the internal fixed columns to <var>_case and compute rates/PAF
    names(grid_df)[names(grid_df) == "ewe_case"]     <- paste0(var, "_case")
    names(grid_df)[names(grid_df) == "ewe_case_pop"] <- paste0(var, "_case_pop")
    names(country_df)[names(country_df) == "ewe_case"] <- paste0(var, "_case")
    names(region_df)[names(region_df) == "ewe_case"]   <- paste0(var, "_case")

    country_df[[paste0(var, "_burden")]] <-
      (country_df[[paste0(var, "_case")]] / country_df$pop_sum) * 100000
    country_df[[paste0(var, "_PAF")]] <-
      (country_df[[paste0(var, "_case")]] / country_df$case) * 100

    region_df[[paste0(var, "_burden")]] <-
      (region_df[[paste0(var, "_case")]] / region_df$pop_sum) * 100000
    region_df[[paste0(var, "_PAF")]] <-
      (region_df[[paste0(var, "_case")]] / region_df$case) * 100

    list(grid = grid_df, country = country_df, region = region_df)
  }

  res_est  <- calc_spatial_burden(frame_year_case)
  res_low  <- calc_spatial_burden(frame_year_case_low)
  res_high <- calc_spatial_burden(frame_year_case_high)

  filter_sk <- function(df) {
    if ("COUNTRY" %in% colnames(df)) df[df$COUNTRY != "South Korea", ] else df
  }
  res_est$grid  <- filter_sk(res_est$grid)
  res_low$grid  <- filter_sk(res_low$grid)
  res_high$grid <- filter_sk(res_high$grid)

  # ---- 4.6 Time trend  ----
  year_val_cols <- (ncol(grid_info) + 1):(ncol(grid_info) + num_years)
  denom_case    <- colSums(frame_mortality_case_LIMCs_month_10_s[, col0:(col0 + num_years - 1)], na.rm = TRUE)
  denom_pop     <- colSums(frame_pop_LIMCs_year_10_s[, col0:(col0 + num_years - 1)], na.rm = TRUE)
  trend_case_cols <- (col0 + 1):(col0 + 1 + num_years - 1)

  trend_case[row_i, trend_case_cols]       <- colSums(frame_year_case[, year_val_cols], na.rm = TRUE)
  trend_paf[row_i, trend_case_cols]        <- colSums(frame_year_case[, year_val_cols], na.rm = TRUE) / denom_case * 100
  trend_burden[row_i, trend_case_cols]     <- colSums(frame_year_case[, year_val_cols], na.rm = TRUE) / denom_pop * 100000

  trend_case_low[row_i, trend_case_cols]   <- colSums(frame_year_case_low[, year_val_cols], na.rm = TRUE)
  trend_paf_low[row_i, trend_case_cols]    <- colSums(frame_year_case_low[, year_val_cols], na.rm = TRUE) / denom_case * 100
  trend_burden_low[row_i, trend_case_cols] <- colSums(frame_year_case_low[, year_val_cols], na.rm = TRUE) / denom_pop * 100000

  trend_case_high[row_i, trend_case_cols]   <- colSums(frame_year_case_high[, year_val_cols], na.rm = TRUE)
  trend_paf_high[row_i, trend_case_cols]    <- colSums(frame_year_case_high[, year_val_cols], na.rm = TRUE) / denom_case * 100
  trend_burden_high[row_i, trend_case_cols] <- colSums(frame_year_case_high[, year_val_cols], na.rm = TRUE) / denom_pop * 100000

  # ---- 4.7 Export ----
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

  write.csv(trend_case,        p(paste0("frame_", key, "_LIMCs_year_10_case_trend.csv")), row.names = FALSE)
  write.csv(trend_paf,         p(paste0("frame_", key, "_LIMCs_year_10_PAF_trend.csv")), row.names = FALSE)
  write.csv(trend_burden,      p(paste0("frame_", key, "_LIMCs_year_10_burden_trend.csv")), row.names = FALSE)
  write.csv(trend_case_low,    p(paste0("frame_", key, "_LIMCs_year_10_case_low_trend.csv")), row.names = FALSE)
  write.csv(trend_paf_low,     p(paste0("frame_", key, "_LIMCs_year_10_PAF_low_trend.csv")), row.names = FALSE)
  write.csv(trend_burden_low,  p(paste0("frame_", key, "_LIMCs_year_10_burden_low_trend.csv")), row.names = FALSE)
  write.csv(trend_case_high,   p(paste0("frame_", key, "_LIMCs_year_10_case_high_trend.csv")), row.names = FALSE)
  write.csv(trend_paf_high,    p(paste0("frame_", key, "_LIMCs_year_10_PAF_high_trend.csv")), row.names = FALSE)
  write.csv(trend_burden_high, p(paste0("frame_", key, "_LIMCs_year_10_burden_high_trend.csv")), row.names = FALSE)

  cat("  done", key, "-> files exported to", out_dir, "\n")

  # Release the large objects of the current exposure to avoid memory accumulation
  rm(list_an, frame, cbTEM, cb_coef, cb_vcov, coef, vcov); gc()
}

cat("\nAll exposures computed and exported.\n")
