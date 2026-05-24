# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE 
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES

# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/

# SCRIPT: Exposure Definition, Spatiotemporal Exposure Matching

# Objective:
#   This comprehensive technical pipeline links  gridded environmental time-series 
#   with individual-level Demographic and Health Surveys Database. 
#   It automates localized climate extreme thresholding, spatial coordinate anchoring,  
#   and the construction of multi-definition retrospective rolling lag matrices.

# NOTE ON DATA PREPARATION:
#   Data source for Temperature, precipitation, and the other four extreme climatic events (floods, droughts, 
#   tropical cyclones, and wildfires) are all gridded datasets. Prior to running this 
#   pipeline, long-term daily time-series data during 1980 to 2023 have already been extracted from the corresponding 
#   grid cells based on the centroid coordinates of each individual's spatial DHS cluster. 
#   This code demonstrates the core subsequent processes: defining multi-dimensional extreme climatic 
#   exposures and performing retrospective spatiotemporal matching based on these pre-extracted 
#   time-series datasets.

# Methodology Overview:
#   1. Localized Metric Thresholding & Run-Length Engine:
#      Parses continuous daily gridded climate inputs cell-by-cell to compute relative thresholds (e.g., 90th/95th 
#      percentiles). A vectorized run-length look-ahead algorithm is executed to isolate discrete extreme weather spells based 
#      on dual constraints of intensity and consecutive duration.
#   2. Spatial Matching & Data Chunking:
#      Extracts unique coordinate centroids from individual health clusters. It implements an optimized 1st-Nearest 
#      Neighbor search algorithm (`FNN::knnx`) to map clusters instantaneously to the closest climate grid node. 
#      Geodesic displacement distances are recorded for quality control. 
#   3. Multi-Definition & Longitudinal Window Allocation:
#      Constructs parallel exposure series for 6 core climate anomalies across varied indices (e.g., event 
#      occurrence vs continuous volume/intensity). For each child observation, the pipeline slices backward 
#      from the precise index date (date of death or matched reference control date) to extract cumulative monthly 
#      dosages across a 13-month retrospective window (Lag 0m through Lag 12m).

# NOTICE: 
# This script utilizes a subset of data from Tanzania as a representative demonstration. 
# Due to the limited sample size and reduced statistical power of this illustrative dataset, 
# modeling estimates (DLNM lag structures) may exhibit instability or non-convergence. 
# Full analytical robustness requires the complete multi-country pooled dataset.

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - 'dplyr' package (v1.1.4): Grammar of Data Manipulation & Vectorized Wrangling
#   - 'lubridate' package (v1.9.3): Functions to Simplify Date-Time Manipulation
#   - 'FNN' package (v1.1.2.2): Fast Nearest Neighbor Search Coordinate Indexing
#=============================================================================================================================
## MODULE 1: HEATWAVE IDENTIFICATION
#=============================================================================================================================
library(foreign)
library(FNN)
library(lubridate)
library(dplyr)

#' Identify Heatwave Events via Vectorized Run-Length Encoding
#' the IDHeatwaves1 function is to automatically identify and construct customized heatwave exposure series 
#' from baseline daily temperature time series.

#' @param threshold Numeric. The temperature threshold defining an extreme day (e.g., localized 90th percentile).
#' @param data Data.frame. Time-series framework containing at least 'date' and 'temp' columns.
#' @param numDays Integer. The minimum consecutive duration constraint (e.g., >= 2 or 3 days).
#' 
#' @return A data.frame aligned to the primary temporal scale containing:
#'   \item{hw1}{Binary indicator profile (1 = Active heatwave day, 0 = Non-event day)}
#'   \item{hw.number1}{Sequential event tracking index (Unique identification number for each discrete spell)}

IDHeatwaves1 <-
function (threshold, data, numDays){
    days <- numDays

    # Convert data frame to a matrix
    data <- cbind(data$date,data$temp)
    colnames(data) <- c("date", "temp")

    # Generate binary indicator vector for threshold exceedance
    tempsExceedingthreshold <- as.numeric(data[, 2] >= threshold)

    # Boundary Fix: Append a 0 at the end so an active heatwave at the very end can terminate properly
    tempsExceedingthreshold <- c(tempsExceedingthreshold, 0)

    # Set up the target pattern consecutive heatwave days
    HeatwaveForm <- rep(1, days) 
    counter <- 1 # ID counter for unique heatwave events
    hwBound <- days - 1 # Look-ahead window size

    # Initialize output data frame with a temporary placeholder row (9, 9)
    hwInfo <- data.frame(hw1 = c(9), hw.number1 = c(9))# Initialize storage matrix with placeholder values

    i <- 1

    # Loop through the entire time series day by day
    while (i <= nrow(data)) {

      # Check if the next N days match our heatwave pattern
        if (identical(tempsExceedingthreshold[i:(i + hwBound)], HeatwaveForm)) {

            # --- CASE A: HEATWAVE FOUND ---
            # Find the next '0' to see exactly how many days this heatwave lasts
            size <- match(0, tempsExceedingthreshold[-(1:i)]) 

            # Record the whole heatwave block (1s for exposure, counter for event ID)
            hwInfo <- data.frame(hw1 = c(hwInfo[, 1], rep(1, size)), 
                hw.number1 = c(hwInfo[, 2], rep(counter, size)))

            counter <- counter + 1 # Move to next event ID
            i <- i + size # Skip the iterator past this entire heatwave block

        } else {

            # --- CASE B: NO HEATWAVE ---
            # Log a single normal day (0, 0) and advance by 1 day
            hwInfo <- rbind(hwInfo, c(0, 0)) 
            i <- i + 1
        }
    }

    # Drop the first placeholder row (9, 9) and return the clean data frame
    return(data.frame(hwInfo[-1, ])) # Drop placeholder row and return time-series matrix
}

# =============================================================================================
# EXPOSURE PROCESSING PIPELINE: MULTI-DEFINITION HEATWAVE EXTRACTION
# =============================================================================================

# Define the relative path to the directory containing the environmental baseline lists
# (Using a generic relative data path instead of a localized system root directory)
input_dir <- "data/ECE_baseline_lists/"

# Retrieve all standardized RData file names within the target directory
file_names <- list.files(path = input_dir, pattern = "\\.RData$", full.names = FALSE)

for (file_idx in file_names) {
    
    # Standardized Console Logging
    message("Processing File: ", file_idx)
    
    # Load the workspace object (expected to contain 'list_temp_country_year')
    load(file.path(input_dir, file_idx))
    
    # Initialize an empty list to store computed heatwave exposure frames for the current file
    list_hw_country_year <- list()
    
    # Iterate through each geographic/administrative unit (e.g., spatial clusters)
    for (spatial_unit in names(list_temp_country_year)) {
        
        message("  -> Processing Spatial Unit: ", spatial_unit)
        
        # Extract the baseline climate data frame for the current spatial unit
        frame_temp <- list_temp_country_year[[spatial_unit]]
        
        # Isolate temporal and thermal attributes for exposure assessment
        frame_temp <- frame_temp[, c(1, 2)]
        colnames(frame_temp) <- c("date", "temp")
        frame_temp$date      <- as.Date(frame_temp$date)
        
        # -------------------------------------------------------------------------------------
        # Methodological Rigor: Calculate 4 combinations of heatwave definitions
        # Matrix: 2 Intensity Thresholds (90th vs 95th %ile) * 2 Duration Constraints (2 vs 3 days)
        # -------------------------------------------------------------------------------------
        # Definition 1: 90th percentile threshold, sustained for >= 2 consecutive days
        frame_hw_90_2 <- IDHeatwaves1(threshold = quantile(frame_temp$temp, probs = 0.90, na.rm = TRUE), 
                                      data = frame_temp, numDays = 2)
        
        # Definition 2: 95th percentile threshold, sustained for >= 2 consecutive days
        frame_hw_95_2 <- IDHeatwaves1(threshold = quantile(frame_temp$temp, probs = 0.95, na.rm = TRUE), 
                                      data = frame_temp, numDays = 2)
        
        # Definition 3: 90th percentile threshold, sustained for >= 3 consecutive days
        frame_hw_90_3 <- IDHeatwaves1(threshold = quantile(frame_temp$temp, probs = 0.90, na.rm = TRUE), 
                                      data = frame_temp, numDays = 3)
        
        # Definition 4: 95th percentile threshold, sustained for >= 3 consecutive days
        frame_hw_95_3 <- IDHeatwaves1(threshold = quantile(frame_temp$temp, probs = 0.95, na.rm = TRUE), 
                                      data = frame_temp, numDays = 3)
        
        # -------------------------------------------------------------------------------------
        # Reconstruct and merge back to spatial clusters
        # -------------------------------------------------------------------------------------
        frame_hw <- list_temp_country_year[[spatial_unit]]
        
        # Append the binary exposure vectors for heatwave (hw1) into the master data frame
        frame_hw$hw_90_2 <- frame_hw_90_2$hw1
        frame_hw$hw_95_2 <- frame_hw_95_2$hw1
        frame_hw$hw_90_3 <- frame_hw_90_3$hw1
        frame_hw$hw_95_3 <- frame_hw_95_3$hw1
        
        # Store the processed exposure grid back into the compiled country-year list
        list_hw_country_year[[spatial_unit]] <- frame_hw
    }
    
    # Save the processed exposure matrix list to a generic output directory
    # output_filename <- paste0("processed_hw_", file_idx)
    # save(list_hw_country_year, file = file.path("data/processed_exposure_lists", output_filename))
}

# =============================================================================================
# KNN MATCHING & RETROSPECTIVE ROLLING LAGS FOR HEATWAVE
#==============================================================================================

library(FNN)   # Fast Nearest Neighbor execution
library(dplyr) # Clean and structured data manipulation

# Initialize master list to store finalized cohort exposure data frames
list_u5mr_hw <- list()

# Define generic relative workspace directories instead of localized root paths
input_dir  <- "data/spatial_merge_frames"
output_dir <- "data/processed_exposure_output"

# Loop through each standardized file identifier (e.g., country-year dynamic lists)
for (file_token in country_year) {
    
    # Standardized Console Logging
    message("------------------------------------------------------------")
    message("Executing Spatiotemporal Exposure Matching for: ", file_token)
    message("------------------------------------------------------------")
    
    # Secure workspace loading via relative paths
    load(file.path(input_dir, paste0(file_token, ".rda")))
    
    # Extract structural components from the loaded database environment
    data     <- list_u5mr_EW[[file_token]]
    exposure <- do.call(rbind, list_hw_country_year)
    
    # -----------------------------------------------------------------------------------------
    # 1. VECTORIZED SPATIAL GRID METADATA GENERATION
    # -----------------------------------------------------------------------------------------
    # Isolate unique grid coordinates to build the environmental anchor profile
    exposure_extraction <- distinct(exposure[, c(3, 4)])
    
    # Replaces the inefficient for-loop with a clean, vectorized framework
    grid_info <- data.frame(
        grid_code = seq_len(nrow(exposure_extraction)),
        lon          = exposure_extraction[[1]],
        lat          = exposure_extraction[[2]]
    )
    
    # Enforce standard column schemas across the long-format exposure grid
    colnames(exposure) <- c("date", "temp", "lon", "lat", "hw_90_2", "hw_95_2", "hw_90_3", "hw_95_3")
    
    # Map calculated grid codes back to the long environmental records
    exposure <- merge(grid_info, exposure, by = c("lon", "lat"))
    
    # Vectorized split: Slice long-format data into sub-dataframes named by grid_code
    data_list <- split(exposure, exposure$grid_code)
    
    # -----------------------------------------------------------------------------------------
    # 2. COORDINATE DISTANCE MATCHING VIA KNN
    # -----------------------------------------------------------------------------------------
    # Execute 1st-Nearest Neighbor lookup to match child clusters to closest weather grids
    nearest.index <- knnx.index(grid_info[, c("lat", "lon")], data[, c("LATNUM", "LONGNUM")], k = 1)
    nearest.dist  <- knnx.dist(grid_info[, c("lat", "lon")], data[, c("LATNUM", "LONGNUM")], k = 1)
    
    # Append the localized pointer indices directly onto individual child profiles
    data$grid  <- grid_info$grid_code[nearest.index]
    data$distance <- nearest.dist[, 1]
    
    # -----------------------------------------------------------------------------------------
    # 3. CHRONOLOGICAL WINDOW EXTRACTION (13 RETROSPECTIVE LAG MONTHS * 4 CO-DEFINITIONS)
    # -----------------------------------------------------------------------------------------
    # Allocate empty memory matrix upfront (Rows: child observations, Columns: 52 exposure indicators)
    output <- matrix(NA, nrow = nrow(data), ncol = 13 * 4)
    
    for (l in seq_len(nrow(data))) {
        if (l %% 100 == 0) {
            message("  -> Processing retrospective matrix: Individual ", l, " / ", nrow(data))
        }
        
        # Enforce character extraction to prevent numeric positional indexing drift
        site <- data_list[[as.character(data$grid[l])]]
        
        output[l, ] <- tryCatch({
            # Columns 6 to 9 hold the binary exposure flags for the 4 distinct heatwave definitions
            c(
                colSums(site[site$date > data$death_date_lag_1m[l]  & site$date <=  data$death_date_lag_1m[l], 6:9], na.rm = TRUE), # lag_0m
                colSums(site[site$date > data$death_date_lag_2m[l]  & site$date <=  data$death_date_lag_1m[l], 6:9], na.rm = TRUE), # lag_1m
                colSums(site[site$date > data$death_date_lag_3m[l]  & site$date <=  data$death_date_lag_2m[l], 6:9], na.rm = TRUE), # lag_2m
                colSums(site[site$date > data$death_date_lag_4m[l]  & site$date <=  data$death_date_lag_3m[l], 6:9], na.rm = TRUE), # lag_3m
                colSums(site[site$date > data$death_date_lag_5m[l]  & site$date <=  data$death_date_lag_4m[l], 6:9], na.rm = TRUE), # lag_4m
                colSums(site[site$date > data$death_date_lag_6m[l]  & site$date <=  data$death_date_lag_5m[l], 6:9], na.rm = TRUE), # lag_5m
                colSums(site[site$date > data$death_date_lag_7m[l]  & site$date <=  data$death_date_lag_6m[l], 6:9], na.rm = TRUE), # lag_6m
                colSums(site[site$date > data$death_date_lag_8m[l]  & site$date <=  data$death_date_lag_7m[l], 6:9], na.rm = TRUE), # lag_7m
                colSums(site[site$date > data$death_date_lag_9m[l]  & site$date <=  data$death_date_lag_8m[l], 6:9], na.rm = TRUE), # lag_8m
                colSums(site[site$date > data$death_date_lag_10m[l] & site$date <=  data$death_date_lag_9m[l], 6:9], na.rm = TRUE), # lag_9m
                colSums(site[site$date > data$death_date_lag_11m[l] & site$date <=  data$death_date_lag_10m[l],6:9], na.rm = TRUE), # lag_10m
                colSums(site[site$date > data$death_date_lag_12m[l] & site$date <=  data$death_date_lag_11m[l],6:9], na.rm = TRUE), # lag_11m
                colSums(site[site$date > data$death_date_lag_12m[l] & site$date <=  data$death_date_lag_12m[l],6:9], na.rm = TRUE)  # lag_12m 
            )
        }, error = function(e) rep(NA, 13 * 4))
    }
    
    # -----------------------------------------------------------------------------------------
    # 4. COLUMNS RE-LABELING & MASTER LIST COMPILATION
    # -----------------------------------------------------------------------------------------
    pollutants <- c("hw_90_2", "hw_95_2", "hw_90_3", "hw_95_3")
    terms      <- paste0("lag_", 0:12, "m") # Aligns directly with tracking window assignments
    
    # Formulate explicit column combinations
    cols             <- paste0(rep(terms, each = 4), "_", rep(pollutants, 13))
    colnames(output) <- cols
    output           <- data.frame(output)
    
    # Combine exposure data back with the full epidemiological file
    DHS_data_u5mr             <- list_u5mr_EW[[file_token]]
    list_u5mr_hw[[file_token]] <- data.frame(DHS_data_u5mr, output)
}

# =============================================================================================
# SPATIAL MATCHING & RETROSPECTIVE ROLLING MONTHLY LAGS FOR EXTREME PRECIPITATION
# =============================================================================================
library(FNN)   # Fast Nearest Neighbor execution
library(dplyr) # Clean and structured data manipulation

# Initialize master list to store finalized cohort exposure data frames
list_u5mr_EP <- list()

# Define generic relative workspace directories instead of localized root paths
input_dir  <- "data/spatial_merge_frames"
output_dir <- "data/processed_exposure_output"

for (file_token in country_year) {
  
  # Standardized Console Logging
  message("------------------------------------------------------------")
  message("Executing Spatiotemporal EP Exposure Matching for: ", file_token)
  message("------------------------------------------------------------")
  
  # Secure workspace loading via relative paths
  load(file.path(input_dir, paste0(file_token, ".rda")))
  data <- list_u5mr_EW[[file_token]]
  
  # -----------------------------------------------------------------------------------------
  # 1. VECTORIZED SPATIAL grid METADATA GENERATION
  # -----------------------------------------------------------------------------------------
  # Isolate unique grid coordinates to build the environmental anchor profile
  exposure_extraction <- distinct(exposure_raw[, c(2, 3)])
  
  # Vectorized generation of the grid metadata data frame
  grid_info <- data.frame(
    grid_code = seq_len(nrow(exposure_extraction)),
    lon          = exposure_extraction[[1]],
    lat          = exposure_extraction[[2]]
  )
  
  # Standardize column naming conventions for incoming long-format precipitation grid
  # Assumes structure: 1: date, 2: lon, 3: lat, 4: daily_tp
  colnames(exposure_raw)[1:6] <- c('date', 'lon', 'lat', 'daily_tp', "extr_tp_day_99", "extr_tp_day_95"))
  
  # -----------------------------------------------------------------------------------------
  # 2. METRIC UNIT CONVERSION & EXPOSURE VARIABLE CONSTRUCTION
  # -----------------------------------------------------------------------------------------
  # Standardize precipitation volume metrics (e.g., converting scale to mm/day)
  exposure_raw$daily_tp <- exposure_raw$daily_tp * 1000
  
  # Build continuous volumetric metrics conditioned on the binary threshold flags
  exposure_raw$extr_tp_day_99 <- as.numeric(exposure_raw$extr_tp_day_99)
  exposure_raw$extr_tp_day_95 <- as.numeric(exposure_raw$extr_tp_day_95)
  
  exposure_raw$extr_tp_99 <- ifelse(exposure_raw$extr_tp_day_99 == 1, exposure_raw$daily_tp, 0)
  exposure_raw$extr_tp_95 <- ifelse(exposure_raw$extr_tp_day_95 == 1, exposure_raw$daily_tp, 0)
  
  # Enforce strict column sequencing to align with the matrix extraction step below
  # Final order: 1.grid_code, 2.lon, 3.lat, 4.date, 5.daily_tp, 6:9. Core EP Indicators
  exposure_tem <- merge(grid_info, exposure_raw, by = c('lon', 'lat'))
  
  exposure <- exposure_tem %>% 
    select(grid_code, lon, lat, date, daily_tp, 
           extr_tp_day_99, extr_tp_day_95, extr_tp_99, extr_tp_95)
  
  # High-speed data splitting: slice grid records into sub-frames keyed by grid code
  data_list <- split(exposure, exposure$grid_code)
  
  # -----------------------------------------------------------------------------------------
  # 3. COORDINATE DISTANCE MATCHING VIA KNN
  # -----------------------------------------------------------------------------------------
  # Execute 1st-Nearest Neighbor lookup to match cohort cluster keys to closest grid
  nearest.index <- knnx.index(grid_info[, c("lat", "lon")], data[, c("LATNUM", "LONGNUM")], k = 1)
  nearest.dist  <- knnx.dist(grid_info[, c("lat", "lon")], data[, c("LATNUM", "LONGNUM")], k = 1)
  
  # Append localized pointer indices and geodesic distance metrics to individual child rows
  data$grid  <- grid_info$grid_code[nearest.index]
  data$distance <- nearest.dist[, 1]
  
  # -----------------------------------------------------------------------------------------
  # 4. CHRONOLOGICAL WINDOW EXTRACTION (13 RETROSPECTIVE LAG MONTHS * 4 CO-DEFINITIONS)
  # -----------------------------------------------------------------------------------------
  # Allocate memory upfront for the matrix (Rows: observations, Columns: 13 lags * 4 metrics = 52)
  total_lags  <- 13
  num_metrics <- 4
  output      <- matrix(NA, nrow = nrow(data), ncol = total_lags * num_metrics)
  
  for (l in seq_len(nrow(data))) {
    if (l %% 100 == 0) {
      message("  -> Processing retrospective matrix: Individual ", l, " / ", nrow(data))
    }
    
    # Force character coercion to safeguard against positional index drifting
    site <- data_list[[as.character(data$grid[l])]]
    
    output[l, ] <- tryCatch({
      # Dynamically compute column sums across designated indices [6:9] for each lag window
      c(
        colSums(site[site$date > data$death_date_lag_1m[l]  & site$date <=  data$death_date_lag_0m[l], 6:9], na.rm = TRUE), # lag_0m
        colSums(site[site$date > data$death_date_lag_2m[l]  & site$date <=  data$death_date_lag_1m[l], 6:9], na.rm = TRUE), # lag_1m
        colSums(site[site$date > data$death_date_lag_3m[l]  & site$date <=  data$death_date_lag_2m[l], 6:9], na.rm = TRUE), # lag_2m
        colSums(site[site$date > data$death_date_lag_4m[l]  & site$date <=  data$death_date_lag_3m[l], 6:9], na.rm = TRUE), # lag_3m
        colSums(site[site$date > data$death_date_lag_5m[l]  & site$date <=  data$death_date_lag_4m[l], 6:9], na.rm = TRUE), # lag_4m
        colSums(site[site$date > data$death_date_lag_6m[l]  & site$date <=  data$death_date_lag_5m[l], 6:9], na.rm = TRUE), # lag_5m
        colSums(site[site$date > data$death_date_lag_7m[l]  & site$date <=  data$death_date_lag_6m[l], 6:9], na.rm = TRUE), # lag_6m
        colSums(site[site$date > data$death_date_lag_8m[l]  & site$date <=  data$death_date_lag_7m[l], 6:9], na.rm = TRUE), # lag_7m
        colSums(site[site$date > data$death_date_lag_9m[l]  & site$date <=  data$death_date_lag_8m[l], 6:9], na.rm = TRUE), # lag_8m
        colSums(site[site$date > data$death_date_lag_10m[l]  & site$date <=  data$death_date_lag_9m[l], 6:9], na.rm = TRUE), # lag_9m
        colSums(site[site$date > data$death_date_lag_11m[l]  & site$date <=  data$death_date_lag_10m[l], 6:9], na.rm = TRUE), # lag_10m
        colSums(site[site$date > data$death_date_lag_12m[l]  & site$date <=  data$death_date_lag_11m[l], 6:9], na.rm = TRUE), # lag_11m
        colSums(site[site$date >= data$death_date_lag_12m[l]  & site$date <=  data$death_date_lag_12m[l],6:9], na.rm = TRUE)  # lag_12m
      )
    }, error = function(e) rep(NA, total_lags * num_metrics))
  }
  
  # -----------------------------------------------------------------------------------------
  # 5. COLUMN RE-LABELING & MASTER LIST COMPILATION
  # -----------------------------------------------------------------------------------------
  pollutants <- c("extr_tp_day_99", "extr_tp_day_95", "extr_tp_99", "extr_tp_95")
  terms      <- paste0("lag_", 0:12, "m")
  
  # Combine indicators and timelines systematically to form 52 structured variable labels
  cols             <- paste0(rep(terms, each = length(pollutants)), "_", rep(pollutants, length(terms)))
  colnames(output) <- cols
  output           <- data.frame(output)
  
  # Merge exposure matrix with epidemiological cohort dataframe and save back to main list
  DHS_data_u5mr         <- list_u5mr_EW[[file_token]]
  list_u5mr_EP[[file_token]] <- data.frame(DHS_data_u5mr, output)
}