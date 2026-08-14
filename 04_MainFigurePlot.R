# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE 
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES

# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/

# ==============================================================================
# SPATIAL VISUALIZATION FRAMEWORK (FIGURE 1)

# Objective: 
# This pipeline generates publication-ready, standardized global maps illustrating 
# country-level study area inclusions (DHS-covered and extrapolated LMICs) alongside 
# the spatial distributions of long-term average exposures to six major Extreme 
# Climatic Events (ECEs).

# DATA NOTICE: 
# Spatial exposure raster inputs (GeoTIFF files) represent grid-level multi-year 
# mean exposure metrics, corresponding to underlying data archived in 'Data_for_Figure_1.csv' 
# and spatial layers processed using ArcGIS Pro. File paths should
# be modified according to local working directory configurations prior to execution.

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - 'tidyverse' package (v2.0.0+): Integrated data manipulation and visualization
#   - 'sf' package (v1.0.0+): Simple features vector spatial data handling
#   - 'terra' package (v1.7.0+): Spatial raster processing and data extraction
#   - 'rnaturalearth' & 'rnaturalearthdata' packages: High-resolution global boundaries

# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

# 2. Define Reusable Custom Theme
# ------------------------------------------------------------------------------
# Centralized theme for consistent map styling across all figures
theme_global_map <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.key.width = unit(2.5, "cm"),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15, face = "bold"),
      
      # Borderless frame with overlaid dashed graticules
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.grid.major = element_line(color = "gray60", linetype = "dashed", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.ontop = TRUE,
      panel.background = element_rect(fill = NA)
    )
}

# 3. Fetch Global Spatial Boundary
# ------------------------------------------------------------------------------
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

# ==============================================================================
# PART 1: Study Area Classification Map
# ==============================================================================

# Country lists
dhs_covered <- c(
  "Albania", "Armenia", "Angola", "Bangladesh", "Bolivia", "Myanmar", "Benin", 
  "Cambodia", "Chad", "Dem. Rep. Congo", "Cameroon", "Comoros", "Colombia", 
  "Central African Rep.", "Dominican Rep.", "Egypt", "Ethiopia", "Gambia", 
  "Gabon", "Ghana", "Guatemala", "Guinea", "Haiti", "Honduras", "India", 
  "Côte d'Ivoire", "Jordan", "Kenya", "Kyrgyzstan", "Liberia", "Lesotho", 
  "Madagascar", "Moldova", "Malawi", "Mali", "Morocco", "Mozambique", "Niger", 
  "Nigeria", "Nepal", "Peru", "Pakistan", "Philippines", "Rwanda", "South Africa", 
  "Senegal", "Sierra Leone", "Thailand", "Tajikistan", "Togo", "Tanzania", 
  "Uganda", "Burkina Faso", "Vietnam", "Namibia", "Eswatini", "Zambia", 
  "Zimbabwe", "Indonesia"
)

dhs_extrapolated <- c(
  "Afghanistan", "Algeria", "Azerbaijan", "Argentina", "Botswana", "Belize",
  "Bosnia and Herz.", "Belarus", "Solomon Is.", "Brazil", "Bhutan", "Burundi",
  "Sri Lanka", "Congo", "China", "Cuba", "Cape Verde", "Djibouti", "Dominica",
  "Ecuador", "Eq. Guinea", "Eritrea", "El Salvador", "Fiji", "Georgia",
  "Grenada", "Palestine", "Iran", "Iraq", "Jamaica", "Korea", "North Korea", "Kiribati",
  "Kazakhstan", "Laos", "Lebanon", "Libya", "Mongolia", "Macedonia", 
  "Marshall Is.", "Maldives", "Mauritius", "Mauritania", "Montenegro", 
  "Mexico", "Malaysia", "Vanuatu", "Suriname", "Nicaragua", "Paraguay", 
  "Papua New Guinea", "Guinea-Bissau", "Serbia", "Sudan", "Syria", "Tonga", 
  "Sao Tome and Principe", "Tunisia", "Turkey", "Taiwan", "Turkmenistan", 
  "Ukraine", "Uzbekistan", "St. Vin. and Gren.", "St. Lucia", "Samoa", 
  "Yemen", "S. Sudan", "Timor-Leste"
)

# Label study area status
world_study_area <- world_sf %>%
  mutate(
    status = case_when(
      name %in% dhs_covered ~ "DHS Covered",
      name %in% dhs_extrapolated ~ "DHS Extrapolated",
      TRUE ~ "Not Included"
    ),
    status = factor(status, levels = c("DHS Covered", "DHS Extrapolated", "Not Included"))
  )

# Plot Study Area Map
p_study_area <- ggplot(data = world_study_area) +
  geom_sf(aes(fill = status), color = "black", linewidth = 0.3) + 
  scale_fill_manual(
    values = c(
      "DHS Covered"      = "#1F78B4",  
      "DHS Extrapolated" = "#A6CEE3", 
      "Not Included"     = "white"     
    ),
    labels = c("DHS Covered (n=59)", "Extrapolated (n=68)", "Not Included"),
    name   = "LMICs category"
  ) +
  coord_sf(expand = FALSE) + 
  theme_global_map() +
  labs(
    title = "(A) Study Area",
    x     = "Longitude",
    y     = "Latitude"
  )

# Save Study Area Map
ggsave(
  filename = file.path("global_average.png"), 
  plot     = p_study_area, 
  width    = 12, 
  height   = 6, 
  dpi      = 300
)

# ==============================================================================
# PART 2: Helper Function for Environmental Exposure Mapping
# ==============================================================================
plot_exposure_raster <- function(tif_filename, 
                                 title_label, 
                                 legend_title, 
                                 color_option = "viridis", 
                                 zero_to_na = FALSE,
                                 output_filename) {
  
  # Load raster and convert to data frame
  r <- terra::rast(file_path)
  r_df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  colnames(r_df) <- c("lon", "lat", "value")
  
  # Optionally convert zero values to NA for cleaner plotting
  if (zero_to_na) {
    r_df$value[r_df$value == 0] <- NA
  }
  
  # Generate plot
  p <- ggplot() +
    geom_raster(data = r_df, aes(x = lon, y = lat, fill = value)) +
    geom_sf(data = world_sf, fill = NA, color = "black", linewidth = 0.3) +
    scale_fill_viridis_c(
      option    = color_option,
      na.value  = "white",
      name      = legend_title,
      direction = 1
    ) +
    coord_sf(expand = FALSE) +
    theme_global_map() +
    labs(
      title = title_label,
      x     = "Longitude",
      y     = "Latitude"
    )
  
  # Save map
  ggsave(
    filename = file.path(output_filename), 
    plot     = p, 
    width    = 12, 
    height   = 6, 
    dpi      = 300
  )
}

# ==============================================================================
# PART 3: Batch Map Generation for Extreme Weather Exposures
# ==============================================================================

# 1. Heatwave
plot_exposure_raster(
  tif_filename    = "frame_hw_LIMCs_month_grid_Raster_mean.tif",
  title_label     = "(B) Heatwave",
  legend_title    = "Number of heatwave day (days)",
  color_option    = "magma",
  zero_to_na      = TRUE,
  output_filename = "global_hw_mean_map.png"
)

# 2. Extreme Precipitation
plot_exposure_raster(
  tif_filename    = "frame_tp_99_LIMCs_month_grid_Raster_mean.tif",
  title_label     = "(C) Extreme precipitation",
  legend_title    = "Cumulative extreme precipitation (mm)",
  color_option    = "mako",
  zero_to_na      = FALSE,
  output_filename = "global_tp_mean_map.png"
)

# 3. Flood
plot_exposure_raster(
  tif_filename    = "frame_flood_LIMCs_month_grid_Raster_mean.tif",
  title_label     = "(D) Flood",
  legend_title    = "Number of flood month (months)",
  color_option    = "viridis",
  zero_to_na      = FALSE,
  output_filename = "global_flood_mean_map.png"
)

# 4. Drought
plot_exposure_raster(
  tif_filename    = "frame_dr_6m_LIMCs_month_grid_revised_Raster_mean.tif",
  title_label     = "(E) Drought",
  legend_title    = "Number of drought months (months)",
  color_option    = "cividis",
  zero_to_na      = FALSE,
  output_filename = "global_dr_6m_mean_map.png"
)

# 5. Tropical Cyclone
plot_exposure_raster(
  tif_filename    = "frame_tc_R34_LIMCs_month_grid_Raster_mean.tif",
  title_label     = "(F) Tropical cyclone",
  legend_title    = "Number of tropical cyclone month (months)",
  color_option    = "mako",
  zero_to_na      = FALSE,
  output_filename = "global_tc_R34_mean_map.png"
)

# 6. Wildfire
plot_exposure_raster(
  tif_filename    = "frame_wf_LIMCs_month_grid_revised_Raster_mean.tif",
  title_label     = "(G) Wildfire",
  legend_title    = "Number of wildfire month (months)",
  color_option    = "magma",
  zero_to_na      = FALSE,
  output_filename = "global_wf_mean_map.png"
)

# ==============================================================================
# FIGURE 2: MULTI-DEFINITION ECE RISK ESTIMATION (FOREST PLOT VISUALIZATION)

# Objective:
# This pipeline generates a standardized multi-panel forest plot visualizing the 
# estimated Odds Ratios (ORs) and 95% Confidence Intervals (CIs) for under-5 child 
# mortality associated with six major Extreme Climatic Events (ECEs) across 
# different operational exposure definitions.

# DATA NOTICE: 
# The input dataset contains point estimates (ORs) and 95% confidence intervals 
# archived in 'Data_for_Figures_2_6.xlsx' (Figure 2 sheet). File paths should be 
# configured locally prior to running the script.

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - 'tidyverse' package (v2.0.0+): Data manipulation and ggplot2 visualization

# ==============================================================================
# Figure 2 Script Execution
# ==============================================================================

# 1. Load Required Packages
# ------------------------------------------------------------------------------
library(tidyverse)

# 2. Ingest and Format Data
# ------------------------------------------------------------------------------
forest_data <- read.csv(data_path, header = TRUE)

# Preserve row order for exposure definitions
forest_data$EWE_define <- factor(forest_data$EWE_define, levels = unique(forest_data$EWE_define))

# Order ECE stressors logically
forest_data$EWE <- factor(
  forest_data$EWE, 
  levels = c("Heatwave", "Extreme Precipitation", "Flood", "Drought", "Tropical Cyclone", "Wildfire")
)

# 3. Define Color Palette for Climate Stressors
# ------------------------------------------------------------------------------
group_colors <- c(
  "Heatwave"              = "darkred", 
  "Extreme Precipitation" = "royalblue", 
  "Flood"                 = "darkgreen",
  "Drought"               = "#D55E00",
  "Tropical Cyclone"      = "#84CAC0",
  "Wildfire"              = "#984EA3"
)

# 4. Construct Forest Plot
# ------------------------------------------------------------------------------
p_forest <- ggplot(forest_data, aes(x = OR_value, y = EWE_define, color = EWE)) +
  # Null effect reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  # Point estimates and horizontal 95% Confidence Intervals
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.3, linewidth = 0.8) +
  geom_point(size = 3.5) +
  
  # Axes scales and limits
  scale_x_continuous(
    limits = c(0.90, 1.30), 
    breaks = seq(0.90, 1.30, 0.05)
  ) +
  scale_color_manual(values = group_colors) +
  
  # Styling and Legends
  theme_classic() +
  labs(
    x = "Odds ratio (95% CI)", 
    y = ""
  ) +
  guides(color = guide_legend(title = "Extreme climatic events")) +
  theme(
    axis.title   = element_text(face = "bold", color = "black", size = 15),
    axis.text    = element_text(face = "bold", color = "black", size = 12),
    legend.title = element_text(face = "bold", color = "black", size = 15),
    legend.text  = element_text(face = "bold", color = "black", size = 15),
    legend.position = "bottom"
  )

# 5. Export High-Resolution Figure
# ------------------------------------------------------------------------------
ggsave(
  filename = output_path,
  plot     = p_forest,
  width    = 11.6,   # ~3500px at 300 DPI
  height   = 8.3,    # ~2500px at 300 DPI
  dpi      = 300,
  device   = "jpeg"
)

# ==============================================================================
# FIGURE 3: REGIONAL HETEROGENEITY OF ECE-ASSOCIATED CHILD MORTALITY RISKS

# Objective:
# This script constructs a grouped forest plot depicting regional stratified Odds 
# Ratios (ORs) and 95% Confidence Intervals (CIs) for under-5 child mortality across 
# six major Extreme Climatic Events (ECEs) compared against global baseline estimates.

# DATA NOTICE: 
# The underlying dataset containing regional stratified point estimates (ORs) and 95% 
# confidence intervals is archived in 'Data_for_Figures_2_6.xlsx' (Figure 3 sheet). File 
# paths should be configured locally prior to running the script.

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - 'tidyverse' package (v2.0.0+): Integrated data manipulation and visualization

# ==============================================================================
# Figure 3 Script Execution
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)

# 2. Ingest and Format Regional Risk Data
# ------------------------------------------------------------------------------
forest_data <- read.csv(data_path, header = TRUE)

# Define explicit ordering for regions (Global pooled + 5 geographic regions)
region_order <- c(
  "Global",
  "Sub-Saharan Africa", 
  "Central & South Asia", 
  "Europe, Middle East & North Africa", 
  "Latin America & Caribbean", 
  "East Asia & Pacific"
)

# Define explicit reverse ordering for climate events (vertical plot stacking)
temp_order <- c(
  "Wildfire", 
  "Tropical cyclone", 
  "Drought", 
  "Flood", 
  "Extreme precipitation", 
  "Heatwave"
)

# Apply factor levels
forest_data$region <- factor(forest_data$region, levels = region_order)
forest_data$temp   <- factor(forest_data$temp, levels = temp_order)

# 3. Define Color Mapping for Geographic Regions
# ------------------------------------------------------------------------------
region_colors <- c(
  "Global"                             = "darkred", 
  "Sub-Saharan Africa"                 = "royalblue", 
  "Central & South Asia"               = "darkgreen",
  "Europe, Middle East & North Africa" = "#D55E00",
  "Latin America & Caribbean"          = "#84CAC0",
  "East Asia & Pacific"                = "#984EA3"
)

# 4. Construct Grouped Forest Plot
# ------------------------------------------------------------------------------
p_region_forest <- ggplot(forest_data, aes(x = OR_value, y = temp, color = region)) +
  # Null reference line (OR = 1.0)
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  # Dodged point estimates and horizontal 95% Confidence Intervals
  geom_errorbarh(
    aes(xmin = OR_low, xmax = OR_high), 
    height   = 0.4, 
    linewidth = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    size     = 3.5, 
    position = position_dodge(width = 0.8)
  ) +
  
  # Axis limits and scale tick marks
  scale_x_continuous(
    limits = c(0.60, 1.80), 
    breaks = seq(0.60, 1.80, 0.20)
  ) +
  scale_color_manual(values = region_colors) +
  
  # Theme and labels
  theme_classic() +
  labs(
    x = "Odds ratio (95% CI)", 
    y = "Extreme Climatic Events"
  ) +
  guides(color = guide_legend(title = "Regions")) +
  theme(
    axis.title   = element_text(face = "bold", color = "black", size = 15),
    axis.text    = element_text(face = "bold", color = "black", size = 12),
    legend.title = element_text(face = "bold", color = "black", size = 13),
    legend.text  = element_text(face = "bold", color = "black", size = 12),
    legend.position = "bottom"
  )

# 5. Export High-Resolution Figure
# ------------------------------------------------------------------------------
ggsave(
  filename = output_path,
  plot     = p_region_forest,
  width    = 11.6,   # ~3500px at 300 DPI
  height   = 8.3,    # ~2500px at 300 DPI
  dpi      = 300,
  device   = "jpeg"
)

# ==============================================================================

# FIGURE 4: GLOBAL AND REGIONAL MORTALITY BURDEN AND TEMPORAL TRENDS (PANELS A, B, C)

# Objective:
# This pipeline constructs Figure 4, which multi-dimensionally quantifies the under-5 
# child mortality burden attributable to six Extreme Climatic Events (ECEs) across LMICs:
#   - Panel 4A: Regional stacked bar plot of attributable mortality rates.
#   - Panel 4B: Spatial raster map showing global distribution of attributable mortality rates.
#   - Panel 4C: Temporal area plot demonstrating Population Attributable Fraction (PAF%) 
#               trends over time (2002–2019).

# DATA NOTICE: 
# Underlying tabular datasets for Panels 4A, 4B, and 4C are archived in 
# 'Data_for_Figures_2_6.xlsx' (Figure 4A and 4C sheet) and "Data_for_Figure_4B.csv".
# File paths should be configured locally prior to execution.

# Methodology Overview:
#   1. Regional Burden Decomposition (4A): Constructs a stacked bar chart displaying 
#      attributable deaths per 100,000 person-years across geographical regions, 
#      utilizing wrapped text labels and custom color palettes.
#   2. High-Resolution Spatial Gridding (4B): Ingests GeoTIFF raster grids via 'terra', 
#      converts spatial coordinates to data frames, overlays Natural Earth national 
#      boundaries via 'sf', and applies a perceptually uniform Magma color scale.
#   3. Temporal PAF Transformation (4C): Dynamically pivots wide-format historical trend 
#      matrixes into tidy long-format structures using 'pivot_longer' to visualize 
#      cumulative PAF percentage shifts across time horizons.
#   4. High-Resolution Graphics Export: Standardizes high-resolution outputs (300 DPI) 
#      using ggsave() across consistent aspect ratios.

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - Packages: 'tidyverse', 'terra', 'rnaturalearth', 'sf', 'RColorBrewer'

# ==============================================================================
# Setup & Library Ingestion
# ==============================================================================

library(tidyverse)
library(terra)
library(rnaturalearth)
library(sf)
library(RColorBrewer)

# Define Shared Palette for ECE Stressors (Diverging RdBu Subset)
ece_palette <- brewer.pal(9, "RdBu")[c(2:4, 6:8)]

# ==============================================================================
# Figure 4A: Regional Attributable Mortality Burden (Stacked Bar Chart)
# ==============================================================================

# Note: Ensure object 'data' contains regional aggregated burden metrics prior to executing
frame_EWE_LIMCs_burden_plot <- data 
write.csv(frame_EWE_LIMCs_burden_plot, file = path_data_a, row.names = FALSE)

# Construct Stacked Bar Chart
p_4a <- ggplot(data, aes(x = region, y = burden, fill = EWE)) +
  geom_bar(
    stat      = "identity", 
    position  = "stack", 
    width     = 0.8, 
    color     = "black", 
    linewidth = 0.3
  ) +
  scale_fill_manual(values = ece_palette) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  theme_bw() +
  labs(
    x    = "Region", 
    y    = "Attributable deaths per 100,000 person-years", 
    fill = "Extreme climatic events"
  ) +
  theme(
    axis.title      = element_text(face = "bold", color = "black", size = 18),
    axis.text       = element_text(face = "bold", color = "black", size = 13),
    legend.title    = element_text(face = "bold", color = "black", size = 18),
    legend.text     = element_text(face = "bold", color = "black", size = 18),
    legend.position = "bottom"
  )

# Export Panel 4A
ggsave(path_out_a, plot = p_4a, width = 12.5, height = 7, dpi = 300)

# ==============================================================================
# Figure 4B: Spatial Distribution of Attributable Mortality (Raster Map)
# ==============================================================================

# Ingest GeoTIFF Spatial Raster
tif_path <- file.path(dir_root, "burden_map_plot/global/frame_EWE_LIMCs_year_10_case_grid_Raster.tif")

if (!file.exists(tif_path)) {
  stop("Target spatial GeoTIFF raster file not found. Check file path: ", tif_path)
}

# Convert Terra SpatRaster to Dataframe
r <- terra::rast(tif_path)
r_df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
colnames(r_df) <- c("lon", "lat", "value")

# Ingest Spatial Country Vectors
world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Construct Spatial Map
p_4b <- ggplot() +
  # Continuous Spatial Grid Layer
  geom_raster(data = r_df, aes(x = lon, y = lat, fill = value)) +
  # Global Administrative Boundaries Layer
  geom_sf(data = world_sf, fill = NA, color = "black", linewidth = 0.3) +
  # Perceptually Uniform Palette Mapping
  scale_fill_viridis_c(
    option   = "magma",
    na.value = "white",
    name     = "Attributable deaths per 100,000 person-years"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude"
  ) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 15, face = "bold"),
    legend.text      = element_text(size = 13),
    legend.key.width = unit(2.5, "cm"),
    axis.text        = element_text(size = 13),
    axis.title       = element_text(size = 15, face = "bold"),
    # Remove outer frames for clean spatial alignment
    panel.border     = element_blank(),
    axis.line        = element_blank(),
    # Top-layered graticule setup
    panel.grid.major = element_line(color = "gray60", linetype = "dashed", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.ontop      = TRUE,
    panel.background = element_rect(fill = NA)
  )

# Export Panel 4B
ggsave(path_out_b, plot = p_4b, width = 12, height = 6, dpi = 300)

# ==============================================================================
# Figure 4C: Temporal Trend of Population Attributable Fraction (Area Plot)
# ==============================================================================

# Ingest Wide-Format Historical Trend Matrix
paf_raw <- read.csv(path_data_c, header = TRUE)

# Tidy Data Transformation: Pivot Wide Matrix to Modern Tidy Dataframe
paf_tidy <- paf_raw %>%
  # Reshape year columns into key-value pairs
  pivot_longer(
    cols      = starts_with("X") | matches("^[0-9]{4}$"), 
    names_to  = "year", 
    values_to = "PAF"
  ) %>%
  # Clean year string formatting (e.g., 'X2002' -> 2002)
  mutate(
    year = as.numeric(gsub("[^0-9]", "", year)),
    PAF  = as.numeric(PAF)
  ) %>%
  # Filter out baseline buffer years (2000-2001) and total aggregated summaries
  filter(
    !year %in% c(2000, 2001),
    EWE != "Total"
  ) %>%
  # Maintain logical factor ordering
  mutate(EWE = factor(EWE, levels = unique(EWE)))

# Construct Stacked Area Plot
p_4c <- ggplot(paf_tidy, aes(x = year, y = PAF, fill = EWE)) +
  geom_area() +
  scale_fill_manual(values = ece_palette) +
  scale_x_continuous(
    name   = "Year",
    breaks = c(2002, 2005, 2010, 2015, 2019),
    labels = c("2002", "2005", "2010", "2015", "2019")
  ) +
  theme_bw() +
  labs(
    x    = "Year", 
    y    = "PAF (%)", 
    fill = "Extreme climatic events"
  ) +
  theme(
    axis.title      = element_text(face = "bold", color = "black", size = 18),
    axis.text       = element_text(face = "bold", color = "black", size = 18),
    legend.title    = element_text(face = "bold", color = "black", size = 18),
    legend.text     = element_text(face = "bold", color = "black", size = 18),
    legend.position = "bottom"
  )

# Export Panel 4C
ggsave(path_out_c, plot = p_4c, width = 12, height = 7, dpi = 300)

# ==============================================================================
# FIGURE 5: GLOBAL SPATIAL DISTRIBUTION OF ATTRIBUTABLE CHILD MORTALITY BURDEN 
# ACROSS SIX EXTREME CLIMATIC EVENTS

# Objective:
# This pipeline processes high-resolution GeoTIFF spatial grids to map attributable 
# under-5 child mortality rates (deaths per 100,000 person-years) globally for 
# six distinct Extreme Climatic Event (ECE) indicators:
#   - Panel 5A: Heatwaves
#   - Panel 5B: Extreme Precipitation
#   - Panel 5C: Floods
#   - Panel 5D: Droughts
#   - Panel 5E: Tropical Cyclones
#   - Panel 5F: Wildfires

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - Packages: 'terra', 'sf', 'rnaturalearth', 'rnaturalearthdata', 'tidyverse'

# ==============================================================================
# Setup & Library Ingestion
# ==============================================================================

library(terra)              # Raster data ingestion and processing
library(sf)                 # Vector spatial boundary support
library(rnaturalearth)      # Natural Earth geographic boundaries
library(rnaturalearthdata)
library(tidyverse)          # Data wrangling and ggplot2 graphics engine

# Ingest Global Country Spatial Boundaries (Shared Across Maps)
world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# ==============================================================================
# Modular Engine Function Definition
# ==============================================================================

#' Plot and Export Spatial Exposure Raster
#'
#' @param input_tif String. Full file path to input GeoTIFF raster.
#' @param output_png String. Full file path for exported PNG image.
#' @param panel_title String. Title overlay for the panel (e.g., "(A) Heatwave").
#' @param color_option String. Viridis color scale option ("magma", "mako", "viridis", "cividis").
#' @param world_bounds sf Object. Spatial polygon data frame of world boundaries.
#'
plot_exposure_raster <- function(input_tif, output_png, panel_title, color_option, world_bounds) {
  
  # 1. Raster Ingestion and Data Frame Conversion
  r <- terra::rast(input_tif)
  r_df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  colnames(r_df) <- c("lon", "lat", "value")
  
  # 2. Construct Map Graphic
  p_map <- ggplot() +
    # Spatial Grid Raster Layer
    geom_raster(data = r_df, aes(x = lon, y = lat, fill = value)) +
    # Administrative Boundaries Overlay
    geom_sf(data = world_bounds, fill = NA, color = "black", linewidth = 0.3) +
    # Perceptually Uniform Palette
    scale_fill_viridis_c(
      option   = color_option,
      na.value = "white",
      name     = "Attributable deaths per 100,000 person-years"
    ) +
    coord_sf(expand = FALSE) +
    theme_minimal() +
    labs(
      title = panel_title,
      x     = "Longitude",
      y     = "Latitude"
    ) +
    theme(
      plot.title       = element_text(size = 22, face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_text(size = 15, lineheight = 0.9),
      legend.text      = element_text(size = 13),
      legend.key.width = unit(2.5, "cm"),
      axis.text        = element_text(size = 13),
      axis.title       = element_text(size = 15, face = "bold"),
      
      # Borderless Frame and Overlaid Graticules
      panel.border     = element_blank(),
      axis.line        = element_blank(),
      panel.grid.major = element_line(color = "gray60", linetype = "dashed", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.ontop      = TRUE,
      panel.background = element_rect(fill = NA)
    )
  
  # 43 Save High-Resolution Map
  ggsave(
    filename = output_png,
    plot     = p_map,
    width    = 12,
    height   = 6,
    dpi      = 300
  )
 
 }

# ==============================================================================
# Configuration Matrix & Batch Execution
# ==============================================================================

# Define Configuration Data Frame for Six ECE Indicators
ece_map_configs <- tibble::tribble(
  ~tif_name,                                  ~out_name,                   ~title,                      ~palette,
  "frame_hw_LIMCs_year_10_case_grid_Raster",    "global_hw_burden_map",      "(A) Heatwave",              "magma",
  "frame_tp_99_LIMCs_year_10_case_grid_Raster", "global_tp_burden_map",      "(B) Extreme precipitation", "mako",
  "frame_flood_LIMCs_year_10_case_grid_Raster", "global_flood_burden_map",   "(C) Flood",                 "viridis",
  "frame_dr_6m_LIMCs_year_10_case_grid_Raster", "global_dr_6m_burden_map",   "(D) Drought",               "cividis",
  "frame_tc_R34_LIMCs_year_10_case_grid_Raster","global_tc_R34_burden_map",  "(E) Tropical cyclone",      "mako",
  "frame_wf_LIMCs_year_10_case_grid_Raster",    "global_wf_burden_map",      "(F) Wildfire",              "magma"
)

# Batch Processing Engine Execution
purrr::pwalk(
  list(
    ece_map_configs$tif_name,
    ece_map_configs$out_name,
    ece_map_configs$title,
    ece_map_configs$palette
  ),
  function(tif, out, title, palette) {
    input_file  <- file.path(base_input_dir, paste0(tif, ".tif"))
    output_file <- file.path(base_output_dir, paste0(out, ".png"))
    
    plot_exposure_raster(
      input_tif    = input_file,
      output_png   = output_file,
      panel_title  = title,
      color_option = palette,
      world_bounds = world_sf
    )
  }
)

# ==============================================================================

# FIGURE 6: TOP 30 COUNTRY-LEVEL ATTRIBUTABLE CHILD MORTALITY NUMBERS AND RATES 
# ACROSS SIX EXTREME CLIMATIC EVENTS

# Objective:
# This script visualizes the top 30 low- and middle-income countries bearing the 
# highest child mortality burden across six extreme climatic event (ECE) types:
#   - Figure 6A: Total annual attributable deaths (1,000s)
#   - Figure 6B: Attributable mortality rate per 100,000 person-years

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - Packages: 'tidyverse', 'RColorBrewer', 'stringr'

# ==============================================================================
# Setup & Library Ingestion
# ==============================================================================

library(tidyverse)          # Core data manipulation and ggplot2 graphing engine
library(RColorBrewer)       # Color palettes for scientific publications
library(stringr)            # String manipulation and text wrapping

# Define Palette for 6 ECE Categories (RColorBrewer 'RdBu' palette subset)
ece_palette <- brewer.pal(9, "RdBu")[c(2:4, 6:8)]

# ==============================================================================
# Modular Engine Function Definition
# ==============================================================================

#' Generate and Save Country-Level Stacked Bar Plot
#'
#' @param data Data frame containing country names, ECE categories, and numeric metric.
#' @param value_var String. Column name of the quantitative metric (e.g., "case_mean", "burden").
#' @param y_label String. Y-axis title for the horizontal bar plot.
#' @param output_path String. File path for exported high-resolution PNG image.
#'
plot_stacked_country_burden <- function(data, value_var, y_label, output_path) {
  
  # 1. Dynamic Sorting: Order countries by total metric sum across all ECEs
  data_ordered <- data %>%
    group_by(country) %>%
    mutate(total_val = sum(.data[[value_var]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(country = reorder(country, total_val))

  # 2. Construct Stacked Horizontal Bar Graphic
  p <- ggplot(data_ordered, aes(x = country, y = .data[[value_var]], fill = EWE)) +
    geom_bar(
      stat      = "identity",
      position  = "stack",
      width     = 0.8,
      color     = "black",
      linewidth = 0.3
    ) +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = ece_palette) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
    labs(
      x    = "Country",
      y    = y_label,
      fill = "Extreme climatic events"
    ) +
    theme(
      axis.title      = element_text(face = "bold", color = "black", size = 25, lineheight = 1),
      axis.text       = element_text(face = "bold", color = "black", size = 19),
      legend.title    = element_text(face = "bold", color = "black", size = 21, lineheight = 1),
      legend.text     = element_text(face = "bold", color = "black", size = 19, lineheight = 1),
      legend.position = "bottom"
    )

  # 3. Save Publication-Grade Graphic (600 DPI)
  ggsave(
    filename = output_path,
    plot     = p,
    width    = 13.8,
    height   = 13.5,
    dpi      = 600
  )
  
}

# ==============================================================================
# Execution & File Ingestion
# ==============================================================================

# 1. Process Absolute Attributable Deaths (Figure 6A)
path_case_data <- file.path(base_dir, "data_case_mean_30.csv")
path_case_out  <- file.path(base_dir, "global_case_mean_plot_country.png")

data_case_mean_30 <- read.csv(path_case_data)
  
plot_stacked_country_burden(
  data        = data_case_mean_30,
  value_var   = "case_mean",
  y_label     = "Annual attributable deaths (1,000)",
  output_path = path_case_out
)


# 2. Process Attributable Mortality Rates (Figure 6B)
path_burden_data <- file.path(base_dir, "data_burden_30.csv")
path_burden_out  <- file.path(base_dir, "global_burden_plot_country.png")

data_burden_30 <- read.csv(path_burden_data)
  
plot_stacked_country_burden(
  data        = data_burden_30,
  value_var   = "burden",
  y_label     = "Attributable deaths per 100,000 person-years",
  output_path = path_burden_out
)


# ==============================================================================
