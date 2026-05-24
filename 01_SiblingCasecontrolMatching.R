# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE 
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES

# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/

# SCRIPT: Sibling-Matched Case-Control Strata & Temporal Alignment Pipeline
# 
# Objective: 
#   This pipeline extracts data from Demographic and Health Surveys (DHS) 
#   to construct sibling-matched case-control strata (clusters grouped by 
#   maternal lines) for conditional risk estimation (conditional logistic regression and DLNM analysis).

# Methodology Overview:
#   1. Case Screening: Isolate index cases that experienced under-5 mortality events (u5mr == 1).
#   2. Temporal Alignment (Lag 0 Anchor): Define the exact calendar date of the death event using survival 
#      duration (death_age_month) and birth dates.
#   3. Informative Strata Pruning (Quality Control): Group records by unique maternal IDs (`caseid`) to match 
#      deceased cases with surviving controls (siblings). Exclude uninformative "dead strata" (families with 
#      only mortality events and no surviving controls) that provide zero statistical contrast to the conditional 
#      likelihood function.

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - 'dplyr' package (v1.1.4): Grammar of Data Manipulation & Vectorized Wrangling
#   - 'lubridate' package (v1.9.3): Functions to Simplify Date-Time Manipulation

library(dplyr) library(lubridate)
#==============================================================================================================================================================================
# STEP 1: Load the Survey Dataset
# Load the pre-structured data containing under-5 child mortality records
load("frame_DHS_health_under5_mortality.rda")

# STEP 2: Identify and Isolate Mortality Cases (Case Selection)
# Subset the cohort to isolate historical individual records that experienced the mortality event (u5mr == 1)
frame_u5mr_1 <- frame_u5mr[which(frame_u5mr$u5mr == 1),]

# Eliminate potential redundant duplicates within the same birth event to ensure strict unique index cases
frame_u5mr_1 <- frame_u5mr_1 %>% distinct(caseid, birth_date, .keep_all = T)

# Extract unique maternal identifiers (caseid) to define the sibling-matched clusters (strata)
caseid_1 <- as.character(frame_u5mr_1$caseid)

# Initialize an empty container to compile the time-aligned case-control strata matrix
frame_u5mr_merge <- c()

# STEP 3: Iterate Within Strata to Establish Temporal Alignment & Imputation
# Loop through each unique family stratum to calculate exact event dates and match controls (siblings)
for (j in caseid_1){

	# Extract all child records sharing the same maternal lineage (matching cases and controls within the family stratum)
	frame_u5mr_1_0 <- frame_u5mr[which(frame_u5mr$caseid == j),]

	# Operationalize the specific age at death within the family stratum (handling multi-child death clusters if present)
	frame_u5mr_1_0$death_date <- frame_u5mr_1_0$birth_date + frame_u5mr_1_0$death_age_month * 30

    # Incrementally bind the structured family stratum back into the master dataset
	frame_u5mr_merge <- rbind(frame_u5mr_merge, frame_u5mr_1_0)

}

# STEP 4: Finalize the Sibling-Matched Case-Control Strata Matrix
# Perform a final deduplication check to ensure that each unique individual child within the 
# family-stratified design occupies a single, unambiguous row for conditional risk estimation.
frame_u5mr_merge <- frame_u5mr_merge %>% distinct(caseid, birth_date, .keep_all = T)

# STEP 5: Optimized Filtering of Uninformative Sibling Clusters (Strata Pruning)
# In sibling-matched or conditional logistic designs, a cluster must possess internal contrast 
# (at least one surviving control, where 'alive == 1') to contribute information to the conditional likelihood. 
frame_u5mr_merge <- frame_u5mr_merge %>%
  group_by(caseid) %>%
  # Keep only informative strata that contain at least one strictly surviving control (alive == 1)
  # This automatically prunes uninformative clusters containing only mortality cases (sum == 0)
  filter(sum(alive == 1, na.rm = TRUE) > 0) %>%
  ungroup()

