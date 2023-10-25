# Documentation -----------------------------------------------------------

# Script to import Excel inputs file, set treatment pahts, and re-export as JSON.

# Libraries ---------------------------------------------------------------

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(readxl)
library (rjson)

# Directories -------------------------------------------------------------
proj_path <- getwd()

file_location <- "inputs" # Inputs directory within project

excel_file_name <- "ivi_mdd_inputs_v7.0.xlsx" # Existing inputs file name, assumed to be Excel

json_file_name <- "ivi_mdd_inputs_v7.0.json" # Inputs file name TO BE OUTPUT IN THIS SCRIPT, assumed to be JSON

# Import ------------------------------------------------------------------

## From Excel, with NAs

import_from_excel <- function(object_name, sheet_name){
  
  importtable <- read_excel(file.path(proj_path,file_location, excel_file_name), sheet = eval(sheet_name)) %>% 
    as.data.table() %>%
    na.omit(., cols="Variable") %>%
    select("Variable", "Base", "Min", "Max", "DSA_category", "Label")
  
  importvector.base <- importtable[, mapply(setNames, Base, Variable)]
  importvector.min <- importtable[, mapply(setNames, Min, Variable)]
  importvector.max <- importtable[, mapply(setNames, Max, Variable)]
  importvector.DSACat <- importtable[, mapply(setNames, DSA_category, Variable)]
  importvector.lbl <- importtable[, setNames(Label, Variable)]
  
  assign(object_name, importvector.base, envir = .GlobalEnv)
  assign(paste0(object_name,"_min"), importvector.min, envir = .GlobalEnv)
  assign(paste0(object_name,"_max"), importvector.max, envir = .GlobalEnv)
  assign(paste0(object_name,"_DSACat"), importvector.DSACat, envir = .GlobalEnv)
  assign(paste0(object_name,"_lbl"), importvector.lbl, envir = .GlobalEnv)
}

object_names <- c("inputs_general","inputs_efficacy","inputs_gaps","inputs_costs","inputs_utilities")
sheet_names <- c("General", "Efficacy", "Gaps", "Costs", "Utilities")

invisible(mapply(import_from_excel, object_names, sheet_names))
  
## From Excel, no NAs

inputs_mortality <- read_excel(file.path(proj_path,file_location, excel_file_name), sheet = "Mortality") %>%
  as.data.table()
#inputs_mortality <- inputs_mortality[, mapply(setNames, Age, Female, p_death)]

## Manual
days_per_yr <- data.table("Variable" = "days_per_yr", "Base" = 365.25)
days_per_yr <- days_per_yr[, mapply(setNames, Base, Variable)]

#  Set Treatment Paths ------------------------------------------------------------------

## Paths

txPaths <- list("pathway1" = rep("No active treatment", 4),
                "pathway2" = rep("SSRI", 4),
                "pathway3" = c("SSRI", rep("SSRI + psychotherapy", 3)),
                "pathway4" = c("SSRI", "SSRI", rep("SSRI + psychotherapy", 2)),
                "pathway5" = c("SSRI", "SNRI", "SNRI + atypical antidepressant", "SNRI + antipsychotic"))

## Line 4 somatic add-on (Yes/No)

ln4.addonsom.YNs <- list(
  pathway1_l4.somatic = "No",
  pathway2_l4.somatic = "No",
  pathway3_l4.somatic = "No",
  pathway4_l4.somatic = "Yes",
  pathway5_l4.somatic = "No")

# Build List ------------------------------------------------------------------

json_list <- list(
  # Base values and mortality data
  "inputs_general" = inputs_general,
  "inputs_efficacy" = inputs_efficacy,
  "inputs_gaps" = inputs_gaps,
  "inputs_costs" = inputs_costs,
  "inputs_utilities" = inputs_utilities,
  "inputs_mortality" = inputs_mortality,
  # Minimum allowable values
  "inputs_general_min" = inputs_general_min,
  "inputs_efficacy_min" = inputs_efficacy_min,
  "inputs_gaps_min" = inputs_gaps_min,
  "inputs_costs_min" = inputs_costs_min,
  "inputs_utilities_min" = inputs_utilities_min,
  # Maximum allowable values
  "inputs_general_max" = inputs_general_max,
  "inputs_efficacy_max" = inputs_efficacy_max,
  "inputs_gaps_max" = inputs_gaps_max,
  "inputs_costs_max" = inputs_costs_max,
  "inputs_utilities_max" = inputs_utilities_max,
  # DSA categories
  "inputs_general_DSACat" = inputs_general_DSACat,
  "inputs_efficacy_DSACat" = inputs_efficacy_DSACat,
  "inputs_gaps_DSACat" = inputs_gaps_DSACat,
  "inputs_costs_DSACat" = inputs_costs_DSACat,
  "inputs_utilities_DSACat" = inputs_utilities_DSACat,
  # Parameter labels
  "inputs_general_lbl" = inputs_general_lbl,
  "inputs_efficacy_lbl" = inputs_efficacy_lbl,
  "inputs_gaps_lbl" = inputs_gaps_lbl,
  "inputs_costs_lbl" = inputs_costs_lbl,
  "inputs_utilities_lbl" = inputs_utilities_lbl,
  # Other
  "days_per_yr" = days_per_yr,
  "txPaths" = txPaths,
  "ln4.addonsom.YNs" = ln4.addonsom.YNs
  )

# Export ------------------------------------------------------------------

write(toJSON(json_list, indent = 1), file.path(proj_path,file_location, json_file_name))

