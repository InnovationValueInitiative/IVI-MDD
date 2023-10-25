# Documentation -----------------------------------------------------------

# Script to import import JSON inputs file, run model, and export results to either Excel or JSON.

# Libraries ---------------------------------------------------------------

library(openxlsx)
library(doParallel)
library(dampack)
library(rjson)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

#setup parallel backend to use 8 processors
cores <- detectCores(all.tests = FALSE, logical = TRUE)
cl<-makeCluster(cores - 1)
registerDoParallel(cl)


# Directories -------------------------------------------------------------
proj_path <- getwd()

file_location <- "inputs" # Inputs directory within project

file_name <- "ivi_mdd_inputs_v7.0.json" # Inputs file name, assumed to be JSON

results_destination <- "results" # Result directory within project

## Other scripts needed to execute simulation
source("process_inputs.R") # inputs cleaning and intermediate calculations based on inputs file
source("ivi_mdd_sim.R") # model engine
source("summarize_results.R") # results compilation
source("dsa.R") # DSA


# Tx name to abbreviation mapping -----------------------------------------

map.tx_name.tx_abbrev <- data.table(
  tx.name = c("SSRI",
              "SSRI + atypical antidepressant",
              "SSRI + psychotherapy",
              "SSRI + antipsychotic",
              "SNRI",
              "SNRI + atypical antidepressant",
              "SNRI + psychotherapy",
              "SNRI + antipsychotic",
              "Atypical antidepressant",
              "Atypical antidepressant + psychotherapy",
              "Atypical antidepressant + antipsychotic",
              "Psychotherapy",
              "Placeholder therapy",
              "No active treatment"),
  tx.abbrev = c("ssri",
                "ssri_atypantidep",
                "ssri_psycther",
                "ssri_antipsych",
                "snri",
                "snri_atypantidep",
                "snri_psycther",
                "snri_antipsych",
                "atypantidep",
                "atypantidep_psycther",
                "atypantidep_antipsych",
                "psycther",
                "phtx",
                "notreat"),
  tx.pharm = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0)
)

# Summary measures
summMeas.tmp <- data.table(
  var1a = c("% completing >=2 or exhausting all available lines of therapy", "any.trd"),
  var1bi = c("% experiencing >=1 serious AE", "any.sae"),
  var1bii = c("% experiencing >=1 placeholder AE 1", "any.ae1"),
  var1biii = c("% experiencing >=1 placeholder AE 2", "any.ae2"),
  var1biv = c("% experiencing >=1 placeholder AE 3", "any.ae3"),
  var1bv = c("% experiencing >=1 placeholder AE 4", "any.ae4"),
  var1ci = c("% ever achieving CR", "any.cr"),
  var1cii = c("% ever achieving PR", "any.pr"),
  var1ciii = c("% ever achieving response (CR or PR)", "any.crpr"),
  var1civ = c("% ever in remission", "any.remission"),
  var1d = c("% deceased prior to end of modeled time horizon", "any.death"),
  var2ai = c("No. serious AEs experienced", "n.sae"),
  var2aii = c("No. placeholder AE 1s experienced", "n.ae1"),
  var2aiii = c("No. placeholder AE 2s experienced", "n.ae2"),
  var2aiv = c("No. placeholder AE 3s experienced", "n.ae3"),
  var2av = c("No. placeholder AE 4s experienced", "n.ae4"),
  var2b = c("No. lines of therapy initiated", "n.tx"),
  var2ci = c("No. CRs", "n.cr"),
  var2cii = c("No. PRs", "n.pr"),
  var2ciii = c("No. responses (CR or PR)", "n.crpr"),
  var2civ = c("No. remissions", "n.remission"),
  var2di = c("No. relapses", "n.relapse"),
  var2dii = c("No. relapses given >=1 relapse", "n.relapse.cond"),
  var2diii = c("No. treatment failures", "n.fail"),
  var3ai = c("Months to first CR given any CR", "mosTo.cr.cond"),
  var3aii = c("Months to first response (CR or PR) given any response", "mosTo.crpr.cond"),
  var3b = c("Months to first relapse given >=1 relapse", "mosTo.relapse.cond"),
  var4a = c("Months on treatment", "mos.onTx"),
  var4bi = c("Months in NR", "mos.nr"),
  var4bii = c("Months in CR", "mos.cr"),
  var4biii = c("Months in PR", "mos.pr"),
  var4biv = c("Months in response (CR or PR)", "mos.crpr"),
  var4bv = c("Months in response (CR or PR) given any response", "mos.crpr.cond"),
  var4bvi = c("Months in remission", "mos.remission"),
  var4bvii = c("Months in remission given any remission", "mos.remission.cond"),
  var5a = c("Total life years, discounted", "lys.disc"),
  var5b = c("Total QALYs, discounted", "qaly.disc"),
  var6a = c("Total costs incurred, discounted", "cost.disc"),
  var6bi = c("Total direct costs incurred, discounted", "cost.d.disc"),
  var6bii = c("Direct treatment costs incurred, discounted", "cost.d.tx.disc"),
  var6biii = c("Other direct disease-related costs incurred, discounted", "cost.d.state.disc"),  
  var6ci = c("Total indirect costs incurred, discounted", "cost.i.disc"),
  var6cii = c("Indirect treatment-related transportation costs incurred, discounted", "cost.i.trnsprt.disc"),
  var6ciii = c("Total indirect work-loss costs incurred, discounted", "cost.i.wl.disc"),
  var6civ = c("Indirect absenteeism costs incurred, discounted", "cost.i.abs.disc"),
  var6cv = c("Indirect presenteeism costs incurred, discounted", "cost.i.pres.disc"),
  var7a = c("Cost per QALY, discounted", "ace.disc"),
  var7b = c("Direct cost per QALY, discounted", "ace.d.disc"),
  var7c = c("Direct treatment cost per QALY, discounted", "ace.d.tx.disc"),
  var8a = c("Net monetary benefit, discounted", "nmb.disc"),
  var8b = c("Net monetary benefit, direct costs only, discounted", "nmb.d.disc"),
  var8c = c("Net monetary benefit, direct treatment costs only, discounted", "nmb.d.tx.disc"),
  varn = c("Measure", "Measure.abbrev"))

summMeas <- transpose(summMeas.tmp, make.names = "varn")


# Import ------------------------------------------------------------------

# Singular JSON file import

json_list <- fromJSON(file = file.path(proj_path,file_location, file_name))

# Extract relevant inputs from JSON file

# Base values
in.base.raw <- unlist(
  c(
    json_list[["inputs_general"]],
    json_list[["inputs_efficacy"]],
    json_list[["inputs_gaps"]],
    json_list[["inputs_costs"]],
    json_list[["inputs_utilities"]],
    json_list[["days_per_yr"]])
  )

# Min, max, and DSA category values to facilitate DSA
f.extractForDSA <- function(suff = ""){
  
  extractedStuff <- unlist(
    c(
      json_list[[paste0("inputs_general", suff)]],
      json_list[[paste0("inputs_efficacy", suff)]],
      json_list[[paste0("inputs_gaps", suff)]],
      json_list[[paste0("inputs_costs", suff)]],
      json_list[[paste0("inputs_utilities", suff)]])
  )
  return(extractedStuff)
  
}
in.min.raw <- f.extractForDSA(suff = "_min")
in.max.raw <- f.extractForDSA(suff = "_max")
in.DSACat.raw <- f.extractForDSA(suff = "_DSACat")
in.lbl.raw <- f.extractForDSA(suff = "_lbl")

inputs_mortality <- as.data.table(json_list[["inputs_mortality"]])

# Prep inputs NOT specific to modeled tx pathway --------------------------

# Convert relapse probabilities and probabilities of spontaneous response to rates
in.base <- f.process_relapse_inputs(in.base.raw)
in.base <- f.process_spontaneous_response(in.base)

# Calculate mortality rate by age and sex
in.mort <- f.process_mortality_inputs(inputs_mortality, in.base[['days_per_yr']]) 

# Designate tx pathways to simulate--------------------------

txPaths <- json_list[["txPaths"]]
ln4.addonsom.YNs <- json_list[["ln4.addonsom.YNs"]] 


# Function to simulate specified tx pathways ------------------------------

# rmIn.txPaths: list of treatment pathways to model; up to 6 pathways permitted
#   per run; each pathway is itself a list of treatments
# rmIn.ln4.addonsom.TNs: list of yeses and nos (one value per pathway) regarding
#   inclusion of add-on somatic therapy in 4th line for the given pathways;
#   may be "Yes" or "No" for each pathway
# rmIn.txMap: data.table mapping treatment names to their aliases
# rmIn.inputs: named list of model input parameter values
# rmIn.mort: data.table of mortality rates by age and sex
# rmIn.summMeas: data.table mapping names of summary measures to their aliases
# rmIn.export.XLS: content to export to Excel; may be "All", "Summary", or "None" (default);
#   if "Summary", simulation and patient-level data will NOT be exported
# rmIn.resDest: file path for exported results
# rmIn.resLabel: label appended to name of exported results Excel file (note:
#   the name of each exported result file includes a time stamp, so an additional
#   label is not necessary to distinguish different model runs)
f.runModel <- function(rmIn.txPaths,
                       rmIn.ln4.addonsom.YNs,
                       rmIn.txMap= map.tx_name.tx_abbrev,
                       rmIn.inputs= in.base,
                       rmIn.mort= in.mort,
                       rmIn.summMeas= summMeas,
                       rmIn.export.JSON= "Yes",
                       rmIn.export.XLS= "None",
                       rmIn.resDest= results_destination,
                       rmIn.resLabel= ""){
  
  # Number of treatment pathways modeled
  n.txPaths <- length(rmIn.txPaths)
  
  # Instantiate list of results for all pathways and create a copy of
  #   the data table containing the names of the summary measures
  rmOut.results <- list()
  resSummAllPaths <- copy(rmIn.summMeas)
  
  # Create data table to hold pathways (for eventual exporting with results)
  dt.txPaths <- data.table(txLine = c("Line 1",
                                      "Line 2",
                                      "Line 3",
                                      "Line 4",
                                      "Include add-on somatic therapy in 4th line?"))
  
  # If electing to export results to Excel, go ahead and create workbook object and
  #   add sheets for treatment pathways and summary results. Also create
  #   styles for formatting exported results.
  # if(rmIn.export.XLS != "None" & rmIn.export.JSON != "Yes"){
  if(rmIn.export.XLS != "None"){
    wb <- createWorkbook()
    sht.path <- addWorksheet(wb, "pathways")
    sht.summ <- addWorksheet(wb, "summStats")
    sht.ICER <- addWorksheet(wb, "ICERs")
    
    # Create styles to format results workbook
    style.bold <- createStyle(textDecoration = "bold")
    style.pct1 <- createStyle(numFmt = "0.0%")
    style.cmm0 <- createStyle(numFmt = "#,##0")
    style.cmm2 <- createStyle(numFmt = "#,##0.00")
    style.hib1 <- createStyle(fgFill = "#DDEBF7")
    style.hib2 <- createStyle(fgFill = "#BDD7EE")
    style.hiy1 <- createStyle(fgFill = "#FFF2CC")
    style.hiy2 <- createStyle(fgFill = "#FFE699")
    style.hig1 <- createStyle(fgFill = "#E2EFDA")
    style.hig2 <- createStyle(fgFill = "#C6E0B4")
  }
  
  # Loop through treatment pathways
  for(p in 1:n.txPaths){
    
    # Add pathway to data table of pathways
    dt.txPaths[, paste0("txPath", p) := c(rmIn.txPaths[[p]], rmIn.ln4.addonsom.YNs[[p]])]
    
    # Calculate efficacy rates by line of therapy for a given treatment path
    in.eff <- f.process_efficacy_inputs(txPath = rmIn.txPaths[[p]],
                                        ln4.addonsom.YN = rmIn.ln4.addonsom.YNs[[p]],
                                        txMap = rmIn.txMap,
                                        inputs = rmIn.inputs)
    
    # Calculate treatment costs (direct + transportation-related) and likelihood
    #   of AEs for treated patients by line of therapy for each treatment path
    in.othTx <- f.process_other_tx_inputs(txPath = rmIn.txPaths[[p]],
                                          ln4.addonsom.YN = rmIn.ln4.addonsom.YNs[[p]],
                                          txMap = rmIn.txMap,
                                          inputs = rmIn.inputs)
    
    # Execute simulation
    simRes <- f.simulate_txPath(in.eff, rmIn.mort, rmIn.inputs)
    print(paste0("Pathway ", p, ": ", nrow(simRes)))
    
    # Summarize results
    rmOut.results[[p]] <- f.summarize(p, in.othTx, rmIn.inputs, simRes)
    
    # Separate simulation data, pt-level data, and overall summary data
    res.sim <- rmOut.results[[p]][[1]]
    res.pt <- rmOut.results[[p]][[2]]
    res.summ <- rmOut.results[[p]][[3]]
    
    # Combine summary results for all modeled pathways
    resSummAllPaths <- resSummAllPaths[res.summ, on = c("Measure.abbrev")]
    
    # Extract total direct costs and QALYs (discounted)
    res.cost.d.disc <- as.vector(as.matrix(resSummAllPaths[Measure.abbrev == "cost.d.disc", !c(1:2)]))
    res.qaly.disc <- as.vector(as.matrix(resSummAllPaths[Measure.abbrev == "qaly.disc", !c(1:2)]))
    txPathIDs <- paste0("txPath", 1:length(res.cost.d.disc))
    
    # Calculate ICERs based on discounted direct costs and discounted QALYs
    ## In export, "Status" indicates non-dominated (ND), extended dominated (ED), or dominated (D) 
    res.ICER <- calculate_icers(cost=res.cost.d.disc, 
                                effect=res.qaly.disc, 
                                strategies=txPathIDs)

    # If exporting simulation and patient-level data...
    if(rmIn.export.XLS == "All"){
      
      # Create two worksheets for exporting simulation data and patient-level data
      sht.sim <- addWorksheet(wb, paste0("sim", p))
      sht.pt <- addWorksheet(wb, paste0("pt", p))
      
      # Add data to the worksheets just created
      writeData(wb, sht.sim, res.sim, colNames = TRUE)
      writeData(wb, sht.pt, res.pt, colNames=TRUE)
      
      # Add styles to the desired cells
      addStyle(wb, sht.sim, style=style.bold, cols=1:100, row=1, gridExpand=TRUE, stack=TRUE)
      setColWidths(wb, sht.sim, cols=1:100, widths="auto")
      freezePane(wb, sht.sim, firstActiveRow = 2, firstActiveCol = 2)
      
      addStyle(wb, sht.pt, style=style.bold, cols=1:100, row=1, gridExpand=TRUE, stack=TRUE)
      setColWidths(wb, sht.pt, cols=1:100, widths="auto")
      freezePane(wb, sht.pt, firstActiveRow = 2, firstActiveCol = 2)
      
    }
    
  }
  
  # If exporting to Excel...
  # if(rmIn.export.XLS != "None" & rmIn.export.JSON != "Yes"){
  if(rmIn.export.XLS != "None"){
    
    # Add data to the worksheets created earlier
    writeData(wb, sht.path, dt.txPaths)
    writeData(wb, sht.summ, resSummAllPaths)
    
    # Map the "Status" abbreviations directly
    res.ICER_mapped <- as.data.table(res.ICER) %>%
      mutate(Status = case_when(
        Status == "D" ~ "Dominated",
        Status == "ED" ~ "Extended Dominated",
        Status == "ND" ~ "Non-dominated",
        TRUE ~ ""
      ))
    writeData(wb, sht.ICER, res.ICER_mapped)
    
    # Add styles to the desired cells
    addStyle(wb, sht.path, style=style.bold, cols=1:(n.txPaths+1), row=1, gridExpand=TRUE, stack=TRUE)
    setColWidths(wb, sht.path, cols=1:(n.txPaths+1), widths="auto")
    freezePane(wb, sht.path, firstActiveRow = 2, firstActiveCol = 2)
    
    addStyle(wb, sht.summ, style=style.bold, cols=1:(n.txPaths+2), row=1, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.pct1, cols=3:(n.txPaths+2), rows=2:12, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.cmm2, cols=3:(n.txPaths+2), rows=13:38, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.cmm0, cols=3:(n.txPaths+2), rows=39:53, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.hib1, cols=1:(n.txPaths+2), rows=c(40,43), gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.hib2, cols=1:(n.txPaths+2), rows=39, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.hiy1, cols=1:(n.txPaths+2), rows=37, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.hiy2, cols=1:(n.txPaths+2), rows=38, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.hig1, cols=1:(n.txPaths+2), rows=52, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.summ, style=style.hig2, cols=1:(n.txPaths+2), rows=51, gridExpand=TRUE, stack=TRUE)
    setColWidths(wb, sht.summ, cols=1:(n.txPaths+2), widths="auto")
    freezePane(wb, sht.summ, firstActiveRow = 2, firstActiveCol = 2)
    
    addStyle(wb, sht.ICER, style=style.bold, cols=1:7, row=1, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.ICER, style=style.cmm0, cols=c(2, 4, 6), rows=2:(n.txPaths+1), gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.ICER, style=style.cmm2, cols=c(3, 5), rows=2:(n.txPaths+1), gridExpand=TRUE, stack=TRUE)
    setColWidths(wb, sht.ICER, cols=1:7, widths="auto")
    
    # Save
    outFl <- file.path(rmIn.resDest, paste0("IVI_MDD_model_results_", rmIn.resLabel, format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".xlsx"))
    saveWorkbook(wb, file = outFl, overwrite = TRUE)
    
  }
  
  # Add compiled table of summary results to final results list
  # rmOut.results[[n.txPaths + 1]] <- list(dt.txPaths, resSummAllPaths, res.ICER)
  rmOut.results[[n.txPaths + 1]] <- dt.txPaths
  rmOut.results[[n.txPaths + 2]] <- resSummAllPaths
  # Map the "Status" abbreviations directly
  res.ICER_mapped <- as.data.table(res.ICER) %>%
    mutate(Status = case_when(
      Status == "D" ~ "Dominated",
      Status == "ED" ~ "Extended Dominated",
      Status == "ND" ~ "Non-dominated",
      TRUE ~ ""
    ))
  rmOut.results[[n.txPaths + 3]] <- res.ICER_mapped
  
  # If exporting to JSON...
  
  if(rmIn.export.JSON == "Yes"){
    
    json_export <- list(
      "txPaths" = rmOut.results[[n.txPaths + 1]],
      "resSummAllPaths" = rmOut.results[[n.txPaths + 2]],
      "resICER" = rmOut.results[[n.txPaths + 3]]
    )
  
    write(toJSON(json_export, indent = 1), file.path(rmIn.resDest, paste0("IVI_MDD_model_results_", rmIn.resLabel, format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".json")))
  }
  
  # Return wad of results
  return(rmOut.results)
  
}

# Model run ----------------------------

hereGoesNothin <- f.runModel(rmIn.txPaths= txPaths,
                             rmIn.ln4.addonsom.YNs= ln4.addonsom.YNs,
                             rmIn.export.JSON= "Yes",
                             rmIn.export.XLS= "Summary",
                             # rmIn.resLabel= "TEST_"
                             rmIn.resLabel= ""
                             )


# DSA run ------------------------------

res.DSA <- f.runDSA(rDSAIn.txPath = txPaths[[in.base["dsa_txpath_n"]]]
                            ,rDSAIn.ln4.addonsom.YN = ln4.addonsom.YNs[[in.base["dsa_txpath_n"]]]
                            ,rDSAIn.txMap = map.tx_name.tx_abbrev
                            ,rDSAIn.rawIn.base= in.base.raw
                            ,rDSAIn.rawIn.min= in.min.raw
                            ,rDSAIn.rawIn.max= in.max.raw
                            ,rDSAIn.rawIn.DSACat= in.DSACat.raw
                            ,rDSAIn.rawIn.lbl= in.lbl.raw
                            ,rDSAIn.SummMeas = summMeas.DSA
                            ,rDSAIn.in.mort= in.mort
                            ,rDSAin.resDest = results_destination
                            ,rmDSA.export.XLS="Yes"
                            ,rmDSA.export.JSON="Yes")
