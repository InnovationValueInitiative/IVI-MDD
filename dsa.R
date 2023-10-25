# Summary measures ----------

summMeas.tmp.DSA <- data.table(
  var1a = c("Input value, low", "in.low"),
  var1b = c("Input value, base", "in.base"),
  var1c = c("Input value, high", "in.high"),
  var2a = c("Total life years, discounted, low input", "lys.disc.low"),
  var2b = c("Total life years, discounted, base input", "lys.disc.base"),
  var2c = c("Total life years, discounted, high input", "lys.disc.high"),
  var2d = c("Total life years, discounted, range", "lys.disc.rng"),
  var2e = c("Total life years, discounted, low minus base result", "lys.disc.lowLessBase"),
  var2f = c("Total life years, discounted, high minus base result", "lys.disc.highLessBase"),
  
  var3a = c("Total QALYs, discounted, low input", "qaly.disc.low"),
  var3b = c("Total QALYs, discounted, base input", "qaly.disc.base"),
  var3c = c("Total QALYs, discounted, high input", "qaly.disc.high"),
  var3d = c("Total QALYs, discounted, range", "qaly.disc.rng"),
  var3e = c("Total QALYs, discounted, low minus base result", "qaly.disc.lowLessBase"),
  var3f = c("Total QALYs, discounted, high minus base result", "qaly.disc.highLessBase"),

  var4a = c("Total costs incurred, discounted, low input", "cost.disc.low"),
  var4b = c("Total costs incurred, discounted, base input", "cost.disc.base"),
  var4c = c("Total costs incurred, discounted, high input", "cost.disc.high"),
  var4d = c("Total costs incurred, discounted, range", "cost.disc.rng"),
  var4e = c("Total costs incurred, discounted, low minus base result", "cost.disc.lowLessBase"),
  var4f = c("Total costs incurred, discounted, high minus base result", "cost.disc.highLessBase"),
  
  var5a = c("Net monetary benefit, discounted, low input", "nmb.disc.low"),
  var5b = c("Net monetary benefit, discounted, base input", "nmb.disc.base"),
  var5c = c("Net monetary benefit, discounted, high input", "nmb.disc.high"),
  var5d = c("Net monetary benefit, discounted, range", "nmb.disc.rng"),
  var5e = c("Net monetary benefit, discounted, low minus base result", "nmb.disc.lowLessBase"),
  var5f = c("Net monetary benefit, discounted, high minus base result", "nmb.disc.highLessBase"),
  varn = c("Measure", "Measure.abbrev"))

summMeas.DSA <- transpose(summMeas.tmp.DSA, make.names = "varn")

# Function to run the model once for a given set of input values ----------

f.runDSA.inner <- function(rDSAinIn.txPath
                           ,rDSAinIn.ln4.addonsom.YN
                           ,rDSAinIn.txMap
                           ,rDSAinIn.rawIn.main
                           ,rDSAinIn.in.mort
                           ,rDSAinIn.simData=NULL
                           ,rDSAinIn.returnSimData=FALSE){
  
  in.list <- rDSAinIn.rawIn.main
  
  # Convert relapse probabilities and probabilities of spontaneous response to rates
  in.list <- f.process_relapse_inputs(in.list)
  in.list <- f.process_spontaneous_response(in.list)
  
  # Calculate efficacy rates by line of therapy for a given treatment path
  in.eff <- f.process_efficacy_inputs(txPath = rDSAinIn.txPath,
                                      ln4.addonsom.YN = rDSAinIn.ln4.addonsom.YN,
                                      txMap = rDSAinIn.txMap,
                                      inputs = in.list)
  
  # Calculate treatment costs (direct + transportation-related) and likelihood
  #   of AEs for treated patients by line of therapy for a given treatment path
  in.othTx <- f.process_other_tx_inputs(txPath = rDSAinIn.txPath,
                                        ln4.addonsom.YN = rDSAinIn.ln4.addonsom.YN,
                                        txMap = rDSAinIn.txMap,
                                        inputs = in.list)
  
  # Execute simulation, if necessary
  if( is.null(rDSAinIn.simData) == TRUE ){
    rDSAinIn.simData <- f.simulate_txPath(in.eff, rDSAinIn.in.mort, in.list)
  }
  
  # Calculate summary measures
  rDSAinOut.summ <- f.summarize(pw=1,
                              inputs.dt.tx=in.othTx,
                              inputs.oth=in.list,
                              simRes=rDSAinIn.simData,
                              summIn.forDSA=TRUE)
  
  if(rDSAinIn.returnSimData == TRUE){
    return( list(rDSAinOut.summ, rDSAinIn.simData) )
  }
  
  return(rDSAinOut.summ)
  
}



f.runDSA <- function(rDSAIn.txPath
                     ,rDSAIn.ln4.addonsom.YN
                     ,rDSAIn.txMap
                     ,rDSAIn.rawIn.base
                     ,rDSAIn.rawIn.min
                     ,rDSAIn.rawIn.max
                     ,rDSAIn.rawIn.DSACat
                     ,rDSAIn.in.mort
                     ,rDSAIn.varPct= 0.2
                     ,rDSAIn.rawIn.lbl
                     ,rDSAIn.SummMeas
                     ,rDSAin.resDest
                     ,rmDSA.export.XLS="No"
                     ,rmDSA.export.JSON="Yes"){
  
  # Function for adding a suffix to the name of a column on a data.table
  f.rnm <- function(currName, suff, dt2rn){
    setnames(dt2rn, currName, paste0(currName, suff))
    }
  
  # Create table shell for eventual pathway-specific results
  res.DSA <- data.table(Parameter=character(),
                        # Name=character(),
                        in.low=double(),
                        in.base=double(),
                        in.high=double(),
                        # Results measures
                        lys.disc.low=double(),
                        lys.disc.base=double(),
                        lys.disc.high=double(),
                        qaly.disc.low=double(),
                        qaly.disc.base=double(),
                        qaly.disc.high=double(),
                        cost.disc.low=double(),
                        cost.disc.base=double(),
                        cost.disc.high=double(),
                        nmb.disc.low=double(),
                        nmb.disc.base=double(),
                        nmb.disc.high=double()
  )
  
  # Create a copy of the inputs list for modifying as we go
  in.DSA.raw <- copy(rDSAIn.rawIn.base)
  
  # Begin by running the model once without changing any input values
  res.base <- f.runDSA.inner(rDSAinIn.txPath= rDSAIn.txPath
                             ,rDSAinIn.ln4.addonsom.YN= rDSAIn.ln4.addonsom.YN
                             ,rDSAinIn.txMap= rDSAIn.txMap
                             ,rDSAinIn.rawIn.main= in.DSA.raw
                             ,rDSAinIn.in.mort= rDSAIn.in.mort
                             ,rDSAinIn.returnSimData=TRUE)
  
  # Extract summary data, save column names, then update column
  #   names to indicate that these results were generated using
  #   the original input values
  res.summ.base <- res.base[[1]]
  outMeas <- names(copy(res.summ.base))
  sapply(names(res.summ.base), f.rnm, suff=".base", dt2rn=res.summ.base)

  # Extract simulation data for use in expediting analysis of
  #   input parameters that do not affect patient trajectories
  res.sim.base <- res.base[[2]]

  # Identify parameters for inclusion in the DSA
  inclTF <- sapply(rDSAIn.rawIn.DSACat, function(x) x %in% c(1, 2) )
  DSAVars <- names(rDSAIn.rawIn.DSACat[inclTF])
  
  # Loop through list of inputs varied in the DSA
  for (i in 1:length(DSAVars)){
  # for (i in 1:3){
    print(i)
    
    # Pull the name, base, min, and maximum values, and the DSA category
    # for the current input parameter
    v.alias <- DSAVars[i]
    v.dt <- data.table("Parameter" = v.alias)
    v.base <- rDSAIn.rawIn.base[[v.alias]]
    v.min <- rDSAIn.rawIn.min[[v.alias]]
    v.max <- rDSAIn.rawIn.max[[v.alias]]
    v.cat <- rDSAIn.rawIn.DSACat[[v.alias]]
    
    # Determine whether this parameter requires a re-run of the simulation
    if( near(v.cat, 1) == TRUE ){
      v.simRes <- res.sim.base
    } else {
      v.simRes <- NULL
    }
    
    # LOW
    
    # Calculate low value of input parameter and set to this value
    #   within the DSA inputs list
    v.low <- min( max( v.base*(1 - rDSAIn.varPct), v.min), v.max)
    in.DSA.raw[[v.alias]] <- v.low
    
    # Run the simulation (if necessary) and calculate select summary
    #   statistics with the current parameter set equal to its low value
    v.summ.low <- f.runDSA.inner(rDSAinIn.txPath= rDSAIn.txPath
                                 ,rDSAinIn.ln4.addonsom.YN= rDSAIn.ln4.addonsom.YN
                                 ,rDSAinIn.txMap= rDSAIn.txMap
                                 ,rDSAinIn.rawIn.main= in.DSA.raw
                                 ,rDSAinIn.in.mort= rDSAIn.in.mort
                                 ,rDSAinIn.simData= v.simRes)
    
    
    # HIGH
    
    # Calculate high value of input parameter and set to this value
    #   within the DSA inputs list
    v.high <- min( max( v.base*(1 + rDSAIn.varPct), v.min), v.max)
    in.DSA.raw[[v.alias]] <- v.high
    
    # Run the simulation (if necessary) and calculate select summary
    #   statistics with the current parameter set equal to its high value
    v.summ.high <- f.runDSA.inner(rDSAinIn.txPath= rDSAIn.txPath
                                 ,rDSAinIn.ln4.addonsom.YN= rDSAIn.ln4.addonsom.YN
                                 ,rDSAinIn.txMap= rDSAIn.txMap
                                 ,rDSAinIn.rawIn.main= in.DSA.raw
                                 ,rDSAinIn.in.mort= rDSAIn.in.mort
                                 ,rDSAinIn.simData= v.simRes)
    
    # Revert parameter value back to its original value in input list
    in.DSA.raw[[v.alias]] <- v.base
    
    # Compile results and append to full results table
    sapply(names(v.summ.low), f.rnm, suff=".low", dt2rn=v.summ.low)
    sapply(names(v.summ.high), f.rnm, suff=".high", dt2rn=v.summ.high)
    
    dt1 <- cbind(data.table(Parameter = v.alias,
                            in.low = v.low,
                            in.base = v.base,
                            in.high = v.high),
                 v.summ.low,
                 res.summ.base,
                 v.summ.high)

    res.DSA <- rbindlist(list(res.DSA, dt1), use.names=TRUE)
    
  }
  
  # Calculate range of results and differences between DSA results and base result
  #   for each variable and outcome measure
  f.calcDiffs <- function(meas, dt2mnplt){
    
    vl <- dt2mnplt[][[paste0(meas,".low")]]
    vb <- dt2mnplt[][[paste0(meas,".base")]]
    vh <- dt2mnplt[][[paste0(meas,".high")]]
    
    dt2mnplt[, paste0(meas,".rng") := pmax(vl, vb, vh) - pmin(vl, vb, vh)]
    dt2mnplt[, paste0(meas,".lowLessBase") := vl - vb]
    dt2mnplt[, paste0(meas,".highLessBase") := vh - vb]
  }
  invisible(sapply(outMeas, f.calcDiffs, dt2mnplt=res.DSA))
  
  # Sort
  # setorderv(res.DSA, cols=paste0(outMeas, ".rng"), order=-1, na.last=TRUE)
  setorderv(res.DSA, cols="nmb.disc.rng", order=-1, na.last=TRUE)
  
  # Prepare exports
  
  ## Set a parameter label column
  
  parameters.labels <- as.data.table(rDSAIn.rawIn.lbl, keep.rownames = "id") %>% 
    rename(Parameter.label = rDSAIn.rawIn.lbl)
  
  res.DSA_w_labels <- left_join(res.DSA, parameters.labels, by = c("Parameter" = "id")) %>%
    relocate(Parameter.label)
  
  ## Excel
  
  if(rmDSA.export.XLS != "None"){
    
    ### Prepare Excel
    
    wb <- createWorkbook()
    sht.map <- addWorksheet(wb, "Measure.definitions") # mapping of measures to definitions
    sht.DSA <- addWorksheet(wb, "DSA") # results
    
    ### Create styles to format results workbook
    
    style.bold <- createStyle(textDecoration = "bold")
    style.cmm0 <- createStyle(numFmt = "#,##0")
    style.cmm2 <- createStyle(numFmt = "#,##0.00")
    style.cmm4 <- createStyle(numFmt = "#,##0.0000")
    
    ### Write to Excel
    
    writeData(wb, sht.map, rDSAIn.SummMeas)
    writeData(wb, sht.DSA, res.DSA_w_labels)
    
    ### Add styles to the desired cells
    addStyle(wb, sht.map, style=style.bold, cols = 1:2, rows = 1,gridExpand=TRUE, stack=TRUE)
    setColWidths(wb, sht.map, cols=1:2, widths="auto")
    freezePane(wb, sht.map, firstRow = TRUE)
    
    addStyle(wb, sht.DSA, style=style.bold, cols = 1:29, rows=1, gridExpand=TRUE, stack=TRUE)
    addStyle(wb, sht.DSA, style=style.cmm2, cols = 3:29, rows=2:39, gridExpand = TRUE, stack = TRUE)
    setColWidths(wb, sht.DSA, cols=1:29, widths="auto")
    freezePane(wb, sht.DSA, firstActiveRow = 2, firstActiveCol = 3)
    
    ### Save
    outFl.DSA <- file.path(rDSAin.resDest, paste0("IVI_MDD_model_results_DSA_", in.base["dsa_txpath_n"], "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".xlsx"))
    saveWorkbook(wb, file = outFl.DSA)
  }
  
  ## JSON
  if(rmDSA.export.JSON != "None"){
    write(toJSON(res.DSA_w_labels, indent = 1), file.path(rDSAin.resDest, paste0("IVI_MDD_model_results_DSA_", in.base["dsa_txpath_n"], "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".json")))
    
  }
    
  return(res.DSA_w_labels)
  
}






