#[.N:.N]

# Mortality rates by age and sex
f.process_mortality_inputs <- function(deathProbs, days_per_yr){
  
  # Calculate mortality rates based on annual probabilities
  # Exponential cdf:    F(t) = 1 - exp(-lambda*t)
  #                 --> lambda = log(1 - F(t)) / -t
  deathRates <- deathProbs[, rate_death := -log(1 - p_death) / days_per_yr] %>%
    select(-c("p_death"))
  
  return(deathRates)
  
}

# Relapse rates by response status and line of therapy (treated patients)
#   and response status (untreated patients)
f.process_relapse_inputs <- function(inputs){
  
  # dpy <- inputs[Variable == 'days_per_yr', Base]
  dpy <- inputs[['days_per_yr']]
  
  # Create empty data table to hold line-of-therapy-specific
  #   relapse rates
  # in.out <- data.table("Variable" = character(), "Base" = double())
  
  # Create variables for storing running maximums of all relapse rates
  #   for treated patients in each state
  maxRate.cr <- 0
  maxRate.pr <- 0
  
  # Pull relapse HRs comparing early v late responders
  # hr.crEarly <- inputs[Variable == "hr_relapse_early_vs_late_cr", Base]
  # hr.prEarly <- inputs[Variable == "hr_relapse_early_vs_late_pr", Base]
  hr.crEarly <- inputs[["hr_relapse_early_vs_late_cr"]]
  hr.prEarly <- inputs[["hr_relapse_early_vs_late_pr"]]
  
  # Pull relapse HR comparing untreated v treated responders
  # hr.notreat <- inputs[Variable == "hr_relapse_notreat_vs_treat", Base]
  hr.notreat <- inputs[["hr_relapse_notreat_vs_treat"]]
  
  # Iterate through lines of therapy
  for(ln in 1:4){
  
    # Pull annual probabilities of relapse among late responders
    # p.crLate.nr <- inputs[Variable == paste0("p_relapse_late_cr_line", ln), Base]
    # p.prLate.nr <- inputs[Variable == paste0("p_relapse_late_pr_line", ln), Base]
    p.crLate.nr <- inputs[[paste0("p_relapse_late_cr_line", ln)]]
    p.prLate.nr <- inputs[[paste0("p_relapse_late_pr_line", ln)]]
    
    # Convert to rates
    # Exponential cdf:    F(t) = 1 - exp(-lambda*t)
    #                 --> lambda = log(1 - F(t)) / -t
    rate.crLate.nr <- -log(1 - p.crLate.nr) / dpy
    rate.prLate.nr <- -log(1 - p.prLate.nr) / dpy
    
    # Add to in.out data.table
    # varsToAdd <- c(paste0("rate_relapse_late_cr_line", ln),
    #                paste0("rate_relapse_late_pr_line", ln))
    # valsToAdd <- c(rate.crLate.nr,
    #                rate.prLate.nr)
    # in.out <- rbindlist(list(in.out,
    #                          data.table("Variable" = varsToAdd,
    #                                     "Base" = valsToAdd)
    # )
    # )
    inputs[[paste0("rate_relapse_late_cr_line", ln)]] <- rate.crLate.nr
    inputs[[paste0("rate_relapse_late_pr_line", ln)]] <- rate.prLate.nr

    # Update running maximum rates
    maxRate.cr <- max(c(maxRate.cr, rate.crLate.nr, rate.crLate.nr * hr.crEarly))
    maxRate.pr <- max(c(maxRate.pr, rate.prLate.nr, rate.prLate.nr * hr.prEarly))

  }
  
  # Calculate relapse rates for untreated patients based on overall
  #   maximum rates for treated patients
  rate.cr.nr.notreat <- maxRate.cr * hr.notreat
  rate.pr.nr.notreat <- maxRate.pr * hr.notreat

  # Add rates (including those for untreated pts) to existing inputs table
  # varsToAdd <- c("rate_relapse_cr_notreat", "rate_relapse_pr_notreat")
  # valsToAdd <- c(rate.cr.nr.notreat, rate.pr.nr.notreat)
  # in.out <- rbindlist(list(inputs,
  #                          in.out,
  #                          data.table("Variable" = varsToAdd,
  #                                     "Base" = valsToAdd)))
  inputs[["rate_relapse_cr_notreat"]] <- rate.cr.nr.notreat
  inputs[["rate_relapse_pr_notreat"]] <- rate.pr.nr.notreat

  # return(in.out)
  return(inputs)
  
}

# Rates of spontaneous response
f.process_spontaneous_response <- function(inputs){
  
  # Days per year
  # dpy <- inputs[Variable == "days_per_yr", Base]
  dpy <- inputs[['days_per_yr']]
  
  # Pull 56-day probabilities of spontaneous response
  # p.cr <- inputs[Variable == "p_cr_notreat", Base]
  # p.pr <- inputs[Variable == "p_pr_notreat", Base]
  p.cr <- inputs[["p_cr_notreat"]]
  p.pr <- inputs[["p_pr_notreat"]]
  
  # Make sure each probability is <=100% and that they sum to <=100%
  p.cr <- min(c(1, p.cr))
  p.pr <- min(c(1-p.cr, p.pr))
  
  # Convert to rates
  # Exponential cdf:    F(t) = 1 - exp(-lambda*t)
  #                 --> lambda = log(1 - F(t)) / -t
  rate.cr <- -log(1 - p.cr) / dpy
  rate.pr <- -log(1 - p.pr) / dpy
  
  # Add to inputs in in.out data.table
  # varsToAdd <- c("rate_cr_notreat", "rate_pr_notreat")
  # valsToAdd <- c(rate.cr, rate.pr)
  # in.out <- rbindlist(list(inputs,
  #                          data.table("Variable" = varsToAdd,
  #                                     "Base" = valsToAdd)
  # )
  # )
  inputs[["rate_cr_notreat"]] <- rate.cr
  inputs[["rate_pr_notreat"]] <- rate.pr

  # return(in.out)
  return(inputs)
}


# Efficacy rates by line of therapy
f.process_efficacy_inputs <- function(txPath,
                                      ln4.addonsom.YN="No",
                                      txMap=map.tx_name.tx_abbrev,
                                      inputs){
  
  # Create empty data table to hold line-of-therapy-specific
  #   efficacy inputs
  # in.out <- data.table("Variable" = character(), "Base" = double())
  in.out <- setNames(numeric(0), character(0))
  
  # Temporarily assume final tx line is 4th line
  # max.ln <- nrow(txPath)
  max.ln <- length(txPath)
  
  # Iterate through lines of therapy in pathway
  for(ln in 1:max.ln){
    
    # Treatment name
    # tx.nm <- txPath[ln,1]
    tx.nm <- txPath[ln]
    # Pull abbreviated tx name for given treatment class
    tx <- txMap[tx.name==tx.nm, tx.abbrev]
    
    # If line of therapy is no active treatment, we are finished
    #   compiling efficacy inputs by line of therapy. Update
    #   number of final (maximum) treatment line and exit loop
    if(tx == "notreat"){
      max.ln <- ln - 1
      break
    }
    
    # Obtain duration of each of initiation phase and initiation-phase extension
    # durPart1 <- inputs[Variable == paste0("phaseDays_init_", tx), Base]
    # durPart2 <- inputs[Variable == paste0("phaseDays_initExt_", tx), Base]
    durPart1 <- inputs[[paste0("phaseDays_init_", tx)]]
    durPart2 <- inputs[[paste0("phaseDays_initExt_", tx)]]
    
    # Obtain maximum likelihood of CR (initiation phase part 1)
    p.nr.cr <- inputs[[paste0("p_cr_init_", tx)]]

    # Obtain maximum likelihood of PR (initiation phase part 1)
    # p.nr.pr.mult <- inputs[Variable == "pRatio_pr_vs_cr_init_nonsomatic", Base]
    # p.nr.pr <- p.nr.cr * p.nr.pr.mult
    p.nr.pr <- p.nr.cr * inputs[["pRatio_pr_vs_cr_init_nonsomatic"]]
    
    # Account for add-on somatic therapy, if applicable (4th line only)
    if(ln==4 & ln4.addonsom.YN=="Yes"){
      # p.nr.cr <- p.nr.cr * inputs[Variable == "rr_cr_init_somatic", Base]
      # p.nr.pr <- p.nr.pr * inputs[Variable == "rr_pr_init_somatic", Base]
      p.nr.cr <- p.nr.cr * inputs[["rr_cr_init_somatic"]]
      p.nr.pr <- p.nr.pr * inputs[["rr_pr_init_somatic"]]
    }
    
    # Obtain maximum likelihood of CR given PR (initiation phase part 2)
    # p.pr.cr <- inputs[Variable == paste0("p_cr_initExt_", tx), Base]
    p.pr.cr <- inputs[[paste0("p_cr_initExt_", tx)]]
    
    # If 1st line, store tx for reference in line 2
    if(ln==1){
      tx.ln1 <- tx
      
    # Otherwise, if line 2+, apply efficacy decrements as needed.
    #   - If line 2, only apply decrements if same tx class as line 1
    #   - If line 3+, apply decrements
    } else if( ln>2 | tx.ln1 == tx ){
      # decr.cr <- inputs[Variable == paste0("rr_cr_line", ln), Base]
      # decr.pr <- inputs[Variable == paste0("rr_pr_line", ln), Base]
      # p.nr.cr <- p.nr.cr * decr.cr
      # p.pr.cr <- p.pr.cr * decr.cr
      # p.nr.pr <- p.nr.pr * decr.pr
      p.nr.cr <- p.nr.cr * inputs[[paste0("rr_cr_line", ln)]]
      p.pr.cr <- p.pr.cr * inputs[[paste0("rr_cr_line", ln)]]
      p.nr.pr <- p.nr.pr * inputs[[paste0("rr_pr_line", ln)]]
    }
    
    # Make sure all probabilities are individually <=100% and that
    #   initiation-phase probabilities sum to <=100%
    if(1 < p.nr.cr){
      p.nr.cr <- 1
    }
    if(1 < p.nr.cr + p.nr.pr){
      p.nr.pr <- 1 - p.nr.cr
    }
    if(1 < p.pr.cr){
      p.pr.cr <- 1
    }
    
    # Calculate efficacy rates based on probabilities & phase durations
    # Exponential cdf:    F(t) = 1 - exp(-lambda*t)
    #                 --> lambda = log(1 - F(t)) / -t
    rt.nr.cr <- -log(1 - p.nr.cr) / durPart1
    rt.nr.pr <- -log(1 - p.nr.pr) / durPart1
    rt.pr.cr <- -log(1 - p.pr.cr) / durPart2
    
    # Add values of interest to results table of inputs
    # varsToAdd <- c(paste0("phaseDays_init_line", ln),
    #                paste0("phaseDays_initExt_line", ln),
    #                paste0("rate_cr_init_line", ln),
    #                paste0("rate_pr_init_line", ln),
    #                paste0("rate_cr_initExt_line", ln))
    # valsToAdd <- c(durPart1,
    #                durPart2,
    #                rt.nr.cr,
    #                rt.nr.pr,
    #                rt.pr.cr)
    # in.out <- rbindlist(list(in.out,
    #                          data.table("Variable" = varsToAdd,
    #                                     "Base" = valsToAdd)
    #                          )
    #                     )
    in.out[[paste0("phaseDays_init_line", ln)]] <- durPart1
    in.out[[paste0("phaseDays_initExt_line", ln)]] <- durPart2
    in.out[[paste0("rate_cr_init_line", ln)]] <- rt.nr.cr
    in.out[[paste0("rate_pr_init_line", ln)]] <- rt.nr.pr
    in.out[[paste0("rate_cr_initExt_line", ln)]] <- rt.pr.cr
    
  }
  
  # Add maximum treatment line to results table of inputs
  # in.out <- rbindlist(list(in.out,
  #                          data.table("Variable" = "max_tx_line",
  #                                     "Base" = max.ln)
  #                          )
  #                     )
  in.out[["max_tx_line"]] <- max.ln
  
  # Drop now obsolete tx-specific inputs from inputs table (taking
  #   care not to drop those for "notreat")
  #!TODO! but for now, just leave 'em
  
  # Return compiled efficacy inputs by line of therapy
  return(in.out)
  
}

# Treatment costs and AE incidence proportions by line of therapy
f.process_other_tx_inputs <- function(txPath,
                                      ln4.addonsom.YN="No",
                                      txMap=map.tx_name.tx_abbrev,
                                      inputs){
  
  # Create data table to hold line-of-therapy-specific cost and AE-incidence
  #   inputs. Start with one row for the pre-treatment period during which
  #   patients cannot sustain treatment-related costs or develop treatment-related
  #   AEs (cost = incidence = 0)
  in.out <- data.table("tx_line" = 0,
                       "onTx" = 0,
                       "p.sae" = 0,
                       "p.ae1" = 0,
                       "p.ae2" = 0,
                       "p.ae3" = 0,
                       "p.ae4" = 0,
                       "cpy.tx" = 0,
                       "c1x.tx" = 0,
                       "cpy.tx.transport" = 0,
                       "c1x.tx.transport" = 0)
  
  # Iterate through lines of therapy in pathway
  for(ln in 1:4){
    
    # Treatment name
    # tx.nm <- txPath[ln,1]
    tx.nm <- txPath[ln]
    # Pull abbreviated tx name for given treatment class
    tx <- txMap[tx.name==tx.nm, tx.abbrev]
    # Pull binary indicator of pharmacotherapy
    tx.rx <- txMap[tx.name==tx.nm, tx.pharm]
    
    # If line of therapy is no active treatment, we are finished compiling
    #   inputs by line of therapy. Exit loop.
    if(tx == "notreat"){
      break
    }
    
    # Split treatment into its component parts
    tx.cs <- unlist(strsplit(tx, "_"))
    
    # Instantiate temporary cost and incidence variables
    p_sae <- 0
    p_ae1 <- 0
    p_ae2 <- 0
    p_ae3 <- 0
    p_ae4 <- 0
    cpy <- 0
    cpy.tp <- 0
    c1x <- 0
    c1x.tp <- 0
    
    # Update AE incidence proportions if treatment is pharmacotherapy or
    #   placeholder therapy
    if(tx.rx == 1 | tx == "phtx"){
      
      if(tx.rx == 1) { suff <- "pharm" } else { suff <- "phtx" }
      
      # p_sae <- inputs[Variable == paste0("p_sae_", suff), Base]
      # p_ae1 <- inputs[Variable == paste0("p_ae1_", suff), Base]
      # p_ae2 <- inputs[Variable == paste0("p_ae2_", suff), Base]
      # p_ae3 <- inputs[Variable == paste0("p_ae3_", suff), Base]
      # p_ae4 <- inputs[Variable == paste0("p_ae4_", suff), Base]
      p_sae <- inputs[[paste0("p_sae_", suff)]]
      p_ae1 <- inputs[[paste0("p_ae1_", suff)]]
      p_ae2 <- inputs[[paste0("p_ae2_", suff)]]
      p_ae3 <- inputs[[paste0("p_ae3_", suff)]]
      p_ae4 <- inputs[[paste0("p_ae4_", suff)]]
      
    }
    
    # Loop through component treatment classes
    for(c in 1:length(tx.cs)){
      
      # Annual direct treatment costs for all but add-on somatic therapy
      # cpy <- cpy + inputs[Variable == paste0("costPerYr_direct_", tx.cs[[c]]), Base]
      cpy <- cpy + inputs[[paste0("costPerYr_direct_", tx.cs[[c]])]]
      
      # Determine annual and one-time transportation costs for all treatments excluding
      #   add-on somatic therapy. If applicable, also determine one-time direct costs
      #   for placeholder therapy. Note: Placeholder therapy is never included
      #   in combination or augmented regimens (apart from placeholder therapy +
      #   add-on somatic therapy), so do not need to take maximums
      if(tx.cs[[c]] %in% c("ssri","snri","atypantidep","antipsych")){
        # cpy.tp <- max(c(cpy.tp,inputs[Variable == "costPerYr_transport_pharm", Base]))
        cpy.tp <- max(c(cpy.tp,inputs[["costPerYr_transport_pharm"]]))
      } else if(tx.cs[[c]] == "psycther"){
        # cpy.tp <- max(c(cpy.tp,inputs[Variable == "costPerYr_transport_psycther", Base]))
        cpy.tp <- max(c(cpy.tp,inputs[["costPerYr_transport_psycther"]]))
      } else if(tx.cs[[c]] == "phtx"){
        # cpy.tp <- inputs[Variable == "costPerYr_transport_phtx", Base]
        # c1x.tp <- inputs[Variable == "costOneX_transport_phtx", Base]
        # c1x <- inputs[Variable == "costOneX_direct_phtx", Base]
        cpy.tp <- inputs[["costPerYr_transport_phtx"]]
        c1x.tp <- inputs[["costOneX_transport_phtx"]]
        c1x <- inputs[["costOneX_direct_phtx"]]
      }
      
    }
    
    # Account for add-on somatic therapy, if applicable (4th line only)
    if(ln==4 & ln4.addonsom.YN=="Yes"){
      
      # One-time costs associated with add-on somatic therapy
      # c1x <- c1x + inputs[Variable == "costOneX_direct_somatic", Base]
      # c1x.tp <- c1x.tp + inputs[Variable == "costOneX_transport_somatic", Base]
      c1x <- c1x + inputs[["costOneX_direct_somatic"]]
      c1x.tp <- c1x.tp + inputs[["costOneX_transport_somatic"]]
      
      # Maximums of current AE incidence props and AE incidence props specific
      #   to somatic therapy
      # p_sae <- max(c(p_sae, inputs[Variable == "p_sae_somatic", Base]))
      # p_ae1 <- max(c(p_ae1, inputs[Variable == "p_ae1_somatic", Base]))
      # p_ae2 <- max(c(p_ae2, inputs[Variable == "p_ae2_somatic", Base]))
      # p_ae3 <- max(c(p_ae3, inputs[Variable == "p_ae3_somatic", Base]))
      # p_ae4 <- max(c(p_ae4, inputs[Variable == "p_ae4_somatic", Base]))
      p_sae <- max(c(p_sae, inputs[["p_sae_somatic"]]))
      p_ae1 <- max(c(p_ae1, inputs[["p_ae1_somatic"]]))
      p_ae2 <- max(c(p_ae2, inputs[["p_ae2_somatic"]]))
      p_ae3 <- max(c(p_ae3, inputs[["p_ae3_somatic"]]))
      p_ae4 <- max(c(p_ae4, inputs[["p_ae4_somatic"]]))
      
    }
    
    # Add values of interest to final inputs table of AE incidence and tx cost
    #   inputs. Include "dummy" rows with proportions and costs all equal to 0
    #   for each treatment line for patients NOT currently receiving treatment
    in.out <- rbindlist(list(in.out,
                             data.table("tx_line" = c(ln, ln),
                                        "onTx" = c(0, 1),
                                        "p.sae" = c(0, p_sae),
                                        "p.ae1" = c(0, p_ae1),
                                        "p.ae2" = c(0, p_ae2),
                                        "p.ae3" = c(0, p_ae3),
                                        "p.ae4" = c(0, p_ae4),
                                        "cpy.tx" = c(0, cpy),
                                        "c1x.tx" = c(0, c1x),
                                        "cpy.tx.transport" = c(0, cpy.tp),
                                        "c1x.tx.transport" = c(0, c1x.tp)
                             )
    )
    )
    
  }
  
  # Return AE incidence and cost inputs by line of therapy
  return(in.out)
  
}


# Function to convert inputs in data table to inputs in a list
# f.collect_inputs <- function(dt, valVar){
#   mylist <- list()
#   for (i in 1:nrow(dt)){
#     if(!is.na(dt$Variable[i])) mylist[[dt$Variable[i]]] <- as.numeric(dt[[valVar]][i])
#   }
#   return(mylist)
# }








