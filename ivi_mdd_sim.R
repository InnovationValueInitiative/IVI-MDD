
#--- HELPER FUNCTION(S) --------------------------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

#----------------------------#
#~~~ New patient creation ~~~#
#----------------------------#

# Function to create new hypothetical patient
f.create_pt <- function(id, female, age.in, SES){
  
  X <- data.table(
    # id: Unique patient identifier
    id = id,
    # female: 1 if patient is female, 0 otherwise
    female = female,
    # SES: socioeconomic status quartile (1 = lowest, 4 = highest)
    SES = SES,
    # age.in: patient age as of time.in
    age.in = age.in,
    # tx_line: current or most recent line of therapy
    tx_line = 0,
    # phase.curr: current phase (i.e., phase patient is in from time.in to
    #   time.out) (onset, init, initExt, maint, gap, relapse, fail, postTx,
    #   or death)
    phase.curr = "onset",
    # phase.next: phase patient will be in after time.out
    phase.next = "onset",
    # time.phase_st: model time at which patient entered the current phase
    time.phase_st = NA_real_,
    # state.curr: current state (i.e., state patient is in from time.in to
    #   time.out) (nr, pr, cr, or x if deceased)
    state.curr = "nr",
    # state.next: state patient will be in after time.out
    state.next = "nr",
    # time.cr_st: model time at which patient's current stint in the cr state
    #   began; equals NA when state.curr does not equal cr
    time.cr_st = NA_real_,
    # time.in: starting model time for the period currently described by the
    #   patient data.table X
    time.in = 0.0,
    # time.out: ending model time for the period currently described by the
    #   patient data.table X
    time.out = 0.0,
    # timeEl.tx_st.prcr: time elapsed from start of current treatment to time
    #   at which patient first achieved pr (cr) if the patient is currently in pr
    #   (cr); if patient is currently off treatment or in the nr or x states,
    #   then timeEl.tx_st.prcr is NA; updated when patient transitions from pr
    #   to cr while on treatment; used to help determine time to response and
    #   therefore relapse rate
    timeEl.tx_st.prcr = NA_real_,
    # remission: 1 when patient is in remission, 0 otherwise
    remission = 0,
    # rate_txReinit: treatment re-initiation rate post gap; equals NA when
    #   patient is actively receiving treatment and during onset, postTx,
    #   and death phases
    rate_txReinit = NA_real_
  )
  
  return(X)
  
}


#------------------------------#
#~~~ Recurring calculations ~~~#
#------------------------------#

# Helper function to randomly sample from exponential distribution
#   given rate parameter lambda and optional scaling factor. If
#   F(x) = scaling factor, then we can randomly sample from conditional
#   exponential distribution with cdf F(t | T<x, t<x)
f.calc.rand_expntl <- function(lambda, cscale=1){
  
  # If rate is 0, event will never happen --> return Inf
  if(lambda == 0){
    return(Inf)
  }
  
  # Randomly sample from Uniform distribution on (0, 1)
  rand.unif <- runif(1)
  
  # Scale to facilitate sampling from conditional exponential
  rand.unif <- rand.unif * cscale
  
  # Calculate t given F(t)=p=rand.unif
  rand.expntl <- -log(1 - rand.unif) / lambda
  
  return(rand.expntl)
}

# Helper function to determine time to death given current age, sex, treatment
#   history, and disease state
f.calc.time_to_death <- function(X, inputs, mort){
  
  # Determine patient's current age in whole years
  age.fl <- floor(X$age.in)
  
  # Pull death rate given patient's age and sex
  rt.mort.pre <- mort[Age == age.fl & Female == X$female, rate_death]
  
  # Pull state-specific HR
  hr.mort <- inputs[[paste0("hr_mort_", X$state.curr)]]
  
  # If applicable, pull HR for patients who have completed at least 2
  #   lines of therapy and multiply by state-specific HR
  if(X$tx_line > 2 | (X$tx_line > 1 & X$phase.curr %in% c("gap","relapse","fail","postTx")) | X$phase.curr == "postTx"){
    if(age.fl < 30){
      hr.mort.trd <- inputs[["hr_mort_trd_v_mdd_age_18_29"]]
    } else if(age.fl < 50){
      hr.mort.trd <- inputs[["hr_mort_trd_v_mdd_age_30_49"]]
    } else if(age.fl < 70){
      hr.mort.trd <- inputs[["hr_mort_trd_v_mdd_age_50_69"]]
    } else {
      hr.mort.trd <- inputs[["hr_mort_trd_v_mdd_age_70_150"]]
    }
    hr.mort <- hr.mort * hr.mort.trd
  }
  
  # Apply hazard ratio
  rt.mort <- rt.mort.pre * hr.mort
  
  # Determine time to death
  tt.death <- f.calc.rand_expntl(rt.mort)
  return(tt.death)
  
}

# Function to calculate number of days to patient's next birthday. Truncate at
#   end of modeled time horizon.
f.calc.time_to_bday <- function(X, inputs){
  
  dpy <- inputs[["days_per_yr"]]
  
  # Determine time to end of modeled time horizon
  tt.stop <- (inputs[["yrs_modeled"]] * dpy) - X$time.in

  # Determine time to next birthday
  tt.bday <- (ceiling(X$age.in) - X$age.in) * dpy
  
  # If birthday is RIGHT NOW, set to yr from now
  if(near(tt.bday, 0)==TRUE){ tt.bday <- dpy }
  
  # Return earlier of time to next birthday and time to end of modeled time horizon
  return(min(c(tt.bday, tt.stop)))
}

# Function to estimate time to relapse
f.calc.time_to_relapse <- function(X, inputs){
  
  # Pull relapse rate and relevant hazard ratio, if applicable, based on
  #   state, treatment status, and time to response (if treated). Multiply
  #   base relapse rate by hazard ratio, if applicable, to calculate rate.
  
  # If time from start of current tx line to response is NA, then the
  #   patient is not currently receiving treatment
  if(is.na(X$timeEl.tx_st.prcr) == TRUE){
    lambda.rel <- inputs[[paste0("rate_relapse_", X$state.curr, "_notreat")]]
    
  # Otherwise, the patient is currently on treatment
  } else {
    lambda.rel <- inputs[[paste0("rate_relapse_late_", X$state.curr, "_line", X$tx_line)]]
    
    #If the patient was an early responder, pull & apply relevant hazard ratio
    if(X$timeEl.tx_st.prcr <= inputs[["maxDays_early_crpr"]]){
      hr.rel <- inputs[[paste0("hr_relapse_early_vs_late_", X$state.curr)]]
      lambda.rel <- lambda.rel * hr.rel
    }    
    
  }
  
  # Estimate time to relapse given rate
  tt.rel <- f.calc.rand_expntl(lambda=lambda.rel)
  
  return(tt.rel)
  
}

# Function to estimate time to remission
f.calc.time_to_remission <- function(X, inputs){
  
  # If patient is not in CR or if patient is already in remission, return
  #   artificially high number for time to remission
  if(X$remission == 1 | X$state.curr != "cr"){
    return(999)
  }
  
  # Otherwise, if in CR but not already in remission, calculate time in CR
  #   as of time.in
  timeInCR <- X$time.in - X$time.cr_st
  
  # Calculate time to remission given time in CR thus far. Currently defining
  #   remission as continuous CR for at least 3 months
  tt.rem <- (3*inputs[["days_per_yr"]]/12) - timeInCR
    
  # Return time to remission
  return(tt.rem)
  
}


#--- ORIGIN PHASE(S)/STATE(S) --------------------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

#---------------------#
#~~~ Disease onset ~~~#
#---------------------#

f.onset <- function(X, inputs, mort){
  
  # If model time 0...
  if( is.na(X$time.phase_st) == TRUE ){
    
    # Update phase start time
    X$time.phase_st <- X$time.in
    
    # If there are no treatments available, update next phase (post tx)
    #   and send back to simulation hub. No other updates needed.
    if( near(inputs[["max_tx_line"]], 0) == TRUE ){
      X$phase.next <- "postTx"
      return(X)
    }
    
    # Otherwise, if there is at least one treatment to be tried and the
    #   patient is NOT delaying its initiation, update next phase (init)
    #   and send back to simulation hub. No other updates needed.
    if( inputs[["p_txDelay_line1"]] < runif(1) ){
      X$phase.next <- "init"
      return(X)
    }
  }
  
  # For patients delaying their first line of therapy...
  
  # Determine time to death given current mortality rate
  tt.death <- f.calc.time_to_death(X=X, inputs=inputs, mort=mort)
  
  # Determine time to next birthday
  tt.bday <- f.calc.time_to_bday(X=X, inputs=inputs)
  
  # Determine time to remission. If not applicable, will be set to artificially
  #   high value
  tt.rem <- f.calc.time_to_remission(X=X, inputs=inputs)
  
  # # If patient is already in remission, update remission flag and set
  # #   time to remission to an artificially high value
  # if(near(tt.rem, 0) == TRUE){
  #   X$remission <- 1
  #   tt.rem <- 999
  # }
  
  # If patient is currently in nr state, determine times to spontaneous
  #   response and treatment initiation
  if(X$state.curr == "nr"){
    
    # Pull rates of spontaneous cr/pr
    lambda.cr <- inputs[['rate_cr_notreat']]
    lambda.pr <- inputs[['rate_pr_notreat']]
    lambda.resp <- lambda.pr + lambda.cr
    
    # Determine time to response (partial or complete)
    tt.resp <- f.calc.rand_expntl(lambda=lambda.resp)
    
    # Determine maximum possible time to treatment initiation assuming
    #   patient does not achieve spontaneous cr/pr in the meantime
    max_phase_days <- inputs[["maxDays_txDelay_line1"]]
    tt.txInit.max <- max_phase_days - X$time.in
    
    # Calculate treatment initiation rate
    # Exponential cdf:    F(t) = 1 - exp(-lambda*t)
    #                 --> lambda = log(1 - F(t)) / -t
    lambda.txInit <- -log(1 - inputs[["p_txDelay_line1"]]) / max_phase_days
    
    # Calculate scaling factor for conditional exponential distribution
    #   (i.e., P(D <= max_phase_days | D > X$time.in) = P(D <= tt.txInit.max)
    #   where D is time to treatment initiation) using treatment initiation
    #   rate calculated above
    pStar <- 1 - exp(-lambda.txInit * tt.txInit.max)
    
    # Determine time to treatment initiation, D, given D <= tt.txInit.max
    tt.txInit <- f.calc.rand_expntl(lambda=lambda.txInit, cscale=pStar)
    
    # If next event is transition to spontaneously-achieved PR/CR...
    if ( tt.resp < min(c(tt.death, tt.txInit, tt.bday)) ){
      
      # Determine whether transitioning to pr or cr and update next state
      if(runif(1) < lambda.pr / (lambda.cr + lambda.pr)){
        X$state.next <- "pr"
      } else {
        X$state.next <- "cr"
      }
      
      # Update time out
      X$time.out <- X$time.in + tt.resp
      
      # Return to simulation hub.
      # Note: Do not need to update next phase because unchanged
      #   (patient will return to onset phase next iteration)
      return(X)
      
      # Otherwise, if next event is initiation of 1st line of therapy...
    } else if ( tt.txInit < min(c(tt.death, tt.bday)) ){
      
      # Update next phase (init). Will be no change in state
      X$phase.next <- "init"
      
      # Update time out
      X$time.out <- X$time.in + tt.txInit
      
      # Return to simulation hub.
      return(X)
      
    }
    
    # Otherwise, if patient has already achieved spontaneous response,
    #   determine time to relapse
  } else {
    
    # Estimate time to relapse
    tt.rel <- f.calc.time_to_relapse(X=X, inputs=inputs)
    
    # If next event is relapse back to NR...
    # Note: Assume untreated patients in CR/PR won't resume treatment
    #   until they relapse, so time to relapse = time to next line of tx
    if ( tt.rel < min(c(tt.death, tt.bday, tt.rem)) ){
      
      # Update next phase (relapse) and state (NR)
      X$phase.next <- "relapse"
      X$state.next <- "nr"
      
      # Update time out
      X$time.out <- X$time.in + tt.rel
      
      # Return to simulation hub.
      return(X)
      
    }
  }
  
  # If patient is not experiencing spontaneous response, relapsing, or
  #   initiating treatment in the immediate future, patient must either be
  #   dying, celebrating a birthday, or entering remission. Update time out
  X$time.out <- X$time.in + min(c(tt.death, tt.bday, tt.rem))
  
  # If dying, update next phase and next state
  if( tt.death < min(c(tt.bday, tt.rem)) ){
    X$phase.next <- "death"
    X$state.next <- "x"
  }
  
  # Return to simulation hub.
  # Note: If next event is patient's birthday or change in remission status,
  #   do not need to update next phase or next state because both are unchanged
  #   (patient will return to onset phase in same state next iteration)  
  return(X)
  
}


#--- TRANSITIONAL PHASE(S)/STATE(S) --------------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

#----------------------------#
#~~~ Treatment initiation ~~~#
#----------------------------#

# Possible outcomes
#   - NR --> CR       >>Return
#   - NR --> PR       >>Return
#   - --> remission   >>Return
#   - Age + 1         >>Return
#   - End of phase    >>Leave
#   - Death           >>Leave
# Variables needed
#   - time.phase_st
#   - tx_line
#   - phase.curr
#   - phase.next
#   - time.in
#   - time.out
#   - timeEl.tx_st.prcr
#   - state.curr
#   - state.next
#   - age.in
#   - rate_txReinit
f.init <- function(X, inputs, mort){
  
  # If patient is entering initiation phase for the first time since failing/
  #   relapsing most recent treatment or since beginning of modeled time
  #   horizon, update phase start time, increment line of therapy, and reset
  #   treatment re-initiation rate variable.
  if( is.na(X$time.phase_st) == TRUE ){
    X$time.phase_st <- X$time.in
    X$tx_line <- X$tx_line + 1
    X$rate_txReinit <- NA
  }
  
  # Determine time to death given current mortality rate
  tt.death <- f.calc.time_to_death(X=X, inputs=inputs, mort=mort)
  
  # Determine time to phase end
  phase_days <- inputs[[paste0("phaseDays_init_line", X$tx_line)]]
  tt.phase_fn <- phase_days - (X$time.in - X$time.phase_st)
  
  # Determine time to next birthday
  tt.bday <- f.calc.time_to_bday(X=X, inputs=inputs)
  
  # Determine time to remission. If not applicable, will be set to
  #   artificially high value
  tt.rem <- f.calc.time_to_remission(X=X, inputs=inputs)

  # If patient is in NR state, determine time to response
  if(X$state.curr == "nr"){
    
    # Pull response rates
    lambda.pr <- inputs[[paste0("rate_pr_init_line", X$tx_line)]]
    lambda.cr <- inputs[[paste0("rate_cr_init_line", X$tx_line)]]
    lambda.resp <- lambda.pr + lambda.cr
    
    # Determine time to response (partial or complete)
    tt.resp <- f.calc.rand_expntl(lambda=lambda.resp)

    # If next event is transition to PR/CR...
    if ( tt.resp < min(c(tt.death, tt.phase_fn, tt.bday)) ){
      
      # Determine whether transitioning to pr or cr and update next state
      if(runif(1) < lambda.pr / (lambda.cr + lambda.pr)){
        X$state.next <- "pr"
      } else {
        X$state.next <- "cr"
      }
      
      # Update time out
      X$time.out <- X$time.in + tt.resp

      # Return to simulation hub.
      # Note: Do not need to update next phase because unchanged
      #   (patient will return to initiation phase next iteration)
      return(X)
      
    }
  }
  
  # For patients for whom next event is NOT a transition to pr/cr...
  
  # If patient is in CR/PR but time from initiation of current line of therapy
  #   to response has not yet been documented, document now.
  if(X$state.curr %in% c("pr","cr") & is.na(X$timeEl.tx_st.prcr)==TRUE){
    X$timeEl.tx_st.prcr <- X$time.in - X$time.phase_st
  }
  
  # Update time out
  X$time.out <- X$time.in + min(c(tt.death, tt.phase_fn, tt.bday, tt.rem))
  
  # If next event is death, set next phase to death
  if( tt.death < min(c(tt.phase_fn, tt.bday, tt.rem)) ){
    X$phase.next <- "death"
    X$state.next <- "x"
    
  # Otherwise, if next event is end of phase, set next phase
  #   to either maintenance, initiation extension, or tx failure
  #   (potential gap) based on current state
  } else if( tt.phase_fn < min(c(tt.bday, tt.rem)) ){
    
    # If in CR, go to maintenance phase
    if( X$state.curr == "cr" ){
      X$phase.next <- "maint"
      
    # If in PR, go to initiation-phase extension
    } else if (X$state.curr == "pr"){
      X$phase.next <- "initExt"

    # If in NR, go to failed tx phase
    } else {
      X$phase.next <- "fail"
    }
    
  }
  
  # Return to simulation hub.
  # Note: If next event is patient's birthday or change in remission status,
  #   do not need to update next phase or state because unchanged (patient
  #   will return to initiation phase in same state next iteration)  
  return(X)
}


#----------------------------#
#~~~ Initiation extension ~~~#
#----------------------------#

# Possible outcomes
#   - PR --> CR       >>Return
#   - --> remission   >>Return
#   - PR/CR --> NR    >>Leave (relapse?)
#   - Age + 1         >>Return
#   - End of phase    >>Leave (failure or maintenance)
#   - Death           >>Leave (death)
f.initExt <- function(X, inputs, mort){
  
  # If patient is entering initiation-phase extension for the first
  #   time while receiving the current line of therapy, update phase
  #   start time
  if( is.na(X$time.phase_st) == TRUE ){
    X$time.phase_st <- X$time.in
  }
  
  # If patient is in CR but time from initiation of current line of therapy to
  #   CR has not yet been documented, document now.
  # Note: Must be updated before calculating time to relapse.
  phaseDays_init <- inputs[[paste0("phaseDays_init_line", X$tx_line)]]
  if(X$state.curr == "cr" & X$timeEl.tx_st.prcr <= phaseDays_init){
    X$timeEl.tx_st.prcr <- X$time.in - X$time.phase_st + phaseDays_init
  }
  
  # Determine time to death given current mortality rate
  tt.death <- f.calc.time_to_death(X=X, inputs=inputs, mort=mort)
  
  # Determine time to phase end
  phase_days <- inputs[[paste0("phaseDays_initExt_line", X$tx_line)]]
  tt.phase_fn <- phase_days - (X$time.in - X$time.phase_st)
  
  # Determine time to next birthday
  tt.bday <- f.calc.time_to_bday(X=X, inputs=inputs)
  
  # Determine time to remission. If not applicable, will be set to
  #   artificially high value
  tt.rem <- f.calc.time_to_remission(X=X, inputs=inputs)
  
  # Estimate time to relapse
  tt.rel <- f.calc.time_to_relapse(X=X, inputs=inputs)
  
  # If patient is in PR state, determine time to CR
  if(X$state.curr == "pr"){
    
    # Pull complete-response rate and determine time to CR
    lambda.cr <- inputs[[paste0("rate_cr_initExt_line", X$tx_line)]]
    tt.cr <- f.calc.rand_expntl(lambda=lambda.cr)
    
    # If next event is transition to CR...
    if ( tt.cr < min(c(tt.death, tt.phase_fn, tt.bday, tt.rel)) ){
      
      # Update next state and time out
      X$state.next <- "cr"
      X$time.out <- X$time.in + tt.cr
      
      # # Update time to CR from time of initiation of current line of therapy
      # phaseDays_init <- inputs[[paste0("phaseDays_init_line", X$tx_line)]]
      # X$timeEl.tx_st.prcr <- X$time.out - X$time.phase_st + phaseDays_init
      
      # Return to simulation hub.
      # Note: Do not need to update next phase because unchanged
      #   (patient will return to initiation-extension phase next iteration)
      return(X)
      
    }
  }
  
  # For patients for whom next event is NOT a transition to CR...
  
  # Update time out
  X$time.out <- X$time.in + min(c(tt.death, tt.phase_fn, tt.bday, tt.rem, tt.rel))
  
  # If next event is death, set next phase to death
  if( tt.death < min(c(tt.phase_fn, tt.bday, tt.rem, tt.rel)) ){
    X$phase.next <- "death"
    X$state.next <- "x"
    
  # Otherwise, if next event is end of phase, set next phase to either
  #   maintenance or tx failure (potential gap) based on current state.
  #   Update next state accordingly.
  } else if( tt.phase_fn < min(c(tt.bday, tt.rem, tt.rel)) ){
    
    # If in CR, go to maintenance phase
    if( X$state.curr == "cr" ){
      X$phase.next <- "maint"
      
    # Otherwise, if still in PR, go to either maintenance or failed tx phase
    } else if(runif(1) < inputs[["p_pr_to_maintenance"]]) {
      X$phase.next <- "maint"
    } else {
      # Patients assumed to discontinue treatment and return to nr state
      X$phase.next <- "fail"
      X$state.next <- "nr"
    }

  # Otherwise, if relapsing, go to relapse phase regardless of time remaining
  #   in initiation-phase extension
  } else if(tt.rel < min(c(tt.bday, tt.rem)) ){
    X$phase.next <- "relapse"
    X$state.next <- "nr"
  }
  
  # Return to simulation hub.
  # Note: If next event is patient's birthday or a change in remission status,
  #   do not need to update next phase or state because unchanged (patient will
  #   return to initiation-extension phase in same state next iteration) 
  return(X)
  
}


#-------------------#
#~~~ Maintenance ~~~#
#-------------------#

# Possible outcomes
#   - --> remission   >>Return
#   - PR/CR --> NR    >>Leave (relapse)
#   - Age + 1         >>Return
#   - Death           >>Leave (death)           >>
f.maint <- function(X, inputs, mort){
  
  # If patient is entering maintenance phase for the first time since
  #   initiating any treatment, update phase start time
  if( is.na(X$time.phase_st) == TRUE ){
    X$time.phase_st <- X$time.in
  }
  
  # Determine time to death given current mortality rate
  tt.death <- f.calc.time_to_death(X=X, inputs=inputs, mort=mort)
  
  # Determine time to next birthday
  tt.bday <- f.calc.time_to_bday(X=X, inputs=inputs)
  
  # Determine time to remission. If not applicable, will be set to artificially
  #   high value
  tt.rem <- f.calc.time_to_remission(X=X, inputs=inputs)
  
  # Estimate time to relapse
  tt.rel <- f.calc.time_to_relapse(X=X, inputs=inputs)
  
  # # If not already in remission, determine time to remission
  # if(X$remission == 0 & X$state.curr == "cr"){
  #   
  #   # Obtain initiation-phase and initiation-phase extension durations
  #   durInit <- inputs[[paste0("phaseDays_init_line", X$tx_line)]]
  #   durInitExt <- inputs[[paste0("phaseDays_initExt_line", X$tx_line)]]
  #   
  #   # Calculate total time in CR (for the current line of therapy) as of
  #   #   time.in. If patient went to initiation-phase extension, account
  #   #   for this.
  #   timeInCR <- (durInit - X$timeEl.tx_st.prcr) + (X$time.in - X$time.phase_st)
  #   if(durInit < X$timeEl.tx_st.prcr){
  #     timeInCR <- timeInCR + durInitExt
  #   }
  #   
  #   # Calculate time to remission given time in CR thus far. Currently defining
  #   #   remission as continuous CR for at least 3 months
  #   tt.rem <- (3*inputs[["days_per_yr"]]/12) - timeInCR
  #   
  #   # If patient is in remission as of time.in, set remission flag to 1
  #   if(near(tt.rem, 0) == TRUE){ X$remission <- 1 }
  #   
  # # Otherwise, set time to remission to an artificially high number
  # } else {
  #   tt.rem <- 999
  # }

  # Update time out
  X$time.out <- X$time.in + min(c(tt.death, tt.bday, tt.rel, tt.rem))
  
  # If next event is death, set next phase to death
  if( tt.death < min(c(tt.bday, tt.rem, tt.rel)) ){
    X$phase.next <- "death"
    X$state.next <- "x"
    
  # Otherwise, if next event is relapse, set next phase to relapse and
  #   next state to nr
  } else if( tt.rel < min(c(tt.bday, tt.rem)) ){
    X$phase.next <- "relapse"
    X$state.next <- "nr"
  }
  
  # Return to simulation hub.
  # Note: If next event is patient's birthday or a change in remission status,
  #   do not need to update next phase because unchanged (patient will return
  #   to maintenance phase in same state next iteration)
  return(X)
  
}


#---------------------#
#~~~ Treatment gap ~~~#
#---------------------#

# Period between consecutive lines of therapy during which patient is untreated

# Possible outcomes
#   - NR --> PR/CR    >>Return
#   - --> remission   >>Return
#   - PR/CR --> NR    >>Leave (relapse --> init)
#   - Resume tx       >>Leave
#   - Age + 1         >>Return
#   - Death           >>Leave (death)
f.gap <- function(X, inputs, mort){
  
  # If patient is entering gap phase for the first time since
  #   most recent line of therapy, update phase start time
  if( is.na(X$time.phase_st) == TRUE ){
    X$time.phase_st <- X$time.in
  }
  
  # Determine time to death given current mortality rate
  tt.death <- f.calc.time_to_death(X=X, inputs=inputs, mort=mort)
  
  # Determine time to next birthday
  tt.bday <- f.calc.time_to_bday(X=X, inputs=inputs)
  
  # Determine time to remission. If not applicable, will be set to artificially
  #   high value
  tt.rem <- f.calc.time_to_remission(X=X, inputs=inputs)
  
  # If patient is currently in nr state, determine times to spontaneous cr/pr
  #   and treatment re-initiation
  if(X$state.curr == "nr"){
    
    # Pull rates of spontaneous cr/pr
    lambda.cr <- inputs[['rate_cr_notreat']]
    lambda.pr <- inputs[['rate_pr_notreat']]
    lambda.resp <- lambda.pr + lambda.cr
    
    # Determine time to response (partial or complete)
    tt.resp <- f.calc.rand_expntl(lambda=lambda.resp)
    
    # Determine time to treatment re-initiation
    tt.txReinit <- f.calc.rand_expntl(lambda=X$rate_txReinit)
    
    # If next event is transition to PR/CR...
    if ( tt.resp < min(c(tt.death, tt.txReinit, tt.bday)) ){
      
      # Determine whether transitioning to pr or cr and update next state
      if(runif(1) < lambda.pr / (lambda.cr + lambda.pr)){
        X$state.next <- "pr"
      } else {
        X$state.next <- "cr"
      }
      
      # Update time out
      X$time.out <- X$time.in + tt.resp

      # Return to simulation hub.
      # Note: Do not need to update next phase because unchanged
      #   (patient will return to gap phase next iteration)
      return(X)
      
    # Otherwise, if next event is treatment resumption...
    } else if ( tt.txReinit < min(c(tt.death, tt.bday)) ){
      
      # Update next phase (init). Will be no change in state
      X$phase.next <- "init"
      
      # Update time out
      X$time.out <- X$time.in + tt.txReinit
      
      # Return to simulation hub.
      return(X)
      
    }
    
  # Otherwise, if patient has already achieved spontaneous response,
  #   determine time to relapse
  } else {
    
    # Estimate time to relapse
    tt.rel <- f.calc.time_to_relapse(X=X, inputs=inputs)
    
    # If next event is relapse back to NR...
    # Note: Assume untreated patients in CR/PR won't resume treatment
    #   until they relapse, so time to relapse = time to next line of tx
    if ( tt.rel < min(c(tt.death, tt.bday, tt.rem)) ){
      
      # Update next phase (relapse) and state (NR)
      X$phase.next <- "relapse"
      X$state.next <- "nr"
      
      # Update time out
      X$time.out <- X$time.in + tt.rel
      
      # Return to simulation hub.
      return(X)
      
    }
  }
  
  # If patient is not experiencing spontaneous response, relapsing, or
  #   resuming treatment in the immediate future, they must be either
  #   dying, celebrating a birthday, or entering remission. Update time out
  X$time.out <- X$time.in + min(c(tt.death, tt.bday, tt.rem))
  
  # If dying, update next phase and next state
  if( tt.death < min(c(tt.bday, tt.rem)) ){
    X$phase.next <- "death"
    X$state.next <- "x"
  }
  
  # Return to simulation hub.
  # Note: If next event is patient's birthday or a change in remission status,
  #   do not need to update next phase or state because both are unchanged
  #   (patient will return to gap phase in same state next iteration)  
  return(X)
  
}



#--- TEMPORARY PSEUDO-PHASE(S)/STATE(S) ----------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

#---------------#
#~~~ Relapse ~~~#
#---------------#

# Applicable only when relapsing precipitates change in phase (i.e., relapse
#   from spontaneous CR/PR during treatment gaps/delays and relapse from
#   treatment-attributed CR/PR)

f.relapse <- function(X, inputs, mort){
  
  # Update phase start time
  X$time.phase_st <- X$time.in

  # Regardless of whether patient is relapsing following spontaneously-
  #   achieved or treatment-attributed CR/PR, if there are no more treatment
  #   options available to the patient, the patient's next stop will be
  #   the post treatment phase. Update next phase accordingly.
  if( inputs[["max_tx_line"]] <= X$tx_line ) {
    X$phase.next <- "postTx"
    
  # If, on the other hand, there are more treatments to be tried...
  } else {
    
    # If variable indicating time from initiation of current line of
    #   therapy to time of response is NA, patient is relapsing from
    #   spontaneously-achieved CR/PR and is now poised to initiate
    #   the next line of therapy. Update next phase accordingly.
    if( is.na(X$timeEl.tx_st.prcr) == TRUE ){
      X$phase.next <- "init"
      
      # If, on the other hand, the variable for time from current tx
      #   initiation to response is NOT NA, patient is relapsing from
      #   treatment-attributed CR/PR...
    } else {
      
      # Pull likelihood of experiencing a gap in treatment following relapse
      p_gap <- inputs[[paste0("p_txGap_line", X$tx_line, "_relapse")]]
      
      # Determine whether the patient will delay initiation of the next line
      #   of therapy. If they will delay...
      if( runif(1) < p_gap ){
        
        # Pull and save treatment re-initiation rate for patients
        #   delaying treatment after relapsing from tx-attributed CR/PR
        #   following successful treatment with the patient's current
        #   (most recent) line of therapy
        X$rate_txReinit <- inputs[[paste0("rate_txReinit_line", X$tx_line, "_relapse")]]
        
        # Prepare to send patient to the treatment-gap phase
        X$phase.next <- "gap"
        
        # Otherwise, if the patient is ready to start the next line of therapy
        #   without delay, prepare to send the patient to the initiation phase
      } else {
        X$phase.next <- "init"
      }
    }
  }
  
  # Reset time from current tx initiation to response
  X$timeEl.tx_st.prcr <- NA
  
  # Return to simulation hub.
  # Note: Temporary state so no time is elapsing --> no need to update X$time.out.
  # Note: Regardless of what happens next, patient returned to nr state upon
  #   relapse and will stay in nr state through entry into the next phase -->
  #   no need to update X$status.next.
  return(X)
  
}


#-------------------------#
#~~~ Treatment failure ~~~#
#-------------------------#

f.fail <- function(X, inputs, mort){
  
  # Update phase start time
  X$time.phase_st <- X$time.in

  # Reset time from current line of therapy initiation to response as
  #   some patients may be considered "treatment failures" despite
  #   having achieved (and sustained) partial response (for a time)
  X$timeEl.tx_st.prcr <- NA
  
  # If there no more treatment options available to the patient, the
  #   patient's next stop will be the post treatment phase. Update next
  #   phase accordingly.
  if( inputs[["max_tx_line"]] <= X$tx_line ) {
    X$phase.next <- "postTx"
    
    # If, on the other hand, there are more treatments to be tried...
  } else {
    
    # Pull likelihood of experiencing a gap in treatment following failure
    p_gap <- inputs[[paste0("p_txGap_line", X$tx_line, "_failure")]]
    
    # Determine whether the patient will delay initiation of the next line
    #   of therapy. If they will delay...
    if( runif(1) < p_gap ){
      
      # Pull and save treatment re-initiation rate for patients
      #   delaying treatment after failing most recent line of therapy
      X$rate_txReinit <- inputs[[paste0("rate_txReinit_line", X$tx_line, "_failure")]]
      
      # Prepare to send patient to the treatment-gap phase
      X$phase.next <- "gap"
      
      # Otherwise, if the patient is ready to start the next line of therapy
      #   without delay, prepare to send the patient back to the initiation phase
    } else {
      X$phase.next <- "init"
    }
  }
  
  # Return to simulation hub.
  # Note: Temporary state so no time is elapsing --> no need to update X$time.out.
  # Note: Regardless of what happens next, patient remained in/returned to nr
  #   state following tx failure and will stay in nr state through entry into the
  #   next phase --> no need to update X$status.next.
  return(X)
  
}



#--- TERMINAL AND PSEUDO-TERMINAL PHASE(S)/STATE(S) ----------------------------
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

#----------------------#
#~~~ POST TREATMENT ~~~#
#----------------------#

# Patient has already tried all available treatments in pathway

# Possible outcomes
#   - NR --> PR/CR    >>Return
#   - PR/CR --> NR    >>Return
#   - Age + 1         >>Return
#   - Death           >>Leave (death)
f.postTx <- function(X, inputs, mort){
  
  # If patient is entering post-treatment phase for the first time, update
  #   phase start time
  if( is.na(X$time.phase_st) == TRUE ){
    X$time.phase_st <- X$time.in
  }
  
  # Determine time to death given current mortality rate
  tt.death <- f.calc.time_to_death(X=X, inputs=inputs, mort=mort)
  
  # Determine time to next birthday
  tt.bday <- f.calc.time_to_bday(X=X, inputs=inputs)
  
  # Determine time to remission. If not applicable, will be set to artificially
  #   high value
  tt.rem <- f.calc.time_to_remission(X=X, inputs=inputs)
  
  # If patient is currently in nr state, determine time to spontaneous cr/pr
  if(X$state.curr == "nr"){
    
    # Pull rates of spontaneous cr/pr
    lambda.cr <- inputs[['rate_cr_notreat']]
    lambda.pr <- inputs[['rate_pr_notreat']]
    lambda.resp <- lambda.pr + lambda.cr
    
    # Determine time to response (partial or complete)
    tt.resp <- f.calc.rand_expntl(lambda=lambda.resp)
    
    # If next event is transition to PR/CR...
    if ( tt.resp < min(c(tt.death, tt.bday)) ){
      
      # Determine whether transitioning to pr or cr and update next state
      if(runif(1) < lambda.pr / (lambda.cr + lambda.pr)){
        X$state.next <- "pr"
      } else {
        X$state.next <- "cr"
      }
      
      # Update time out
      X$time.out <- X$time.in + tt.resp
      
      # Return to simulation hub.
      # Note: Do not need to update next phase because unchanged
      #   (patient will always return to post-tx phase unless deceased)
      return(X)
      
    }
    
  # Otherwise, if patient has already achieved spontaneous response,
  #   determine time to relapse
  } else {
    
    # Estimate time to relapse
    tt.rel <- f.calc.time_to_relapse(X=X, inputs=inputs)
    
    # If next event is relapse back to NR...
    if ( tt.rel < min(c(tt.death, tt.bday, tt.rem)) ){
      
      # Update next state (NR).
      # Note: No need to send patient to relapse "phase" because the only *true*
      #   phases available to patients who have exhausted all treatment options
      #   are this one (post tx) and death, and only thing that changes when a
      #   patient relapses while in post-tx phase is the patient's state
      X$state.next <- "nr"
      
      # Update time out
      X$time.out <- X$time.in + tt.rel
      
      # Return to simulation hub.
      return(X)
      
    }
  }
  
  # If patient is not experiencing spontaneous response or relapsing, the
  #   patient must either be dying, celebrating a birthday, or entering
  #   remission. Update time out.
  X$time.out <- X$time.in + min(c(tt.death, tt.bday, tt.rem))
  
  # If dying, update next phase and next state
  if( tt.death < min(c(tt.bday, tt.rem)) ){
    X$phase.next <- "death"
    X$state.next <- "x"
  }
  
  # Return to simulation hub.
  # Note: If next event is patient's birthday or a change in remission status,
  #   do not need to update next phase or state because both are unchanged
  #   (patient will return to post-tx phase in same state next iteration)  
  return(X)
  
}


#-------------#
#~~~ Death ~~~#
#-------------#

f.death <- function(X, inputs, mort){
  
  # Update phase start time
  X$time.phase_st <- X$time.in
  
  # Fast forward to end of modeled time horizon
  X$time.out <- inputs[['yrs_modeled']] * inputs[['days_per_yr']]
  
  # Return to simulation hub.
  return(X)
  
}



#--- SIMULATION FUNCTIONS ------------------------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

#------------------------------#
#~~~ Patient simulation hub ~~~#
#------------------------------#

f.simulate_pt <- function(X, inputs, mort){
  
  t <- 0

  stopTime <- inputs[['yrs_modeled']] * inputs[['days_per_yr']]
  
  n_iter <- 0
  
  # Create empty trace data set to which to append patient data
  trc <- copy(X)[0L]
  
  while (t < stopTime & n_iter < 5000){
    
    # Set phase start time to NA if switching phases
    if(X$phase.curr != X$phase.next){
      X$time.phase_st <- NA
    }
    
    # Update starting phase and state to prior ending phase and state.
    # Also update age.in and time.in
    X$phase.curr <- X$phase.next
    X$state.curr <- X$state.next
    X$age.in <- X$age.in + (X$time.out - X$time.in) / inputs[['days_per_yr']]
    X$time.in <- X$time.out
    
    # If patient is NOT in CR, set time at which patient entered CR to NA and
    #   set remission flag to 0
    if(X$state.curr != "cr"){
      X$time.cr_st <- NA
      X$remission <- 0
      
    # Otherwise, if patient is in CR...
    } else {
      
      # If time at which patient entered CR has not yet been documented, do so now
      if(is.na(X$time.cr_st) == TRUE){
        X$time.cr_st <- X$time.in
      }
      
      # Calculate time to remission
      tt.rem <- f.calc.time_to_remission(X=X, inputs=inputs)
      
      # If time to remission is 0, flag patient as in remission
      if(near(tt.rem, 0) == TRUE){
        X$remission <- 1
      }
    }

    # Call appropriate transition function using current phase name
    X <- eval(call(paste0("f.", X$phase.curr), X, inputs, mort))

    # Append results to output data table
    trc <- rbindlist(list(trc, X))
    
    # Increment time
    t <- X$time.out
    
    n_iter <- n_iter + 1

  }
  
  return(trc)
  
}


#----------------------------------------#
#~~~ Treatment pathway simulation hub ~~~#
#----------------------------------------#

# Generate a random patient, run through simulation, save results
f.simulate_txPath <- function(inputs.eff, inputs.mort, inputs.oth, seed=112358){
  
  # Start the clock!
  ptm <- proc.time()
  
  # Set seed
  set.seed(seed)

  # Combine efficacy and other inputs into single data table
  #inputs <- rbindlist(list(inputs.eff, inputs.oth))
  # Combine efficacy and other inputs into single list
  inputs <- c(inputs.eff, inputs.oth)
  
  # Number of patients to simulate
  n_pts <- inputs[['n_pts']]
  
  
  #--- Determine starting age of each simulated patient ---#
  
  # Create probability distribution for age bands
  in.p.age <- c(inputs[['pct_age_18_24']],
                inputs[['pct_age_25_34']],
                inputs[['pct_age_35_44']],
                inputs[['pct_age_45_54']],
                inputs[['pct_age_55_64']])
  
  # Determine the lower bound of the age band into which each patient falls
  ages.pre <- sample(x=c(18, 25, 35, 45, 55),
                     size=n_pts,
                     replace=T,
                     prob=in.p.age)
  
  # Assume age is uniformly distributed within each band. Note, each band
  #   is 10 years wide save for the first which is only 7 years wide [18, 25)
  ages <- mapply(function(x, y){if(x==18){x + (7*y)} else{x + (10*y)}},
                 x=ages.pre,
                 y=runif(n_pts))
  
  
  #--- Determine sex of each simulated patient ---#
  
  # Create probability distribution for sex
  in.p.sex <- c(1 - inputs[['pct_female']], inputs[['pct_female']])
  
  # Randomly sample from sex distribution
  sexes <- sample(x=0:1, size=n_pts, replace=T, prob=in.p.sex)
  
  
  #--- Determine SES quartile of each simulated patient ---#
  
  # Create probability distribution for SES quartiles
  in.p.SES <- c(inputs[['pct_SES_qrtl1']],
                inputs[['pct_SES_qrtl2']],
                inputs[['pct_SES_qrtl3']],
                inputs[['pct_SES_qrtl4']])
  
  # Randomly sample from SES distribution
  SESes <- sample(x=1:4, size=n_pts, replace=T, prob=in.p.SES)
  
  
  #--- Conduct individual patient simulation ---#
  
  # for (i in 1:n_pts){
  # 
  #   # Set seed for reproducibility
  #   set.seed(seed+i)
  # 
  #   X <- f.create_pt(id=i, female=sexes[[i]], age.in=ages[[i]], SES=SESes[[i]])
  # 
  #   trace <- f.simulate_pt(X=X, inputs=inputs, mort=inputs.mort)
  # 
  #   if(i==1){results <- copy(trace)} else {results <- rbindlist(list(results, trace))}
  # }
  
  # parallel processing loop using foreach
  results <- foreach(
    i=1:n_pts,
    # .combine = 'bind_rows',
    .combine = function(x,y)rbindlist(list(x,y)),
    .packages = c('data.table','dtplyr','dplyr'),
    .export = c('f.create_pt', 'f.calc.rand_expntl',
                'f.calc.time_to_death', 'f.calc.time_to_bday',
                'f.calc.time_to_relapse', 'f.calc.time_to_remission',
                'f.onset', 'f.init', 'f.initExt', 'f.maint', 'f.gap',
                'f.relapse', 'f.fail', 'f.postTx', 'f.death',
                'f.simulate_pt')
  ) %dopar% {

    # Set seed for reproducibility
    set.seed(seed+i)

    X <- f.create_pt(id=i, female=sexes[[i]], age.in=ages[[i]], SES=SESes[[i]])

    trace <- f.simulate_pt(X=X, inputs=inputs, mort=inputs.mort)

    return(trace)

  }
  
  # Stop the clock
  print(proc.time() - ptm)
  
  return(results) 
}  



