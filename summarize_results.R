
f.summarize <- function(pw, inputs.dt.tx, inputs.oth, simRes, summIn.forDSA=FALSE, seed=581321){
  
  set.seed(seed)
  
  dpy <- inputs.oth[['days_per_yr']]
  disc.cost <- inputs.oth[['discount_cost']]
  disc.health <- inputs.oth[['discount_health']]
  maxWTP <- inputs.oth[['maxWTP']]
  
  res.sim <- copy(simRes)
  
  # Add the following binary variables and variables indicating model year(s):
  #   --> onTx: binary variable indicating treatment status (1 = treated)
  #   --> cr/pr/nr: binary variables indicating current disease status
  #   --> trd: binary variable indicating treatment experience (1 = completed 2+ lines of therapy)
  #   --> init_st: binary variable indicating start of each initiation phase
  #   --> yr.in: model year corresponding to time.in (1st yr is yr 0)
  #   --> yr.out: model year corresponding to time.out (1st yr is yr 0)
  res.sim[, `:=`(onTx = ifelse(phase.curr %in% c("init", "initExt", "maint"), 1, 0),
                 cr = ifelse(state.curr == "cr", 1, 0),
                 pr = ifelse(state.curr == "pr", 1, 0),
                 nr = ifelse(state.curr == "nr", 1, 0),
                 trd = ifelse(tx_line >2 | (tx_line >1 & phase.curr %in% c("gap","relapse","fail","postTx")) | phase.curr == "postTx", 1, 0),
                 init_st = 1*(phase.curr == "init" & near(time.in, time.phase_st) == TRUE),
                 yr.in = floor(time.in / dpy),
                 yr.out = floor(time.out / dpy))]
  
  # Create a lagged version of the remission flag to help count number of unique
  #   instances of remission
  res.sim[,`:=` (lag.remission = c(0, remission[-.N])), by=id]
  
  # Add the following duration and discount-factor variables:
  #   --> dur: time (days) elapsed
  #   --> dur1: time (days) elapsed in model year corresponding to time.in (yr.in)
  #   --> disc.cost.in: discount factor for costs in model year corresponding to
  #       time.in (yr.in)
  #   --> disc.cost.out: discount factor for costs in model year corresponding to
  #       time.out (yr.out)
  #   --> disc.hlth.in: discount factor for health benefits in model year
  #       corresponding to time.in (yr.in)
  #   --> disc.hlth.out: discount factor for health benefits in model year
  #       corresponding to time.out (yr.out)
  #   --> disc.cost.wa: weighted average of the two discount factors for costs
  #   --> disc.hlth.wa: weighted average of the two discount factors for health benefits
  res.sim[,`:=`(dur = time.out - time.in,
                dur1 = ifelse(yr.in < yr.out, (yr.out*dpy) - time.in, time.out - time.in),
                disc.cost.in = 1 / (1+disc.cost)^yr.in,
                disc.cost.out = 1 / (1+disc.cost)^yr.out,
                disc.hlth.in = 1 / (1+disc.health)^yr.in,
                disc.hlth.out = 1 / (1+disc.health)^yr.out)]
  res.sim[,`:=`(disc.cost.wa = ifelse(dur==0, disc.cost.in, (dur1*disc.cost.in + (dur-dur1)*disc.cost.out) / dur),
                disc.hlth.wa = ifelse(dur==0, disc.hlth.in, (dur1*disc.hlth.in + (dur-dur1)*disc.hlth.out) / dur))]
  
  # Add treatment-related unit costs and AE incidence rates by treatment status
  #   and line of therapy, then determine who experienced which AE(s) when
  # res.tmp2 <- res.sim[inputs.dt.tx, on = c("tx_line","onTx")
  res.tmp2 <- inputs.dt.tx[res.sim, on = c("tx_line","onTx")
  ][,`:=`(sae = init_st*(runif(.N) <= p.sae),
          ae1 = init_st*(runif(.N) <= p.ae1),
          ae2 = init_st*(runif(.N) <= p.ae2),
          ae3 = init_st*(runif(.N) <= p.ae3),
          ae4 = init_st*(runif(.N) <= p.ae4)
  )]
  
  # Replace NAs with 0s
  # res.tmp2[is.na(res.tmp2),] <- 0
  
  # Add the following additional unit cost(s), annual wages, work-loss rates,
  # utility multipliers, and disutilities:
  #   --> cpy.state: state-specific cost per year (also dependent on trd status)
  #   --> wage: annual earnings
  #   --> pct_abs: excess absenteeism (as percent of earnings)
  #   --> pct_pres: excess presenteeism (as percent of earnings)
  #   --> util.state: state-specific utility multiplier
  #   --> disutil.ae: disutility due to adverse event(s)
  res.tmp2[,`:=`(cpy.state = (inputs.oth[["costPerYr_direct_mdd_cr"]]*cr*(1-trd)) +
                   (inputs.oth[["costPerYr_direct_mdd_pr"]]*pr*(1-trd)) +
                   (inputs.oth[["costPerYr_direct_mdd_nr"]]*nr*(1-trd)) +
                   (inputs.oth[["costPerYr_direct_trd_cr"]]*cr*trd) +
                   (inputs.oth[["costPerYr_direct_trd_pr"]]*pr*trd) +
                   (inputs.oth[["costPerYr_direct_trd_nr"]]*nr*trd),
                 wage = (inputs.oth[["incomePerYr_age_18_24_male"]]*(female==0 & age.in>=18 & age.in<25)) +
                   (inputs.oth[["incomePerYr_age_25_34_male"]]*(female==0 & age.in>=25 & age.in<35)) +
                   (inputs.oth[["incomePerYr_age_35_44_male"]]*(female==0 & age.in>=35 & age.in<45)) +
                   (inputs.oth[["incomePerYr_age_45_54_male"]]*(female==0 & age.in>=45 & age.in<55)) +
                   (inputs.oth[["incomePerYr_age_55_64_male"]]*(female==0 & age.in>=55 & age.in<65)) +
                   (inputs.oth[["incomePerYr_age_65_150_male"]]*(female==0 & age.in>=65)) +
                   (inputs.oth[["incomePerYr_age_18_24_female"]]*(female==1 & age.in>=18 & age.in<25)) +
                   (inputs.oth[["incomePerYr_age_25_34_female"]]*(female==1 & age.in>=25 & age.in<35)) +
                   (inputs.oth[["incomePerYr_age_35_44_female"]]*(female==1 & age.in>=35 & age.in<45)) +
                   (inputs.oth[["incomePerYr_age_45_54_female"]]*(female==1 & age.in>=45 & age.in<55)) +
                   (inputs.oth[["incomePerYr_age_55_64_female"]]*(female==1 & age.in>=55 & age.in<65)) +
                   (inputs.oth[["incomePerYr_age_65_150_female"]]*(female==1 & age.in>=65)),
                 pct_abs = inputs.oth[["pct_absenteeism_cr"]]*cr +
                   inputs.oth[["pct_absenteeism_pr"]]*pr +
                   inputs.oth[["pct_absenteeism_nr"]]*nr,
                 pct_pres = inputs.oth[["pct_presenteeism_cr"]]*cr +
                   inputs.oth[["pct_presenteeism_pr"]]*pr +
                   inputs.oth[["pct_presenteeism_nr"]]*nr,
                 util.state = (inputs.oth[["util_mdd_cr"]]*cr*(1-trd)) +
                   (inputs.oth[["util_mdd_pr"]]*pr*(1-trd)) +
                   (inputs.oth[["util_mdd_nr"]]*nr*(1-trd)) +
                   (inputs.oth[["util_trd_cr"]]*cr*trd) +
                   (inputs.oth[["util_trd_pr"]]*pr*trd) +
                   (inputs.oth[["util_trd_nr"]]*nr*trd),
                 disutil.ae = (inputs.oth[["disutil_sae"]]*sae) +
                   (inputs.oth[["disutil_ae1"]]*ae1) +
                   (inputs.oth[["disutil_ae2"]]*ae2) +
                   (inputs.oth[["disutil_ae3"]]*ae3) +
                   (inputs.oth[["disutil_ae4"]]*ae4) )]
  
  # Calculate costs per row (i.e., apply unit costs to applicable periods).
  res.tmp2[, `:=`(cost.d.tx.disc = (cpy.tx/dpy)*dur*disc.cost.wa + c1x.tx*disc.cost.in,
                  cost.i.trnsprt.disc = (cpy.tx.transport/dpy)*dur*disc.cost.wa + c1x.tx.transport*disc.cost.in,
                  cost.d.state.disc = (cpy.state/dpy)*dur*disc.cost.wa,
                  cost.i.abs.disc = (wage/dpy)*dur*pct_abs*disc.cost.wa,
                  cost.i.pres.disc = (wage/dpy)*dur*pct_pres*disc.cost.wa
  )]
  
  # Summarize by patient
  #   --> relapse: transition from PR/CR to NR, except when at the end of the 
  #           initiation-phase extension due to treatment discontinuation
  #           precipitated by lack of improvement from PR to CR
  #   --> fail: failure to achieve CR during the initiation phase and
  #           initiation-phase extension
  #             >> Note: Achievement of sustained PR and subsequent procession to
  #                 maintenance (while in PR), only possible when
  #                 p_pr_to_maintenance > 0, is NOT considered treatment failure
  res.pt <- res.tmp2[,
                     # Direct costs
                     .(cost.d.tx.disc = sum(cost.d.tx.disc),
                       cost.d.state.disc = sum(cost.d.state.disc),
                       # Indirect costs
                       cost.i.trnsprt.disc = sum(cost.i.trnsprt.disc),
                       cost.i.abs.disc = sum(cost.i.abs.disc),
                       cost.i.pres.disc = sum(cost.i.pres.disc),
                       # Number of lines of therapy initiated
                       n.tx = sum(init_st),
                       # Number of AEs experienced
                       n.sae = sum(sae),
                       n.ae1 = sum(ae1),
                       n.ae2 = sum(ae2),
                       n.ae3 = sum(ae3),
                       n.ae4 = sum(ae4),
                       # Numbers of relapses and failures
                       n.relapse = sum(1 * (state.curr !="nr" & state.next =="nr" & phase.next !="fail")),
                       n.fail = sum(1 * (phase.curr == "fail")),
                       # Numbers of CRs/PRs/remissions
                       n.cr = sum(1 * (state.curr !="cr" & state.next =="cr")),
                       n.pr = sum(1 * (state.curr == "nr" & state.next == "pr")),
                       n.crpr = sum(1 * (state.curr =="nr" & state.next %in% c("cr","pr"))),
                       n.remission = sum(remission * (lag.remission != 1)),
                       # Time to CR/PR/relapse conditional on any CR/PR/relapse
                       mosTo.cr.cond = min(ifelse(state.curr == "cr", time.in*12/dpy, 99999), na.rm = TRUE),
                       mosTo.crpr.cond = min(ifelse(state.curr %in% c("cr", "pr"), time.in*12/dpy, 99999), na.rm = TRUE),
                       mosTo.relapse.cond = min(ifelse(state.curr !="nr" & state.next =="nr" & phase.next !="fail", time.in*12/dpy, 99999), na.rm = TRUE),
                       # Months on treatment
                       mos.onTx = sum(dur*onTx*12/dpy),
                       # Months in CR/PR/NR/remission
                       mos.cr = sum(cr*dur*12/dpy),
                       mos.pr = sum(pr*dur*12/dpy),
                       mos.crpr = sum((cr+pr)*dur*12/dpy),
                       mos.remission = sum(remission*dur*12/dpy),
                       mos.nr = sum(nr*dur*12/dpy),
                       # Life years and QALYs
                       lys.disc = sum(dur*(phase.curr !="death")*disc.hlth.wa/dpy),
                       qaly.disc = sum((dur*util.state*disc.hlth.wa/dpy) + (disutil.ae*disc.hlth.in)),
                       # Completed at least two lines of therapy
                       any.trd = max(trd),
                       # Died during follow-up
                       any.death = max(1 * (phase.curr == "death"))
                     ),
                     by = .(id)]
  
  # Replace 99999 and 0 with NA for conditional measures, create binary flags for
  #   any sae/ae1/.../ae4/remission/CR/PR/maintenance, and calculate total costs
  res.pt[,`:=` (n.relapse.cond = ifelse(near(n.relapse, 0), NA, n.relapse),
                mosTo.cr.cond = ifelse(mosTo.cr.cond > 99998, NA, mosTo.cr.cond),
                mosTo.crpr.cond = ifelse(mosTo.crpr.cond > 99998, NA, mosTo.crpr.cond),
                mosTo.relapse.cond = ifelse(mosTo.relapse.cond > 99998, NA, mosTo.relapse.cond),
                mos.crpr.cond = ifelse(near(mos.crpr, 0), NA, mos.crpr),
                mos.remission.cond = ifelse(near(mos.remission, 0), NA, mos.remission),
                cost.d.disc = cost.d.tx.disc + cost.d.state.disc,
                cost.i.wl.disc = cost.i.abs.disc + cost.i.pres.disc,
                cost.i.disc = cost.i.trnsprt.disc + cost.i.abs.disc + cost.i.pres.disc)]
  res.pt[,`:=` (cost.disc = cost.d.disc + cost.i.disc,
                any.sae = 1 * (n.sae > 0),
                any.ae1 = 1 * (n.ae1 > 0),
                any.ae2 = 1 * (n.ae2 > 0),
                any.ae3 = 1 * (n.ae3 > 0),
                any.ae4 = 1 * (n.ae4 > 0),
                any.cr = 1 * (n.cr > 0),
                any.pr = 1 * (n.pr > 0),
                any.crpr = 1 * (n.crpr > 0),
                any.remission = 1 * (n.remission > 0))]
  
  # Calculate means
  varsToSummarize <- c("any.trd",
                       "any.sae", "any.ae1", "any.ae2", "any.ae3", "any.ae4",
                       "any.cr", "any.pr", "any.crpr", "any.remission",
                       "any.death",
                       "n.sae", "n.ae1", "n.ae2", "n.ae3", "n.ae4",
                       "n.tx",
                       "n.cr", "n.pr", "n.crpr", "n.remission",
                       "n.relapse", "n.relapse.cond", "n.fail",
                       "mosTo.cr.cond", "mosTo.crpr.cond",
                       "mosTo.relapse.cond",
                       "mos.onTx",
                       "mos.nr", "mos.cr", "mos.pr", "mos.crpr", "mos.crpr.cond",
                       "mos.remission", "mos.remission.cond",
                       "lys.disc", "qaly.disc",
                       "cost.disc", "cost.d.disc", "cost.d.tx.disc", "cost.d.state.disc",
                       "cost.i.disc", "cost.i.trnsprt.disc", "cost.i.wl.disc",
                       "cost.i.abs.disc", "cost.i.pres.disc")
  res.tmp3 <- res.pt[, lapply(.SD, mean, na.rm=TRUE), .SDcols=varsToSummarize ]
  
  # Add average cost per QALY and net monetary benefit
  res.tmp3[, `:=`(ace.disc = cost.disc / qaly.disc,
                  ace.d.disc = cost.d.disc / qaly.disc,
                  ace.d.tx.disc = cost.d.tx.disc / qaly.disc,
                  nmb.disc = (qaly.disc * maxWTP) - cost.disc,
                  nmb.d.disc = (qaly.disc * maxWTP) - cost.d.disc,
                  nmb.d.tx.disc = (qaly.disc * maxWTP) - cost.d.tx.disc)]
  
  # If for DSA, limit to the output measures of interest and return
  if(summIn.forDSA == TRUE){
    
    keepMeas = c("lys.disc", "qaly.disc", "cost.disc", "nmb.disc")
    results.forDSA <- res.tmp3[, keepMeas, with=FALSE]
    return(results.forDSA)

  }
  
  # Otherwise...
  
  # Add variable to indicate that these are means for treatment pathway pw
  res.tmp3[, statLabel := paste0("means.txPath", pw)]
  
  #Transpose to have one row per measure
  # res.final <- res.tmp3[, data.table(t(.SD), keep.rownames=TRUE)]
  res.final <- transpose(res.tmp3,
                         keep.names = "Measure.abbrev",
                         make.names = "statLabel")
  
  # Combine simulation data, patient-level data, and overall summary data
  #   into single list for export
  results <- list(res.sim, res.pt, res.final)
  
  # Return patient-level data and overall summary measures
  return(results)
}


