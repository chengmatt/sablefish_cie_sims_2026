#' Threshold Harvest Control Rule
#'
#' Implements a threshold HCR where F is reduced linearly as biomass declines
#' below the target biomass reference point.
#'
#' @param x Current biomass
#' @param frp Fishing mortality reference point (target F)
#' @param brp Biomass reference point (target biomass)
#' @param alpha Minimum biomass threshold (as fraction of brp)
#'
#' @return Fishing mortality rate
HCR_threshold <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp

  if (stock_status >= 1) {
    f <- frp
  } else if (stock_status > alpha) {
    f <- frp * (stock_status - alpha) / (1 - alpha)
  } else {
    f <- 0
  }

  return(f)
}

#' Compute Single-Region Survey Apportionment
#'
#' This function calculates survey-based regional apportionment for a given year
#' by imputing missing survey observations (based on alternating survey timing)
#' and then computing rolling-average apportionment across regions.
#'
#' @param feedback_start_yr Integer. First year of the feedback simulation.
#' @param n_yrs Integer. Total number of years in the simulation.
#' @param y Integer. The current year for which apportionment is computed.
#' @param srv_idx Numeric matrix. Survey index values by region and year.
#' @param rolling_avg_yrs Integer. Column index of the rolling-average year.
#'
#' @return A numeric vector of apportionment proportions across regions.
get_single_region_survey_apportionment <- function(feedback_start_yr, n_yrs, y, srv_idx, rolling_avg_yrs) {
  # get years in simulation
  x <- (feedback_start_yr + 1):n_yrs
  # get goa years (odd)
  odd_yrs <- x[x %% 2 == 1]
  # get bsai years (even)
  even_yrs <- x[x %% 2 == 0]
  # get survey index
  if(y %in% even_yrs) srv_idx[1:2, rolling_avg_yrs] <- srv_idx[1:2, rolling_avg_yrs - 1] # even years, impute bsai
  if(y %in% odd_yrs) srv_idx[3:5, rolling_avg_yrs] <- srv_idx[3:5, rolling_avg_yrs - 1] # odd years, impute goa
  # get rolling average apportionment
  apportionment <- colMeans(srv_idx / rowSums(srv_idx))
  return(apportionment)
}

#' Generate a Forward Population Projection
#'
#' This function prepares terminal-year population, fishery, and biological inputs
#' and then performs a forward population projection over the specified number of
#' years. The projection uses stored numbers-at-age, weights, maturity, selectivity,
#' mortality, movement, and recruitment information to produce projected stock
#' dynamics and fishing mortality outcomes.
#'
#' @param data A named list containing biological and model dimensions, including
#'   regions, ages, sexes, fleets, and corresponding weight, maturity, mortality,
#'   and selectivity arrays.
#' @param rep A model report object containing historical values such as
#'   numbers-at-age, selectivity, movement, fishing mortality, and recruitment.
#' @param y Integer. Terminal model year used as the starting point for projections.
#' @param n_proj_yrs Integer. Number of projection years. Defaults to 2.
#' @param f_ref_pt Numeric. Fishing mortality reference point used in the projection.
#'
#' @return A list or object produced by `Do_Population_Projection`, containing
#'   projected numbers-at-age, fishing mortality, and associated derived quantities.
get_population_projection <- function(data, rep, y, n_proj_yrs = 2, f_ref_pt) {

  # Get inputs for projection
  tmp_terminal_NAA <- array(rep$NAA[,y,,], dim = c(data$n_regions, length(data$ages), data$n_sexes)) # terminal numbers at age
  tmp_terminal_NAA0 <- array(rep$NAA0[,y,,], dim = c(data$n_regions, length(data$ages), data$n_sexes)) # terminal unfished numbers at age
  tmp_WAA <- array(rep(data$WAA[,y,,], each = n_proj_yrs), dim = c(data$n_regions, n_proj_yrs, length(data$ages), data$n_sexes)) # weight at age
  tmp_WAA_fish <- array(rep(data$WAA_fish[,y,,,], each = n_proj_yrs), dim = c(data$n_regions, n_proj_yrs, length(data$ages), data$n_sexes, data$n_fish_fleets)) # weight at age fishery
  tmp_MatAA <- array(rep(data$MatAA[,y,,], each = n_proj_yrs), dim = c(data$n_regions, n_proj_yrs, length(data$ages), data$n_sexes)) # maturity at age
  tmp_fish_sel <- array(rep(rep$fish_sel[,y,,,], each = n_proj_yrs), dim = c(data$n_regions, n_proj_yrs, length(data$ages), data$n_sexes, data$n_fish_fleets)) # selectivity
  tmp_terminal_F <- array(rep$Fmort[,y,], dim = c(data$n_regions, data$n_fish_fleets)) # terminal fishing mortality
  tmp_natmort <- array(rep(rep$natmort[,y,,], each = n_proj_yrs), dim = c(data$n_regions, n_proj_yrs, length(data$ages), data$n_sexes)) # natural mortality
  tmp_recruitment <- array(rep$Rec[,1:y], dim = c(data$n_regions, length(1:y))) # recruitment to use for projections
  tmp_sexratio <- array(replicate(n = n_proj_yrs, rep$sexratio[,y,]), dim = c(data$n_regions, n_proj_yrs, data$n_sexes)) # recruitment sex ratio
  tmp_Movement <- array(dim = c(data$n_regions, data$n_regions, n_proj_yrs, length(data$ages), data$n_sexes))
  for(proj_yr in 1:n_proj_yrs) tmp_Movement[,,proj_yr,,] <- rep$Movement[,,y,,] # Movement projections

  # Do projection to get TAC
  proj <- SPoRC::Do_Population_Projection(
    n_proj_yrs = n_proj_yrs,
    n_regions = data$n_regions,
    n_ages = length(data$ages),
    n_sexes = data$n_sexes,
    sexratio = tmp_sexratio,
    n_fish_fleets = data$n_fish_fleets,
    do_recruits_move = 0,
    recruitment = tmp_recruitment,
    terminal_NAA = tmp_terminal_NAA,
    terminal_NAA0 = tmp_terminal_NAA0,
    terminal_F = tmp_terminal_F,
    natmort = tmp_natmort,
    WAA = tmp_WAA,
    WAA_fish = tmp_WAA_fish,
    MatAA = tmp_MatAA,
    fish_sel = tmp_fish_sel,
    Movement = tmp_Movement,
    f_ref_pt = f_ref_pt,
    b_ref_pt = NULL,
    HCR_function = NULL,
    recruitment_opt = "mean_rec",
    fmort_opt = 'Input',
    t_spawn = sim_env$t_spawn,
    bh_rec_opt = NULL
  )

  return(proj)

}


#' Aggregate Multi-Region Simulation Data Into a Single-Region Dataset
#'
#' Aggregates catches, survey indices, and composition data across regions to
#' create single-region inputs for SPoRC estimation models. The function applies
#' region-weighted compositions, simulates observation error, imputes survey
#' indices for unsampled years under different longline survey (LLS) design
#' types, and writes aggregated values back into `sim_env`.
#'
#' @param sim_data A list of SPoRC-formatted simulated data for year `y` and
#'   simulation `sim`, typically produced by `simulation_data_to_SPoRC()`.
#' @param sim_env A simulation environment list containing operating-model
#'   outputs and containers for aggregated observations.
#' @param y Integer. The terminal year to aggregate (1 to `n_yrs`).
#' @param sim Integer simulation index.
#' @param lls_design_type Character. Longline survey sampling design:
#'   `"current"` for alternating BSAI/GOA sampling, `"historical"` for
#'   alternating BS/AI with annual GOA sampling, or `"all"` to sample all
#'   regions each year.
#'
#' @return Updates `sim_env` in place with aggregated observations:
#'   catches, survey indices, age compositions, length compositions, and ISS
#'   samples.
agg_data_to_single_rg <- function(sim_data, sim_env, y, sim, lls_design_type) {

  # get years in simulation
  x <- (sim_env$feedback_start_yr + 1):sim_env$n_yrs
  odd_yrs <- x[x %% 2 == 1] # bsai years
  even_yrs <- x[x %% 2 == 0] # goa years

  # dimensions
  n_regions <- 1
  n_yrs <- length(1:y)
  n_sexes <- sim_env$n_sexes
  n_ages <- sim_env$n_ages
  n_lens <- sim_env$n_lens
  n_fish_fleets <- sim_env$n_fish_fleets
  n_srv_fleets <- sim_env$n_srv_fleets

  # Catches (using true catches, and then applying error to aggregated true catches)
  if(y == sim_env$feedback_start_yr) {
    aggregated_true_catch <- apply(sim_env$TrueCatch[,1:y,,sim], c(2,3), sum) # get true aggregated catch
    sim_env$Agg_ObsCatch[1,1:y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0, mean(exp(sim_data$ln_sigmaC))))
  } else{
    aggregated_true_catch <- apply(sim_env$TrueCatch[,y,,sim, drop = FALSE], c(2,3), sum) # get true aggregated catch
    sim_env$Agg_ObsCatch[1,1:y,,sim] <- rbind(sim_env$Agg_ObsCatch[1,1:(y-1),,sim], # bind previous catches with current catch
                                              aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0, mean(exp(sim_data$ln_sigmaC)))))
  }

  # Fishery Compositions
  # Get catch weighting
  total_catches <- apply(sim_env$TrueCatch[,1:y,,sim], c(2, 3), sum)  # Sum across regions for each year-fleet
  total_catches <- array(total_catches, dim = c(n_regions, dim(total_catches))) # coerce into correct format
  catch_prop <- sweep(sim_env$TrueCatch[,1:y,,sim], c(2, 3), total_catches, "/") # get catch proportion by region

  # loop through to get catch weighted compositions
  for(yr in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      sim_env$Agg_ObsFishAgeComps[,yr,,,f,sim] <- apply(sim_data$ObsFishAgeComps[,yr,,,f] * catch_prop[,yr,f], c(2,3), sum)
      sim_env$Agg_ObsFishLenComps[,yr,,,f,sim] <- apply(sim_data$ObsFishLenComps[,yr,,,f] * catch_prop[,yr,f], c(2,3), sum)
      sim_env$Agg_ISS_FishAgeComps[,yr,1,f,sim] <- sum(sim_data$ObsFishAgeComps[,yr,,,f])
      sim_env$Agg_ISS_FishLenComps[,yr,1,f,sim] <- sum(sim_data$ObsFishLenComps[,yr,,,f])
    } # end f loop
  } # end yr loop

  # Survey Index
  if(y == sim_env$feedback_start_yr) {
    agg_true_srvidx <- sim_env$TrueSrvIdx[,1:y,,sim] # get true survey index
    # Impute historical years
    agg_true_srvidx[1,31:37,1] <- agg_true_srvidx[1,38,1] # impute most recent sampled point to historical years for bs
    agg_true_srvidx[1,seq(39,63, 2),1] <- agg_true_srvidx[1,seq(39,63, 2)-1,1] # impute most recent sampled point to historical years for bs
    agg_true_srvidx[2,31:36,1] <- agg_true_srvidx[2,37,1] # impute most recent sampled point to historical years for ai
    agg_true_srvidx[2,seq(38,64,2),1] <- agg_true_srvidx[2,seq(38,64,2)-1,1] # impute most recent sampled point to historical years for ai
    agg_imputed_srvidx <- apply(agg_true_srvidx, c(2,3), sum) # aggregate imputed index
    sim_env$Agg_TrueSrvIdx[1,1:y,,sim] <- agg_imputed_srvidx # save "true imputed index"
    sim_env$Agg_ObsSrvIdx[1,1:y,,sim] <- sim_env$Agg_TrueSrvIdx[1,1:y,,sim] * exp(rnorm(length(agg_imputed_srvidx), 0, as.vector(sim_data$ObsSrvIdx_SE) )) # simulate devs
  } else{

    # LLS Design Type = Current Design, Alternate BSAI and GOA
    if(lls_design_type == 'current') {
      if(y == sim_env$feedback_start_yr + 1) {
        bs_true_srvidx <- sim_env$TrueSrvIdx[1,y-2,1,sim] # get 2024 true survey index
        ai_true_srvidx <- sim_env$TrueSrvIdx[2,y-3,1,sim] # get 2023 true survey index
        agg_true_srvidx <- sim_env$TrueSrvIdx[,y,,sim] # get true survey index this year (2025)
        agg_true_srvidx[1:2,1] <- c(bs_true_srvidx, ai_true_srvidx) # impute with most recent years bsai index
      } else{
        prev_true_srvidx <- sim_env$TrueSrvIdx[,y-1,,sim] # get previous true survey index
        agg_true_srvidx <- sim_env$TrueSrvIdx[,y,,sim] # get true survey index
        if(y %in% even_yrs) agg_true_srvidx[1:2,1] <- prev_true_srvidx[1:2,1] # impute previous years value (goa year, impute bsai)
        if(y %in% odd_yrs) agg_true_srvidx[3:5,1] <- prev_true_srvidx[3:5,1] # impute previous years value (bsai year, impute goa)
      }
    }

    # LLS Design Type = Historical Design, Alternate BS and AI, GOA sampled every year
    if(lls_design_type == 'historical') {
      if(y == sim_env$feedback_start_yr + 1) {
        ai_true_srvidx <- sim_env$TrueSrvIdx[2,y-3,1,sim] # get 2023 true survey index
        agg_true_srvidx <- sim_env$TrueSrvIdx[,y,,sim] # get true survey index this year (2025)
        agg_true_srvidx[2,1] <- c(ai_true_srvidx) # impute with most recent years ai index (even year)
      } else{
        prev_true_srvidx <- sim_env$TrueSrvIdx[,y-1,,sim] # get previous true survey index
        agg_true_srvidx <- sim_env$TrueSrvIdx[,y,,sim] # get true survey index
        if(y %in% even_yrs) agg_true_srvidx[2,1] <- prev_true_srvidx[2,1] # impute previous years value (even year, sample bs)
        if(y %in% odd_yrs) agg_true_srvidx[1,1] <- prev_true_srvidx[1,1] # impute previous years value (odd year, sample ai)
      }
    }

    # LLS Design Type = Sample Every Region
    if(lls_design_type == 'all') {
      agg_true_srvidx <- sim_env$TrueSrvIdx[,y,,sim] # get true survey index
    }

    agg_imputed_srvidx <- colSums(agg_true_srvidx) # get aggregated imputed survey (no data in 2025)
    sim_env$Agg_TrueSrvIdx[1,y,,sim] <- agg_imputed_srvidx # save "true imputed index"
    sim_env$Agg_ObsSrvIdx[1,1:y,,sim] <- rbind(sim_env$Agg_ObsSrvIdx[1,1:(y-1),,sim],
                                               sim_env$Agg_TrueSrvIdx[1,y,,sim] * exp(rnorm(length(agg_imputed_srvidx), 0, colMeans(sim_data$ObsSrvIdx_SE[,y,]))))
  }

  # Survey Compositions
  # Survey index weighting
  true_srvidx <- sim_env$TrueSrvIdx[,1:y,,sim] # get true survey index
  total_srv <- apply(true_srvidx, c(2, 3), sum)  # Sum across regions for each year-fleet
  srv_prop <- sweep(true_srvidx, c(2, 3), total_srv, "/") # get catch proportion by region

  # Set historical years at 0s in years not sampled
  om_values$data$UseSrvIdx[,,1]
  srv_prop[1,c(31:37, seq(39, 65, 2)),1] <- 0 # set historical bs years at 0
  srv_prop[2,c(31:36, seq(38, 64, 2)),1] <- 0 # set historical ai years a 0
  srv_prop[,65,1] <- 0 # no data in 2024

  # Loop through to get survey weighted compositions
  for(yr in 1:n_yrs) {

    # make survey comp weighting 0s in years not sampled
    if(yr >= sim_env$feedback_start_yr + 1) {
      if(lls_design_type == 'current') { # current survey design
        if(yr %in% even_yrs) srv_prop[1:2,yr,1] <- 0 # bsai not sampled
        if(yr %in% odd_yrs) srv_prop[3:5,yr,1] <- 0 # goa not sampled
      }
      if(lls_design_type == 'historical') {
        if(yr %in% even_yrs) srv_prop[2,yr,1] <- 0 # ai not sampled
        if(yr %in% odd_yrs) srv_prop[1,yr,1] <- 0 # bs not sampled
      }
      # 'all' design: no zeroing needed
    }

    srv_prop[,yr,1] <- srv_prop[,yr,1] / sum(srv_prop[,yr,1]) # renormalize survey weighting

    for(f in 1:n_srv_fleets) {
      sim_env$Agg_ObsSrvAgeComps[,yr,,,f,sim] <- apply(sim_data$ObsSrvAgeComps[,yr,,,f] * srv_prop[,yr,f], c(2,3), sum)
      sim_env$Agg_ObsSrvLenComps[,yr,,,f,sim] <- apply(sim_data$ObsSrvLenComps[,yr,,,f] * srv_prop[,yr,f], c(2,3), sum)
      # Determine regions to sum (for longline survey)
      if(f == 1 && yr > sim_env$feedback_start_yr + 1) {
        regions <- if(yr %in% odd_yrs) 1:2 else 3:5
      } else {
        regions <- 1:n_regions
      }
      sim_env$Agg_ISS_SrvAgeComps[,yr,1,f,sim] <- sum(sim_data$ObsSrvAgeComps[regions,yr,,,f])
      sim_env$Agg_ISS_SrvLenComps[,yr,1,f,sim] <- sum(sim_data$ObsSrvLenComps[regions,yr,,,f])
    }
  }

}

#' Aggregate Multi-Region Operating Model Data to a Single Region
#'
#' This function aggregates multi-region operating-model outputs (catch,
#' compositions, survey indices, and survey compositions) into a single
#' combined region for use in single-region assessment or feedback simulations.
#' The function applies region-specific weighting, observation-error simulation,
#' and survey-design imputation rules (current, historical, or all-regions-sampled)
#' to produce internally consistent single-region inputs.
#'
#' @param sim_data A list containing simulated observation data, including
#'   fishery and survey age- and length-compositions, observation-error
#'   variances, and survey index standard errors.
#' @param sim_env An environment or list storing truth-level quantities,
#'   including true catches, true survey indices, simulated compositions,
#'   model dimensions, and storage arrays for aggregated outputs.
#' @param y Integer. Current simulation year being aggregated.
#' @param sim Integer. Simulation replicate index.
#' @param lls_design_type Character string specifying the longline survey
#'   design to emulate: `"current"` (alternating BSAI and GOA), `"historical"`
#'   (alternating BS and AI, GOA sampled annually), or `"all"` (all regions
#'   sampled each year).
#'
#' @return The function updates `sim_env` in place with aggregated single-region
#'   catch, composition, and survey quantities. No value is returned.
single_region_use_indicators <- function(y, n_fish_fleets, n_srv_fleets, age_lag = 1) {

  # Data indicators
  usecatch <- array(1, dim = c(1, y, n_fish_fleets)) # use every year
  usesrvidx <- array(0, dim = c(1, y, n_srv_fleets))
  usefishage <- array(0, c(1, y, n_fish_fleets))
  usefishlen <- array(0, c(1, y, n_fish_fleets))
  usesrvage <- array(0, c(1, y, n_srv_fleets))
  usesrvlen <- array(0, c(1, y, n_srv_fleets)) # not used at all

  # Fixed-gear fleet
  usefishage[,40:(y-age_lag),1] <- 1 # use data starting year 40, with age lag
  usefishlen[,31:39,1] <- 1 # use data only from 31 - 39

  # Trawl gear fleet (no ages)
  usefishlen[,31:y,2] <- 1 # use data starting year 31, without any lags

  # Longline Survey (Domestic)
  usesrvage[,37:(y-age_lag),1] <- 1 # use data starting year 37, with age lag
  usesrvage[,65,1] <- 0 # no survey in 2024
  usesrvidx[,31:y,1] <- 1 # use data starting year 31
  usesrvidx[,65,1] <- 0 # no survey in 2024

  # Longline Survey (JP)
  usesrvage[,seq(26,34,2),3] <- 1 # use data in select years
  usesrvidx[,20:35,3] <- 1  # use data in select years

  return(list(usefishage = usefishage, usefishlen = usefishlen,
              usesrvage = usesrvage, usesrvlen = usesrvlen,
              usecatch = usecatch, usesrvidx = usesrvidx))

}

#' Run a Single-Region Closed-Loop Simulation Iteration
#'
#' This function executes one full closed-loop simulation for a single replicate
#' under a single-region assessment framework. For each simulation year, the
#' function updates population dynamics, performs stock assessment (when within
#' the feedback period), computes reference points, applies a harvest control
#' rule, conducts short-term projections to derive a TAC, apportions that TAC
#' across fleets, and solves for the fishing mortality rates needed to achieve
#' those catches in the operating model.
#'
#' @param sim_env An environment or list containing all operating-model state
#'   variables, true quantities, storage arrays, and simulation settings.
#' @param sim Integer. Simulation replicate index to run.
#' @param fleet_allocation Numeric vector giving the proportion of regional TAC
#'   assigned to each fishery fleet.
#' @param lls_design_type Character string specifying the longline survey design:
#'   `"current"`, `"historical"`, or `"all"`. Passed to the assessment model
#'   builder for appropriate survey imputation logic.
#'
#' @return The function updates `sim_env` in place (population states, catches,
#'   F rates, and assessment outputs). No value is returned.
run_single_rg_closedloop_i <- function(sim_env, sim, fleet_allocation, lls_design_type) {

  # Run Closed Loop ---------------------------------------------------------
  for(y in 1:sim_env$n_yrs) {

    # Execute annual population dynamics
    run_annual_cycle(y, sim, sim_env)

    # Run feedback
    if(y >= sim_env$feedback_start_yr) {

      ### Assessment --------------------------------------------------------------
      asmt_list <- single_region_em(sim_env = sim_env, y = y, sim = sim, srv_idx_se = 0.2, age_lag = 1, lls_design_type) # get assessment data
      model <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = TRUE) # get model

      ### Reference Points --------------------------------------------------------
      ref_pts <- get_closed_loop_reference_points(
        use_true_values = FALSE,
        sim_env = sim_env,
        asmt_data = asmt_list$data,
        asmt_rep = model$rep,
        y = y,
        sim = sim,

        # single region reference points
        reference_points_opt = list(
          n_avg_yrs = 1,
          SPR_x = 0.4,
          calc_rec_st_yr = 3,
          rec_age = 4,
          type = 'single_region',
          what = "SPR",
          B_x = 0.4
        ),
        n_proj_yrs = 2
      )

      ### HCR ---------------------------------------------------------------------
      agg_ssb <- model$rep$SSB[,y] # get ssb
      exp_frate <- HCR_threshold(x = agg_ssb, frp = ref_pts$f_ref_pt[,2], brp = ref_pts$b_ref_pt[,2], alpha = 0.05) # get prescribed f

      ### Projection --------------------------------------------------------------
      projection <- get_population_projection(data = asmt_list$data, rep = model$rep, y = y, n_proj_yrs = 2, f_ref_pt = ref_pts$f_ref_pt) # get 1 year projection
      tac <- sum(projection$proj_Catch[,2,]) # get tac

      ### Apportionment -----------------------------------------------------------
      # get survey apportionment (survey biomass)
      tmp_Srv_BiomIA <- sim_env$SrvIAA[,((y - 4): y),,,1,sim] * sim_env$WAA_srv[,((y - 4): y),,,1,sim]
      tmp_srvidx_biom <- apply(tmp_Srv_BiomIA, 1:2, sum) # get biomass
      apportionment <- get_single_region_survey_apportionment(feedback_start_yr = sim_env$feedback_start_yr,
                                                              n_yrs = sim_env$n_yrs, y = y,
                                                              srv_idx = tmp_srvidx_biom, # using true survey biomass values to do apportionment
                                                              rolling_avg_yrs = 5)
      tac_r <- tac * apportionment # get regional tac
      tac_rf <- tac_r * fleet_allocation # allocate regional tac by fleet

      ### TAC to F ----------------------------------------------------------------
      if (y < sim_env$n_yrs) {

        rf_grid <- expand.grid(r = seq_len(sim_env$n_regions), f = seq_len(sim_env$n_fish_fleets)) # set up region, fleet grid to bisection across
        tmp_f <- mapply(function(r, f) { # do bisection to go from region and fleet specific catch to region and fleet specific F rates
          bisection_F(
            f_guess = 0.05, # guess for fishing mortality rate
            catch = tac_rf[r,f], # catch values to use
            NAA = sim_env$NAA[r, y+1, , , sim], # numbers at age in simulation (truth)
            WAA = sim_env$WAA[r, y+1, , , sim], # weight-at-age in simulation (truth)
            natmort  = sim_env$natmort[r, y+1, , , sim], # natural mortality in simulation (truth)
            fish_sel = sim_env$fish_sel[r, y+1, , , f, sim] # fishery selectivity in simulation (truth)
          )
        }, r = rf_grid$r, f = rf_grid$f)

        # assign bisection values back into simulation
        sim_env$Fmort[,y+1,,sim] <- array(tmp_f, dim = c(sim_env$n_regions, sim_env$n_fish_fleets))
      }

    } # end feedback

  } # end y loop

}

#' Run Single-Region Closed-Loop Simulations in Parallel
#'
#' This function executes multiple single-region closed-loop simulations in
#' parallel using `future` and `furrr`. Each simulation replicate is run with
#' its own index and the results are merged back into the shared simulation
#' environment. Arrays whose final dimension corresponds to simulation number
#' are updated appropriately.
#'
#' @param sim_env An environment or list containing all operating-model truth,
#'   state variables, and storage arrays modified during simulation.
#' @param n_sims Integer. Number of simulation replicates to run.
#' @param fleet_allocation Numeric vector giving TAC allocation fractions across
#'   fleets for each region.
#' @param lls_design_type Character string specifying the longline survey design:
#'   `"current"`, `"historical"`, or `"all"`. Passed to the inner simulation.
#' @param n_cores Integer. Number of parallel workers to use.
#'
#' @return The updated `sim_env` object containing results across all simulation
#'   replicates. The function modifies `sim_env` in place and also returns it.
run_single_rg_closedloop_parallel <- function(sim_env, n_sims, fleet_allocation, lls_design_type, n_cores) {

  plan(multisession, workers = n_cores)
  options(future.globals.maxSize = 5e9)

  # run in parrallel and return simulation environment
  with_progress({
    env_list <- future_map(
      1:n_sims,
      ~{
        run_single_rg_closedloop_i(sim_env, .x, fleet_allocation, lls_design_type)
        sim_env
      },
      .progress = TRUE
    )
  }, handlers = progressr::handler_progress(format = "[:bar] :percent"))

  # Merge results back in
  for(i in 1:n_sims) {
    for(var_name in ls(sim_env, all.names = TRUE)) { # loop through sim_env to get variable names
      arr <- env_list[[i]][[var_name]] # get array
      if(is.array(arr)) { # check if array
        ndim <- length(dim(arr)) # check array dimensions
        last_dim <- dim(arr)[ndim] # get last dimension
        if(last_dim == n_sims) { # if last dimension matches number of sims
          comma_str <- paste(rep(",", ndim - 1), collapse = "") # build comma structure, e.g., 3d array gives ,,,
          expr <- paste0("sim_env[[\"", var_name, "\"]][", comma_str, "i] <- env_list[[i]][[\"", var_name, "\"]][", comma_str, "i]") # write expression for array
          eval(parse(text = expr)) # parse expression
        }
      }
    }
  }
  return(sim_env)
}

#' Add Aggregated Observation Objects to the Simulation Environment
#'
#' This function initializes storage arrays within the simulation environment
#' for aggregated (single-region) catches, survey indices, and age- and
#' length-composition data. These arrays are used when collapsing multi-region
#' quantities into a single aggregated region during closed-loop simulations.
#'
#' @param sim_env An environment or list containing simulation dimensions
#'   (years, fleets, regions, ages, sexes, lengths) and operating-model
#'   state/observation objects. The object is modified in place.
#'
#' @return The updated `sim_env` object with newly allocated aggregated arrays.
#'   The function modifies `sim_env` and also returns it.
add_aggregated_obj_to_simenv <- function(sim_env) {
  # allow for aggregated catches, indices, and compositions
  sim_env$Agg_ObsCatch <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_fish_fleets, n_sims))
  sim_env$Agg_ObsSrvIdx <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_srv_fleets, n_sims))
  sim_env$Agg_TrueSrvIdx <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_srv_fleets, n_sims))
  sim_env$Agg_ObsFishAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
  sim_env$Agg_ObsFishLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
  sim_env$Agg_ObsSrvAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_ages, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
  sim_env$Agg_ObsSrvLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
  sim_env$Agg_ISS_FishAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
  sim_env$Agg_ISS_FishLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
  sim_env$Agg_ISS_SrvAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
  sim_env$Agg_ISS_SrvLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
  return(sim_env)
}
