#' Solve for fishing mortality rates that achieve target catches for multiple fleets
#'
#' @param target_catch Numeric vector of target catch values for each fleet
#' @param NAA Matrix of numbers-at-age (ages x sexes)
#' @param WAA Matrix of weight-at-age (ages x sexes)
#' @param natmort Matrix of natural mortality (ages x sexes)
#' @param fish_sel 3D array of fishery selectivity (ages x sexes x fleets)
#' @param f_init Initial guess for F values (scalar or vector)
#' @param control List of control parameters for nleqslv
#' @return Numeric vector of F values for each fleet
solve_multifleet_F <- function(target_catch, NAA, WAA, natmort, fish_sel,
                               f_init = 0.05, control = list(btol = 1e-6)) {

  n_fleets <- length(target_catch)

  # Expand f_init if scalar
  if(length(f_init) == 1) f_init <- rep(f_init, n_fleets)

  # Function to minimize: difference between predicted and target catch for all fleets
  catch_diff <- function(f_vec) {
    pred_catches <- numeric(n_fleets)

    for(f in 1:n_fleets) {
      # F-at-age for this fleet
      FAA <- f_vec[f] * fish_sel[, , f]

      # Total Z includes F from ALL fleets
      ZAA_total <- natmort
      for(ff in 1:n_fleets) {
        ZAA_total <- ZAA_total + f_vec[ff] * fish_sel[, , ff]
      }

      # Predicted catch for this fleet (Baranov catch equation)
      pred_catches[f] <- sum((FAA / ZAA_total * NAA * (1 - exp(-ZAA_total))) * WAA)
    }

    return(pred_catches - target_catch)  # Difference from target
  }

  # Solve for F vector
  result <- nleqslv::nleqslv(f_init, catch_diff, control = control)

  return(result$x)
}

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
#' @param lls_design_type Either historical or current design or all (if all, doesn't do any imputing)
#'
#' @return A numeric vector of apportionment proportions across regions.
get_single_region_survey_apportionment <- function(feedback_start_yr, n_yrs, y, srv_idx, rolling_avg_yrs, lls_design_type) {
  # get years in simulation
  x <- (feedback_start_yr + 1):n_yrs
  # get goa years (odd)
  odd_yrs <- x[x %% 2 == 1]
  # get bsai years (even)
  even_yrs <- x[x %% 2 == 0]

  # get survey index
  if(lls_design_type == 'current') {
    if(y %in% even_yrs) srv_idx[1:2, rolling_avg_yrs] <- srv_idx[1:2, rolling_avg_yrs - 1] # even years, impute bsai
    if(y %in% odd_yrs) srv_idx[3:5, rolling_avg_yrs] <- srv_idx[3:5, rolling_avg_yrs - 1] # odd years, impute goa
  }

  if(lls_design_type == 'historical') {
    if(y %in% even_yrs) srv_idx[2, rolling_avg_yrs] <- srv_idx[2, rolling_avg_yrs - 1] # even years, impute ai
    if(y %in% odd_yrs) srv_idx[1, rolling_avg_yrs] <- srv_idx[1, rolling_avg_yrs - 1] # odd years, impute bs
  }

  # get rolling average apportionment
  prop_yr <- sweep(srv_idx, 2, colSums(srv_idx, na.rm = TRUE), FUN = "/")
  apportionment <- rowMeans(prop_yr, na.rm = TRUE)
  return(apportionment)
}

#' Compute Three-Region Model Survey Apportionment
#'
#' This function calculates survey-based regional apportionment for a given year
#' by imputing missing survey observations (based on alternating survey timing)
#' and then computing rolling-average apportionment across the GOA regions.
#'
#' @param feedback_start_yr Integer. First year of the feedback simulation.
#' @param n_yrs Integer. Total number of years in the simulation.
#' @param y Integer. The current year for which apportionment is computed.
#' @param srv_idx Numeric matrix. Survey index values by region and year.
#' @param rolling_avg_yrs Integer. Column index of the rolling-average year.
#' @param lls_design_type Either historical or current design or all (if all or historical, doesn't do any imputing for the GOA since sampled every year)
#'
#' @return A numeric vector of apportionment proportions across the GOA regions.
get_three_rg_survey_apportionment <- function(feedback_start_yr, n_yrs, y, srv_idx, rolling_avg_yrs, lls_design_type) {

  # get years in simulation
  x <- (feedback_start_yr + 1):n_yrs
  # get goa years (odd)
  odd_yrs <- x[x %% 2 == 1]
  # get bsai years (even)
  even_yrs <- x[x %% 2 == 0]

  srv_idx_3 <- srv_idx[3:5, , drop = FALSE]

  # get survey index
  if(lls_design_type == 'current') {
    if(y %in% odd_yrs) srv_idx_3[,rolling_avg_yrs] <- srv_idx_3[,rolling_avg_yrs - 1] # odd years, impute goa
  }

  # get rolling average apportionment
  prop_yr <- sweep(srv_idx_3, 2, colSums(srv_idx_3, na.rm = TRUE), FUN = "/")
  apportionment <- rowMeans(prop_yr, na.rm = TRUE)
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
#' @param srv_wgt 'numbers', 'biomass', or 'eqwt', weighting for survey comps (default, numbers)
#' @param fish_wgt numbers', 'biomass', or 'eqwt' weighting for fishery comps (default, biomass)
#' @param srv_idx_se Survey index se to simulate from
#'
#' @return Updates `sim_env` in place with aggregated observations:
#'   catches, survey indices, age compositions, length compositions, and ISS
#'   samples.
agg_data_to_single_rg <- function(sim_data, sim_env, y, sim, lls_design_type, srv_wgt = 'numbers', fish_wgt = 'biomass', srv_idx_se) {

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
    sim_env$Agg_ObsCatch[1,y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0,
                                                                                mean(exp(sim_data$ln_sigmaC))))
  }

  # Fishery Compositions
  # Get catch weighting
  if(fish_wgt == 'biomass') {
    total_catches <- apply(sim_env$TrueCatch[,1:y,,sim], c(2, 3), sum)  # Sum across regions for each year-fleet
    total_catches <- array(total_catches, dim = c(n_regions, dim(total_catches))) # coerce into correct format
    catch_prop <- sweep(sim_env$TrueCatch[,1:y,,sim], c(2, 3), total_catches, "/") # get catch proportion by region
  }

  if(fish_wgt == 'numbers') {
    catch_prop <- apply(sim_env$CAA[,1:y,,,,sim], c(1,2,5), sum)
    total_by_year <- apply(catch_prop, c(2,3), sum)
    catch_prop <- sweep(catch_prop, c(2, 3), total_by_year, "/") # get catch proportion by region
  }

  if(fish_wgt == 'eqwt') {
    # get catches to look for 0s
    catch_array <- sim_env$TrueCatch[, 1:y, , sim]
    # fxn to assign equal proportions to non-zero entries
    equal_nonzero <- function(x) {
      nz <- x != 0
      n_nz <- sum(nz)
      if (n_nz == 0) return(x)  # all zeros, leave as is
      x[nz] <- 1 / n_nz
      x
    }
    catch_prop <- apply(catch_array, c(2,3), equal_nonzero)
  }

  # loop through to get catch weighted compositions
  for(yr in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      # figure out which regions have catches
      regions <- (1:sim_env$n_regions)[-which(catch_prop[,yr,f] == 0)]
      if(length(regions) == 0) regions <- 1:sim_env$n_regions

      # age comps
      regional_props_age <- sweep(sim_data$ObsFishAgeComps[,yr,,,f], 1,
                                  apply(sim_data$ObsFishAgeComps[,yr,,,f], 1, sum), "/")
      agg_age <- apply(regional_props_age * catch_prop[,yr,f], c(2,3), sum)
      sim_env$Agg_ObsFishAgeComps[,yr,,,f,sim] <- agg_age / sum(agg_age)

      # Calculate n_eff for age comps
      regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,f], 1, sum)
      regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,f], 1, regional_n_age, "/")
      pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions,yr,f], "*")) / sum(catch_prop[,yr,f])
      numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
      denominator_age <- sum(catch_prop[regions,yr,f]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
      sim_env$Agg_ISS_FishAgeComps[,yr,1,f,sim] <- numerator_age / denominator_age

      # length comps
      regional_props_len <- sweep(sim_data$ObsFishLenComps[,yr,,,f], 1,
                                  apply(sim_data$ObsFishLenComps[,yr,,,f], 1, sum), "/")
      agg_len <- apply(regional_props_len * catch_prop[,yr,f], c(2,3), sum)
      sim_env$Agg_ObsFishLenComps[,yr,,,f,sim] <- agg_len / sum(agg_len)

      # Calculate n_eff for length comps
      regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,f], 1, sum)
      regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,f], 1, regional_n_len, "/")
      pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions,yr,f], "*")) / sum(catch_prop[,yr,f])
      numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
      denominator_len <- sum(catch_prop[regions,yr,f]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
      sim_env$Agg_ISS_FishLenComps[,yr,1,f,sim] <- numerator_len / denominator_len

    } # end f loop
  } # end yr loop

  # Survey Index
  if(y == sim_env$feedback_start_yr) {
    agg_true_srvidx <- sim_env$TrueSrvIdx[,1:y,,sim] # get true survey index
    # Impute historical years
    if(lls_design_type %in% c("current", "historical")) {
      agg_true_srvidx[1,31:37,1] <- agg_true_srvidx[1,38,1] # impute most recent sampled point to historical years for bs
      agg_true_srvidx[1,seq(39,63, 2),1] <- agg_true_srvidx[1,seq(39,63, 2)-1,1] # impute most recent sampled point to historical years for bs
      agg_true_srvidx[2,31:36,1] <- agg_true_srvidx[2,37,1] # impute most recent sampled point to historical years for ai
      agg_true_srvidx[2,seq(38,64,2),1] <- agg_true_srvidx[2,seq(38,64,2)-1,1] # impute most recent sampled point to historical years for ai
    }
    agg_imputed_srvidx <- apply(agg_true_srvidx, c(2,3), sum) # aggregate imputed index
    sim_env$Agg_TrueSrvIdx[1,1:y,,sim] <- agg_imputed_srvidx # save "true imputed index"
    sim_env$Agg_ObsSrvIdx[1,1:y,,sim] <- sim_env$Agg_TrueSrvIdx[1,1:y,,sim] * exp(rnorm(length(agg_imputed_srvidx), 0, srv_idx_se )) # simulate devs
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
    sim_env$Agg_ObsSrvIdx[1,y,,sim] <- sim_env$Agg_TrueSrvIdx[1,y,,sim] * exp(rnorm(length(agg_imputed_srvidx), 0, srv_idx_se))
  }

  # Survey Compositions
  # Survey index weighting
  if(srv_wgt == 'numbers') {
    true_srvidx <- sim_env$TrueSrvIdx[,1:y,,sim] # get true survey index
    total_srv <- apply(true_srvidx, c(2, 3), sum)  # Sum across regions for each year-fleet
    srv_prop <- sweep(true_srvidx, c(2, 3), total_srv, "/") # get catch proportion by region
  }

  if(srv_wgt == 'biomass') {
    true_srvidx <- apply(sim_env$SrvIAA[,1:y,,,,sim] * sim_env$WAA_srv[,1:y,,,,sim], c(1,2,5), sum) # get true survey index
    total_srv <- apply(true_srvidx, c(2, 3), sum)  # Sum across regions for each year-fleet
    srv_prop <- sweep(true_srvidx, c(2, 3), total_srv, "/") # get catch proportion by region
  }

  if(srv_wgt == 'eqwt') {
    srv_prop <- array(1/sim_env$n_regions, dim = c(5, length(1:y), n_srv_fleets)) # replace all with 1s
  }

  # Set historical years at 0s in years not sampled
  if(lls_design_type %in% c("current", 'historical')){
    srv_prop[1,c(31:37, seq(39, 65, 2)),1] <- 0 # set historical bs years at 0
    srv_prop[2,c(31:36, seq(38, 64, 2)),1] <- 0 # set historical ai years a 0
  }
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
      regional_props_srv <- sweep(sim_data$ObsSrvAgeComps[,yr,,,f], 1,
                                  apply(sim_data$ObsSrvAgeComps[,yr,,,f], 1, sum), "/")
      agg_srv <- apply(regional_props_srv * srv_prop[,yr,f], c(2,3), sum)
      sim_env$Agg_ObsSrvAgeComps[,yr,,,f,sim] <- agg_srv / sum(agg_srv)

      regions <- (1:sim_env$n_regions)[-which(srv_prop[,yr,f] == 0)]
      if(length(regions) == 0) regions <- 1:sim_env$n_regions

      # get neff
      regional_n_age <- apply(sim_data$ObsSrvAgeComps[regions,yr,,,f], 1, sum)
      regional_p_age <- sweep(sim_data$ObsSrvAgeComps[regions,yr,,,f], 1, regional_n_age, "/")
      pooled_comp_age <- colSums(sweep(regional_p_age, 1, srv_prop[regions,yr,f], "*")) / sum(srv_prop[regions,yr,f])
      numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
      denominator_age <- sum(srv_prop[regions,yr,f]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
      sim_env$Agg_ISS_SrvAgeComps[,yr,1,f,sim] <- numerator_age / denominator_age

      # note: no survey lengths used
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
#' @param srv_idx_se Survey idnex SE
#' @param age_lag Age lag
#' @param srv_wgt numbers or biomass comp weighting
#' @param fish_wgt numbers or biomass comp weighting
#'
#' @return The updated `sim_env` object containing results across all simulation
#'   replicates. The function modifies `sim_env` in place and also returns it.
run_single_rg_closedloop_parallel <- function(sim_env, n_sims, fleet_allocation, srv_idx_se = 0.2, lls_design_type,
                                              age_lag = 1, srv_wgt = 'numbers', fish_wgt = 'numbers', n_cores) {


  plan(multisession, workers = n_cores)
  options(future.globals.maxSize = 8e9)
  handlers(handler_progress(format = "[:bar] :percent"))

  # storage containers for models
  n_feedback_yrs <- sim_env$n_yrs - sim_env$feedback_start_yr + 1
  sim_env$models <- vector("list", n_feedback_yrs)
  for(i in seq_len(n_feedback_yrs)) sim_env$models[[i]] <- vector("list", n_sims)

  # run in parrallel and return simulation environment
  with_progress({
    env_list <- future_map(
      1:n_sims,
      ~{
        run_single_rg_closedloop_i(sim_env, .x, fleet_allocation,
                                   lls_design_type, srv_idx_se,
                                   age_lag, srv_wgt, fish_wgt)
        sim_env
      },
      .progress = TRUE,
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

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

    # merge model results back in
    for(yr_idx in 1:n_feedback_yrs) {
      sim_env$models[[yr_idx]][[i]] <- env_list[[i]]$models[[yr_idx]][[i]]
    }

  }


  return(sim_env)
}

#' Run Three-Region Closed-Loop Simulations in Parallel
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
#' @param srv_idx_se Survey idnex SE
#' @param age_lag Age lag
#' @param srv_wgt numbers or biomass comp weighting
#' @param fish_wgt numbers or biomass comp weighting
#'
#' @return The updated `sim_env` object containing results across all simulation
#'   replicates. The function modifies `sim_env` in place and also returns it.
run_three_rg_closedloop_parallel <- function(sim_env, n_sims, fleet_allocation, lls_design_type, srv_idx_se,
                                             age_lag, srv_wgt, fish_wgt, n_cores) {


  plan(multisession, workers = n_cores)
  options(future.globals.maxSize = 8e9)
  handlers(handler_progress(format = "[:bar] :percent"))

  # storage containers for models
  n_feedback_yrs <- sim_env$n_yrs - sim_env$feedback_start_yr + 1
  sim_env$models <- vector("list", n_feedback_yrs)
  for(i in seq_len(n_feedback_yrs)) sim_env$models[[i]] <- vector("list", n_sims)

  # run in parallel and return simulation environment
  with_progress({
    env_list <- future_map(
      1:n_sims,
      ~{
        tryCatch({
          run_three_rg_closedloop_i(sim_env, .x, fleet_allocation,
                                    lls_design_type, srv_idx_se,
                                    age_lag, srv_wgt, fish_wgt)
          sim_env  # return on success
        }, error = function(e) {
          warning(paste0("Simulation ", .x, " failed: ", e$message))
          NULL  # return NULL on failure
        })
      },
      .progress = TRUE
    )
  })

  # Merge results back in (skip failed sims)
  for(i in 1:n_sims) {
    # Skip if this sim failed
    if(is.null(env_list[[i]])) next

    # Merge arrays
    for(var_name in ls(sim_env, all.names = TRUE)) {
      arr <- env_list[[i]][[var_name]]
      if(is.array(arr)) {
        ndim <- length(dim(arr))
        last_dim <- dim(arr)[ndim]
        if(last_dim == n_sims) {
          comma_str <- paste(rep(",", ndim - 1), collapse = "")
          expr <- paste0("sim_env[[\"", var_name, "\"]][", comma_str, "i] <- env_list[[i]][[\"", var_name, "\"]][", comma_str, "i]")
          eval(parse(text = expr))
        }
      }
    }

    # Merge model results back in
    for(yr_idx in 1:n_feedback_yrs) {
      sim_env$models[[yr_idx]][[i]] <- env_list[[i]]$models[[yr_idx]][[i]]
    }

  }

  return(sim_env)
}


#' Run Three-Region Closed-Loop Simulations in Parallel
#'
#' This function executes multiple single-region closed-loop simulations in
#' parallel using `future` and `furrr`. Each simulation replicate is run with
#' its own index and the results are merged back into the shared simulation
#' environment. Arrays whose final dimension corresponds to simulation number
#' are updated appropriately.
#'
#' @param sim_env A list or environment containing the simulation environment, including population dynamics,
#'   survey data, movement matrices, natural mortality, and other model parameters.
#' @param fleet_allocation Numeric vector. Proportions used to allocate regional TAC to fleets.
#' @param lls_design_type Character. Type of longline survey design used for apportionment.
#' @param srv_idx_se Numeric. Standard error of survey indices.
#' @param age_lag Integer. Age lag between recruitment and survey observation.
#' @param srv_wgt Character Weighting applied to survey data (e.g., "numbers" or "biomass").
#' @param fish_wgt Character Weighting applied to fishery data (e.g., "numbers" or "biomass").
#' @param faa_n_fish_fleets Integer. Number of fishing fleets in the feedback assessment.
#' @param faa_n_srv_fleets Integer. Number of survey fleets in the feedback assessment. Typically 4.
#' @param fish_sel_model Character vector. Fleet-specific selectivity models for the fishing fleets.
#' @param srv_sel_model Character vector. Fleet-specific selectivity models for survey fleets.
#' @param fish_selex_prior LPriors for fishing fleet selectivity parameters. Constructed based on
#'   fleet blocks, sex, and region.
#' @param srv_selex_prior  Priors for survey fleet selectivity parameters.
#' @param n_cores Number of cores
#' @param n_sims Number of simulations
#' @param fish_sel_blocks Fishery blocks
#' @param srv_sel_blocks Survey blocks
#'
#' @return The updated `sim_env` object containing results across all simulation
#'   replicates. The function modifies `sim_env` in place and also returns it.
run_faa_closedloop_parallel <- function(sim_env, n_sims, fleet_allocation,
                                        lls_design_type, srv_idx_se,
                                        age_lag, srv_wgt,
                                        fish_wgt, faa_n_fish_fleets,
                                        faa_n_srv_fleets,
                                        fish_sel_model,
                                        srv_sel_model,
                                        fish_selex_prior,
                                        srv_selex_prior,
                                        fish_sel_blocks,
                                        srv_sel_blocks,
                                        n_cores) {


  plan(multisession, workers = n_cores)
  options(future.globals.maxSize = 8e9)
  handlers(handler_progress(format = "[:bar] :percent"))

  # storage containers for models
  n_feedback_yrs <- sim_env$n_yrs - sim_env$feedback_start_yr + 1
  sim_env$models <- vector("list", n_feedback_yrs)
  for(i in seq_len(n_feedback_yrs)) sim_env$models[[i]] <- vector("list", n_sims)

  # run in parallel and return simulation environment
  with_progress({
    env_list <- future_map(
      1:n_sims,
      ~{
        tryCatch({
          run_faa_closedloop_i(sim_env, .x,
                               fleet_allocation,
                               lls_design_type,
                               srv_idx_se,
                               age_lag,
                               srv_wgt,
                               fish_wgt,
                               faa_n_fish_fleets,
                               faa_n_srv_fleets,
                               fish_sel_model,
                               srv_sel_model,
                               fish_selex_prior,
                               srv_selex_prior,
                               fish_sel_blocks,
                               srv_sel_blocks
                               )

          sim_env  # return on success
        }, error = function(e) {
          warning(paste0("Simulation ", .x, " failed: ", e$message))
          NULL  # return NULL on failure
        })
      },
      .progress = TRUE,
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

  # Merge results back in (skip failed sims)
  for(i in 1:n_sims) {
    # Skip if this sim failed
    if(is.null(env_list[[i]])) next

    # Merge arrays
    for(var_name in ls(sim_env, all.names = TRUE)) {
      arr <- env_list[[i]][[var_name]]
      if(is.array(arr)) {
        ndim <- length(dim(arr))
        last_dim <- dim(arr)[ndim]
        if(last_dim == n_sims) {
          comma_str <- paste(rep(",", ndim - 1), collapse = "")
          expr <- paste0("sim_env[[\"", var_name, "\"]][", comma_str, "i] <- env_list[[i]][[\"", var_name, "\"]][", comma_str, "i]")
          eval(parse(text = expr))
        }
      }
    }


    # Merge model results back in
    for(yr_idx in 1:n_feedback_yrs) {
      sim_env$models[[yr_idx]][[i]] <- env_list[[i]]$models[[yr_idx]][[i]]
    }

  }


  return(sim_env)
}


#' Run Five-Region Closed-Loop Simulations in Parallel
#'
#' This function executes multiple single-region closed-loop simulations in
#' parallel using `future` and `furrr`. Each simulation replicate is run with
#' its own index and the results are merged back into the shared simulation
#' environment. Arrays whose final dimension corresponds to simulation number
#' are updated appropriately.
#'
#' @param sim_env A list or environment containing the simulation environment, including population dynamics,
#'   survey data, movement matrices, natural mortality, and other model parameters.
#' @param fleet_allocation Numeric vector. Proportions used to allocate regional TAC to fleets.
#' @param n_cores Number of cores
#' @param n_sims Number of simulations
#' @param hcr_type How HCR is implemented: global_ssb_global_b40, local_ssb_local_b40, local_ssb_global_b40
#'
#' @return The updated `sim_env` object containing results across all simulation
#'   replicates. The function modifies `sim_env` in place and also returns it.
run_five_rg_closedloop_parallel <- function(sim_env, n_sims, fleet_allocation, hcr_type, n_cores) {


  plan(multisession, workers = n_cores)
  options(future.globals.maxSize = 8e9)
  handlers(handler_progress(format = "[:bar] :percent"))

  # run in parrallel and return simulation environment
  with_progress({
    env_list <- future_map(
      1:n_sims,
      ~{
        run_five_rg_closedloop_i(sim_env, .x,
                                 fleet_allocation,
                                 hcr_type)
        sim_env
      },
      .progress = TRUE,
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

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

  # Merge model results back in
  for(i in 1:n_sims) sim_env$models[[i]] <- env_list[[i]]$models[[i]]

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
#' @param type
#' @param faa_n_fish_fleets
#' @param faa_n_srv_fleets
#'
#' @return The updated `sim_env` object with newly allocated aggregated arrays.
#'   The function modifies `sim_env` and also returns it.
add_aggregated_obj_to_simenv <- function(sim_env, type = 'sgl_rg', faa_n_fish_fleets, faa_n_srv_fleets) {

  if(type == 'sgl_rg') {
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
  }

  if(type == 'faa') {
    # allow for aggregated catches, indices, and compositions
    sim_env$Agg_ObsCatch <- array(NA, dim = c(1, sim_env$n_yrs, faa_n_fish_fleets, n_sims))
    sim_env$Agg_ObsSrvIdx <- array(NA, dim = c(1, sim_env$n_yrs, faa_n_srv_fleets, n_sims))
    sim_env$Agg_TrueSrvIdx <- array(NA, dim = c(1, sim_env$n_yrs, faa_n_srv_fleets, n_sims))
    sim_env$Agg_ObsFishAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_ages, sim_env$n_sexes, faa_n_fish_fleets, n_sims))
    sim_env$Agg_ObsFishLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_lens, sim_env$n_sexes, faa_n_fish_fleets, n_sims))
    sim_env$Agg_ObsSrvAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_ages, sim_env$n_sexes, faa_n_srv_fleets, n_sims))
    sim_env$Agg_ObsSrvLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_lens, sim_env$n_sexes, faa_n_srv_fleets, n_sims))
    sim_env$Agg_ISS_FishAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, faa_n_fish_fleets, n_sims))
    sim_env$Agg_ISS_FishLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, faa_n_fish_fleets, n_sims))
    sim_env$Agg_ISS_SrvAgeComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, faa_n_srv_fleets, n_sims))
    sim_env$Agg_ISS_SrvLenComps <- array(NA, dim = c(1, sim_env$n_yrs, sim_env$n_sexes, faa_n_srv_fleets, n_sims))
  }

  if(type == 'three_rg') {
    # allow for aggregated catches, indices, and compositions
    sim_env$Agg_ObsCatch <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_fish_fleets, n_sims))
    sim_env$Agg_ObsSrvIdx <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_srv_fleets, n_sims))
    sim_env$Agg_TrueSrvIdx <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_srv_fleets, n_sims))
    sim_env$Agg_ObsFishAgeComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
    sim_env$Agg_ObsFishLenComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
    sim_env$Agg_ObsSrvAgeComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_ages, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
    sim_env$Agg_ObsSrvLenComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
    sim_env$Agg_ISS_FishAgeComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
    sim_env$Agg_ISS_FishLenComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_fish_fleets, n_sims))
    sim_env$Agg_ISS_SrvAgeComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
    sim_env$Agg_ISS_SrvLenComps <- array(NA, dim = c(3, sim_env$n_yrs, sim_env$n_sexes, sim_env$n_srv_fleets, n_sims))
  }

  return(sim_env)
}

#' Aggregate Simulation Data to Form FAA Inputs
#'
#' This function aggregates simulated fishery and survey data to produce
#' catch, age–composition, and length–composition inputs for a
#' fully–aggregated assessment (FAA). It constructs aggregated catches,
#' catch–weighted compositions, and corresponding effective sample sizes
#' across regions and fleets, depending on the number of fishery and
#' survey fleets used in the operating model.
#'
#' @param sim_data
#'   A list containing simulated observed data, including components
#'   such as `ObsFishAgeComps`, `ObsFishLenComps`, and associated
#'   sampling variances.
#'
#' @param sim_env
#'   A list containing the simulation environment, including
#'   region-specific population states, true catches, true age
#'   compositions, survey observations, and arrays for storing aggregated
#'   results.
#'
#' @param y
#'   The current simulation year index (integer), relative to the start
#'   of the simulation.
#'
#' @param sim
#'   The simulation replicate index (integer).
#'
#' @param faa_n_fish_fleets
#'   The number of fishery fleets to be represented in the FAA model
#'   (typically \code{2}, \code{4}, \code{5}, or \code{6}).
#'
#' @param faa_n_srv_fleets
#'   The number of survey fleets to be represented in the FAA model.
#'
#' @param fish_wgt
#'   The weighting scheme used when aggregating fishery age and length
#'   compositions. Accepts \code{"biomass"} or \code{"numbers"}.
#'   Under \code{"biomass"}, regional catches are used. Under
#'   \code{"numbers"}, regional age composition totals are used.
#'
#' @param srv_wgt
#'   The weighting scheme used when aggregating survey data. Accepts
#'   \code{"numbers"} or \code{"biomass"}; usage is analogous to
#'   \code{fish_wgt}.
#'
#' @param srv_idx_se
#'   A vector or array of standard errors associated with survey indices.
#'   These values are used when adding observation error to aggregated
#'   survey indices.
#'
#' @returns
#'   The function returns \emph{invisibly}, but updates elements inside
#'   \code{sim_env} in place. The following components are populated or
#'   extended:
#'   \itemize{
#'     \item \code{Agg_ObsCatch}: Aggregated observed catches.
#'     \item \code{Agg_ObsFishAgeComps}: Aggregated fishery age compositions.
#'     \item \code{Agg_ObsFishLenComps}: Aggregated fishery length compositions.
#'     \item \code{Agg_ISS_FishAgeComps}: Effective sample sizes for age compositions.
#'     \item \code{Agg_ISS_FishLenComps}: Effective sample sizes for length compositions.
#'     \item Corresponding aggregated survey quantities.
#'   }
#'
agg_data_to_faa <- function(sim_data,
                            sim_env,
                            y,
                            sim,
                            faa_n_fish_fleets,
                            faa_n_srv_fleets,
                            srv_wgt = 'numbers',
                            fish_wgt = 'biomass',
                            srv_idx_se
                            ) {

  # dimensions
  n_regions <- 1
  n_yrs <- length(1:y)
  n_sexes <- sim_env$n_sexes
  n_ages <- sim_env$n_ages
  n_lens <- sim_env$n_lens
  n_fish_fleets <- faa_n_fish_fleets
  n_srv_fleets <- faa_n_srv_fleets

  # Catches (using true catches, and then applying error to aggregated true catches)
  if(faa_n_fish_fleets == 2) { # Two Fishery Fleets (i.e., fishery is not FAA)
    if(y == sim_env$feedback_start_yr) {
      aggregated_true_catch <- apply(sim_env$TrueCatch[,1:y,,sim], c(2,3), sum) # get true aggregated catch
      sim_env$Agg_ObsCatch[1,1:y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0, mean(exp(sim_data$ln_sigmaC))))
    } else{
      aggregated_true_catch <- apply(sim_env$TrueCatch[,y,,sim, drop = FALSE], c(2,3), sum) # get true aggregated catch
      sim_env$Agg_ObsCatch[1,y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0, mean(exp(sim_data$ln_sigmaC))))
    }
  }

  if(faa_n_fish_fleets == 4) { # Four Fishery Fleets (i.e., fixed-gear is FAA, trawl is not)
    if(y == sim_env$feedback_start_yr) {
      aggregated_true_catch <- array(0, dim = c(1, length(1:y), faa_n_fish_fleets))
      aggregated_true_catch[1,,1] <- sim_env$TrueCatch[1,1:y,1,sim] # BS Fixed Gear
      aggregated_true_catch[1,,2] <- sim_env$TrueCatch[2,1:y,1,sim] # AI Fixed Gear
      aggregated_true_catch[1,,3] <- colSums(sim_env$TrueCatch[3:5,1:y,1,sim]) # GOA Fixed Gear
      aggregated_true_catch[1,,4] <- colSums(sim_env$TrueCatch[,1:y,2,sim]) # Trawl Gear
      sim_env$Agg_ObsCatch[1,1:y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0,   mean(exp(sim_data$ln_sigmaC)))) # add error
    } else{
      aggregated_true_catch <- array(0, dim = c(n_fish_fleets))
      aggregated_true_catch[1] <- sim_env$TrueCatch[1,y,1,sim] # BS Fixed Gear
      aggregated_true_catch[2] <- sim_env$TrueCatch[2,y,1,sim] # AI Fixed Gear
      aggregated_true_catch[3] <- sum(sim_env$TrueCatch[3:5,y,1,sim]) # GOA Fixed Gear
      aggregated_true_catch[4] <- sum(sim_env$TrueCatch[,y,2,sim]) # Trawl Gear
      sim_env$Agg_ObsCatch[1,y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0,
                                                                                  mean(exp(sim_data$ln_sigmaC))))
    }
  }

  if(faa_n_fish_fleets == 5) { # Five Fishery Fleets (i.e., fixed gear FAA, trawl gear BS is FAA, and AI and GOA are combined)
    if(y == sim_env$feedback_start_yr) {
      aggregated_true_catch <- array(0, dim = c(1, length(1:y), faa_n_fish_fleets))
      aggregated_true_catch[1,,1] <- sim_env$TrueCatch[1,1:y,1,sim] # BS Fixed Gear
      aggregated_true_catch[1,,2] <- sim_env$TrueCatch[2,1:y,1,sim] # AI Fixed Gear
      aggregated_true_catch[1,,4] <- sim_env$TrueCatch[1,1:y,2,sim] # BS Trawl Gear

      # Do GOA Fleets
      aggregated_true_catch[1,,3] <- colSums(sim_env$TrueCatch[3:5,1:y,1,sim]) # GOA Fixed Gear
      aggregated_true_catch[1,,5] <- colSums(sim_env$TrueCatch[2:5,1:y,2,sim]) # AI + GOA Trawl Gear
      sim_env$Agg_ObsCatch[1,1:y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0,   mean(exp(sim_data$ln_sigmaC)))) # add error

    } else{
      aggregated_true_catch <- array(0, dim = c(n_fish_fleets))
      aggregated_true_catch[1] <- sim_env$TrueCatch[1,y,1,sim] # BS Fixed Gear
      aggregated_true_catch[2] <- sim_env$TrueCatch[2,y,1,sim] # AI Fixed Gear
      aggregated_true_catch[4] <- sim_env$TrueCatch[1,y,2,sim] # BS Trawl Gear

      # Do GOA Fleets
      aggregated_true_catch[3] <- sum(sim_env$TrueCatch[3:5,y,1,sim]) # GOA Fixed Gear
      aggregated_true_catch[5] <- sum(sim_env$TrueCatch[2:5,y,2,sim]) # AI + GOA Trawl Gear
      sim_env$Agg_ObsCatch[1,y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0,
                                                                                  mean(exp(sim_data$ln_sigmaC))))
    }
  }

  if(faa_n_fish_fleets == 6) { # Six Fishery Fleets (i.e., fixed gear and trawl gear are FAA)
    if(y == sim_env$feedback_start_yr) {
      aggregated_true_catch <- array(0, dim = c(1, length(1:y), faa_n_fish_fleets))
      aggregated_true_catch[1,,1] <- sim_env$TrueCatch[1,1:y,1,sim] # BS Fixed Gear
      aggregated_true_catch[1,,2] <- sim_env$TrueCatch[2,1:y,1,sim] # AI Fixed Gear
      aggregated_true_catch[1,,4] <- sim_env$TrueCatch[1,1:y,2,sim] # BS Trawl Gear
      aggregated_true_catch[1,,5] <- sim_env$TrueCatch[2,1:y,2,sim] # AI Trawl Gear

      # Do GOA Fleets
      aggregated_true_catch[1,,3] <- colSums(sim_env$TrueCatch[3:5,1:y,1,sim]) # GOA Fixed Gear
      aggregated_true_catch[1,,6] <- colSums(sim_env$TrueCatch[3:5,1:y,2,sim]) # GOA Trawl Gear
      sim_env$Agg_ObsCatch[1,1:y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0,   mean(exp(sim_data$ln_sigmaC)))) # add error

    } else{
      aggregated_true_catch <- array(0, dim = c(n_fish_fleets))
      aggregated_true_catch[1] <- sim_env$TrueCatch[1,y,1,sim] # BS Fixed Gear
      aggregated_true_catch[2] <- sim_env$TrueCatch[2,y,1,sim] # AI Fixed Gear
      aggregated_true_catch[4] <- sim_env$TrueCatch[1,y,2,sim] # BS Trawl Gear
      aggregated_true_catch[5] <- sim_env$TrueCatch[2,y,2,sim] # AI Trawl Gear

      # Do GOA Fleets
      aggregated_true_catch[3] <- sum(sim_env$TrueCatch[3:5,y,1,sim]) # GOA Fixed Gear
      aggregated_true_catch[6] <- sum(sim_env$TrueCatch[3:5,y,2,sim]) # GOA Trawl Gear
      sim_env$Agg_ObsCatch[1,y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0,
                                                                                  mean(exp(sim_data$ln_sigmaC))))
    }
  }

  # Fishery Age Compositions
  # Figure out weighting schemes
  if(faa_n_fish_fleets == 2) {
    if(fish_wgt == 'biomass') {
      total_catches <- apply(sim_env$TrueCatch[,1:y,,sim], c(2, 3), sum)  # Sum across regions for each year-fleet
      total_catches <- array(total_catches, dim = c(n_regions, dim(total_catches))) # coerce into correct format
      catch_prop <- sweep(sim_env$TrueCatch[,1:y,,sim], c(2, 3), total_catches, "/") # get catch proportion by region
      catch_fixed_gear_prop <- catch_prop[,,1]
      catch_trawl_gear_prop <- catch_prop[,,2]
    }

    if(fish_wgt == 'numbers') {
      catch_prop <- apply(sim_env$CAA[,1:y,,,,sim], c(1,2,5), sum)
      total_by_year <- apply(catch_prop, c(2,3), sum)
      catch_prop <- sweep(catch_prop, c(2, 3), total_by_year, "/") # get catch proportion by region
      catch_fixed_gear_prop <- catch_prop[,,1]
      catch_trawl_gear_prop <- catch_prop[,,2]
    }
  }

  if(faa_n_fish_fleets == 4) {
    if(fish_wgt == 'biomass') {

      # Get fixed gear weights
      fixed_gear_wts <- sim_env$TrueCatch[3:5,1:y,1,sim] # only weighting GOA for fixed-gear weights
      total_fixed_gear_catches <- colSums(fixed_gear_wts)
      catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

      # get trawl gear weights
      trawl_gear_wts <- sim_env$TrueCatch[,1:y,2,sim]
      total_trawl_gear_catches <- colSums(trawl_gear_wts)
      catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
    }

    if(fish_wgt == 'numbers') {
      # Get fixed gear weights
      fixed_gear_wts <- sim_env$CAA[3:5,1:y,,,1,sim] # only weighting GOA for fixed-gear weights
      fixed_gear_wts <- apply(fixed_gear_wts, c(1,2), sum) # sum numbers up
      total_fixed_gear_catches <- colSums(fixed_gear_wts)
      catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

      # get trawl gear weights
      trawl_gear_wts <- sim_env$CAA[,1:y,,,2,sim]
      trawl_gear_wts <- apply(trawl_gear_wts, c(1,2), sum) # sum numbers up
      total_trawl_gear_catches <- colSums(trawl_gear_wts)
      catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
    }
  }

  if(faa_n_fish_fleets == 5) {
    if(fish_wgt == 'biomass') {
      # Get fixed gear weights
      fixed_gear_wts <- sim_env$TrueCatch[3:5,1:y,1,sim] # only weighting GOA for fixed-gear weights
      total_fixed_gear_catches <- colSums(fixed_gear_wts)
      catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

      # get trawl gear weights
      trawl_gear_wts <- sim_env$TrueCatch[2:5,1:y,2,sim]
      total_trawl_gear_catches <- colSums(trawl_gear_wts)
      catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
    }

    if(fish_wgt == 'numbers') {
      # Get fixed gear weights
      fixed_gear_wts <- sim_env$CAA[3:5,1:y,,,1,sim] # only weighting GOA for fixed-gear weights
      fixed_gear_wts <- apply(fixed_gear_wts, c(1,2), sum) # sum numbers up
      total_fixed_gear_catches <- colSums(fixed_gear_wts)
      catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

      # get trawl gear weights
      trawl_gear_wts <- sim_env$CAA[2:5,1:y,,,2,sim]
      trawl_gear_wts <- apply(trawl_gear_wts, c(1,2), sum) # sum numbers up
      total_trawl_gear_catches <- colSums(trawl_gear_wts)
      catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
    }
  }

  if(faa_n_fish_fleets == 6) {
    if(fish_wgt == 'biomass') {
      # Get fixed gear weights
      fixed_gear_wts <- sim_env$TrueCatch[3:5,1:y,1,sim] # only weighting GOA for fixed-gear weights
      total_fixed_gear_catches <- colSums(fixed_gear_wts)
      catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

      # get trawl gear weights
      trawl_gear_wts <- sim_env$TrueCatch[3:5,1:y,2,sim]
      total_trawl_gear_catches <- colSums(trawl_gear_wts)
      catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
    }

    if(fish_wgt == 'numbers') {
      # Get fixed gear weights
      fixed_gear_wts <- sim_env$CAA[3:5,1:y,,,1,sim] # only weighting GOA for fixed-gear weights
      fixed_gear_wts <- apply(fixed_gear_wts, c(1,2), sum) # sum numbers up
      total_fixed_gear_catches <- colSums(fixed_gear_wts)
      catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

      # get trawl gear weights
      trawl_gear_wts <- sim_env$CAA[3:5,1:y,,,2,sim]
      trawl_gear_wts <- apply(trawl_gear_wts, c(1,2), sum) # sum numbers up
      total_trawl_gear_catches <- colSums(trawl_gear_wts)
      catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
    }
  }

  # loop through to get catch weighted compositions
  if(faa_n_fish_fleets == 2) {
    for(yr in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        # figure out which regions have catches
        if(f == 1) catch_prop <- catch_fixed_gear_prop
        if(f == 2) catch_prop <- catch_trawl_gear_prop
        regions <- (1:sim_env$n_regions)[-which(catch_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 1:sim_env$n_regions

        # age comps
        sim_env$Agg_ObsFishAgeComps[,yr,,,f,sim] <- apply(sim_data$ObsFishAgeComps[,yr,,,f] * catch_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsFishAgeComps[,yr,,,f,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,f,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,f,sim])

        # Calculate n_eff for age comps
        regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,f], 1, sum)
        regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,f], 1, regional_n_age, "/")
        pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions,yr], "*")) / sum(catch_prop[,yr])
        numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
        denominator_age <- sum(catch_prop[regions,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
        sim_env$Agg_ISS_FishAgeComps[,yr,1,f,sim] <- numerator_age / denominator_age

        # length comps
        sim_env$Agg_ObsFishLenComps[,yr,,,f,sim] <- apply(sim_data$ObsFishLenComps[,yr,,,f] * catch_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsFishLenComps[,yr,,,f,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,f,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,f,sim])

        # Calculate n_eff for length comps
        regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,f], 1, sum)
        regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,f], 1, regional_n_len, "/")
        pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions,yr], "*")) / sum(catch_prop[,yr])
        numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
        denominator_len <- sum(catch_prop[regions,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
        sim_env$Agg_ISS_FishLenComps[,yr,1,f,sim] <- numerator_len / denominator_len

      } # end f loop
    } # end yr loop
  }

  if(faa_n_fish_fleets == 4) {
    for(yr in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        # figure out which regions have catches
        if(f == 1) catch_prop <- catch_fixed_gear_prop
        if(f == 2) catch_prop <- catch_trawl_gear_prop

        # age comps
        if(f == 1) { # fixed gear

          regions <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
          sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] <- sim_data$ObsFishAgeComps[1,yr,,,1] # BS
          sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] <- sim_data$ObsFishAgeComps[2,yr,,,1] # AI
          sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim])
          sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] <- apply(sim_data$ObsFishAgeComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_FishAgeComps[,yr,1,1,sim] <- sum(sim_data$ObsFishAgeComps[1,yr,,,1]) # sum ISS for BS fixed gear
          sim_env$Agg_ISS_FishAgeComps[,yr,1,2,sim] <- sum(sim_data$ObsFishAgeComps[2,yr,,,1]) # sum ISS for AI fixed gear

          # get n_eff for GOA fixed gear
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, sum)
            regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_FishAgeComps[,yr,1,3,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_FishAgeComps[,yr,1,3,sim] <- sum(sim_data$ObsFishAgeComps[regions,yr,,,1])
          }
        }

        # trawl age comps (not used)
        if(f == 2) {
          regions <- (1:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 1:5
          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] <- apply(sim_data$ObsFishAgeComps[,yr,,,2] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim])
          # get age compositions
          regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, sum)
          regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, regional_n_age, "/")
          pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions,yr], "*")) / sum(catch_prop[,yr])
          numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
          denominator_age <- sum(catch_prop[regions,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
          sim_env$Agg_ISS_FishAgeComps[,yr,1,4,sim] <- numerator_age / denominator_age
        }

        # len comps
        if(f == 1) { # fixed gear

          regions <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
          sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] <- sim_data$ObsFishLenComps[1,yr,,,1] # BS
          sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] <- sim_data$ObsFishLenComps[2,yr,,,1] # AI
          sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,1,sim])
          sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] <- apply(sim_data$ObsFishLenComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,3,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_FishLenComps[,yr,1,1,sim] <- sum(sim_data$ObsFishLenComps[1,yr,,,1]) # sum ISS for BS fixed gear
          sim_env$Agg_ISS_FishLenComps[,yr,1,2,sim] <- sum(sim_data$ObsFishLenComps[2,yr,,,1]) # sum ISS for AI fixed gear

          # get n_eff for GOA fixed gear
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,1], 1, sum)
            regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,1], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_FishLenComps[,yr,1,3,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_FishLenComps[,yr,1,3,sim] <- sum(sim_data$ObsFishLenComps[regions,yr,,,1])
          }
        }

        # trawl len comps
        if(f == 2) {
          regions <- (1:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 1:5
          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] <- apply(sim_data$ObsFishLenComps[,yr,,,2] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,4,sim])
          # get len compositions
          regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,2], 1, sum)
          regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,2], 1, regional_n_len, "/")
          pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions,yr], "*")) / sum(catch_prop[,yr])
          numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
          denominator_len <- sum(catch_prop[regions,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
          sim_env$Agg_ISS_FishLenComps[,yr,1,4,sim] <- numerator_len / denominator_len
        }

      } # end f loop
    } # end yr loop
  }

  if(faa_n_fish_fleets == 5) {
    for(yr in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        # figure out which regions have catches
        if(f == 1) catch_prop <- catch_fixed_gear_prop
        if(f == 2) catch_prop <- catch_trawl_gear_prop

        # age comps
        if(f == 1) { # fixed gear

          regions <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
          sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] <- sim_data$ObsFishAgeComps[1,yr,,,1] # BS
          sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] <- sim_data$ObsFishAgeComps[2,yr,,,1] # AI
          sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim])
          sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] <- apply(sim_data$ObsFishAgeComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_FishAgeComps[,yr,1,1,sim] <- sum(sim_data$ObsFishAgeComps[1,yr,,,1]) # sum ISS for BS fixed gear
          sim_env$Agg_ISS_FishAgeComps[,yr,1,2,sim] <- sum(sim_data$ObsFishAgeComps[2,yr,,,1]) # sum ISS for AI fixed gear

          # get n_eff for GOA fixed gear
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, sum)
            regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_FishAgeComps[,yr,1,3,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_FishAgeComps[,yr,1,3,sim] <- sum(sim_data$ObsFishAgeComps[regions,yr,,,1])
          }
        }

        # trawl age comps (not used)
        if(f == 2) {
          regions <- (2:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 2:5
          regions_offset <- regions - 1 # get offset for indexing catch_prop

          # input BS, into appropriate fleets for trawl-gear fishery (no weighting)
          sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] <- sim_data$ObsFishAgeComps[1,yr,,,2] # BS
          sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim])

          # input AI + GOA fleet weighting by proportions
          sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim] <- apply(sim_data$ObsFishAgeComps[2:5,yr,,,2] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_FishAgeComps[,yr,1,4,sim] <- sum(sim_data$ObsFishAgeComps[1,yr,,,2]) # sum ISS for BS trawl gear

          # get n_eff for AI + GOA trawl gear
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, sum)
            regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_FishAgeComps[,yr,1,5,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_FishAgeComps[,yr,1,5,sim] <- sum(sim_data$ObsFishAgeComps[regions,yr,,,2])
          }
        }

        # len comps
        if(f == 1) { # fixed gear

          regionf <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
          sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] <- sim_data$ObsFishLenComps[1,yr,,,1] # BS
          sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] <- sim_data$ObsFishLenComps[2,yr,,,1] # AI
          sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,1,sim])
          sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] <- apply(sim_data$ObsFishLenComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,3,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_FishLenComps[,yr,1,1,sim] <- sum(sim_data$ObsFishLenComps[1,yr,,,1]) # sum ISS for BS fixed gear
          sim_env$Agg_ISS_FishLenComps[,yr,1,2,sim] <- sum(sim_data$ObsFishLenComps[2,yr,,,1]) # sum ISS for AI fixed gear

          # get n_eff for GOA fixed gear
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,1], 1, sum)
            regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,1], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_FishLenComps[,yr,1,3,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_FishLenComps[,yr,1,3,sim] <- sum(sim_data$ObsFishLenComps[regions,yr,,,1])
          }
        }

        # trawl len comps
        if(f == 2) {
          regions <- (2:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 2:5
          regions_offset <- regions - 1 # get offset for indexing catch_prop

          # input BS, into appropriate fleets for trawl-gear fishery (no weighting)
          sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] <- sim_data$ObsFishLenComps[1,yr,,,2] # BS
          sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,4,sim])

          # input AI + GOA fleet weighting by proportions
          sim_env$Agg_ObsFishLenComps[,yr,,,5,sim] <- apply(sim_data$ObsFishLenComps[2:5,yr,,,2] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishLenComps[,yr,,,5,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,5,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,5,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_FishLenComps[,yr,1,4,sim] <- sum(sim_data$ObsFishLenComps[1,yr,,,2]) # sum ISS for BS trawl gear

          # get n_eff for GOA trawl gear
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,2], 1, sum)
            regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,2], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_FishLenComps[,yr,1,5,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_FishLenComps[,yr,1,5,sim] <- sum(sim_data$ObsFishLenComps[regions,yr,,,2])
          }
        }

      } # end f loop
    } # end yr loop
  }

  if(faa_n_fish_fleets == 6) {
    for(yr in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        # figure out which regions have catches
        if(f == 1) catch_prop <- catch_fixed_gear_prop
        if(f == 2) catch_prop <- catch_trawl_gear_prop

        # age comps
        if(f == 1) { # fixed gear

          regions <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
          sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] <- sim_data$ObsFishAgeComps[1,yr,,,1] # BS
          sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] <- sim_data$ObsFishAgeComps[2,yr,,,1] # AI
          sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,1,sim])
          sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] <- apply(sim_data$ObsFishAgeComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,3,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_FishAgeComps[,yr,1,1,sim] <- sum(sim_data$ObsFishAgeComps[1,yr,,,1]) # sum ISS for BS fixed gear
          sim_env$Agg_ISS_FishAgeComps[,yr,1,2,sim] <- sum(sim_data$ObsFishAgeComps[2,yr,,,1]) # sum ISS for AI fixed gear

          # get n_eff for GOA fixed gear
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, sum)
            regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_FishAgeComps[,yr,1,3,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_FishAgeComps[,yr,1,3,sim] <- sum(sim_data$ObsFishAgeComps[regions,yr,,,1])
          }
        }

        # trawl age comps (not used)
        if(f == 2) {
          regions <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for trawl-gear fishery (no weighting)
          sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] <- sim_data$ObsFishAgeComps[1,yr,,,2] # BS
          sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim] <- sim_data$ObsFishAgeComps[2,yr,,,2] # AI
          sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,4,sim])
          sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,5,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishAgeComps[,yr,,,6,sim] <- apply(sim_data$ObsFishAgeComps[3:5,yr,,,2] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishAgeComps[,yr,,,6,sim] <- sim_env$Agg_ObsFishAgeComps[,yr,,,6,sim] / sum(sim_env$Agg_ObsFishAgeComps[,yr,,,6,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_FishAgeComps[,yr,1,4,sim] <- sum(sim_data$ObsFishAgeComps[1,yr,,,2]) # sum ISS for BS trawl gear
          sim_env$Agg_ISS_FishAgeComps[,yr,1,5,sim] <- sum(sim_data$ObsFishAgeComps[2,yr,,,2]) # sum ISS for AI trawl gear

          # get n_eff for GOA trawl gear
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, sum)
            regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_FishAgeComps[,yr,1,6,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_FishAgeComps[,yr,1,6,sim] <- sum(sim_data$ObsFishAgeComps[regions,yr,,,2])
          }
        }

        # len comps
        if(f == 1) { # fixed gear

          regionf <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
          sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] <- sim_data$ObsFishLenComps[1,yr,,,1] # BS
          sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] <- sim_data$ObsFishLenComps[2,yr,,,1] # AI
          sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,1,sim])
          sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] <- apply(sim_data$ObsFishLenComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,3,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_FishLenComps[,yr,1,1,sim] <- sum(sim_data$ObsFishLenComps[1,yr,,,1]) # sum ISS for BS fixed gear
          sim_env$Agg_ISS_FishLenComps[,yr,1,2,sim] <- sum(sim_data$ObsFishLenComps[2,yr,,,1]) # sum ISS for AI fixed gear

          # get n_eff for GOA fixed gear
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,1], 1, sum)
            regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,1], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_FishLenComps[,yr,1,3,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_FishLenComps[,yr,1,3,sim] <- sum(sim_data$ObsFishLenComps[regions,yr,,,1])
          }
        }

        # trawl len comps
        if(f == 2) {
          regions <- (3:5)[-which(catch_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing catch_prop

          # input BS, AI, into appropriate fleets for trawl-gear fishery (no weighting)
          sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] <- sim_data$ObsFishLenComps[1,yr,,,2] # BS
          sim_env$Agg_ObsFishLenComps[,yr,,,5,sim] <- sim_data$ObsFishLenComps[2,yr,,,2] # AI
          sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,4,sim])
          sim_env$Agg_ObsFishLenComps[,yr,,,5,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,5,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,5,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsFishLenComps[,yr,,,6,sim] <- apply(sim_data$ObsFishLenComps[3:5,yr,,,2] * catch_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsFishLenComps[,yr,,,6,sim] <- sim_env$Agg_ObsFishLenComps[,yr,,,6,sim] / sum(sim_env$Agg_ObsFishLenComps[,yr,,,6,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_FishLenComps[,yr,1,4,sim] <- sum(sim_data$ObsFishLenComps[1,yr,,,2]) # sum ISS for BS trawl gear
          sim_env$Agg_ISS_FishLenComps[,yr,1,5,sim] <- sum(sim_data$ObsFishLenComps[2,yr,,,2]) # sum ISS for AI trawl gear

          # get n_eff for GOA trawl gear
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,2], 1, sum)
            regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,2], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_FishLenComps[,yr,1,6,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_FishLenComps[,yr,1,6,sim] <- sum(sim_data$ObsFishLenComps[regions,yr,,,2])
          }
        }

      } # end f loop
    } # end yr loop
  }

  if(y == sim_env$feedback_start_yr) {
    true_srvidx <- sim_env$TrueSrvIdx[,1:y,,sim] # get true survey index
    sim_env$Agg_TrueSrvIdx[,1:y,1,sim] <- true_srvidx[1,1:y,1] # BS Domestic
    sim_env$Agg_TrueSrvIdx[,1:y,2,sim] <- true_srvidx[2,1:y,1] # AI Domestic
    sim_env$Agg_TrueSrvIdx[,1:y,3,sim] <- colSums(true_srvidx[3:5,1:y,1]) # GOA Domestic
    if(faa_n_srv_fleets == 4) sim_env$Agg_TrueSrvIdx[,1:y,4,sim] <- colSums(true_srvidx[,1:y,3]) # BS + AI + GOA JP
    if(faa_n_srv_fleets == 6) {
      sim_env$Agg_TrueSrvIdx[,1:y,4,sim] <- true_srvidx[1,1:y,3] # BS Domestic
      sim_env$Agg_TrueSrvIdx[,1:y,5,sim] <- true_srvidx[2,1:y,3] # AI Domestic
      sim_env$Agg_TrueSrvIdx[,1:y,6,sim] <- colSums(true_srvidx[3:5,1:y,3]) # GOA Domestic
    }
    # add devs
    sim_env$Agg_ObsSrvIdx[1,1:y,,sim] <- sim_env$Agg_TrueSrvIdx[1,1:y,,sim] * exp(rnorm(length(sim_env$Agg_TrueSrvIdx[1,1:y,,sim]), 0, srv_idx_se)) # simulate devs
  } else {

    true_srvidx <- sim_env$TrueSrvIdx[,y,,sim] # get true survey index
    if(faa_n_srv_fleets == 4) {
      agg_true_srvidx <- array(0, dim = c(4))
      agg_true_srvidx[1] <- true_srvidx[1,1] # BS Domestic
      agg_true_srvidx[2] <- true_srvidx[2,1] # AI Domestic
      agg_true_srvidx[3] <- sum(true_srvidx[3:5,1]) # GOA Domestic
      agg_true_srvidx[4] <- sum(true_srvidx[,3]) # BS, AI, GOA JP
    }

    if(faa_n_srv_fleets == 6) {
      agg_true_srvidx <- array(0, dim = c(6))
      agg_true_srvidx[1] <- true_srvidx[1,1] # BS Domestic
      agg_true_srvidx[2] <- true_srvidx[2,1] # AI Domestic
      agg_true_srvidx[3] <- sum(true_srvidx[3:5,1]) # GOA Domestic
      agg_true_srvidx[4] <- true_srvidx[1,3] # BS JP
      agg_true_srvidx[5] <- true_srvidx[2,3] # AI JP
      agg_true_srvidx[6] <- sum(true_srvidx[3:5,3]) # GOA JP
    }

    # save true aggregated index
    sim_env$Agg_TrueSrvIdx[1,y,,sim] <- agg_true_srvidx # save "true index"
    sim_env$Agg_ObsSrvIdx[1,y,,sim] <-
      agg_true_srvidx * exp(rnorm(length(agg_true_srvidx), 0, srv_idx_se))

  }

  if(faa_n_srv_fleets == 4) {
    if(srv_wgt == 'biomass') {

      # Get domestic weights
      domestic_wts <- apply(sim_env$SrvIAA[3:5,1:y,,,1,sim] * sim_env$WAA_srv[3:5,1:y,,,1,sim], c(1,2), sum)  # only weighting GOA for domestic weights
      total_domestic_idx <- colSums(domestic_wts)
      idx_domestic_prop <- sweep(domestic_wts, 2, total_domestic_idx, '/')

      # get jp weights
      jp_wts <- apply(sim_env$SrvIAA[3:5,1:y,,,3,sim] * sim_env$WAA_srv[3:5,1:y,,,3,sim], c(1,2), sum)
      total_jp_idx <- colSums(jp_wts)
      idx_jp_prop <- sweep(jp_wts, 2, total_jp_idx, '/')
    }

    if(srv_wgt == 'numbers') {
      # Get domestic weights
      domestic_wts <- sim_env$SrvIAA[3:5,1:y,,,1,sim] # only weighting GOA for domestic weights
      domestic_wts <- apply(domestic_wts, c(1,2), sum) # sum numbers up
      total_domestic_idx <- colSums(domestic_wts)
      idx_domestic_prop <- sweep(domestic_wts, 2, total_domestic_idx, '/')

      # get jp weights
      jp_wts <- sim_env$SrvIAA[,1:y,,,3,sim]
      jp_wts <- apply(jp_wts, c(1,2), sum) # sum numbers up
      total_jp_idx <- colSums(jp_wts)
      idx_jp_prop <- sweep(jp_wts, 2, total_jp_idx, '/')
    }
  }

  if(faa_n_srv_fleets == 6) {
    if(srv_wgt == 'biomass') {
      # Get domestic weights
      domestic_wts <- apply(sim_env$SrvIAA[3:5,1:y,,,1,sim] * sim_env$WAA_srv[3:5,1:y,,,1,sim], c(1,2), sum)  # only weighting GOA for domestic weights
      total_domestic_idx <- colSums(domestic_wts)
      idx_domestic_prop <- sweep(domestic_wts, 2, total_domestic_idx, '/')

      # get jp weights
      jp_wts <- apply(sim_env$SrvIAA[3:5,1:y,,,3,sim] * sim_env$WAA_srv[3:5,1:y,,,3,sim], c(1,2), sum)
      total_jp_idx <- colSums(jp_wts)
      idx_jp_prop <- sweep(jp_wts, 2, total_jp_idx, '/')
    }

    if(srv_wgt == 'numbers') {
      # Get domestic weights
      domestic_wts <- sim_env$SrvIAA[3:5,1:y,,,1,sim] # only weighting GOA for domestic weights
      domestic_wts <- apply(domestic_wts, c(1,2), sum) # sum numbers up
      total_domestic_idx <- colSums(domestic_wts)
      idx_domestic_prop <- sweep(domestic_wts, 2, total_domestic_idx, '/')

      # get jp weights
      jp_wts <- sim_env$SrvIAA[3:5,1:y,,,3,sim]
      jp_wts <- apply(jp_wts, c(1,2), sum) # sum numbers up
      total_jp_idx <- colSums(jp_wts)
      idx_jp_prop <- sweep(jp_wts, 2, total_jp_idx, '/')
    }
  }

  if(faa_n_srv_fleets == 4) {
    for(yr in 1:n_yrs) {
      for(f in c(1,3)) { # looping through only domestic and jp fleet

        # figure out which regions have catches
        if(f == 1) idx_prop <- idx_domestic_prop
        if(f == 3) idx_prop <- idx_jp_prop

        # age comps
        if(f == 1) { # domestic

          regions <- (3:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing idx_prop

          # input BS, AI, into appropriate fleets for domestic srvery (no weighting)
          sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim] <- sim_data$ObsSrvAgeComps[1,yr,,,1] # BS
          sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim] <- sim_data$ObsSrvAgeComps[2,yr,,,1] # AI
          sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim])
          sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim] <- apply(sim_data$ObsSrvAgeComps[3:5,yr,,,1] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_SrvAgeComps[,yr,1,1,sim] <- sum(sim_data$ObsSrvAgeComps[1,yr,,,1]) # sum ISS for BS domestic
          sim_env$Agg_ISS_SrvAgeComps[,yr,1,2,sim] <- sum(sim_data$ObsSrvAgeComps[2,yr,,,1]) # sum ISS for AI domestic

          # get n_eff for GOA domestic
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsSrvAgeComps[regions,yr,,,1], 1, sum)
            regional_p_age <- sweep(sim_data$ObsSrvAgeComps[regions,yr,,,1], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_SrvAgeComps[,yr,1,3,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_SrvAgeComps[,yr,1,3,sim] <- sum(sim_data$ObsSrvAgeComps[regions,yr,,,1])
          }
        }

        # jp age comps
        if(f == 3) {
          regions <- (1:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 1:5
          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim] <- apply(sim_data$ObsSrvAgeComps[,yr,,,3] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim])
          # get age compositions
          regional_n_age <- apply(sim_data$ObsSrvAgeComps[regions,yr,,,3], 1, sum)
          regional_p_age <- sweep(sim_data$ObsSrvAgeComps[regions,yr,,,3], 1, regional_n_age, "/")
          pooled_comp_age <- colSums(sweep(regional_p_age, 1, idx_prop[regions,yr], "*")) / sum(idx_prop[,yr])
          numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
          denominator_age <- sum(idx_prop[regions,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
          sim_env$Agg_ISS_SrvAgeComps[,yr,1,4,sim] <- numerator_age / denominator_age
        }

        # len comps (not used)
        if(f == 1) { # domestic

          regions <- (3:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing idx_prop

          # input BS, AI, into appropriate fleets for domestic srvery (no weighting)
          sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim] <- sim_data$ObsSrvLenComps[1,yr,,,1] # BS
          sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim] <- sim_data$ObsSrvLenComps[2,yr,,,1] # AI
          sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim])
          sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim] <- apply(sim_data$ObsSrvLenComps[3:5,yr,,,1] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_SrvLenComps[,yr,1,1,sim] <- sum(sim_data$ObsSrvLenComps[1,yr,,,1]) # sum ISS for BS domestic
          sim_env$Agg_ISS_SrvLenComps[,yr,1,2,sim] <- sum(sim_data$ObsSrvLenComps[2,yr,,,1]) # sum ISS for AI domestic

          # get n_eff for GOA domestic
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsSrvLenComps[regions,yr,,,1], 1, sum)
            regional_p_len <- sweep(sim_data$ObsSrvLenComps[regions,yr,,,1], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_SrvLenComps[,yr,1,3,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_SrvLenComps[,yr,1,3,sim] <- sum(sim_data$ObsSrvLenComps[regions,yr,,,1])
          }
        }

        # jp len comps
        if(f == 3) {
          regions <- (1:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 1:5
          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim] <- apply(sim_data$ObsSrvLenComps[,yr,,,3] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim])
          # get len compositions
          regional_n_len <- apply(sim_data$ObsSrvLenComps[regions,yr,,,3], 1, sum)
          regional_p_len <- sweep(sim_data$ObsSrvLenComps[regions,yr,,,3], 1, regional_n_len, "/")
          pooled_comp_len <- colSums(sweep(regional_p_len, 1, idx_prop[regions,yr], "*")) / sum(idx_prop[,yr])
          numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
          denominator_len <- sum(idx_prop[regions,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
          sim_env$Agg_ISS_SrvLenComps[,yr,1,4,sim] <- numerator_len / denominator_len
        }

      } # end f loop
    } # end yr loop
  }

  if(faa_n_srv_fleets == 6) {
    for(yr in 1:n_yrs) {
      for(f in c(1,3)) { # loop through only domestic and jp

        # figure out which regions have catches
        if(f == 1) idx_prop <- idx_domestic_prop
        if(f == 3) idx_prop <- idx_jp_prop

        # age comps
        if(f == 1) { # domestic

          regions <- (3:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing idx_prop

          # input BS, AI, into appropriate fleets for domestic srvery (no weighting)
          sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim] <- sim_data$ObsSrvAgeComps[1,yr,,,1] # BS
          sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim] <- sim_data$ObsSrvAgeComps[2,yr,,,1] # AI
          sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,1,sim])
          sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim] <- apply(sim_data$ObsSrvAgeComps[3:5,yr,,,1] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_SrvAgeComps[,yr,1,1,sim] <- sum(sim_data$ObsSrvAgeComps[1,yr,,,1]) # sum ISS for BS domestic
          sim_env$Agg_ISS_SrvAgeComps[,yr,1,2,sim] <- sum(sim_data$ObsSrvAgeComps[2,yr,,,1]) # sum ISS for AI domestic

          # get n_eff for GOA domestic
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsSrvAgeComps[regions,yr,,,1], 1, sum)
            regional_p_age <- sweep(sim_data$ObsSrvAgeComps[regions,yr,,,1], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_SrvAgeComps[,yr,1,3,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_SrvAgeComps[,yr,1,3,sim] <- sum(sim_data$ObsSrvAgeComps[regions,yr,,,1])
          }
        }

        # jp age comps
        if(f == 3) {
          regions <- (3:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing idx_prop

          # input BS, AI, into appropriate fleets for jp-gear srvery (no weighting)
          sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim] <- sim_data$ObsSrvAgeComps[1,yr,,,3] # BS
          sim_env$Agg_ObsSrvAgeComps[,yr,,,5,sim] <- sim_data$ObsSrvAgeComps[2,yr,,,3] # AI
          sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,4,sim])
          sim_env$Agg_ObsSrvAgeComps[,yr,,,5,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,5,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,5,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvAgeComps[,yr,,,6,sim] <- apply(sim_data$ObsSrvAgeComps[3:5,yr,,,3] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvAgeComps[,yr,,,6,sim] <- sim_env$Agg_ObsSrvAgeComps[,yr,,,6,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,6,sim])

          # Calculate n_eff for age comps
          sim_env$Agg_ISS_SrvAgeComps[,yr,1,4,sim] <- sum(sim_data$ObsSrvAgeComps[1,yr,,,3]) # sum ISS for BS jp
          sim_env$Agg_ISS_SrvAgeComps[,yr,1,5,sim] <- sum(sim_data$ObsSrvAgeComps[2,yr,,,3]) # sum ISS for AI jp

          # get n_eff for GOA jp
          if(length(regions) > 1) {
            regional_n_age <- apply(sim_data$ObsSrvAgeComps[regions,yr,,,3], 1, sum)
            regional_p_age <- sweep(sim_data$ObsSrvAgeComps[regions,yr,,,3], 1, regional_n_age, "/")
            pooled_comp_age <- colSums(sweep(regional_p_age, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
            numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
            denominator_age <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
            sim_env$Agg_ISS_SrvAgeComps[,yr,1,6,sim] <- numerator_age / denominator_age
          } else {
            sim_env$Agg_ISS_SrvAgeComps[,yr,1,6,sim] <- sum(sim_data$ObsSrvAgeComps[regions,yr,,,3])
          }
        }

        # len comps (not used)
        if(f == 1) { # domestic

          regions <- (3:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing idx_prop

          # input BS, AI, into appropriate fleets for domestic srvery (no weighting)
          sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim] <- sim_data$ObsSrvLenComps[1,yr,,,1] # BS
          sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim] <- sim_data$ObsSrvLenComps[2,yr,,,1] # AI
          sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,1,sim])
          sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,2,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim] <- apply(sim_data$ObsSrvLenComps[3:5,yr,,,1] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,3,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_SrvLenComps[,yr,1,1,sim] <- sum(sim_data$ObsSrvLenComps[1,yr,,,1]) # sum ISS for BS domestic
          sim_env$Agg_ISS_SrvLenComps[,yr,1,2,sim] <- sum(sim_data$ObsSrvLenComps[2,yr,,,1]) # sum ISS for AI domestic

          # get n_eff for GOA domestic
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsSrvLenComps[regions,yr,,,1], 1, sum)
            regional_p_len <- sweep(sim_data$ObsSrvLenComps[regions,yr,,,1], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_SrvLenComps[,yr,1,3,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_SrvLenComps[,yr,1,3,sim] <- sum(sim_data$ObsSrvLenComps[regions,yr,,,1])
          }
        }

        # jp len comps
        if(f == 3) {
          regions <- (3:5)[-which(idx_prop[,yr] == 0)]
          if(length(regions) == 0) regions <- 3:5
          regions_offset <- regions - 2 # get offset for indexing idx_prop

          # input BS, AI, into appropriate fleets for jp-gear srvery (no weighting)
          sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim] <- sim_data$ObsSrvLenComps[1,yr,,,3] # BS
          sim_env$Agg_ObsSrvLenComps[,yr,,,5,sim] <- sim_data$ObsSrvLenComps[2,yr,,,3] # AI
          sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,4,sim])
          sim_env$Agg_ObsSrvLenComps[,yr,,,5,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,5,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,5,sim])

          # input GOA fleet weighting by proportions
          sim_env$Agg_ObsSrvLenComps[,yr,,,6,sim] <- apply(sim_data$ObsSrvLenComps[3:5,yr,,,3] * idx_prop[,yr], c(2,3), sum)
          sim_env$Agg_ObsSrvLenComps[,yr,,,6,sim] <- sim_env$Agg_ObsSrvLenComps[,yr,,,6,sim] / sum(sim_env$Agg_ObsSrvLenComps[,yr,,,6,sim])

          # Calculate n_eff for len comps
          sim_env$Agg_ISS_SrvLenComps[,yr,1,4,sim] <- sum(sim_data$ObsSrvLenComps[1,yr,,,3]) # sum ISS for BS jp
          sim_env$Agg_ISS_SrvLenComps[,yr,1,5,sim] <- sum(sim_data$ObsSrvLenComps[2,yr,,,3]) # sum ISS for AI jp

          # get n_eff for GOA jp
          if(length(regions) > 1) {
            regional_n_len <- apply(sim_data$ObsSrvLenComps[regions,yr,,,3], 1, sum)
            regional_p_len <- sweep(sim_data$ObsSrvLenComps[regions,yr,,,3], 1, regional_n_len, "/")
            pooled_comp_len <- colSums(sweep(regional_p_len, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
            numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
            denominator_len <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
            sim_env$Agg_ISS_SrvLenComps[,yr,1,6,sim] <- numerator_len / denominator_len
          } else {
            sim_env$Agg_ISS_SrvLenComps[,yr,1,6,sim] <- sum(sim_data$ObsSrvLenComps[regions,yr,,,3])
          }
        }

      } # end f loop
    } # end yr loop
  }


}

#' Construct Data-Use Indicators for an Aggregated Assessment (FAA)
#'
#' This function generates logical indicator arrays that determine which
#' fishery catches, fishery compositions, and survey data are used in a
#' fully aggregated assessment (FAA) model. Indicators depend on the
#' number of fishery and survey fleets, the current year within the
#' simulation, the age-composition lag, and the assumed survey design
#' type.
#'
#' @param sim_env
#'   A list containing the simulation environment. It must include
#'   information such as \code{feedback_start_yr} and \code{n_yrs}, which
#'   define the simulation timeline.
#'
#' @param y
#'   The current year index within the simulation (integer). Indicators
#'   are constructed only up to this year.
#'
#' @param faa_n_fish_fleets
#'   The number of fishery fleets represented in the aggregated
#'   assessment (for example \code{2}, \code{4}, \code{5}, or \code{6}).
#'
#' @param faa_n_srv_fleets
#'   The number of survey fleets represented in the aggregated
#'   assessment.
#'
#' @param age_lag
#'   The number of years by which age–composition data are lagged
#'   relative to the year they inform (default \code{1}).
#'
#' @param lls_design_type
#'   A character string specifying the survey design used to determine
#'   which survey data are included. Accepts \code{"all"},
#'   \code{"historical"}, or \code{"current"}.
#'
#' @returns
#'   A list containing six indicator arrays, each having dimensions
#'   \code{c(1, y, n_fleet)}:
#'   \describe{
#'     \item{\code{usefishage}}{Indicators for fishery age–composition data.}
#'     \item{\code{usefishlen}}{Indicators for fishery length–composition data.}
#'     \item{\code{usesrvage}}{Indicators for survey age–composition data.}
#'     \item{\code{usesrvlen}}{Indicators for survey length–composition data.}
#'     \item{\code{usecatch}}{Indicators for fishery catch observations.}
#'     \item{\code{usesrvidx}}{Indicators for survey biomass or abundance indices.}
#'   }
#'
faa_use_indicators <- function(sim_env,
                               y,
                               faa_n_fish_fleets,
                               faa_n_srv_fleets,
                               age_lag = 1,
                               lls_design_type) {

  # get years in simulation
  x <- (sim_env$feedback_start_yr + 1):sim_env$n_yrs
  odd_yrs <- x[x %% 2 == 1] # bsai years
  even_yrs <- x[x %% 2 == 0] # goa years

  # Data indicators
  usecatch <- array(0, dim = c(1, y, faa_n_fish_fleets))
  usesrvidx <- array(0, dim = c(1, y, faa_n_srv_fleets))
  usefishage <- array(0, c(1, y, faa_n_fish_fleets))
  usefishlen <- array(0, c(1, y, faa_n_fish_fleets))
  usesrvage <- array(0, c(1, y, faa_n_srv_fleets))
  usesrvlen <- array(0, c(1, y, faa_n_srv_fleets)) # not used at all

  # Catches
  if(faa_n_fish_fleets == 2) { # availiable every year
    usecatch[] <- 1
  }

  if(faa_n_fish_fleets == 4) {
    usecatch[,,1] <- 1 # BS fixed gear catch availiable every year
    usecatch[,-c(1:3),2] <- 1 # AI fixed gear catch not availiable in first three years
    usecatch[,,3] <- 1 # GOA fixed gear catch availaible in every year
    usecatch[,,4] <- 1 # aggregated trawl catch availiable every year
  }

  if(faa_n_fish_fleets == 5) {
    usecatch[,,1] <- 1 # BS fixed gear catch availiable every year
    usecatch[,-c(1:3),2] <- 1 # AI fixed gear catch not availiable in first three years
    usecatch[,,3] <- 1 # GOA fixed gear catch availaible in every year
    usecatch[,,4] <- 1 # BS trawl gear catch availiable every year
    usecatch[,,5] <- 1 # aggregated AI + GOA trawl gear catch availaible in every year
  }

  if(faa_n_fish_fleets == 6) {
    usecatch[,,1] <- 1 # BS fixed gear catch availiable every year
    usecatch[,-c(1:3),2] <- 1 # AI fixed gear catch not availiable in first three years
    usecatch[,,3] <- 1 # GOA fixed gear catch availaible in every year
    usecatch[,,4] <- 1 # BS trawl gear catch availiable every year
    usecatch[,-c(1:3),5] <- 1 # AI trawl gear catch not availiable in first three years
    usecatch[,,6] <- 1 # aggregated GOA trawl gear catch availaible in every year
  }

  # Fishery Compositions
  if(faa_n_fish_fleets == 2) { # availiable every year
    usefishage[,40:(y-age_lag),1] <- 1 # use data starting year 40, with age lag for all fixed gear
    usefishlen[,31:39,1] <- 1 # use data only from 31 - 39 for all fixed gear
    usefishlen[,c(31,32,35:37,39:y),2] <- 1 # trawl gear length comps (aggregated)
  }

  if(faa_n_fish_fleets == 4) {
    usefishage[,40:(y-age_lag),1:3] <- 1 # use data starting year 40, with age lag for all fixed gear
    usefishage[,55,1] <- 0 # no samples in BS fleet in year 55
    usefishlen[,31:39,1:3] <- 1 # use data only from 31 - 39 for all fixed gear
    usefishlen[,c(31,32,35:37,39:y),4] <- 1 # trawl gera length comps (aggregated)
  }

  if(faa_n_fish_fleets == 5) {
    usefishage[,40:(y-age_lag),1:3] <- 1 # use data starting year 40, with age lag for all fixed gear
    usefishage[,55,1] <- 0 # no samples in BS fleet in year 55
    usefishlen[,31:39,1:3] <- 1 # use data only from 31 - 39 for all fixed gear
    usefishlen[,c(31,32,35,36,40:43, 45:51, 53, 54, 57:y),4] <- 1 # BS trawl
    usefishlen[,c(31,32,35:37, 39:y),5] <- 1 # AI + GOA trawl
  }

  if(faa_n_fish_fleets == 6) {
    usefishage[,40:(y-age_lag),1:3] <- 1 # use data starting year 40, with age lag for all fixed gear
    usefishage[,55,1] <- 0 # no samples in BS fleet in year 55
    usefishlen[,31:39,1:3] <- 1 # use data only from 31 - 39 for all fixed gear
    usefishlen[,c(31,32,35,36,40:43, 45:51, 53, 54, 57:y),4] <- 1 # BS trawl
    usefishlen[,c(31, 35, 36, 44, 46, 58, 61:y),5] <- 1 # AI trawl
    usefishlen[,c(31,32,35:37, 39:y),6] <- 1 # GOA trawl
  }

  # Survey Index
  if(lls_design_type == 'all') { # if design = all, then fitting data for all fleets
    usesrvage[,31:(y-age_lag),1:3] <- usesrvidx[,31:y,1:3] <- 1
    usesrvage[,65,1:3] <- usesrvidx[,65,1:3] <- 0 # no data in year 65
  }

  if(lls_design_type == 'historical') {
    usesrvage[,seq(38, 64, 2),1] <- usesrvidx[,seq(38, 64, 2),1] <- 1 # use data in bs domestic fleet in even years
    usesrvage[,seq(37, 64, 2),2] <- usesrvidx[,seq(37, 64, 2),2] <- 1 # use data in ai domestic fleet in odd years
    usesrvidx[,31:y,3] <- 1 # use data in goa domestic fleet every single year
    usesrvage[,37:(y-age_lag),3] <- 1 # use data in goa domestic fleet every single year
    usesrvidx[,65,] <- usesrvage[,65,] <- 0 # no survey in 65

    # use indicators in subsequent years
    bs_idx_indicator <- even_yrs[even_yrs <= y]
    ai_idx_indicator <- odd_yrs[odd_yrs <= y]
    bs_age_indicator <- even_yrs[even_yrs <= y - age_lag]
    ai_age_indicator <- odd_yrs[odd_yrs <= y - age_lag]

    if(sum(bs_idx_indicator) != 0) usesrvidx[,bs_idx_indicator,1] <- 1 # use bs domestic fleet in even years
    if(sum(ai_idx_indicator) != 0) usesrvidx[,ai_idx_indicator,2] <- 1 # use ai domestic fleet in odd years
    if(sum(bs_age_indicator) != 0) usesrvage[,bs_age_indicator,1] <- 1 # use bs domestic fleet in even years
    if(sum(ai_age_indicator) != 0) usesrvage[,ai_age_indicator,2] <- 1 # use ai domestic fleet in odd years
  }

  if(lls_design_type == 'current') {
    # set up the historical indicators first
    usesrvage[,seq(38, 64, 2),1] <- usesrvidx[,seq(38, 64, 2),1] <- 1 # use data in bs domestic fleet in even years
    usesrvage[,seq(37, 64, 2),2] <- usesrvidx[,seq(37, 64, 2),2] <- 1 # use data in ai domestic fleet in odd years
    usesrvidx[,31:64,3] <- 1 # use data in goa domestic fleet every single year
    usesrvage[,37:64,3] <- 1 # use data in goa domestic fleet every single year
    usesrvidx[,65,] <- usesrvage[,65,] <- 0 # no survey in 65

    # setup current indicators
    goa_idx_indicator <- even_yrs[even_yrs <= y]
    bsai_idx_indicator <- odd_yrs[odd_yrs <= y]
    goa_age_indicator <- even_yrs[even_yrs <= y - age_lag]
    bsai_age_indicator <- odd_yrs[odd_yrs <= y - age_lag]

    if(sum(goa_idx_indicator) != 0) usesrvidx[,goa_idx_indicator,3] <- 1 # use goa domestic fleet in even years
    if(sum(bsai_idx_indicator) != 0) usesrvidx[,bsai_age_indicator,c(1,2)] <- 1 # use bsai domestic fleet in odd years
    if(sum(goa_age_indicator) != 0) usesrvage[,goa_age_indicator,3] <- 1 # use goa domestic fleet in even years
    if(sum(bsai_age_indicator) != 0) usesrvage[,bsai_age_indicator,c(1,2)] <- 1 # use bsai domestic fleet in odd years
  }

  # Figure out Japanese LLS
  if(faa_n_srv_fleets == 4) {
    usesrvidx[,20:35,4:4] <- 1 # survey index
    usesrvage[,c(22, seq(26, 34, 2)),4] <- 1 # age compositions
  }

  if(faa_n_srv_fleets == 6) {
    # compositions
    usesrvage[,c(26, 30, 34),4] <- 1 # BS
    usesrvage[,c(26, 28, 45),5] <- 1 # AI
    usesrvage[,c(22, seq(26, 34, 2)),6] <- 1 # GOA
    usesrvidx[,20:35,4:6] <- 1 # survey index
  }

  return(list(usefishage = usefishage, usefishlen = usefishlen,
              usesrvage = usesrvage, usesrvlen = usesrvlen,
              usecatch = usecatch, usesrvidx = usesrvidx))

}

#' Construct Selectivity Prior Structure
#'
#' This function constructs a prior structure for selectivity parameters by
#' expanding all combinations of fleets and their corresponding blocks, then
#' merging the result with a table of sex-specific parameters. The function
#' returns a data frame containing all fleet–block–sex combinations along with
#' region identifiers and prior parameters.
#'
#' @param fleets_blocks A named list in which each element corresponds to a
#'   fleet and contains a numeric vector of block identifiers. For example:
#'   \code{list("1" = 1:3, "2" = 1:3, "4" = 1)}.
#' @param sex_par A data frame containing the sex-specific parameters that are
#'   merged with the fleet–block structure.
#' @param region A numeric region identifier to be included in the output.
#'   Defaults to 1.
#' @param mu A numeric prior mean applied uniformly across all rows. Defaults to 1.
#' @param sd A numeric prior standard deviation applied uniformly across all rows.
#'   Defaults to 2.
#'
#' @returns A data frame containing all combinations of fleets, blocks, and
#'   sex parameters, augmented with the region identifier and prior parameters.
#'
build_selex_prior <- function(fleets_blocks, sex_par, region = 1, mu = 1, sd = 2) {

  # convert list to long data.frame of fleet/block combinations
  fleet_block_df <- do.call(rbind, lapply(names(fleets_blocks), function(fleet_id) {
    data.frame(
      fleet = as.numeric(fleet_id),
      block = fleets_blocks[[fleet_id]],
      stringsAsFactors = FALSE
    )
  }))

  # merge with sex_par to create all fleet x block x sex combinations
  selex_structure <- merge(fleet_block_df, sex_par)

  # attach region and priors
  out <- cbind(
    region = region,
    selex_structure,
    mu = mu,
    sd = sd
  )

  rownames(out) <- NULL
  return(out)
}



#' Aggregate simulated data to three assessment regions
#'
#' Aggregates simulated fishery, survey, and tagging data from the full spatial
#' resolution to three assessment regions: Bering Sea (BS), Aleutian Islands (AI),
#' and Gulf of Alaska (GOA). GOA data are formed by aggregating regions 3–5 using
#' either biomass- or numbers-based weighting, depending on user specification.
#'
#' The function updates the \code{sim_env} object in place by creating aggregated
#' observations, true values, effective sample sizes, and tag-release structures.
#' Observation error is applied to aggregated catches and survey indices.
#'
#' @param sim_data A list containing simulated observation-level data, including
#'   age and length compositions, survey indices, and tagging information.
#'
#' @param sim_env A list containing simulation settings and true underlying values.
#'   This object is modified in place to include aggregated quantities.
#'
#' @param y Integer. The terminal year of data to be aggregated.
#'
#' @param sim Integer. Simulation replicate index.
#'
#' @param srv_wgt Character string specifying survey aggregation weights.
#'   Either \code{"numbers"} or \code{"biomass"}.
#'
#' @param fish_wgt Character string specifying fishery aggregation weights.
#'   Either \code{"numbers"} or \code{"biomass"}.
#'
#' @param srv_idx_se Numeric. Log-scale standard error used to generate observation
#'   error for aggregated survey indices.
#'
agg_data_to_three_rg <- function(sim_data,
                                 sim_env,
                                 y,
                                 sim,
                                 srv_wgt = 'numbers',
                                 fish_wgt = 'biomass',
                                 srv_idx_se) {


  # get years in simulation
  x <- (sim_env$feedback_start_yr + 1):sim_env$n_yrs
  odd_yrs <- x[x %% 2 == 1] # bsai years
  even_yrs <- x[x %% 2 == 0] # goa years

  # dimensions
  n_regions <- 3
  n_yrs <- length(1:y)
  n_sexes <- sim_env$n_sexes
  n_ages <- sim_env$n_ages
  n_lens <- sim_env$n_lens
  n_fish_fleets <- sim_env$n_fish_fleets
  n_srv_fleets <- sim_env$n_srv_fleets

  # Catches (using true catches, and then applying error to aggregated true catches)
  if(y == sim_env$feedback_start_yr) {
    aggregated_true_catch <- array(0, dim = c(3, length(1:y), n_fish_fleets))
    aggregated_true_catch[1,,] <- sim_env$TrueCatch[1,1:y,,sim] # BS
    aggregated_true_catch[2,,] <- sim_env$TrueCatch[2,1:y,,sim] # AI
    aggregated_true_catch[3,,] <- apply(sim_env$TrueCatch[3:5,1:y,,sim], 2:3, sum) # GOA
    sim_env$Agg_ObsCatch[,1:y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0, mean(exp(sim_data$ln_sigmaC)))) # add error
  } else{
    aggregated_true_catch <- array(0, dim = c(n_regions, n_fish_fleets))
    aggregated_true_catch[1,] <- sim_env$TrueCatch[1,y,,sim] # BS Gear
    aggregated_true_catch[2,] <- sim_env$TrueCatch[2,y,,sim] # AI
    aggregated_true_catch[3,] <- colSums(sim_env$TrueCatch[3:5,y,,sim]) # GOA
    sim_env$Agg_ObsCatch[,y,,sim] <- aggregated_true_catch * exp(rnorm(length(aggregated_true_catch), 0, mean(exp(sim_data$ln_sigmaC))))
  }

  if(fish_wgt == 'biomass') {
    # Get fixed gear weights
    fixed_gear_wts <- sim_env$TrueCatch[3:5,1:y,1,sim] # only weighting GOA for fixed-gear weights
    total_fixed_gear_catches <- colSums(fixed_gear_wts)
    catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

    # get trawl gear weights
    trawl_gear_wts <- sim_env$TrueCatch[3:5,1:y,2,sim] # only weighting GOA
    total_trawl_gear_catches <- colSums(trawl_gear_wts)
    catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
  }

  if(fish_wgt == 'numbers') {
    # Get fixed gear weights
    fixed_gear_wts <- sim_env$CAA[3:5,1:y,,,1,sim] # only weighting GOA for fixed-gear weights
    fixed_gear_wts <- apply(fixed_gear_wts, c(1,2), sum) # sum numbers up
    total_fixed_gear_catches <- colSums(fixed_gear_wts)
    catch_fixed_gear_prop <- sweep(fixed_gear_wts, 2, total_fixed_gear_catches, '/')

    # get trawl gear weights
    trawl_gear_wts <- sim_env$CAA[3:5,1:y,,,2,sim] # only weighting GOA
    trawl_gear_wts <- apply(trawl_gear_wts, c(1,2), sum) # sum numbers up
    total_trawl_gear_catches <- colSums(trawl_gear_wts)
    catch_trawl_gear_prop <- sweep(trawl_gear_wts, 2, total_trawl_gear_catches, '/')
  }

  for(yr in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      # figure out which regions have catches
      if(f == 1) catch_prop <- catch_fixed_gear_prop
      if(f == 2) catch_prop <- catch_trawl_gear_prop

      # age comps
      if(f == 1) { # fixed gear

        regions <- (3:5)[-which(catch_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing catch_prop

        # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
        sim_env$Agg_ObsFishAgeComps[1,yr,,,1,sim] <- sim_data$ObsFishAgeComps[1,yr,,,1] # BS
        sim_env$Agg_ObsFishAgeComps[2,yr,,,1,sim] <- sim_data$ObsFishAgeComps[2,yr,,,1] # AI
        sim_env$Agg_ObsFishAgeComps[1,yr,,,1,sim] <- sim_env$Agg_ObsFishAgeComps[1,yr,,,1,sim] / sum(sim_env$Agg_ObsFishAgeComps[1,yr,,,1,sim])
        sim_env$Agg_ObsFishAgeComps[2,yr,,,1,sim] <- sim_env$Agg_ObsFishAgeComps[2,yr,,,1,sim] / sum(sim_env$Agg_ObsFishAgeComps[2,yr,,,1,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsFishAgeComps[3,yr,,,1,sim] <- apply(sim_data$ObsFishAgeComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsFishAgeComps[3,yr,,,1,sim] <- sim_env$Agg_ObsFishAgeComps[3,yr,,,1,sim] / sum(sim_env$Agg_ObsFishAgeComps[3,yr,,,1,sim])

        # Calculate n_eff for age comps
        sim_env$Agg_ISS_FishAgeComps[1,yr,1,1,sim] <- sum(sim_data$ObsFishAgeComps[1,yr,,,1]) # sum ISS for BS fixed gear
        sim_env$Agg_ISS_FishAgeComps[2,yr,1,1,sim] <- sum(sim_data$ObsFishAgeComps[2,yr,,,1]) # sum ISS for AI fixed gear

        # get n_eff for GOA fixed gear
        if(length(regions) > 1) {
          regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, sum)
          regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,1], 1, regional_n_age, "/")
          pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
          numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
          denominator_age <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
          sim_env$Agg_ISS_FishAgeComps[3,yr,1,1,sim] <- numerator_age / denominator_age
        } else {
          sim_env$Agg_ISS_FishAgeComps[3,yr,1,1,sim] <- sum(sim_data$ObsFishAgeComps[regions,yr,,,1])
        }
      }

      # trawl age comps (not used)
      if(f == 2) {
        regions <- (3:5)[-which(catch_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing catch_prop

        # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
        sim_env$Agg_ObsFishAgeComps[1,yr,,,2,sim] <- sim_data$ObsFishAgeComps[1,yr,,,2] # BS
        sim_env$Agg_ObsFishAgeComps[2,yr,,,2,sim] <- sim_data$ObsFishAgeComps[2,yr,,,2] # AI
        sim_env$Agg_ObsFishAgeComps[1,yr,,,2,sim] <- sim_env$Agg_ObsFishAgeComps[1,yr,,,2,sim] / sum(sim_env$Agg_ObsFishAgeComps[1,yr,,,2,sim])
        sim_env$Agg_ObsFishAgeComps[2,yr,,,2,sim] <- sim_env$Agg_ObsFishAgeComps[2,yr,,,2,sim] / sum(sim_env$Agg_ObsFishAgeComps[2,yr,,,2,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsFishAgeComps[3,yr,,,2,sim] <- apply(sim_data$ObsFishAgeComps[3:5,yr,,,2] * catch_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsFishAgeComps[3,yr,,,2,sim] <- sim_env$Agg_ObsFishAgeComps[3,yr,,,2,sim] / sum(sim_env$Agg_ObsFishAgeComps[3,yr,,,2,sim])

        # Calculate n_eff for age comps
        sim_env$Agg_ISS_FishAgeComps[1,yr,1,2,sim] <- sum(sim_data$ObsFishAgeComps[1,yr,,,2]) # sum ISS for BS fixed gear
        sim_env$Agg_ISS_FishAgeComps[2,yr,1,2,sim] <- sum(sim_data$ObsFishAgeComps[2,yr,,,2]) # sum ISS for AI fixed gear

        # get n_eff for GOA fixed gear
        if(length(regions) > 1) {
          regional_n_age <- apply(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, sum)
          regional_p_age <- sweep(sim_data$ObsFishAgeComps[regions,yr,,,2], 1, regional_n_age, "/")
          pooled_comp_age <- colSums(sweep(regional_p_age, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
          numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
          denominator_age <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
          sim_env$Agg_ISS_FishAgeComps[3,yr,1,2,sim] <- numerator_age / denominator_age
        } else {
          sim_env$Agg_ISS_FishAgeComps[3,yr,1,2,sim] <- sum(sim_data$ObsFishAgeComps[regions,yr,,,2])
        }
      }

      # len comps
      if(f == 1) { # fixed gear

        regions <- (3:5)[-which(catch_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing catch_prop

        # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
        sim_env$Agg_ObsFishLenComps[1,yr,,,1,sim] <- sim_data$ObsFishLenComps[1,yr,,,1] # BS
        sim_env$Agg_ObsFishLenComps[2,yr,,,1,sim] <- sim_data$ObsFishLenComps[2,yr,,,1] # AI
        sim_env$Agg_ObsFishLenComps[1,yr,,,1,sim] <- sim_env$Agg_ObsFishLenComps[1,yr,,,1,sim] / sum(sim_env$Agg_ObsFishLenComps[1,yr,,,1,sim])
        sim_env$Agg_ObsFishLenComps[2,yr,,,1,sim] <- sim_env$Agg_ObsFishLenComps[2,yr,,,1,sim] / sum(sim_env$Agg_ObsFishLenComps[2,yr,,,1,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsFishLenComps[3,yr,,,1,sim] <- apply(sim_data$ObsFishLenComps[3:5,yr,,,1] * catch_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsFishLenComps[3,yr,,,1,sim] <- sim_env$Agg_ObsFishLenComps[3,yr,,,1,sim] / sum(sim_env$Agg_ObsFishLenComps[3,yr,,,1,sim])

        # Calculate n_eff for len comps
        sim_env$Agg_ISS_FishLenComps[1,yr,1,1,sim] <- sum(sim_data$ObsFishLenComps[1,yr,,,1]) # sum ISS for BS fixed gear
        sim_env$Agg_ISS_FishLenComps[2,yr,1,1,sim] <- sum(sim_data$ObsFishLenComps[2,yr,,,1]) # sum ISS for AI fixed gear

        # get n_eff for GOA fixed gear
        if(length(regions) > 1) {
          regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,1], 1, sum)
          regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,1], 1, regional_n_len, "/")
          pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
          numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
          denominator_len <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
          sim_env$Agg_ISS_FishLenComps[3,yr,1,1,sim] <- numerator_len / denominator_len
        } else {
          sim_env$Agg_ISS_FishLenComps[3,yr,1,1,sim] <- sum(sim_data$ObsFishLenComps[regions,yr,,,1])
        }
      }

      # trawl len comps (not used)
      if(f == 2) {
        regions <- (3:5)[-which(catch_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing catch_prop

        # input BS, AI, into appropriate fleets for fixed-gear fishery (no weighting)
        sim_env$Agg_ObsFishLenComps[1,yr,,,2,sim] <- sim_data$ObsFishLenComps[1,yr,,,2] # BS
        sim_env$Agg_ObsFishLenComps[2,yr,,,2,sim] <- sim_data$ObsFishLenComps[2,yr,,,2] # AI
        sim_env$Agg_ObsFishLenComps[1,yr,,,2,sim] <- sim_env$Agg_ObsFishLenComps[1,yr,,,2,sim] / sum(sim_env$Agg_ObsFishLenComps[1,yr,,,2,sim])
        sim_env$Agg_ObsFishLenComps[2,yr,,,2,sim] <- sim_env$Agg_ObsFishLenComps[2,yr,,,2,sim] / sum(sim_env$Agg_ObsFishLenComps[2,yr,,,2,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsFishLenComps[3,yr,,,2,sim] <- apply(sim_data$ObsFishLenComps[3:5,yr,,,2] * catch_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsFishLenComps[3,yr,,,2,sim] <- sim_env$Agg_ObsFishLenComps[3,yr,,,2,sim] / sum(sim_env$Agg_ObsFishLenComps[3,yr,,,2,sim])

        # Calculate n_eff for len comps
        sim_env$Agg_ISS_FishLenComps[1,yr,1,2,sim] <- sum(sim_data$ObsFishLenComps[1,yr,,,2]) # sum ISS for BS fixed gear
        sim_env$Agg_ISS_FishLenComps[2,yr,1,2,sim] <- sum(sim_data$ObsFishLenComps[2,yr,,,2]) # sum ISS for AI fixed gear

        # get n_eff for GOA fixed gear
        if(length(regions) > 1) {
          regional_n_len <- apply(sim_data$ObsFishLenComps[regions,yr,,,2], 1, sum)
          regional_p_len <- sweep(sim_data$ObsFishLenComps[regions,yr,,,2], 1, regional_n_len, "/")
          pooled_comp_len <- colSums(sweep(regional_p_len, 1, catch_prop[regions_offset,yr], "*")) / sum(catch_prop[,yr])
          numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
          denominator_len <- sum(catch_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
          sim_env$Agg_ISS_FishLenComps[3,yr,1,2,sim] <- numerator_len / denominator_len
        } else {
          sim_env$Agg_ISS_FishLenComps[3,yr,1,2,sim] <- sum(sim_data$ObsFishLenComps[regions,yr,,,2])
        }
      }
    } # end f loop
  } # end yr loop

  # Survey Index
  if(y == sim_env$feedback_start_yr) {
    true_srvidx <- sim_env$TrueSrvIdx[,1:y,,sim] # get true survey index
    sim_env$Agg_TrueSrvIdx[1,1:y,1,sim] <- true_srvidx[1,1:y,1] # BS Domestic
    sim_env$Agg_TrueSrvIdx[2,1:y,1,sim] <- true_srvidx[2,1:y,1] # AI Domestic
    sim_env$Agg_TrueSrvIdx[3,1:y,1,sim] <- colSums(true_srvidx[3:5,1:y,1]) # GOA Domestic
    sim_env$Agg_TrueSrvIdx[1,1:y,3,sim] <- true_srvidx[1,1:y,3] # BS JP
    sim_env$Agg_TrueSrvIdx[2,1:y,3,sim] <- true_srvidx[2,1:y,3] # AI JP
    sim_env$Agg_TrueSrvIdx[3,1:y,3,sim] <- colSums(true_srvidx[3:5,1:y,3]) # GOA JP
    # add devs
    sim_env$Agg_ObsSrvIdx[,1:y,,sim] <- sim_env$Agg_TrueSrvIdx[,1:y,,sim] * exp(rnorm(length(sim_env$Agg_TrueSrvIdx[,1:y,,sim]), 0, srv_idx_se)) # simulate devs
  } else {
    true_srvidx <- sim_env$TrueSrvIdx[,y,,sim] # get true survey index
    agg_true_srvidx <- array(NA, dim = c(n_regions, n_srv_fleets))
    agg_true_srvidx[1,1] <- true_srvidx[1,1] # BS Domestic
    agg_true_srvidx[2,1] <- true_srvidx[2,1] # AI Domestic
    agg_true_srvidx[3,1] <- sum(true_srvidx[3:5,1]) # GOA Domestic
    agg_true_srvidx[1,3] <- true_srvidx[1,3] # BS JP
    agg_true_srvidx[2,3] <- true_srvidx[2,3] # AI JP
    agg_true_srvidx[3,3] <- sum(true_srvidx[3:5,3]) # GOA JP
    # save true aggregated index
    sim_env$Agg_TrueSrvIdx[,y,,sim] <- agg_true_srvidx # save "true index"
    sim_env$Agg_ObsSrvIdx[,y,,sim] <- agg_true_srvidx * exp(rnorm(length(agg_true_srvidx), 0, srv_idx_se))
  }

  if(srv_wgt == 'biomass') {
    # Get domestic weights
    domestic_wts <- apply(sim_env$SrvIAA[3:5,1:y,,,1,sim] * sim_env$WAA_srv[3:5,1:y,,,1,sim], c(1,2), sum)  # only weighting GOA for domestic weights
    total_domestic_idx <- colSums(domestic_wts)
    idx_domestic_prop <- sweep(domestic_wts, 2, total_domestic_idx, '/')

    # get jp weights
    jp_wts <- apply(sim_env$SrvIAA[3:5,1:y,,,3,sim] * sim_env$WAA_srv[3:5,1:y,,,3,sim], c(1,2), sum)
    total_jp_idx <- colSums(jp_wts)
    idx_jp_prop <- sweep(jp_wts, 2, total_jp_idx, '/')
  }

  if(srv_wgt == 'numbers') {
    # Get domestic weights
    domestic_wts <- sim_env$SrvIAA[3:5,1:y,,,1,sim] # only weighting GOA for domestic weights
    domestic_wts <- apply(domestic_wts, c(1,2), sum) # sum numbers up
    total_domestic_idx <- colSums(domestic_wts)
    idx_domestic_prop <- sweep(domestic_wts, 2, total_domestic_idx, '/')

    # get jp weights
    jp_wts <- sim_env$SrvIAA[3:5,1:y,,,3,sim]
    jp_wts <- apply(jp_wts, c(1,2), sum) # sum numbers up
    total_jp_idx <- colSums(jp_wts)
    idx_jp_prop <- sweep(jp_wts, 2, total_jp_idx, '/')
  }

  for(yr in 1:n_yrs) {
    for(f in c(1,3)) { # loop through only domestic and jp

      # figure out which regions have catches
      if(f == 1) idx_prop <- idx_domestic_prop
      if(f == 2) idx_prop <- idx_jp_prop

      # age comps
      if(f == 1) { # domestic

        regions <- (3:5)[-which(idx_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing idx_prop

        # input BS, AI, into appropriate fleets for domestic srvery (no weighting)
        sim_env$Agg_ObsSrvAgeComps[1,yr,,,1,sim] <- sim_data$ObsSrvAgeComps[1,yr,,,1] # BS
        sim_env$Agg_ObsSrvAgeComps[2,yr,,,1,sim] <- sim_data$ObsSrvAgeComps[2,yr,,,1] # AI
        sim_env$Agg_ObsSrvAgeComps[1,yr,,,1,sim] <- sim_env$Agg_ObsSrvAgeComps[1,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvAgeComps[1,yr,,,1,sim])
        sim_env$Agg_ObsSrvAgeComps[2,yr,,,1,sim] <- sim_env$Agg_ObsSrvAgeComps[2,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvAgeComps[2,yr,,,1,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsSrvAgeComps[3,yr,,,1,sim] <- apply(sim_data$ObsSrvAgeComps[3:5,yr,,,1] * idx_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsSrvAgeComps[3,yr,,,1,sim] <- sim_env$Agg_ObsSrvAgeComps[3,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvAgeComps[3,yr,,,1,sim])

        # Calculate n_eff for age comps
        sim_env$Agg_ISS_SrvAgeComps[1,yr,1,1,sim] <- sum(sim_data$ObsSrvAgeComps[1,yr,,,1]) # sum ISS for BS domestic
        sim_env$Agg_ISS_SrvAgeComps[2,yr,1,1,sim] <- sum(sim_data$ObsSrvAgeComps[2,yr,,,1]) # sum ISS for AI domestic

        # get n_eff for GOA domestic
        if(length(regions) > 1) {
          regional_n_age <- apply(sim_data$ObsSrvAgeComps[regions,yr,,,1], 1, sum)
          regional_p_age <- sweep(sim_data$ObsSrvAgeComps[regions,yr,,,1], 1, regional_n_age, "/")
          pooled_comp_age <- colSums(sweep(regional_p_age, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
          numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
          denominator_age <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
          sim_env$Agg_ISS_SrvAgeComps[3,yr,1,1,sim] <- numerator_age / denominator_age
        } else {
          sim_env$Agg_ISS_SrvAgeComps[3,yr,1,1,sim] <- sum(sim_data$ObsSrvAgeComps[regions,yr,,,1])
        }
      }

      # jp age comps
      if(f == 3) {
        regions <- (3:5)[-which(idx_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing idx_prop

        # input BS, AI, into appropriate fleets for jp-gear srvery (no weighting)
        sim_env$Agg_ObsSrvAgeComps[1,yr,,,3,sim] <- sim_data$ObsSrvAgeComps[1,yr,,,3] # BS
        sim_env$Agg_ObsSrvAgeComps[2,yr,,,3,sim] <- sim_data$ObsSrvAgeComps[2,yr,,,3] # AI
        sim_env$Agg_ObsSrvAgeComps[1,yr,,,3,sim] <- sim_env$Agg_ObsSrvAgeComps[1,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvAgeComps[1,yr,,,3,sim])
        sim_env$Agg_ObsSrvAgeComps[2,yr,,,3,sim] <- sim_env$Agg_ObsSrvAgeComps[2,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvAgeComps[2,yr,,,3,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsSrvAgeComps[3,yr,,,3,sim] <- apply(sim_data$ObsSrvAgeComps[3:5,yr,,,3] * idx_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsSrvAgeComps[3,yr,,,3,sim] <- sim_env$Agg_ObsSrvAgeComps[3,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvAgeComps[,yr,,,3,sim])

        # Calculate n_eff for age comps
        sim_env$Agg_ISS_SrvAgeComps[1,yr,1,3,sim] <- sum(sim_data$ObsSrvAgeComps[1,yr,,,3]) # sum ISS for BS jp
        sim_env$Agg_ISS_SrvAgeComps[2,yr,1,3,sim] <- sum(sim_data$ObsSrvAgeComps[2,yr,,,3]) # sum ISS for AI jp

        # get n_eff for GOA jp
        if(length(regions) > 1) {
          regional_n_age <- apply(sim_data$ObsSrvAgeComps[regions,yr,,,3], 1, sum)
          regional_p_age <- sweep(sim_data$ObsSrvAgeComps[regions,yr,,,3], 1, regional_n_age, "/")
          pooled_comp_age <- colSums(sweep(regional_p_age, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
          numerator_age <- sum(pooled_comp_age * (1 - pooled_comp_age))
          denominator_age <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_age) * apply(regional_p_age * (1 - regional_p_age), 1, sum))
          sim_env$Agg_ISS_SrvAgeComps[3,yr,1,3,sim] <- numerator_age / denominator_age
        } else {
          sim_env$Agg_ISS_SrvAgeComps[3,yr,1,3,sim] <- sum(sim_data$ObsSrvAgeComps[regions,yr,,,3])
        }
      }

      # len comps (not used)
      if(f == 1) { # domestic

        regions <- (3:5)[-which(idx_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing idx_prop

        # input BS, AI, into appropriate fleets for domestic srvery (no weighting)
        sim_env$Agg_ObsSrvLenComps[1,yr,,,1,sim] <- sim_data$ObsSrvLenComps[1,yr,,,1] # BS
        sim_env$Agg_ObsSrvLenComps[2,yr,,,1,sim] <- sim_data$ObsSrvLenComps[2,yr,,,1] # AI
        sim_env$Agg_ObsSrvLenComps[1,yr,,,1,sim] <- sim_env$Agg_ObsSrvLenComps[1,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvLenComps[1,yr,,,1,sim])
        sim_env$Agg_ObsSrvLenComps[2,yr,,,1,sim] <- sim_env$Agg_ObsSrvLenComps[2,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvLenComps[2,yr,,,1,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsSrvLenComps[3,yr,,,1,sim] <- apply(sim_data$ObsSrvLenComps[3:5,yr,,,1] * idx_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsSrvLenComps[3,yr,,,1,sim] <- sim_env$Agg_ObsSrvLenComps[3,yr,,,1,sim] / sum(sim_env$Agg_ObsSrvLenComps[3,yr,,,1,sim])

        # Calculate n_eff for len comps
        sim_env$Agg_ISS_SrvLenComps[1,yr,1,1,sim] <- sum(sim_data$ObsSrvLenComps[1,yr,,,1]) # sum ISS for BS domestic
        sim_env$Agg_ISS_SrvLenComps[2,yr,1,1,sim] <- sum(sim_data$ObsSrvLenComps[2,yr,,,1]) # sum ISS for AI domestic

        # get n_eff for GOA domestic
        if(length(regions) > 1) {
          regional_n_len <- apply(sim_data$ObsSrvLenComps[regions,yr,,,1], 1, sum)
          regional_p_len <- sweep(sim_data$ObsSrvLenComps[regions,yr,,,1], 1, regional_n_len, "/")
          pooled_comp_len <- colSums(sweep(regional_p_len, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
          numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
          denominator_len <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
          sim_env$Agg_ISS_SrvLenComps[3,yr,1,1,sim] <- numerator_len / denominator_len
        } else {
          sim_env$Agg_ISS_SrvLenComps[3,yr,1,1,sim] <- sum(sim_data$ObsSrvLenComps[regions,yr,,,1])
        }
      }

      # jp len comps
      if(f == 3) {
        regions <- (3:5)[-which(idx_prop[,yr] == 0)]
        if(length(regions) == 0) regions <- 3:5
        regions_offset <- regions - 2 # get offset for indexing idx_prop

        # input BS, AI, into appropriate fleets for jp-gear srvery (no weighting)
        sim_env$Agg_ObsSrvLenComps[1,yr,,,3,sim] <- sim_data$ObsSrvLenComps[1,yr,,,3] # BS
        sim_env$Agg_ObsSrvLenComps[2,yr,,,3,sim] <- sim_data$ObsSrvLenComps[2,yr,,,3] # AI
        sim_env$Agg_ObsSrvLenComps[1,yr,,,3,sim] <- sim_env$Agg_ObsSrvLenComps[1,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvLenComps[1,yr,,,3,sim])
        sim_env$Agg_ObsSrvLenComps[2,yr,,,3,sim] <- sim_env$Agg_ObsSrvLenComps[2,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvLenComps[2,yr,,,3,sim])

        # input GOA fleet weighting by proportions
        sim_env$Agg_ObsSrvLenComps[3,yr,,,3,sim] <- apply(sim_data$ObsSrvLenComps[3:5,yr,,,3] * idx_prop[,yr], c(2,3), sum)
        sim_env$Agg_ObsSrvLenComps[3,yr,,,3,sim] <- sim_env$Agg_ObsSrvLenComps[3,yr,,,3,sim] / sum(sim_env$Agg_ObsSrvLenComps[3,yr,,,1,sim])

        # Calculate n_eff for len comps
        sim_env$Agg_ISS_SrvLenComps[1,yr,1,3,sim] <- sum(sim_data$ObsSrvLenComps[1,yr,,,3]) # sum ISS for BS jp
        sim_env$Agg_ISS_SrvLenComps[2,yr,1,3,sim] <- sum(sim_data$ObsSrvLenComps[2,yr,,,3]) # sum ISS for AI jp

        # get n_eff for GOA jp
        if(length(regions) > 1) {
          regional_n_len <- apply(sim_data$ObsSrvLenComps[regions,yr,,,3], 1, sum)
          regional_p_len <- sweep(sim_data$ObsSrvLenComps[regions,yr,,,3], 1, regional_n_len, "/")
          pooled_comp_len <- colSums(sweep(regional_p_len, 1, idx_prop[regions_offset,yr], "*")) / sum(idx_prop[,yr])
          numerator_len <- sum(pooled_comp_len * (1 - pooled_comp_len))
          denominator_len <- sum(idx_prop[regions_offset,yr]^2 * (1 / regional_n_len) * apply(regional_p_len * (1 - regional_p_len), 1, sum))
          sim_env$Agg_ISS_SrvLenComps[3,yr,1,3,sim] <- numerator_len / denominator_len
        } else {
          sim_env$Agg_ISS_SrvLenComps[3,yr,1,3,sim] <- sum(sim_data$ObsSrvLenComps[regions,yr,,,3])
        }
      }

    } # end f loop
  } # end yr loop

  # Tag Releases (modifieds sim_data)
  # GOA processing
  goa_releases <- sim_data$tag_release_indicator[sim_data$tag_release_indicator[,1] %in% c(3:5),]
  goa_cohorts <- which(sim_data$tag_release_indicator[,1] %in% c(3:5))
  unique_goa_yrs <- unique(goa_releases[,2])

  # Aggregated GOA Tagged Fish and indicators
  agg_goa_tagged_fish <- array(NA, dim = c(length(unique_goa_yrs), 30, 2))
  agg_goa_recap <- array(NA, dim = c(15, length(unique_goa_yrs), 3, 30, 2))
  agg_goa_indicator <- matrix(NA, nrow = length(unique_goa_yrs), ncol = 2)

  for(i in 1:length(unique_goa_yrs)) {
    tmp_goa_cohort <- goa_releases[which(goa_releases[,2] == unique_goa_yrs[i]),]
    tag_cohort_match <- which(sim_data$tag_release_indicator[,1] %in% tmp_goa_cohort[,1] &
                                sim_data$tag_release_indicator[,2] %in% tmp_goa_cohort[,2])
    agg_goa_tagged_fish[i,,] <- apply(sim_data$Tagged_Fish[tag_cohort_match,,], c(2:3), sum)
    # aggregate recaptures
    agg_goa_recap[,i,1,,] <- apply(sim_data$Obs_Tag_Recap[,tag_cohort_match,1,,], c(1,3,4), sum)
    agg_goa_recap[,i,2,,] <- apply(sim_data$Obs_Tag_Recap[,tag_cohort_match,2,,], c(1,3,4), sum)
    agg_goa_recap[,i,3,,] <- apply(sim_data$Obs_Tag_Recap[,tag_cohort_match,3:5,,], c(1,4,5), sum)
    # Create consolidated GOA indicator row
    agg_goa_indicator[i,] <- c(3, unique_goa_yrs[i])
  }

  # Non-GOA cohorts (keep as-is)
  non_goa_cohorts <- which(!sim_data$tag_release_indicator[,1] %in% c(3:5))
  non_goa_indicator <- sim_data$tag_release_indicator[non_goa_cohorts,]
  non_goa_tagged_fish <- sim_data$Tagged_Fish[non_goa_cohorts,,]
  non_goa_recap <- sim_data$Obs_Tag_Recap[,non_goa_cohorts,,,]

  # Consolidate non-GOA regions to match GOA structure
  non_goa_recap_consolidated <- array(NA, dim = c(15, length(non_goa_cohorts), 3, 30, 2))
  non_goa_recap_consolidated[,,1,,] <- non_goa_recap[,,1,,]
  non_goa_recap_consolidated[,,2,,] <- non_goa_recap[,,2,,]
  non_goa_recap_consolidated[,,3,,] <- apply(non_goa_recap[,,3:5,,], c(1,2,4,5), sum)

  # Combine aggregated GOA + non-GOA
  colnames(agg_goa_indicator) <- c("regions", "tag_yrs")
  sim_env$Agg_tag_release_indicator <- rbind(agg_goa_indicator, non_goa_indicator)
  sim_env$Agg_Tagged_Fish <- abind::abind(agg_goa_tagged_fish, non_goa_tagged_fish, along = 1)
  sim_env$Agg_Obs_Tag_Recap <- abind::abind(agg_goa_recap, non_goa_recap_consolidated, along = 2)

}

#' Generate data-use indicators for three-region assessment models
#'
#' Constructs binary indicators specifying which observations are available for
#' fitting in a three-region (BS, AI, GOA) assessment model. Indicators are created
#' for catches, survey indices, and fishery and survey composition data.
#'
#' The survey data availability pattern depends on the specified longline survey
#' design type and accounts for alternating survey schedules and age-composition
#' lags.
#'
#' @param sim_env A list containing simulation settings, including the feedback
#'   start year and total number of years.
#'
#' @param y Integer. Terminal assessment year.
#'
#' @param age_lag Integer. Number of years between observation year and the most
#'   recent age available for composition data. Default is 1.
#'
#' @param lls_design_type Character string specifying the survey design.
#'   Must be one of \code{"all"}, \code{"historical"}, or \code{"current"}.
#'
#' @return A named list of logical indicator arrays:
#' \itemize{
#'   \item \code{usefishage}: Fishery age-composition indicators
#'   \item \code{usefishlen}: Fishery length-composition indicators
#'   \item \code{usesrvage}: Survey age-composition indicators
#'   \item \code{usesrvlen}: Survey length-composition indicators
#'   \item \code{usecatch}: Catch indicators
#'   \item \code{usesrvidx}: Survey index indicators
#' }
#'
threerg_use_indicators <- function(sim_env,
                               y,
                               age_lag = 1,
                               lls_design_type) {

  # get years in simulation
  x <- (sim_env$feedback_start_yr + 1):sim_env$n_yrs
  odd_yrs <- x[x %% 2 == 1] # bsai years
  even_yrs <- x[x %% 2 == 0] # goa years

  # Data indicators
  usecatch <- array(0, dim = c(3, y, 2))
  usesrvidx <- array(0, dim = c(3, y, 3))
  usefishage <- array(0, c(3, y, 2))
  usefishlen <- array(0, c(3, y, 2))
  usesrvage <- array(0, c(3, y, 3))
  usesrvlen <- array(0, c(3, y, 3)) # not used at all

  usecatch[1,,1] <- 1 # BS fixed gear catch availiable every year
  usecatch[2,-c(1:3),1] <- 1 # AI fixed gear catch not availiable in first three years
  usecatch[3,,1] <- 1 # GOA fixed gear catch availaible in every year
  usecatch[1,,2] <- 1 # BS trawl gear catch availiable every year
  usecatch[2,-c(1:3),2] <- 1 # AI trawl gear catch not availiable in first three years
  usecatch[3,,2] <- 1 # aggregated GOA trawl gear catch availaible in every year

  # Fish Comps
  usefishage[1:3,40:(y-age_lag),1] <- 1 # use data starting year 40, with age lag for all fixed gear
  usefishage[1,55,1] <- 0 # no samples in BS fleet in year 55
  usefishlen[1:3,31:39,1] <- 1 # use data only from 31 - 39 for all fixed gear
  usefishlen[1,c(31,32,35,36,40:43, 45:51, 53, 54, 57:y),2] <- 1 # BS trawl
  usefishlen[2,c(31, 35, 36, 44, 46, 58, 61:y),2] <- 1 # AI trawl
  usefishlen[3,c(31,32,35:37, 39:y),2] <- 1 # GOA trawl

  # Survey Index
  if(lls_design_type == 'all') { # if design = all, then fitting data for all fleets
    usesrvage[1:3,31:(y-age_lag),1] <- usesrvidx[1:3,31:y,1] <- 1
    usesrvage[1:3,65,1] <- usesrvidx[1:3,65,1] <- 0 # no data in year 65
  }

  if(lls_design_type == 'historical') {
    usesrvage[1,seq(38, 64, 2),1] <- usesrvidx[1,seq(38, 64, 2),1] <- 1 # use data in bs domestic fleet in even years
    usesrvage[2,seq(37, 64, 2),1] <- usesrvidx[2,seq(37, 64, 2),1] <- 1 # use data in ai domestic fleet in odd years
    usesrvidx[3,31:y,1] <- 1 # use data in goa domestic fleet every single year
    usesrvage[3,37:(y-age_lag),1] <- 1 # use data in goa domestic fleet every single year
    usesrvidx[,65,] <- usesrvage[,65,] <- 0 # no survey in 65

    # use indicators in subsequent years
    bs_idx_indicator <- even_yrs[even_yrs <= y]
    ai_idx_indicator <- odd_yrs[odd_yrs <= y]
    bs_age_indicator <- even_yrs[even_yrs <= y - age_lag]
    ai_age_indicator <- odd_yrs[odd_yrs <= y - age_lag]

    if(sum(bs_idx_indicator) != 0) usesrvidx[1,bs_idx_indicator,1] <- 1 # use bs domestic fleet in even years
    if(sum(ai_idx_indicator) != 0) usesrvidx[2,ai_idx_indicator,1] <- 1 # use ai domestic fleet in odd years
    if(sum(bs_age_indicator) != 0) usesrvage[1,bs_age_indicator,1] <- 1 # use bs domestic fleet in even years
    if(sum(ai_age_indicator) != 0) usesrvage[2,ai_age_indicator,1] <- 1 # use ai domestic fleet in odd years
  }

  if(lls_design_type == 'current') {
    # set up the historical indicators first
    usesrvage[1,seq(38, 64, 2),1] <- usesrvidx[1,seq(38, 64, 2),1] <- 1 # use data in bs domestic fleet in even years
    usesrvage[2,seq(37, 64, 2),1] <- usesrvidx[2,seq(37, 64, 2),1] <- 1 # use data in ai domestic fleet in odd years
    usesrvidx[3,31:64,1] <- 1 # use data in goa domestic fleet every single year
    usesrvage[3,37:64,1] <- 1 # use data in goa domestic fleet every single year
    usesrvidx[,65,] <- usesrvage[,65,] <- 0 # no survey in 65

    # setup current indicators
    goa_idx_indicator <- even_yrs[even_yrs <= y]
    bsai_idx_indicator <- odd_yrs[odd_yrs <= y]
    goa_age_indicator <- even_yrs[even_yrs <= y - age_lag]
    bsai_age_indicator <- odd_yrs[odd_yrs <= y - age_lag]

    if(sum(goa_idx_indicator) != 0) usesrvidx[3,goa_idx_indicator,1] <- 1 # use goa domestic fleet in even years
    if(sum(bsai_idx_indicator) != 0) usesrvidx[c(1,2),bsai_age_indicator,1] <- 1 # use bsai domestic fleet in odd years
    if(sum(goa_age_indicator) != 0) usesrvage[3,goa_age_indicator,1] <- 1 # use goa domestic fleet in even years
    if(sum(bsai_age_indicator) != 0) usesrvage[c(1,2),bsai_age_indicator,1] <- 1 # use bsai domestic fleet in odd years
  }

  # compositions
  usesrvage[1,c(26, 30, 34),3] <- 1 # BS
  usesrvage[2,c(26, 28, 45),3] <- 1 # AI
  usesrvage[3,c(22, seq(26, 34, 2)),3] <- 1 # GOA
  usesrvidx[1:3,20:35,3] <- 1 # survey index

  return(list(usefishage = usefishage, usefishlen = usefishlen,
              usesrvage = usesrvage, usesrvlen = usesrvlen,
              usecatch = usecatch, usesrvidx = usesrvidx))

}

#' Construct data-availability indicators for the five-region model
#'
#' Creates logical indicator arrays specifying which data streams are used
#' by region, year, and fleet for the five-region estimation model. Indicators
#' are constructed for catch, fishery age and length compositions, survey
#' indices, and survey age and length compositions. Availability depends on
#' region-specific historical sampling patterns, an optional age lag, and the
#' specified longline survey design.
#'
#' The five regions are assumed to follow the ordering:
#' Bering Sea (BS), Aleutian Islands (AI), Western GOA (WGOA),
#' Central GOA (CGOA), and Eastern GOA (EGOA).
#'
#' @param sim_env Simulation environment containing the number of years,
#'   feedback start year, and other operating model settings.
#' @param y Terminal assessment year.
#' @param age_lag Integer specifying the lag (in years) between observation and
#'   recruitment year for age-composition data.
#' @param lls_design_type Character string specifying the longline survey design.
#'   Supported values are:
#'   \itemize{
#'     \item \code{"all"}: all regions sampled every year,
#'     \item \code{"historical"}: alternating BS/AI sampling with continuous GOA,
#'     \item \code{"current"}: current alternating BS/AI versus GOA design.
#'   }
#'
fiverg_use_indicators <- function(sim_env,
                                  y,
                                  age_lag = 1,
                                  lls_design_type
) {

  # get years in simulation
  x <- (sim_env$feedback_start_yr + 1):sim_env$n_yrs
  odd_yrs <- x[x %% 2 == 1] # bsai years
  even_yrs <- x[x %% 2 == 0] # goa years

  # Data indicators
  usecatch <- array(0, dim = c(5, y, 2))
  usesrvidx <- array(0, dim = c(5, y, 3))
  usefishage <- array(0, c(5, y, 2))
  usefishlen <- array(0, c(5, y, 2))
  usesrvage <- array(0, c(5, y, 3))
  usesrvlen <- array(0, c(5, y, 3)) # not used at all

  usecatch[1,,1] <- 1 # BS fixed gear catch availiable every year
  usecatch[2,-c(1:3),1] <- 1 # AI fixed gear catch not availiable in first three years
  usecatch[3,-c(1:3),1] <- 1 # WGOA fixed gear catch
  usecatch[4,-c(1:3),1] <- 1 # CGOA fixed gear catch
  usecatch[5,,1] <- 1 # EGOA fixed gear catch availaible in every year
  usecatch[1,,2] <- 1 # BS trawl gear catch availiable every year
  usecatch[2,-c(1:3),2] <- 1 # AI trawl gear catch not availiable in first three years
  usecatch[3,-c(1:3),2] <- 1 # WGOA trawl gear catch
  usecatch[4,-c(1:3),2] <- 1 # CGOA trawl gear catch
  usecatch[5,-c(23:25),2] <- 1 # EGOA trawl gear catch

  # Fish Comps
  usefishage[1:5,40:(y-age_lag),1] <- 1 # use data starting year 40, with age lag for all fixed gear
  usefishage[1,55,1] <- 0 # no samples in BS fleet in year 55
  usefishlen[1:5,31:39,1] <- 1 # use data only from 31 - 39 for all fixed gear
  usefishlen[1,c(31,32,35,36,40:43, 45:51, 53, 54, 57:y),2] <- 1 # BS trawl
  usefishlen[2,c(31, 35, 36, 44, 46, 58, 61:y),2] <- 1 # AI trawl
  usefishlen[3,c(31,42,43,44,45,46,47,48,50,54,56,59:y),2] <- 1 # WGOA trawl
  usefishlen[4,c(31,32,36:y),2] <- 1 # CGOA trawl
  usefishlen[5,c(31,35,40,42:y),2] <- 1 # CGOA trawl

  # Survey Index
  if(lls_design_type == 'all') { # if design = all, then fitting data for all fleets
    usesrvage[1:5,31:(y-age_lag),1] <- usesrvidx[1:5,31:y,1] <- 1
    usesrvage[1:5,65,1] <- usesrvidx[1:5,65,1] <- 0 # no data in year 65
  }

  if(lls_design_type == 'historical') {
    usesrvage[1,seq(38, 64, 2),1] <- usesrvidx[1,seq(38, 64, 2),1] <- 1 # use data in bs domestic fleet in even years
    usesrvage[2,seq(37, 64, 2),1] <- usesrvidx[2,seq(37, 64, 2),1] <- 1 # use data in ai domestic fleet in odd years
    usesrvidx[3:5,31:y,1] <- 1 # use data in goa domestic fleet every single year
    usesrvage[3:5,37:(y-age_lag),1] <- 1 # use data in goa domestic fleet every single year
    usesrvidx[,65,] <- usesrvage[,65,] <- 0 # no survey in 65

    # use indicators in subsequent years
    bs_idx_indicator <- even_yrs[even_yrs <= y]
    ai_idx_indicator <- odd_yrs[odd_yrs <= y]
    bs_age_indicator <- even_yrs[even_yrs <= y - age_lag]
    ai_age_indicator <- odd_yrs[odd_yrs <= y - age_lag]

    if(sum(bs_idx_indicator) != 0) usesrvidx[1,bs_idx_indicator,1] <- 1 # use bs domestic fleet in even years
    if(sum(ai_idx_indicator) != 0) usesrvidx[2,ai_idx_indicator,1] <- 1 # use ai domestic fleet in odd years
    if(sum(bs_age_indicator) != 0) usesrvage[1,bs_age_indicator,1] <- 1 # use bs domestic fleet in even years
    if(sum(ai_age_indicator) != 0) usesrvage[2,ai_age_indicator,1] <- 1 # use ai domestic fleet in odd years
  }

  if(lls_design_type == 'current') {
    # set up the historical indicators first
    usesrvage[1,seq(38, 64, 2),1] <- usesrvidx[1,seq(38, 64, 2),1] <- 1 # use data in bs domestic fleet in even years
    usesrvage[2,seq(37, 64, 2),1] <- usesrvidx[2,seq(37, 64, 2),1] <- 1 # use data in ai domestic fleet in odd years
    usesrvidx[3:5,31:64,1] <- 1 # use data in goa domestic fleet every single year
    usesrvage[3:5,37:64,1] <- 1 # use data in goa domestic fleet every single year
    usesrvidx[,65,] <- usesrvage[,65,] <- 0 # no survey in 65

    # setup current indicators
    goa_idx_indicator <- even_yrs[even_yrs <= y]
    bsai_idx_indicator <- odd_yrs[odd_yrs <= y]
    goa_age_indicator <- even_yrs[even_yrs <= y - age_lag]
    bsai_age_indicator <- odd_yrs[odd_yrs <= y - age_lag]

    if(sum(goa_idx_indicator) != 0) usesrvidx[3:5,goa_idx_indicator,1] <- 1 # use goa domestic fleet in even years
    if(sum(bsai_idx_indicator) != 0) usesrvidx[c(1,2),bsai_age_indicator,1] <- 1 # use bsai domestic fleet in odd years
    if(sum(goa_age_indicator) != 0) usesrvage[3:5,goa_age_indicator,1] <- 1 # use goa domestic fleet in even years
    if(sum(bsai_age_indicator) != 0) usesrvage[c(1,2),bsai_age_indicator,1] <- 1 # use bsai domestic fleet in odd years
  }

  # compositions
  usesrvage[1,c(26, 30, 34),3] <- 1 # BS
  usesrvage[2,c(26, 28, 45),3] <- 1 # AI
  usesrvage[3:5,c(22, seq(26, 34, 2)),3] <- 1 # GOA
  usesrvidx[1:5,20:35,3] <- 1 # survey index

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
#' @param srv_idx_se Survey idnex SE
#' @param age_lag Age lag
#' @param srv_wgt numbers or biomass comp weighting
#' @param fish_wgt numbers or biomass comp weighting
#'
#' @return The function updates `sim_env` in place (population states, catches,
#'   F rates, and assessment outputs). No value is returned.
run_single_rg_closedloop_i <- function(sim_env,
                                       sim,
                                       fleet_allocation,
                                       lls_design_type,
                                       srv_idx_se,
                                       age_lag,
                                       srv_wgt,
                                       fish_wgt
                                       ) {

  # Run Closed Loop ---------------------------------------------------------
  for(y in 1:sim_env$n_yrs) {

    # Execute annual population dynamics
    run_annual_cycle(y, sim, sim_env)

    # Run feedback
    if(y >= sim_env$feedback_start_yr) {

      ### Assessment --------------------------------------------------------------
      asmt_list <- single_region_em(sim_env = sim_env, y = y, sim = sim, srv_idx_se = srv_idx_se,
                                    age_lag = age_lag, lls_design_type = lls_design_type,
                                    srv_wgt = srv_wgt, fish_wgt = fish_wgt) # get assessment data
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
          calc_rec_st_yr = 20,
          rec_age = 2,
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
      projection <- get_population_projection(data = asmt_list$data, rep = model$rep,
                                              y = y, n_proj_yrs = 2,
                                              f_ref_pt = array(exp_frate, c(asmt_list$data$n_regions, 2))) # get 1 year projection
      tac <- sum(projection$proj_Catch[,2,]) # get tac

      ### Apportionment -----------------------------------------------------------
      # get survey apportionment (survey biomass)
      tmp_Srv_BiomIA <- sim_env$SrvIAA[,((y - 4): y),,,1,sim] * sim_env$WAA_srv[,((y - 4): y),,,1,sim]
      tmp_srvidx_biom <- apply(tmp_Srv_BiomIA, 1:2, sum) # get biomass
      apportionment <- get_single_region_survey_apportionment(feedback_start_yr = sim_env$feedback_start_yr,
                                                              n_yrs = sim_env$n_yrs, y = y,
                                                              srv_idx = tmp_srvidx_biom, # using true survey biomass values to do apportionment
                                                              rolling_avg_yrs = 5,
                                                              lls_design_type = lls_design_type)
      tac_r <- tac * apportionment # get regional tac
      tac_rf <- tac_r * fleet_allocation # allocate regional tac by fleet

      ### TAC to F ----------------------------------------------------------------
      if (y < sim_env$n_yrs) {
        move_age <- if(sim_env$do_recruits_move == 0) 2 else 1
        tmp_naa_moved <- array(NA, dim = dim(sim_env$NAA[, y+1, , , sim]))
        for(a in move_age:sim_env$n_ages) for(s in 1:sim_env$n_sexes)
          tmp_naa_moved[,a,s] <- t(sim_env$NAA[,y+1,a,s,sim]) %*% sim_env$Movement[,,y+1,a,s,sim]
        if(move_age == 2) tmp_naa_moved[,1,] <- sim_env$NAA[, y+1,1, , sim]

        # Solve for all F simultaneously for each region
        for(r in 1:sim_env$n_regions) {
          sim_env$Fmort[r, y+1, , sim] <- solve_multifleet_F(
            target_catch = tac_rf[r, ],
            NAA = tmp_naa_moved[r, , ],
            WAA = sim_env$WAA[r, y+1, , , sim],
            natmort = sim_env$natmort[r, y+1, , , sim],
            fish_sel = sim_env$fish_sel[r, y+1, , , , sim],
            f_init = 0.05
          )
        } # end r loop
      } # end if

      # save models
      sim_env$models[[y - sim_env$feedback_start_yr + 1]][[sim]] <- model
      if(y == sim_env$n_yrs) sim_env$models[[y - sim_env$feedback_start_yr + 1]][[sim]]$sd_rep <- RTMB::sdreport(model)

    } # end feedback

  } # end y loop

}

#' Run a Closed-Loop Feedback FAA Assessment
#'
#' This function simulates a fishery over multiple years using a closed-loop
#' framework. It performs annual population dynamics, feedback assessments,
#' calculates reference points, applies harvest control rules, projects
#' population and catch, and allocates total allowable catch (TAC) across regions
#' and fleets.
#'
#' @param sim_env A list or environment containing the simulation environment, including population dynamics,
#'   survey data, movement matrices, natural mortality, and other model parameters.
#' @param sim Integer. Index of the simulation replicate.
#' @param fleet_allocation Numeric vector. Proportions used to allocate regional TAC to fleets.
#' @param lls_design_type Character. Type of longline survey design used for apportionment.
#' @param srv_idx_se Numeric. Standard error of survey indices.
#' @param age_lag Integer. Age lag between recruitment and survey observation.
#' @param srv_wgt Character Weighting applied to survey data (e.g., "numbers" or "biomass").
#' @param fish_wgt Character Weighting applied to fishery data (e.g., "numbers" or "biomass").
#' @param faa_n_fish_fleets Integer. Number of fishing fleets in the feedback assessment.
#' @param faa_n_srv_fleets Integer. Number of survey fleets in the feedback assessment. Typically 4.
#' @param fish_sel_model Character vector. Fleet-specific selectivity models for the fishing fleets.
#' @param srv_sel_model Character vector. Fleet-specific selectivity models for survey fleets.
#' @param fish_selex_prior LPriors for fishing fleet selectivity parameters. Constructed based on
#'   fleet blocks, sex, and region.
#' @param srv_selex_prior  Priors for survey fleet selectivity parameters.
#' @param fish_sel_blocks Fishery blocks
#' @param srv_sel_blocks Survey blocks
#'
#' @returns None. The function updates `sim_env` in-place with population dynamics,
#'   fishing mortality, and assessment model outputs. The last year's model object
#'   is stored in `sim_env$models[[sim]]`.
#'
run_faa_closedloop_i <- function(sim_env,
                                 sim,
                                 fleet_allocation,
                                 lls_design_type,
                                 srv_idx_se,
                                 age_lag,
                                 srv_wgt,
                                 fish_wgt,
                                 faa_n_fish_fleets,
                                 faa_n_srv_fleets,
                                 fish_sel_model,
                                 srv_sel_model,
                                 fish_selex_prior,
                                 srv_selex_prior,
                                 fish_sel_blocks,
                                 srv_sel_blocks
                                 ) {

  # Run Closed Loop ---------------------------------------------------------
  for(y in 1:sim_env$n_yrs) {

    # Execute annual population dynamics
    run_annual_cycle(y, sim, sim_env)

    # Run feedback
    if(y >= sim_env$feedback_start_yr) {

      ### Assessment --------------------------------------------------------------
      # get faa data
      asmt_list <- faa_em(sim_env = sim_env,
                          y = y,
                          sim = sim,
                          srv_idx_se = srv_idx_se,
                          age_lag = age_lag,
                          lls_design_type = lls_design_type,
                          faa_n_fish_fleets = faa_n_fish_fleets,
                          faa_n_srv_fleets = faa_n_srv_fleets,
                          srv_wgt = srv_wgt,
                          fish_wgt = fish_wgt,
                          fish_sel_blocks = fish_sel_blocks,
                          fish_sel_model = fish_sel_model,
                          fish_fixed_sel_pars_spec = rep("est_all", faa_n_fish_fleets),
                          fish_selex_prior = fish_selex_prior,
                          srv_sel_blocks =  srv_sel_blocks,
                          srv_sel_model = srv_sel_model,
                          srv_fixed_sel_pars_spec = rep("est_all", faa_n_srv_fleets),
                          srv_selex_prior = srv_selex_prior,
                          cross_testing = FALSE
      )


      # some starting values
      asmt_list$par$ln_fish_fixed_sel_pars[,1,,,1] <- log(4)
      asmt_list$par$ln_fish_fixed_sel_pars[,2,,,1] <- log(3)
      # asmt_list$par$ln_srv_fixed_sel_pars[] <- log(5)

      model <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = F) # get model

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
          calc_rec_st_yr = 20,
          rec_age = 2,
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
      projection <- get_population_projection(data = asmt_list$data, rep = model$rep,
                                              y = y, n_proj_yrs = 2,
                                              f_ref_pt = array(exp_frate, c(asmt_list$data$n_regions, 2))) # get 1 year projection
      tac <- sum(projection$proj_Catch[,2,]) # get tac

      ### Apportionment -----------------------------------------------------------
      # get survey apportionment (survey biomass)
      tmp_Srv_BiomIA <- sim_env$SrvIAA[,((y - 4): y),,,1,sim] * sim_env$WAA_srv[,((y - 4): y),,,1,sim]
      tmp_srvidx_biom <- apply(tmp_Srv_BiomIA, 1:2, sum) # get biomass
      apportionment <- get_single_region_survey_apportionment(feedback_start_yr = sim_env$feedback_start_yr,
                                                              n_yrs = sim_env$n_yrs, y = y,
                                                              srv_idx = tmp_srvidx_biom, # using true survey biomass values to do apportionment
                                                              rolling_avg_yrs = 5,
                                                              lls_design_type = lls_design_type)
      tac_r <- tac * apportionment # get regional tac
      tac_rf <- tac_r * fleet_allocation # allocate regional tac by fleet

      ### TAC to F ----------------------------------------------------------------
      if (y < sim_env$n_yrs) {
        move_age <- if(sim_env$do_recruits_move == 0) 2 else 1
        tmp_naa_moved <- array(NA, dim = dim(sim_env$NAA[, y+1, , , sim]))
        for(a in move_age:sim_env$n_ages) for(s in 1:sim_env$n_sexes)
          tmp_naa_moved[,a,s] <- t(sim_env$NAA[,y+1,a,s,sim]) %*% sim_env$Movement[,,y+1,a,s,sim]
        if(move_age == 2) tmp_naa_moved[,1,] <- sim_env$NAA[, y+1,1, , sim]

        # Solve for all F simultaneously for each region
        for(r in 1:sim_env$n_regions) {
          sim_env$Fmort[r, y+1, , sim] <- solve_multifleet_F(
            target_catch = tac_rf[r, ],
            NAA = tmp_naa_moved[r, , ],
            WAA = sim_env$WAA[r, y+1, , , sim],
            natmort = sim_env$natmort[r, y+1, , , sim],
            fish_sel = sim_env$fish_sel[r, y+1, , , , sim],
            f_init = 0.05
          )
        } # end r loop
      } # end if

      # save models
      sim_env$models[[y - sim_env$feedback_start_yr + 1]][[sim]] <- model
      if(y == sim_env$n_yrs) sim_env$models[[y - sim_env$feedback_start_yr + 1]][[sim]]$sd_rep <- RTMB::sdreport(model)

    } # end feedback

  } # end y loop

}

#' Run a Closed-Loop Feedback Three Region Assessment
#'
#' This function simulates a fishery over multiple years using a closed-loop
#' framework. It performs annual population dynamics, feedback assessments,
#' calculates reference points, applies harvest control rules, projects
#' population and catch, and allocates total allowable catch (TAC) across regions
#' and fleets.
#'
#' @param sim_env A list or environment containing the simulation environment, including population dynamics,
#'   survey data, movement matrices, natural mortality, and other model parameters.
#' @param sim Integer. Index of the simulation replicate.
#' @param fleet_allocation Numeric vector. Proportions used to allocate regional TAC to fleets.
#' @param lls_design_type Character. Type of longline survey design used for apportionment - only for the GOA
#' @param srv_idx_se Numeric. Standard error of survey indices.
#' @param age_lag Integer. Age lag between recruitment and survey observation.
#' @param srv_wgt Character Weighting applied to survey data for the GOA (e.g., "numbers" or "biomass").
#' @param fish_wgt Character Weighting applied to fishery data for the GOA (e.g., "numbers" or "biomass").
#'
#' @returns None. The function updates `sim_env` in-place with population dynamics,
#'   fishing mortality, and assessment model outputs. The last year's model object
#'   is stored in `sim_env$models[[sim]]`.
#'
run_three_rg_closedloop_i <- function(sim_env,
                                      sim,
                                      fleet_allocation,
                                      lls_design_type,
                                      srv_idx_se,
                                      age_lag,
                                      srv_wgt,
                                      fish_wgt
                                      ) {

  # Run Closed Loop ---------------------------------------------------------
  for(y in 1:sim_env$n_yrs) {

    # Execute annual population dynamics
    run_annual_cycle(y, sim, sim_env)

    # Run feedback
    if(y >= sim_env$feedback_start_yr) {

      ### Assessment --------------------------------------------------------------
      # get three region data
      asmt_list <- three_rg_em(
        sim_env = sim_env,
        y = y,
        sim = sim,
        srv_idx_se = srv_idx_se,
        age_lag = age_lag,
        lls_design_type = lls_design_type,
        srv_wgt = srv_wgt,
        fish_wgt = fish_wgt,
        UseTagging = 1,
        cross_testing = FALSE
      )

      model <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = F) # get model

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
          calc_rec_st_yr = 20,
          rec_age = 2,
          type = 'multi_region',
          what = "global_SPR",
          B_x = 0.4
        ),
        n_proj_yrs = 2
      )

      ### HCR ---------------------------------------------------------------------
      agg_ssb <- sum(model$rep$SSB[,y]) # get ssb
      exp_frate <- HCR_threshold(x = agg_ssb, frp = unique(ref_pts$f_ref_pt[,2]), brp = sum(ref_pts$b_ref_pt[,2]), alpha = 0.05) # get prescribed f

      ### Projection --------------------------------------------------------------
      projection <- get_population_projection(data = asmt_list$data, rep = model$rep,
                                              y = y, n_proj_yrs = 2,
                                              f_ref_pt = array(exp_frate, c(asmt_list$data$n_regions, 2))) # get 1 year projection
      tac <- rowSums(projection$proj_Catch[,2,]) # get tac by region

      ### Apportionment -----------------------------------------------------------
      # get survey apportionment only for the GOA (survey biomass)
      tmp_Srv_BiomIA <- sim_env$SrvIAA[,((y - 4): y),,,1,sim] * sim_env$WAA_srv[,((y - 4): y),,,1,sim]
      tmp_srvidx_biom <- apply(tmp_Srv_BiomIA, 1:2, sum) # get biomass
      apportionment <- get_three_rg_survey_apportionment(
        feedback_start_yr = sim_env$feedback_start_yr,
        n_yrs = sim_env$n_yrs, y = y,
        srv_idx = tmp_srvidx_biom, # using true survey biomass values to do apportionment
        rolling_avg_yrs = 5,
        lls_design_type = lls_design_type
      )
      tac_r <- c(tac[1:2], tac[3] * apportionment) # get regional tac (use f40 regional abc for bs and ai, and then apportion goa)
      tac_rf <- tac_r * fleet_allocation # allocate regional tac by fleet

      ### TAC to F ----------------------------------------------------------------
      if (y < sim_env$n_yrs) {
        move_age <- if(sim_env$do_recruits_move == 0) 2 else 1
        tmp_naa_moved <- array(NA, dim = dim(sim_env$NAA[, y+1, , , sim]))
        for(a in move_age:sim_env$n_ages) for(s in 1:sim_env$n_sexes)
          tmp_naa_moved[,a,s] <- t(sim_env$NAA[,y+1,a,s,sim]) %*% sim_env$Movement[,,y+1,a,s,sim]
        if(move_age == 2) tmp_naa_moved[,1,] <- sim_env$NAA[, y+1,1, , sim]

        # Solve for all F simultaneously for each region
        for(r in 1:sim_env$n_regions) {
          sim_env$Fmort[r, y+1, , sim] <- solve_multifleet_F(
            target_catch = tac_rf[r, ],
            NAA = tmp_naa_moved[r, , ],
            WAA = sim_env$WAA[r, y+1, , , sim],
            natmort = sim_env$natmort[r, y+1, , , sim],
            fish_sel = sim_env$fish_sel[r, y+1, , , , sim],
            f_init = 0.05
          )
        } # end r loop
      } # end if

      # save models
      sim_env$models[[y - sim_env$feedback_start_yr + 1]][[sim]] <- model
      if(y == sim_env$n_yrs) sim_env$models[[y - sim_env$feedback_start_yr + 1]][[sim]]$sd_rep <- RTMB::sdreport(model)

      } # end feedback

  } # end y loop

}

#' Run a Closed-Loop Feedback Five Region Model using true values
#'
#' This function simulates a fishery over multiple years using a closed-loop
#' framework. It performs annual population dynamics, feedback assessments,
#' calculates reference points, applies harvest control rules, projects
#' population and catch, and allocates total allowable catch (TAC) across regions
#' and fleets.
#'
#' @param sim_env A list or environment containing the simulation environment, including population dynamics,
#'   survey data, movement matrices, natural mortality, and other model parameters.
#' @param sim Integer. Index of the simulation replicate.
#' @param fleet_allocation Numeric vector. Proportions used to allocate regional TAC to fleets.
#' @param hcr_type How HCR is implemented: global_ssb_global_b40, local_ssb_local_b40, local_ssb_global_b40
#'
#' @returns None. The function updates `sim_env` in-place with population dynamics,
#'   fishing mortality, and assessment model outputs. The last year's model object
#'   is stored in `sim_env$models[[sim]]`.
#'
run_five_rg_closedloop_i <- function(sim_env,
                                    sim,
                                    fleet_allocation,
                                    hcr_type
                                    ) {

  # Run Closed Loop ---------------------------------------------------------
  for(y in 1:sim_env$n_yrs) {

    # Execute annual population dynamics
    run_annual_cycle(y, sim, sim_env)

    # Run feedback
    if(y >= sim_env$feedback_start_yr) {

      ### Assessment --------------------------------------------------------------
      # Not using assessment

      ### Reference Points --------------------------------------------------------
      ref_pts <- get_closed_loop_reference_points(
        use_true_values = TRUE,
        sim_env = sim_env,
        asmt_data = NULL,
        asmt_rep = NULL,
        y = y,
        sim = sim,

        # single region reference points
        reference_points_opt = list(
          n_avg_yrs = 1,
          SPR_x = 0.4,
          calc_rec_st_yr = 20,
          rec_age = 2,
          type = 'multi_region',
          what = "global_SPR",
          B_x = 0.4
        ),
        n_proj_yrs = 2
      )

      ### HCR ---------------------------------------------------------------------
      f_ref_pt <- array(NA, c(sim_env$n_regions, 2))

      # compare global ssb to global reference points
      if(hcr_type == 'global_ssb_global_b40') {
        agg_ssb <- sum(sim_env$SSB[,y,sim]) # get ssb
        exp_frate <- HCR_threshold(x = agg_ssb, frp = unique(ref_pts$f_ref_pt[,2]), brp = sum(ref_pts$b_ref_pt[,2]), alpha = 0.05) # get prescribed f
        f_ref_pt[] <- exp_frate
      }

      # compare local ssb to local reference points
      if(hcr_type == 'local_ssb_local_b40') {
        ssb_r <- sim_env$SSB[,y,sim] # get ssb
        for(rr in 1:sim_env$n_regions) {
          f_ref_pt[rr,] <- HCR_threshold(x = ssb_r[rr], frp = unique(ref_pts$f_ref_pt[,2]),
                                        brp = ref_pts$b_ref_pt[rr,2], alpha = 0.05) # get prescribed f
        }
      }

      # compare local ssb to global b40
      if(hcr_type == 'local_ssb_global_b40') {
        ssb_r <- sim_env$SSB[,y,sim] # get ssb
        for(rr in 1:sim_env$n_regions) {
          f_ref_pt[rr,] <- HCR_threshold(x = ssb_r[rr], frp = unique(ref_pts$f_ref_pt[,2]),
                                         brp = sum(ref_pts$b_ref_pt[,2]), alpha = 0.05) # get prescribed f
        }
      }

      ### Projection --------------------------------------------------------------
      # create asmt_list and model to store true values in
      asmt_list <- list()
      model <- list()
      asmt_list$data <- list()
      model$rep <- list()

      # make inputs for projection
      model$rep$NAA <- sim_env$NAA[,,,,sim]
      model$rep$NAA0 <- sim_env$NAA0[,,,,sim]
      asmt_list$data$WAA <- sim_env$WAA[,,,,sim]
      asmt_list$data$WAA_fish <- sim_env$WAA_fish[,,,,,sim]
      asmt_list$data$MatAA <- sim_env$MatAA[,,,,sim]
      model$rep$fish_sel <- sim_env$fish_sel[,,,,,sim]
      model$rep$Fmort <- sim_env$Fmort[,,,sim]
      model$rep$natmort <- sim_env$natmort[,,,,sim]
      model$rep$Rec <- sim_env$Rec[,,sim]
      model$rep$sexratio <- sim_env$sexratio[,,,sim]
      model$rep$Movement <- sim_env$Movement[,,,,,sim]
      asmt_list$data$n_regions <- sim_env$n_regions
      asmt_list$data$ages <- 1:sim_env$n_ages
      asmt_list$data$n_sexes <- sim_env$n_sexes
      asmt_list$data$n_fish_fleets <- sim_env$n_fish_fleets

      projection <- get_population_projection(data = asmt_list$data, rep = model$rep, y = y, n_proj_yrs = 2,
                                              f_ref_pt = f_ref_pt) # get 1 year projection
      tac_r <- rowSums(projection$proj_Catch[,2,]) # get tac by region

      ### Apportionment -----------------------------------------------------------
      tac_rf <- tac_r * fleet_allocation # allocate regional tac by fleet

      ### TAC to F ----------------------------------------------------------------
      if (y < sim_env$n_yrs) {
        move_age <- if(sim_env$do_recruits_move == 0) 2 else 1
        tmp_naa_moved <- array(NA, dim = dim(sim_env$NAA[, y+1, , , sim]))
        for(a in move_age:sim_env$n_ages) for(s in 1:sim_env$n_sexes)
          tmp_naa_moved[,a,s] <- t(sim_env$NAA[,y+1,a,s,sim]) %*% sim_env$Movement[,,y+1,a,s,sim]
        if(move_age == 2) tmp_naa_moved[,1,] <- sim_env$NAA[, y+1,1, , sim]

        # Solve for all F simultaneously for each region
        for(r in 1:sim_env$n_regions) {
          sim_env$Fmort[r, y+1, , sim] <- solve_multifleet_F(
            target_catch = tac_rf[r, ],
            NAA = tmp_naa_moved[r, , ],
            WAA = sim_env$WAA[r, y+1, , , sim],
            natmort = sim_env$natmort[r, y+1, , , sim],
            fish_sel = sim_env$fish_sel[r, y+1, , , , sim],
            f_init = 0.05
          )
        } # end r loop
      } # end if

    } # end feedback

    # print(sum(sim_env$SSB[,y,sim]))
    # par(mfrow = c(1,2))
    # plot(colSums(sim_env$SSB[,,sim]))
    # abline(h = sum(ref_pts$b_ref_pt[,2]))
    # plot(rowSums(sim_env$Fmort[4,,,sim]))
    # abline(h = unique(ref_pts$f_ref_pt[,2]))
    # # median(colSums(sim_env$SSB[,65:y,sim]))
  } # end y loop

}
