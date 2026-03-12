#' Set up a single-region FAA estimation model
#'
#' Constructs and returns an SPoRC input list for a single-region,
#' fleet-aggregated age-structured (FAA) estimation model using simulated data.
#' The function aggregates simulated observations to the specified number of
#' fishery and survey fleets, configures data availability indicators, and
#' initializes all major model components including recruitment, biology,
#' fishing mortality, indices, composition data, and selectivity. Fishery and
#' survey selectivity structures are fully user-specified via function
#' arguments.
#'
#' @param sim_env Simulation environment containing operating model outputs and
#'   aggregated observation objects.
#' @param y Terminal assessment year.
#' @param sim Simulation replicate index.
#' @param srv_idx_se Observation error (CV) applied to all survey indices.
#' @param age_lag Integer specifying the lag (in years) between observation and
#'   recruitment year used when constructing data availability indicators.
#' @param lls_design_type Character string specifying the longline survey design
#'   used to determine which data streams are available.
#' @param faa_n_fish_fleets Number of aggregated fishery fleets in the FAA model.
#' @param faa_n_srv_fleets Number of aggregated survey fleets in the FAA model.
#' @param srv_wgt Weighting type used when aggregating survey observations
#'   (e.g., \code{"numbers"} or \code{"biomass"}).
#' @param fish_wgt Weighting type used when aggregating fishery observations
#'   (e.g., \code{"numbers"} or \code{"biomass"}).
#' @param fish_sel_blocks Character vector defining fishery selectivity blocks
#'   passed to \code{Setup_Mod_Fishsel_and_Q}.
#' @param fish_sel_model Character vector specifying the functional form of
#'   fishery selectivity by fleet.
#' @param fish_fixed_sel_pars_spec Character vector specifying parameter sharing
#'   for fishery selectivity (e.g., shared or fleet-specific).
#' @param fish_selex_prior Data frame defining prior distributions for fishery
#'   selectivity parameters.
#' @param srv_sel_blocks Character vector defining survey selectivity blocks
#'   passed to \code{Setup_Mod_Srvsel_and_Q}.
#' @param srv_sel_model Character vector specifying the functional form of survey
#'   selectivity by fleet.
#' @param srv_fixed_sel_pars_spec Character vector specifying parameter sharing
#'   for survey selectivity.
#' @param srv_selex_prior Data frame defining prior distributions for survey
#'   selectivity parameters.
#' @param cross_testing Boolean on whether cross testing to add aggregated objects to simulation environment
#'
faa_em <- function(sim_env,
                   y,
                   sim,
                   srv_idx_se = 0.2,
                   age_lag = 1,
                   lls_design_type,
                   faa_n_fish_fleets,
                   faa_n_srv_fleets,
                   srv_wgt = 'numbers',
                   fish_wgt = 'biomass',
                   fish_sel_blocks,
                   fish_sel_model,
                   fish_fixed_sel_pars_spec,
                   fish_selex_prior,
                   srv_sel_blocks,
                   srv_sel_model,
                   srv_fixed_sel_pars_spec,
                   srv_selex_prior,
                   cross_testing = TRUE
                   ) {

  # get simulated data
  if(cross_testing && y == sim_env$feedback_start_yr) add_aggregated_obj_to_simenv(sim_env, 'faa', faa_n_fish_fleets, faa_n_srv_fleets)
  sim_data <- simulation_data_to_SPoRC(sim_env, y, sim)
  sim_agg_data <- agg_data_to_faa(sim_data, sim_env, y, sim, faa_n_fish_fleets,
                                  faa_n_srv_fleets, srv_wgt, fish_wgt, srv_idx_se)
  use_indcators <- faa_use_indicators(sim_env, y, faa_n_fish_fleets, faa_n_srv_fleets, age_lag, lls_design_type)

  if(faa_n_fish_fleets == 2) {
    fixed_gear_rep <- 1
    trawl_gear_rep <- 1
    fixed_gear_st <- 1
    trawl_gear_st <- 2
  }

  if(faa_n_fish_fleets == 4) {
    fixed_gear_rep <- 3
    trawl_gear_rep <- 1
    fixed_gear_st <- 1
    trawl_gear_st <- 4
  }

  if(faa_n_fish_fleets == 5) {
    fixed_gear_rep <- 3
    trawl_gear_rep <- 2
    fixed_gear_st <- 1
    trawl_gear_st <- 4
  }

  if(faa_n_fish_fleets == 6) {
    fixed_gear_rep <- 3
    trawl_gear_rep <- 3
    fixed_gear_st <- 1
    trawl_gear_st <- 4
  }

  ### Setup Model -------------------------------------------------------------
  # Model dimensions
  input_list <- Setup_Mod_Dim(years = 1:y, # vector of years
                              ages = 1:sim_env$n_ages, # vector of ages
                              lens = 1:sim_env$n_lens, # number of lengths
                              n_regions = 1, # number of regions
                              n_sexes = sim_env$n_sexes, # number of sexes
                              n_fish_fleets = faa_n_fish_fleets, # number of fishery fleet
                              n_srv_fleets = faa_n_srv_fleets, # number of survey fleets
                              verbose = F
  )

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                              do_rec_bias_ramp = 0,
                              sigmaR_switch = 16, # switch sigmaR
                              dont_est_recdev_last = 1, # don't estimate last rec dev
                              # Model options
                              rec_model = "mean_rec", # recruitment model
                              sigmaR_spec = "fix", # fixing
                              ln_sigmaR = log(c(0.4, 0.9)),
                              # values to fix sigmaR at, or starting values
                              ln_global_R0 = log(18)
  )

  # Setup biological stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      WAA = sim_data$WAA[1,,,,drop = FALSE], # weight at age
                                      MatAA = sim_data$MatAA[1,,,,drop = FALSE], # maturity at age
                                      AgeingError = sim_data$AgeingError,
                                      # ageing error matrix
                                      fit_lengths = 1, # fitting lengths
                                      SizeAgeTrans = sim_data$SizeAgeTrans[1,,,,,drop = FALSE],
                                      # size age transition matrix
                                      M_spec = "fix", # fix natural mortality
                                      # values to fix natural mortality at
                                      Fixed_natmort = array(0.0988975, dim = c(input_list$data$n_regions,
                                                                               length(input_list$data$years),
                                                                               length(input_list$data$ages),
                                                                               input_list$data$n_sexes))
  )


  # Configure movement and tagging (no tagging used)
  input_list <- Setup_Mod_Tagging(input_list = input_list, UseTagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  # setup catches
  # setup single region model
  input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                      # Data inputs
                                      ObsCatch = array(sim_env$Agg_ObsCatch[,1:y,,sim], dim = c(1, length(1:y), input_list$data$n_fish_fleets)),
                                      UseCatch = use_indcators$usecatch,
                                      # Model options
                                      Use_F_pen = 1,
                                      # whether to use f penalty, == 0 don't use, == 1 use
                                      sigmaC_spec = 'fix',
                                      ln_sigmaC = array(log(0.05), dim = c(1, length(1:y), input_list$data$n_fish_fleets))
  )

  # get fishery comps
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            # data inputs
                                            ObsFishIdx = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            ObsFishIdx_SE = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            UseFishIdx =  array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            ObsFishAgeComps = array(sim_env$Agg_ObsFishAgeComps[,1:y,,,,sim], dim = c(1, length(1:y), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                            UseFishAgeComps = use_indcators$usefishage,
                                            ISS_FishAgeComps = array(sim_env$Agg_ISS_FishAgeComps[,1:y,,,sim], dim = c(1, length(1:y), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                            ObsFishLenComps = array(sim_env$Agg_ObsFishLenComps[,1:y,,,,sim], dim = c(1, length(1:y), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                            UseFishLenComps = use_indcators$usefishlen,
                                            ISS_FishLenComps = array(sim_env$Agg_ISS_FishLenComps[,1:y,,,sim], dim = c(1, length(1:y), input_list$data$n_sexes, input_list$data$n_fish_fleets)),

                                            # Model options
                                            fish_idx_type = rep("none", input_list$data$n_fish_fleets),
                                            FishAgeComps_LikeType =
                                              c(rep("Multinomial", fixed_gear_rep), rep("none", trawl_gear_rep)),
                                            FishLenComps_LikeType =
                                              rep("Multinomial", input_list$data$n_fish_fleets),
                                            FishAgeComps_Type =
                                              c(paste("spltRjntS_Year_1-terminal_Fleet_", fixed_gear_st:fixed_gear_rep, sep = ''),
                                                paste("none_Year_1-terminal_Fleet_", trawl_gear_st:faa_n_fish_fleets, sep = '')),
                                            FishLenComps_Type =
                                              paste("spltRjntS_Year_1-terminal_Fleet_", 1:faa_n_fish_fleets, sep = '')
  )

  # survey index
  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           # data inputs
                                           ObsSrvIdx = array(sim_env$Agg_ObsSrvIdx[,1:y,,sim], dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                           ObsSrvIdx_SE = array(srv_idx_se, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                           UseSrvIdx =  use_indcators$usesrvidx,
                                           ObsSrvAgeComps = array(sim_env$Agg_ObsSrvAgeComps[,1:y,,,,sim], dim = c(1, length(1:y), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                           UseSrvAgeComps = use_indcators$usesrvage,
                                           ISS_SrvAgeComps = array(sim_env$Agg_ISS_SrvAgeComps[,1:y,,,sim], dim = c(1, length(1:y), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                           ObsSrvLenComps = array(sim_env$Agg_ObsSrvLenComps[,1:y,,,,sim], dim = c(1, length(1:y), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                           UseSrvLenComps = use_indcators$usesrvlen,
                                           ISS_SrvLenComps = array(sim_env$Agg_ISS_SrvLenComps[,1:y,,,sim], dim = c(1, length(1:y), input_list$data$n_sexes, input_list$data$n_srv_fleets)),

                                           # Model options
                                           srv_idx_type = rep("abd", input_list$data$n_srv_fleets),
                                           SrvAgeComps_LikeType =
                                             c(rep("Multinomial", 3), rep("Multinomial", length(4:faa_n_srv_fleets))),
                                           SrvLenComps_LikeType =
                                             rep("none", input_list$data$n_srv_fleets),
                                           SrvAgeComps_Type =
                                             c(paste("spltRjntS_Year_1-terminal_Fleet_", 1:3, sep = ''),
                                               paste("spltRjntS_Year_1-terminal_Fleet_", 4:faa_n_srv_fleets, sep = '')),
                                           SrvLenComps_Type =
                                             paste("none_Year_1-terminal_Fleet_", 1:faa_n_srv_fleets, sep = '')
  )

  # Fishery selectivity model
  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                        cont_tv_fish_sel = paste('none_Fleet_', 1:faa_n_fish_fleets, sep = ""),
                                        fish_sel_blocks = fish_sel_blocks,
                                        fish_sel_model = fish_sel_model,
                                        fish_q_blocks = paste('none_Fleet_', 1:faa_n_fish_fleets, sep = ""),
                                        fish_q_spec = rep('fix', faa_n_fish_fleets),
                                        fish_fixed_sel_pars_spec = fish_fixed_sel_pars_spec,
                                        Use_fish_selex_prior = 1,
                                        fish_selex_prior = fish_selex_prior
  )

  # Survey selectivity model
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                       cont_tv_srv_sel = paste('none_Fleet_', 1:faa_n_srv_fleets, sep = ""),
                                       srv_sel_blocks = srv_sel_blocks,
                                       srv_sel_model = srv_sel_model,
                                       srv_q_blocks = paste('none_Fleet_', 1:faa_n_srv_fleets, sep = ""),
                                       srv_q_spec =  rep('est_all', faa_n_srv_fleets),
                                       srv_fixed_sel_pars_spec = srv_fixed_sel_pars_spec,
                                       Use_srv_selex_prior = 1,
                                       srv_selex_prior = srv_selex_prior

  )

  # set up model weighting stuff
  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_F = 1,
                                    Wt_Tagging = 1,
                                    Wt_FishAgeComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_fish_fleets)),
                                    Wt_FishLenComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_fish_fleets)),
                                    Wt_SrvAgeComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_srv_fleets)),
                                    Wt_SrvLenComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_srv_fleets))
  )

  return(input_list)

}
