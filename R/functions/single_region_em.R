#' Fit a Single-Region Estimation Model to Simulated Data
#'
#' Constructs and configures a single-region SPoRC estimation model using
#' simulated data from a multi-region operating model. The function aggregates
#' inputs, applies indicator settings, prepares biological and fishery
#' structures, and returns a fully assembled `input_list` ready for model
#' estimation.
#'
#' @param sim_env A simulation environment list containing operating-model
#'   outputs, dimensions, and aggregated observation containers.
#' @param y Integer. Terminal year to fit the estimation model to (1 to `n_yrs`).
#' @param sim Integer simulation index.
#' @param srv_idx_se Numeric. Standard error applied to survey indices.
#'   Default is `0.2`.
#' @param age_lag Integer. Lag determining which age compositions are used.
#'   Default is `1`.
#' @param lls_design_type Character. Sampling design type used when aggregating
#'   length-limited surveys.
#'
#' @return A fully configured SPoRC estimation-model `input_list` for a
#'   single-region model.
single_region_em <- function(sim_env, y, sim, srv_idx_se = 0.2, age_lag = 1, lls_design_type) {

  # Get simulated data
  sim_data <- simulation_data_to_SPoRC(sim_env, y, sim)
  agg_data_to_single_rg(sim_data, sim_env, y, sim, lls_design_type) # aggregate data; updates sim_env
  use_indcators <- single_region_use_indicators(y, sim_env$n_fish_fleets, sim_env$n_srv_fleets, age_lag = age_lag)

  ### Setup Model -------------------------------------------------------------
  # Model dimensions
  input_list <- Setup_Mod_Dim(years = 1:y, # vector of years
                              ages = 1:sim_env$n_ages, # vector of ages
                              lens = 1:sim_env$n_lens, # number of lengths
                              n_regions = 1, # number of regions
                              n_sexes = sim_env$n_sexes, # number of sexes
                              n_fish_fleets = sim_env$n_fish_fleets, # number of fishery fleet
                              n_srv_fleets = sim_env$n_srv_fleets, # number of survey fleets
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
                              ln_global_R0 = log(25)
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
                                      ObsCatch = array(sim_env$Agg_ObsCatch[,1:y,,sim],
                                                       dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                      UseCatch = use_indcators$usecatch,
                                      # Model options
                                      Use_F_pen = 1,
                                      # whether to use f penalty, == 0 don't use, == 1 use
                                      sigmaC_spec = 'fix',
                                      ln_sigmaC = array(log(0.05), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  )

  # get fishery comps
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            # data inputs
                                            ObsFishIdx = sim_data$ObsFishIdx[1,,,drop = FALSE],
                                            ObsFishIdx_SE = sim_data$ObsFishIdx[1,,,drop = FALSE],
                                            UseFishIdx =  array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            ObsFishAgeComps = array(sim_env$Agg_ObsFishAgeComps[,1:y,,,,sim],
                                                                    dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes,
                                                                            input_list$data$n_fish_fleets)),
                                            UseFishAgeComps = use_indcators$usefishage,
                                            ISS_FishAgeComps = array(sim_env$Agg_ISS_FishAgeComps[,1:y,,,sim],
                                                                     dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes,
                                                                             input_list$data$n_fish_fleets)),
                                            ObsFishLenComps = array(sim_env$Agg_ObsFishLenComps[,1:y,,,,sim],
                                                                    dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes,
                                                                            input_list$data$n_fish_fleets)),
                                            UseFishLenComps = use_indcators$usefishlen,
                                            ISS_FishLenComps = array(sim_env$Agg_ISS_FishLenComps[,1:y,,,sim],
                                                                     dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes,
                                                                             input_list$data$n_fish_fleets)),

                                            # Model options
                                            fish_idx_type = c("none", "none"),
                                            FishAgeComps_LikeType =
                                              c("Multinomial", "none"),
                                            FishLenComps_LikeType =
                                              c("Multinomial", "Multinomial"),
                                            FishAgeComps_Type =
                                              c("spltRjntS_Year_1-terminal_Fleet_1",
                                                "none_Year_1-terminal_Fleet_2"),
                                            FishLenComps_Type =
                                              c("spltRjntS_Year_1-terminal_Fleet_1",
                                                "spltRjntS_Year_1-terminal_Fleet_2")
  )


  # survey index
  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           # data inputs
                                           ObsSrvIdx = array(sim_env$Agg_ObsSrvIdx[,1:y,,sim], dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                           ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE[1,,,drop = FALSE],
                                           UseSrvIdx =  use_indcators$usesrvidx,
                                           ObsSrvAgeComps = array(sim_env$Agg_ObsSrvAgeComps[,1:y,,,,sim],
                                                                  dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes,
                                                                          input_list$data$n_srv_fleets)),
                                           UseSrvAgeComps = use_indcators$usesrvage,
                                           ISS_SrvAgeComps = array(sim_env$Agg_ISS_SrvAgeComps[,1:y,,,sim],
                                                                   dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes,
                                                                           input_list$data$n_srv_fleets)),
                                           ObsSrvLenComps = array(sim_env$Agg_ObsSrvLenComps[,1:y,,,,sim],
                                                                  dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes,
                                                                          input_list$data$n_srv_fleets)),
                                           UseSrvLenComps = use_indcators$usesrvlen,
                                           ISS_SrvLenComps = array(sim_env$Agg_ISS_SrvLenComps[,1:y,,,sim],
                                                                   dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes,
                                                                           input_list$data$n_srv_fleets)),

                                           # Model options
                                           srv_idx_type = c("abd", 'biom', "abd"),
                                           SrvAgeComps_LikeType =
                                             c("Multinomial", 'none', "Multinomial"),
                                           SrvLenComps_LikeType =
                                             c("none", "none", "none"),
                                           SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                                                                "none_Year_1-terminal_Fleet_2",
                                                                "spltRjntS_Year_1-terminal_Fleet_3"),
                                           SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                                                                "none_Year_1-terminal_Fleet_2",
                                                                "none_Year_1-terminal_Fleet_3")
  )


  # Fishery Selex
  # defining priors
  sex_par <- expand.grid(sex = 1:2, par = 1:2)
  fleet_blocks <- data.frame(
    fleet = c(1, 1, 1, 2),
    block = c(1, 2, 3, 1)
  )

  # Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
  fish_selex_prior <- cbind(
    region = 1,
    fish_selex_structure,
    mu = 1,                                                                      # All selex means = 1 (means should be defined in normal space)
    sd = 5                                                                       # All selex sd = 5
  )

  fish_selex_prior_tf <- fish_selex_prior %>%                                    # set tighter selex prior for TF
    dplyr::filter((fleet == 2 & par == 1)) %>%
    dplyr::mutate(mu = 2, sd = 1) %>%
    dplyr::full_join(fish_selex_prior %>%  dplyr::filter(!(fleet == 2 & par == 1 )))

  fish_selex_prior_tf <- fish_selex_prior_tf %>%                                    # set tighter selex prior for TF
    dplyr::filter((fleet == 2 & par == 2)) %>%
    dplyr::mutate(mu = 1, sd = 2) %>%
    dplyr::full_join(fish_selex_prior_tf %>%  dplyr::filter(!(fleet == 2 & par == 2)))

  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                        cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                        fish_sel_blocks =                        # fishery selectivity time blocks if not TV specified above for a given fleet
                                          c("Block_1_Year_1-35_Fleet_1",         # pre-IFQ time block for fixed gear fishery 1994 and before
                                            "Block_2_Year_36-56_Fleet_1",        # IFQ time block for fixed gear fishery-- 1995 to 2015
                                            "Block_3_Year_57-terminal_Fleet_1",  # Recent time block for fixed gear fishery--2016 to terminal year
                                            "none_Fleet_2"),
                                        fish_sel_model =
                                          c("logist1_Fleet_1", "gamma_Fleet_2"),
                                        fish_q_blocks =
                                          c("none_Fleet_1", "none_Fleet_2"),
                                        fish_q_spec =
                                          c("fix", "fix"),
                                        fish_fixed_sel_pars_spec = c("est_all", "est_all"),
                                        Use_fish_selex_prior = 1,
                                        fish_selex_prior = fish_selex_prior
  )

  # setup survey selectivity
  # Define sex and parameter combinations
  sex_par <- expand.grid(sex = 1:2, par = 1:2)

  # Define valid fleet-block combinations (only estimating domestic and jp LLS)
  fleet_blocks <- data.frame(
    fleet = c(1, 1, 3),
    block = c(1, 2, 1)
  )

  # Merge to get all valid combinations
  srv_selex_structure <- merge(fleet_blocks, sex_par)

  # Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
  srv_selex_prior <- cbind(
    region = 1,
    srv_selex_structure,
    mu = 1,
    sd = 5
  ) %>%
  mutate(mu = ifelse(fleet == 3, 2, mu),
         sd = ifelse(fleet == 3, 3, sd))

  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                       # Model options
                                       cont_tv_srv_sel =
                                         c("none_Fleet_1",
                                           "none_Fleet_2",
                                           "none_Fleet_3"),
                                       srv_sel_blocks =
                                         c("Block_1_Year_1-56_Fleet_1",
                                           "Block_2_Year_57-terminal_Fleet_1",
                                           "none_Fleet_2",
                                           "none_Fleet_3"
                                         ),
                                       srv_sel_model =
                                         c("logist1_Fleet_1",
                                           "exponential_Fleet_2",
                                           "logist1_Fleet_3"
                                         ),
                                       srv_q_blocks =
                                         c("none_Fleet_1",
                                           "none_Fleet_2",
                                           "none_Fleet_3"),
                                       srv_fixed_sel_pars_spec =
                                         c("est_all",
                                           "fix",
                                           "est_all"),
                                       srv_q_spec =
                                         c("est_all",
                                           "fix",
                                           "est_all"),
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
