#' Construct a three-region estimation model from simulated data
#'
#' Builds a complete three-region (Bering Sea, Aleutian Islands, Gulf of Alaska)
#' spatial assessment model using simulated data. The function aggregates
#' fine-scale simulation output to three assessment regions, constructs
#' data-availability indicators, and sequentially configures all major model
#' components including recruitment, biology, movement, tagging, fisheries,
#' surveys, selectivity, and likelihood weighting.
#'
#' This function serves as a high-level wrapper that prepares the full
#' \code{input_list} object required for estimation, using a consistent set of
#' assumptions for aggregation, survey design, and parameter sharing across
#' regions.
#'
#' @param sim_env A list containing the simulation environment, including true
#'   population dynamics, observation processes, and storage objects for
#'   aggregated data. This object is modified in place to store three-region
#'   aggregated observations.
#'
#' @param y Integer. Terminal assessment year.
#'
#' @param sim Integer. Simulation replicate index.
#'
#' @param srv_idx_se Numeric. Log-scale standard error applied to aggregated survey
#'   indices. Default is 0.2.
#'
#' @param age_lag Integer. Number of years between observation year and the most
#'   recent age used in composition data. Default is 1.
#'
#' @param lls_design_type Character string specifying the longline survey design.
#'   Must be one of \code{"all"}, \code{"historical"}, or \code{"current"}.
#'
#' @param srv_wgt Character string specifying how surveys are aggregated across
#'   regions. Either \code{"numbers"} or \code{"biomass"}.
#'
#' @param fish_wgt Character string specifying how fisheries are aggregated across
#'   regions. Either \code{"numbers"} or \code{"biomass"}.
#'
#' @param UseTagging Integer (0 or 1). Indicator for whether tagging data are
#'   included in the model.
#' @param cross_testing Boolean on whether this is cross testing
#'
#' @return A fully populated \code{input_list} object suitable for fitting a
#'   three-region spatial assessment model.
#'
three_rg_em <- function(sim_env,
                        y,
                        sim,
                        srv_idx_se = 0.2,
                        age_lag = 1,
                        lls_design_type,
                        srv_wgt = 'numbers',
                        fish_wgt = 'biomass',
                        UseTagging = 1,
                        cross_testing = TRUE
) {

  # Get simulated data
  if(cross_testing && y == sim_env$feedback_start_yr) add_aggregated_obj_to_simenv(sim_env, 'three_rg', sim_env$n_fish_fleets, sim_env$n_srv_fleets)
  sim_data <- simulation_data_to_SPoRC(sim_env, y, sim)

  # get years in simulation
  x <- (sim_env$feedback_start_yr + 1):sim_env$n_yrs
  odd_yrs <- x[x %% 2 == 1]
  even_yrs <- x[x %% 2 == 0]

  # Deal with LLS Design and Tagging
  if(lls_design_type == 'historical') {
    # GOA released every year, AI released in odd years, BS released in even years
    bs_remove_cohort <- which(sim_data$tag_release_indicator[,1] == 1 & sim_data$tag_release_indicator[,2] %in% odd_yrs)
    ai_remove_cohort <- which(sim_data$tag_release_indicator[,1] == 2 & sim_data$tag_release_indicator[,2] %in% even_yrs)
    remove_cohorts <- c(bs_remove_cohort, ai_remove_cohort) # cohorts to remove
    if(length(remove_cohorts) > 1) {
      sim_data$tag_release_indicator <- sim_data$tag_release_indicator[-remove_cohorts, ]
      sim_data$Tagged_Fish <- sim_data$Tagged_Fish[-remove_cohorts,,]
      sim_data$Obs_Tag_Recap <- sim_data$Obs_Tag_Recap[,-remove_cohorts,,,]
    }
  }

  if(lls_design_type == 'current') {
    # GOA released every year, AI released in odd years, BS released in even years
    bsai_remove_cohort <- which(sim_data$tag_release_indicator[,1] %in% c(1:2) & sim_data$tag_release_indicator[,2] %in% even_yrs)
    goa_remove_cohort <- which(sim_data$tag_release_indicator[,1] %in% c(3:5) & sim_data$tag_release_indicator[,2] %in% odd_yrs)
    remove_cohorts <- c(bsai_remove_cohort, goa_remove_cohort) # cohorts to remove
    if(length(remove_cohorts) > 1) {
      sim_data$tag_release_indicator <- sim_data$tag_release_indicator[-remove_cohorts, ]
      sim_data$Tagged_Fish <- sim_data$Tagged_Fish[-remove_cohorts,,]
      sim_data$Obs_Tag_Recap <- sim_data$Obs_Tag_Recap[,-remove_cohorts,,,]
    }
  }

  # Aggregate up data
  agg_data_to_three_rg(sim_data, sim_env, y, sim, srv_wgt, fish_wgt, srv_idx_se)
  use_indcators <- threerg_use_indicators(sim_env, y,  age_lag, lls_design_type)

  ### Setup Model -------------------------------------------------------------
  # Model dimensions
  input_list <- Setup_Mod_Dim(years = 1:y, # vector of years
                              ages = 1:sim_env$n_ages, # vector of ages
                              lens = 1:sim_env$n_lens, # number of lengths
                              n_regions = 3, # number of regions
                              n_sexes = sim_env$n_sexes, # number of sexes
                              n_fish_fleets = sim_env$n_fish_fleets, # number of fishery fleet
                              n_srv_fleets = sim_env$n_srv_fleets, # number of survey fleets
                              verbose = F
  )

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                              do_rec_bias_ramp = 0,
                              sigmaR_switch = as.integer(length(1960:1975)),
                              dont_est_recdev_last = 1, # don't estimate last rec dev

                              # Model options
                              rec_model = "mean_rec", # recruitment model
                              sigmaR_spec = "fix", # fixing
                              InitDevs_spec = "est_shared_r",
                              # initial deviations are shared across regions,
                              # but recruitment deviations are region specific
                              ln_sigmaR = log(c(0.4, 0.9)),
                              # values to fix sigmaR at, or starting values
                              ln_global_R0 = log(20),
                              Use_Rec_prop_Prior = 1,
                              Rec_prop_prior = 1.5
  )

  # Setup biological stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      WAA = sim_data$WAA[1:3,,,], # weight at age
                                      MatAA = sim_data$MatAA[1:3,,,], # maturity at age
                                      AgeingError = sim_data$AgeingError,
                                      # ageing error matrix
                                      fit_lengths = 1, # fitting lengths
                                      SizeAgeTrans = sim_data$SizeAgeTrans[1:3,,,,],
                                      # size age transition matrix
                                      M_spec = "fix", # fix natural mortality
                                      # values to fix natural mortality at
                                      Fixed_natmort = array(0.0988975, dim = c(input_list$data$n_regions,
                                                                               length(input_list$data$years),
                                                                               length(input_list$data$ages),
                                                                               input_list$data$n_sexes))
  )

  # setting up movement parameterization
  Movement_prior <- expand.grid(
    region_from = 1:3, # regions
    age = c(6,7,16), # age blocks
    sex = 1, # sex
    alpha = I(list(rep(5, 3))) # prior alpha to each row
  )

  input_list <- Setup_Mod_Movement(input_list = input_list,
                                   # Model options
                                   Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)),
                                   # estimating movement in 3 age blocks
                                   # (ages 1-6, ages 7-15, ages 16-30)
                                   Movement_yearblk_spec = "constant", # time-invariant movement
                                   Movement_sexblk_spec = "constant", # sex-invariant movement
                                   do_recruits_move = 0, # recruits do not move
                                   use_fixed_movement = 0, # estimating movement
                                   Use_Movement_Prior = 1, # priors used for movement
                                   Movement_prior = Movement_prior # vague prior to penalize movement away from the extremes
  )

  # setup tagging priors
  tag_prior <- data.frame(
    region = 1,
    block = c(1,2,3),
    mu = NA, # no mean, since symmetric beta
    sd = 5, # sd = 5
    type = 0 # symmetric beta
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list,
                                  UseTagging = UseTagging, # using tagging data
                                  max_tag_liberty = 15, # maximum number of years to track a cohort

                                  # Data Inputs
                                  tag_release_indicator = sim_env$Agg_tag_release_indicator,
                                  # tag release indicator (first col = tag region,
                                  # second col = tag year),
                                  # total number of rows = number of tagged cohorts
                                  Tagged_Fish = sim_env$Agg_Tagged_Fish, # Released fish
                                  # dimensioned by total number of tagged cohorts, (implicitly
                                  # tracks the release year and region), age, and sex
                                  Obs_Tag_Recap = sim_env$Agg_Obs_Tag_Recap,
                                  # dimensioned by max tag liberty, tagged cohorts, regions,
                                  # ages, and sexes

                                  # Model options
                                  Tag_LikeType = "Multinomial_Release", # Multinomial Release Conditioned
                                  mixing_period = 2, # Don't fit tagging until release year + 1
                                  t_tagging = 0.5, # tagging happens midway through the year,
                                  # movement does not occur within that year
                                  tag_selex = "SexSp_DomFleet", # tagging recapture selectivity is the dominant fleet (fixed-gear)
                                  tag_natmort = "AgeSp_SexSp", # tagging natural mortality is age and sex-specific
                                  Use_TagRep_Prior = 1, # tag reporting rate priors are used
                                  TagRep_Prior = tag_prior,
                                  move_age_tag_pool = as.list(1:30), # whether or not to pool tagging data when fitting (for computational cost)
                                  move_sex_tag_pool = as.list(1:2), # whether or not to pool sex-specific data whezn fitting
                                  Init_Tag_Mort_spec = "fix", # fixing initial tag mortality
                                  Tag_Shed_spec = "fix", # fixing chronic shedding
                                  TagRep_spec = "est_shared_r", # tag reporting rates are not region specific
                                  Tag_Reporting_blocks = c(
                                    paste("Block_1_Year_1-36_Region_", c(1:input_list$data$n_regions), sep = ''),
                                    paste("Block_2_Year_36-56_Region_", c(1:input_list$data$n_regions), sep = ''),
                                    paste("Block_3_Year_57-terminal_Region_", c(1:input_list$data$n_regions), sep = '')
                                  ),
                                  # Specify starting values or fixing values
                                  ln_Init_Tag_Mort = log(0.1), # fixing initial tag mortality
                                  ln_Tag_Shed = log(0.02),  # fixing tag shedding
                                  ln_tag_theta = log(0.5), # starting value for tagging overdispersion
                                  Tag_Reporting_Pars = array(log(0.2 / (1-0.2)), dim = c(input_list$data$n_regions, 3))  # starting values for tag reporting pars
  )

  # setup catches
  # setup single region model
  input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                      # Data inputs
                                      ObsCatch = array(sim_env$Agg_ObsCatch[,1:y,,sim],
                                                       dim = c(3, length(1:y), input_list$data$n_fish_fleets)),
                                      UseCatch = use_indcators$usecatch,
                                      # Model options
                                      Use_F_pen = 1,
                                      # whether to use f penalty, == 0 don't use, == 1 use
                                      sigmaC_spec = 'fix',
                                      ln_sigmaC = array(log(0.05), dim = c(3, length(1:y), input_list$data$n_fish_fleets))
  )

  # get fishery comps
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            # data inputs
                                            ObsFishIdx = sim_data$ObsFishIdx[1:3,,,drop = FALSE],
                                            ObsFishIdx_SE = sim_data$ObsFishIdx[1:3,,,drop = FALSE],
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
                                           ObsSrvIdx_SE = array(srv_idx_se, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
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


  # Fishery Selectivity and Catchability

  # defining priors
  sex_par <- expand.grid(sex = 1:2, par = 1:2)
  fleet_blocks <- data.frame(
    fleet = c(1, 2),
    block = 1
  )

  # merge together (note that unlike the operational assessment, selectivity
  # blocks are reduced from 3 to 2)
  fish_selex_structure <- merge(fleet_blocks, sex_par)

  # Merge to get all valid combinations
  fish_selex_structure <- merge(fleet_blocks, sex_par) %>%
    dplyr::filter(!(fleet == 1 & block == 1 & sex == 2 & par == 2)) %>%              # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER input_list$map
    dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 1))                  # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER input_list$map

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

                                        # Model options
                                        cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                        # fishery selectivity, whether continuous time-varying

                                        # fishery selectivity blocks
                                        fish_sel_blocks =
                                          c("none_Fleet_1",
                                            "none_Fleet_2"),
                                        # no blocks for trawl fishery

                                        # fishery selectivity form
                                        fish_sel_model =
                                          c("logist1_Fleet_1", "gamma_Fleet_2"),

                                        # fishery catchability blocks
                                        fish_q_blocks =
                                          c("none_Fleet_1", "none_Fleet_2"),
                                        # no blocks since q is not estimated

                                        # sharing fishery selex input_list$par
                                        fish_fixed_sel_pars =
                                          c("est_shared_r", "est_shared_r"),

                                        # whether to estimate all fixed effects
                                        # for fishery catchability
                                        fish_q_spec =
                                          c("fix", "fix"),
                                        Use_fish_selex_prior = 1,
                                        fish_selex_prior = fish_selex_prior
  )

  # setup survey selectivity
  # Define sex and parameter combinations
  sex_par <- expand.grid(sex = 1:2, par = 1:2)

  # Define valid fleet-block combinations (only estimating domestic and jp LLS)
  fleet_blocks <- data.frame(
    fleet = c(1, 3),
    block = c(1, 1)
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
    filter(!(fleet == 3 & par == 2 & sex == 2)) %>%
    mutate(mu = ifelse(fleet == 3, 2, mu),
           sd = ifelse(fleet == 3, 3, sd))

  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                       # Model options
                                       # survey selectivity, whether continuous time-varying
                                       cont_tv_srv_sel =
                                         c("none_Fleet_1",
                                           "none_Fleet_2",
                                           "none_Fleet_3"),

                                       # survey selectivity blocks
                                       srv_sel_blocks =                          # survey selectivity time blocks if not TV specified above for a given fleet
                                         c("none_Fleet_1",
                                           "none_Fleet_2",                       # No blocks for trawl survey
                                           "none_Fleet_3"                        # No blocks for JPN LLS
                                         ),

                                       # survey selectivity form
                                       srv_sel_model =
                                         c("logist1_Fleet_1",
                                           "exponential_Fleet_2",
                                           "logist1_Fleet_3"
                                         ),

                                       # survey catchability blocks
                                       srv_q_blocks =
                                         c("none_Fleet_1",
                                           "none_Fleet_2",
                                           "none_Fleet_3"),

                                       # whether to estiamte all fixed effects
                                       # for survey selectivity and later
                                       # modify to fix/share input_list$par
                                       srv_fixed_sel_pars_spec =
                                         c("est_shared_r",
                                           "fix",
                                           "est_shared_r"),

                                       # whether to estiamte all
                                       # fixed effects for survey catchability
                                       # spatially-invariant q
                                       srv_q_spec =
                                         c("est_shared_r",
                                           "fix",
                                           "est_shared_r"),
                                       Use_srv_selex_prior = 1,
                                       srv_selex_prior = srv_selex_prior
  )

  # Map off early delta for fishery
  map_fish_fixed <- array(input_list$map$ln_fish_fixed_sel_pars, dim = dim(input_list$par$ln_fish_fixed_sel_pars))
  map_fish_fixed[,2,1,2,1]  <- map_fish_fixed[,2,1,1,1] # share deltas

  # Map off bmax for trawl females
  map_fish_fixed[,1,1,2,2]  <- map_fish_fixed[,1,1,1,2] # share deltas
  input_list$map$ln_fish_fixed_sel_pars <- factor(map_fish_fixed)

  # Map off delta for JP LLS
  map_srv_fixed <- array(input_list$map$ln_srv_fixed_sel_pars, dim = dim(input_list$par$ln_srv_fixed_sel_pars))
  map_srv_fixed[,2,1,2,3]  <- map_srv_fixed[,2,1,1,3] # share deltas
  input_list$map$ln_srv_fixed_sel_pars <- factor(map_srv_fixed)

  # set up model weighting stuff
  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_F = 1,
                                    Wt_Tagging = 0.1,
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
