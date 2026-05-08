# Purpose: To run a five region OM, with a FAA EM (base recruitment)
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 11/26/25

# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)
library(furrr)
library(progressr)

# Read in model output
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "faa_em.R"))

# Setup EMs ---------------------------------------------------------------
ems_grid <- expand.grid(
  fishery_structure = c("1F_1T", "3F_1T", "3F_2T", "3F_3T"),
  fishery_selex     = c("All_Logist", "BS_Gamma", "BS_AI_Gamma", "All_Gamma"),
  survey_selex      = c("All_Logist", "BS_Gamma", "BS_AI_Gamma", "All_Gamma"),
  fish_blocks_type  = c("none", "blocks"),
  srv_blocks_type   = c("none", "blocks"),
  stringsAsFactors  = FALSE
)

# Remove invalid combinations:
# 1F_1T can only use All_Logist or All_Gamma for fishery selex
ems_grid <- ems_grid[!(ems_grid$fishery_structure == "1F_1T" &
                         ems_grid$fishery_selex %in% c("BS_Gamma", "BS_AI_Gamma")), ]
ems_grid$model_id <- seq_len(nrow(ems_grid))


# Helper: base fleet-structure specs
get_fish_structure <- function(struct_name) {
  specs <- list(
    "1F_1T" = list(n_fish_fleets = 2, n_srv_fleets = 4),
    "3F_1T" = list(n_fish_fleets = 4, n_srv_fleets = 4),
    "3F_2T" = list(n_fish_fleets = 5, n_srv_fleets = 4),
    "3F_3T" = list(n_fish_fleets = 6, n_srv_fleets = 4)
  )
  specs[[struct_name]]
}

# Helper: fishery selectivity models
get_fish_selex <- function(selex_name, n_fleets) {
  specs <- list(
    "All_Logist" = list(
      "2" = c("logist1_Fleet_1", "gamma_Fleet_2"),
      "4" = c(paste0("logist1_Fleet_", 1:3), "gamma_Fleet_4"),
      "5" = c(paste0("logist1_Fleet_", 1:3), paste0("gamma_Fleet_", 4:5)),
      "6" = c(paste0("logist1_Fleet_", 1:3), paste0("gamma_Fleet_", 4:6))
    ),
    "All_Gamma" = list(
      "2" = c("gamma_Fleet_1", "gamma_Fleet_2"),
      "4" = c(paste0("gamma_Fleet_", 1:3), "gamma_Fleet_4"),
      "5" = c(paste0("gamma_Fleet_", 1:3), paste0("gamma_Fleet_", 4:5)),
      "6" = c(paste0("gamma_Fleet_", 1:3), paste0("gamma_Fleet_", 4:6))
    ),
    "BS_Gamma" = list(
      "4" = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:3), "gamma_Fleet_4"),
      "5" = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:3), paste0("gamma_Fleet_", 4:5)),
      "6" = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:3), paste0("gamma_Fleet_", 4:6))
    ),
    "BS_AI_Gamma" = list(
      "4" = c(paste0("gamma_Fleet_",   1:2), "logist1_Fleet_3", "gamma_Fleet_4"),
      "5" = c(paste0("gamma_Fleet_",   1:2), "logist1_Fleet_3", paste0("gamma_Fleet_", 4:5)),
      "6" = c(paste0("gamma_Fleet_",   1:2), "logist1_Fleet_3", paste0("gamma_Fleet_", 4:6))
    )
  )
  specs[[selex_name]][[as.character(n_fleets)]]
}

# Helper: survey selectivity models
get_survey_selex <- function(selex_name) {
  specs <- list(
    "All_Logist"  = paste0("logist1_Fleet_", 1:4),
    "BS_Gamma"    = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:4)),
    "BS_AI_Gamma" = c(paste0("gamma_Fleet_", 1:2), paste0("logist1_Fleet_", 3:4)),
    "All_Gamma"   = c(paste0("gamma_Fleet_", 1:3), "logist1_Fleet_4")
  )
  specs[[selex_name]]
}

# Helper: fishery selectivity blocks
# Fixed-gear fleets (1:n_fixed) get two time blocks; trawl fleets get "none".
# Returns fish_sel_blocks vector and matching prior_blocks list.
get_fish_blocks <- function(struct_name, blocks_type) {
  n_total <- get_fish_structure(struct_name)$n_fish_fleets
  n_fixed <- switch(struct_name, "1F_1T" = 1, 3)
  n_trawl <- n_total - n_fixed

  if (blocks_type == "none") {
    blk_vec   <- paste0("none_Fleet_", seq_len(n_total))
    prior_blk <- setNames(as.list(rep(1, n_total)), as.character(seq_len(n_total)))
  } else {
    blk_vec <- c(
      as.vector(sapply(seq_len(n_fixed), function(f) c(
        paste0("Block_1_Year_1-55_Fleet_",       f),
        paste0("Block_2_Year_56-terminal_Fleet_", f)
      ))),
      paste0("none_Fleet_", (n_fixed + 1):n_total)
    )
    # fixed-gear fleets -> 2 blocks each, trawl fleets -> 1 block each
    prior_blk <- setNames(
      as.list(c(rep(2, n_fixed), rep(1, n_trawl))),
      as.character(seq_len(n_total))
    )
  }
  list(fish_sel_blocks = blk_vec, fish_prior_blocks = prior_blk)
}

# Helper: survey selectivity blocks ─────────────────────────────────────────
# Surveys fleets 1-3 get two time blocks; fleet 4 gets none.
get_srv_blocks <- function(blocks_type) {
  if (blocks_type == "none") {
    blk_vec   <- paste0("none_Fleet_", 1:4)
    prior_blk <- setNames(as.list(rep(1, 4)), as.character(1:4))
  } else {
    blk_vec <- c(
      "Block_1_Year_1-55_Fleet_1",       "Block_2_Year_56-terminal_Fleet_1",
      "Block_1_Year_1-55_Fleet_2",       "Block_2_Year_56-terminal_Fleet_2",
      "Block_1_Year_1-55_Fleet_3",       "Block_2_Year_56-terminal_Fleet_3",
      "none_Fleet_4"
    )
    prior_blk <- list("1" = 2, "2" = 2, "3" = 2, "4" = 1)
  }
  list(srv_sel_blocks = blk_vec, srv_prior_blocks = prior_blk)
}

# Build EM list
ems <- vector("list", nrow(ems_grid))

for (i in seq_len(nrow(ems_grid))) {

  fish_struct  <- get_fish_structure(ems_grid$fishery_structure[i])
  fish_selex   <- get_fish_selex(ems_grid$fishery_selex[i], fish_struct$n_fish_fleets)
  surv_selex   <- get_survey_selex(ems_grid$survey_selex[i])
  fish_blk     <- get_fish_blocks(ems_grid$fishery_structure[i], ems_grid$fish_blocks_type[i])
  srv_blk      <- get_srv_blocks(ems_grid$srv_blocks_type[i])

  ems[[i]] <- list(
    faa_n_fish_fleets = fish_struct$n_fish_fleets,
    faa_n_srv_fleets  = fish_struct$n_srv_fleets,

    # Fishery
    fish_sel_blocks          = fish_blk$fish_sel_blocks,
    fish_sel_model           = fish_selex,
    fish_fixed_sel_pars_spec = rep("est_all", fish_struct$n_fish_fleets),
    fish_selex_prior         = build_selex_prior(
      fleets_blocks = fish_blk$fish_prior_blocks,
      sex_par       = expand.grid(sex = 1:2, par = 1:2),
      region        = 1,
      mu            = 3.5,
      sd            = 2
    ),

    # Survey
    srv_sel_blocks           = srv_blk$srv_sel_blocks,
    srv_sel_model            = surv_selex,
    srv_fixed_sel_pars_spec  = rep("est_all", 4),
    srv_selex_prior          = build_selex_prior(
      fleets_blocks = srv_blk$srv_prior_blocks,
      sex_par       = expand.grid(sex = 1:2, par = 1:2),
      region        = 1,
      mu            = 3.5,
      sd            = 2
    )
  )
}

# Condition OM ------------------------------------------------------------

# Extract model components from sablefish model
data <- om_values$data
rep <- om_values$rep
parameters <- om_values$parameters
mapping <- om_values$mapping
sd_rep <- om_values$sd_rep

# Define operating model parameters
closed_loop_yrs <- 30      # Years to project forward
n_years <- length(data$years)  # number of years
burnin_years <- 1:n_years  # Historical conditioning period
n_sims <- 50              # Number of replicate simulations
n_regions <- 5             # number of regions
n_fish_fleets <- 2         # number of fishery fleets
n_srv_fleets <- 3          # number of survey fleets

# get fleet allocation by region
fleet_allocation <- array(NA, dim = c(n_regions, n_fish_fleets))
fleet_allocation[,1] <- c(0.5,0.75,0.8,0.8,0.95) # from fmp
fleet_allocation[,2] <- 1 - fleet_allocation[,1] # trawl gear allocation

for(i in 1:3) {
  # Condition closed-loop simulations, random recruitment
  sim_list <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = n_sims,
    data = data,
    parameters = parameters,
    mapping = mapping,
    sd_rep = sd_rep,
    rep = rep,
    random = NULL,
    recruitment_opt = 'resample_from_input',
    # setup variances
    ISS_FishAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
    ISS_FishLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
    ISS_SrvAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
    ISS_SrvLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
    ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
    ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets))
  )

  # Get constant B40 for comparison
  global_spr <- Get_Reference_Points(data = om_values$data,
                                     rep = om_values$rep,
                                     SPR_x = 0.4,
                                     type = 'multi_region',
                                     what = 'global_SPR',
                                     calc_rec_st_yr = 20,
                                     rec_age = 2
  )

  # Run MSEs ----------------------------------------------------------------
  sim_env_current <- Setup_sim_env(sim_list = sim_list)
  sim_env_current <- add_aggregated_obj_to_simenv(sim_env = sim_env_current, type = 'faa', faa_n_fish_fleets = 4, faa_n_srv_fleets = 4)
  sim_env_current <- run_faa_closedloop_parallel(sim_env = sim_env_current,
                                                 n_sims = n_sims,
                                                 fleet_allocation = fleet_allocation,
                                                 srv_idx_se =  0.2,
                                                 age_lag = 1,
                                                 lls_design_type = 'current',
                                                 srv_wgt = 'numbers',
                                                 fish_wgt = 'numbers',
                                                 faa_n_fish_fleets = ems[[19]]$faa_n_fish_fleets,
                                                 faa_n_srv_fleets = ems[[19]]$faa_n_srv_fleets,
                                                 fish_sel_model = ems[[19]]$fish_sel_model,
                                                 srv_sel_model = ems[[19]]$srv_sel_model,
                                                 fish_selex_prior = ems[[19]]$fish_selex_prior,
                                                 srv_selex_prior = ems[[19]]$srv_selex_prior,
                                                 fish_sel_blocks = ems[[19]]$fish_sel_blocks,
                                                 srv_sel_blocks = ems[[19]]$srv_sel_blocks,
                                                 n_cores = 13)

  saveRDS(sim_env_current, here("outputs", "mse_results", "spatial_noblock_scenarios", paste("faa_base", "_", i, ".RDS", sep = "")))

}

# Time block survey selex
for(i in 1:3) {
  # Condition closed-loop simulations, random recruitment
  sim_list <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = n_sims,
    data = data,
    parameters = parameters,
    mapping = mapping,
    sd_rep = sd_rep,
    rep = rep,
    random = NULL,
    recruitment_opt = 'resample_from_input',
    # setup variances
    ISS_FishAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
    ISS_FishLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
    ISS_SrvAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
    ISS_SrvLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
    ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
    ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets))
  )

  # Get constant B40 for comparison
  global_spr <- Get_Reference_Points(data = om_values$data,
                                     rep = om_values$rep,
                                     SPR_x = 0.4,
                                     type = 'multi_region',
                                     what = 'global_SPR',
                                     calc_rec_st_yr = 20,
                                     rec_age = 2
  )

  # Run MSEs ----------------------------------------------------------------
  sim_env_current <- Setup_sim_env(sim_list = sim_list)
  sim_env_current <- add_aggregated_obj_to_simenv(sim_env = sim_env_current, type = 'faa', faa_n_fish_fleets = 4, faa_n_srv_fleets = 4)
  sim_env_current <- run_faa_closedloop_parallel(sim_env = sim_env_current,
                                                 n_sims = n_sims,
                                                 fleet_allocation = fleet_allocation,
                                                 srv_idx_se =  0.2,
                                                 age_lag = 1,
                                                 lls_design_type = 'current',
                                                 srv_wgt = 'numbers',
                                                 fish_wgt = 'numbers',
                                                 faa_n_fish_fleets = ems[[187]]$faa_n_fish_fleets,
                                                 faa_n_srv_fleets = ems[[187]]$faa_n_srv_fleets,
                                                 fish_sel_model = ems[[187]]$fish_sel_model,
                                                 srv_sel_model = ems[[187]]$srv_sel_model,
                                                 fish_selex_prior = ems[[187]]$fish_selex_prior,
                                                 srv_selex_prior = ems[[187]]$srv_selex_prior,
                                                 fish_sel_blocks = ems[[187]]$fish_sel_blocks,
                                                 srv_sel_blocks = ems[[187]]$srv_sel_blocks,
                                                 n_cores = 13)

  saveRDS(sim_env_current, here("outputs", "mse_results", "spatial_noblock_scenarios", paste("faa_base_tb", "_", i, ".RDS", sep = "")))

}
