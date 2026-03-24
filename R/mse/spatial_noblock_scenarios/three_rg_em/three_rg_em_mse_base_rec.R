# Purpose: To run a five region OM, with a three region EM (base recruitment)
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
source(here("R", "functions", "three_rg_em.R"))

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
n_sims <- 10              # Number of replicate simulations
n_regions <- 5             # number of regions
n_fish_fleets <- 2         # number of fishery fleets
n_srv_fleets <- 3          # number of survey fleets

# get fleet allocation by region
fleet_allocation <- array(NA, dim = c(n_regions, n_fish_fleets))
fleet_allocation[,1] <- c(0.5,0.75,0.8,0.8,0.95) # from fmp
fleet_allocation[,2] <- 1 - fleet_allocation[,1] # trawl gear allocation

# To simulate more tags or resimulate tags with new sample sizes
historical_tags <- rbind(data$tag_release_indicator,
                         unname(as.matrix(data.frame(regions = c(2:5), tag_yr = 63))),
                         unname(as.matrix(data.frame(regions = c(1,3:5), tag_yr = 64)))
)

# Get new simulated tags in
new_goa_tags <- expand.grid(regions = 3:5, tag_yr = seq((n_years + 1), (n_years + closed_loop_yrs), 2))
new_bsai_tags <- expand.grid(regions = 1:2, tag_yr = seq((n_years + 2), (n_years + closed_loop_yrs), 2))
data$tag_release_indicator <- rbind(historical_tags, unname(as.matrix(new_goa_tags)), unname(as.matrix(new_bsai_tags)))

# Condition closed-loop simulations, random recruitment
for(i in 1:25) {
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
    ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets)),
    n_tags_rel_input = rep(2e3, nrow(data$tag_release_indicator))
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
  # Single-region, current design
  sim_env_current <- Setup_sim_env(sim_list = sim_list)
  sim_env_current <- add_aggregated_obj_to_simenv(sim_env = sim_env_current, type = "three_rg")
  sim_env_current <- run_three_rg_closedloop_parallel(sim_env = sim_env_current,
                                                      n_sims = n_sims,
                                                      fleet_allocation = fleet_allocation,
                                                      lls_design_type = "current",
                                                      srv_idx_se = 0.2,
                                                      age_lag = 1,
                                                      srv_wgt = 'numbers',
                                                      fish_wgt = 'numbers',
                                                      n_cores = 14)

  saveRDS(sim_env_current, here("outputs", "mse_results", "spatial_noblock_scenarios", paste("three_region_base", "_", i, ".RDS", sep = "")))
}
