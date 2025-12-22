# Purpose: To simulate data for cross-testing models
# Notes: Runs an MSE loop for one year, but anything prior to n_years can essentially
# be used as a self / cross test without any feedback. A bit slower, but easy enough for our purposes here.
# Creator: Matthew LH. Cheng
# Date: 12/9/25

# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)
library(future.apply)
library(progressr)
library(furrr)

# Load in functions
source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "single_region_em.R"))

# Read in model output
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_model_results.RDS"))

# Condition OM ------------------------------------------------------------

# Extract model components from sablefish model
data <- om_values$data
rep <- om_values$rep
parameters <- om_values$parameters
mapping <- om_values$mapping
sd_rep <- om_values$sd_rep

# Define operating model parameters
closed_loop_yrs <- 1      # Years to project forward
n_years <- length(data$years)  # number of years
burnin_years <- 1:n_years  # Historical conditioning period
n_sims <- 100              # Number of replicate simulations
n_regions <- 5             # number of regions
n_fish_fleets <- 2         # number of fishery fleets
n_srv_fleets <- 3          # number of survey fleets

# get fleet allocation by region
fleet_allocation <- array(NA, dim = c(n_regions, n_fish_fleets))
fleet_allocation[,1] <- c(0.5,0.75,0.8,0.8,0.95) # from fmp
fleet_allocation[,2] <- 1 - fleet_allocation[,1] # trawl gear allocation

# Condition closed-loop simulations, mean recruitment
sim_list_highsamp <- condition_closed_loop_simulations(
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
  ISS_FishAgeComps = array(500, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(500, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(500, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(500, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets))
)

# Condition closed-loop simulations, mean recruitment
sim_list_lowsamp <- condition_closed_loop_simulations(
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

# Condition closed-loop simulations, mean recruitment w/ 30 years closed loop
closed_loop_yrs <- 30
sim_list_lowsamp_long <- condition_closed_loop_simulations(
  closed_loop_yrs = closed_loop_yrs, # use 30 closed loop years instead of 1,
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

# Regime (Low recruitment, long)
closed_loop_yrs <- 30
R0_input <- array(0, dim = c(n_regions, n_years + closed_loop_yrs))
R0_input[] <- om_values$rep$R0 * om_values$rep$Rec_trans_prop # mean recruitment from EM in all years
R0_input[,((n_years + 1):(n_years + 15))] <- om_values$rep$R0 * om_values$rep$Rec_trans_prop * 0.15
R0_input <- replicate(n = n_sims, R0_input)

# Condition closed-loop simulations, regime recruitment
sim_list_lowregime <- condition_closed_loop_simulations(
  closed_loop_yrs = closed_loop_yrs,
  n_sims = n_sims,
  data = data,
  parameters = parameters,
  mapping = mapping,
  sd_rep = sd_rep,
  rep = rep,
  random = NULL,
  recruitment_opt = 'mean_rec',
  # setup variances
  ISS_FishAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets)),
  R0_input = R0_input
)

# Regime (High recruitment, long)
closed_loop_yrs <- 30
R0_input <- array(0, dim = c(n_regions, n_years + closed_loop_yrs))
R0_input[] <- om_values$rep$R0 * om_values$rep$Rec_trans_prop # mean recruitment from EM in all years
R0_input[,((n_years + 1):(n_years + 15))] <- om_values$rep$R0 * om_values$rep$Rec_trans_prop * 3
R0_input <- replicate(n = n_sims, R0_input)

# Condition closed-loop simulations, regime recruitment
sim_list_highregime <- condition_closed_loop_simulations(
  closed_loop_yrs = closed_loop_yrs,
  n_sims = n_sims,
  data = data,
  parameters = parameters,
  mapping = mapping,
  sd_rep = sd_rep,
  rep = rep,
  random = NULL,
  recruitment_opt = 'mean_rec',
  # setup variances
  ISS_FishAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets)),
  R0_input = R0_input
)


# Run MSEs to get "data" ----------------------------------------------------------------
# A bit inefficent, but easy enough
# Single-region high sample size
sim_env_highsamp <- Setup_sim_env(sim_list = sim_list_highsamp)
sim_env_highsamp <- add_aggregated_obj_to_simenv(sim_env = sim_env_highsamp)
sim_env_highsamp <- run_single_rg_closedloop_parallel(sim_env = sim_env_highsamp, n_sims = n_sims,
                                                      fleet_allocation = fleet_allocation,
                                                      lls_design_type = "all", n_cores = 7)
saveRDS(sim_env_highsamp, here("outputs", "cross_test", "spatial_blocked_scenarios", "spt_rand_OM_highsamp.RDS"))

# Single-region low sample size
sim_env_lowsamp <- Setup_sim_env(sim_list = sim_list_lowsamp)
sim_env_lowsamp <- add_aggregated_obj_to_simenv(sim_env = sim_env_lowsamp)
sim_env_lowsamp <- run_single_rg_closedloop_parallel(sim_env = sim_env_lowsamp, n_sims = n_sims,
                                                     fleet_allocation = fleet_allocation,
                                                     lls_design_type = "all", n_cores = 7)
saveRDS(sim_env_lowsamp, here("outputs", "cross_test", "spatial_blocked_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Single-region low sample size (long closed loop year)
sim_env_lowsamp_long <- Setup_sim_env(sim_list = sim_list_lowsamp_long)
sim_env_lowsamp_long <- add_aggregated_obj_to_simenv(sim_env = sim_env_lowsamp_long)
sim_env_lowsamp_long <- run_single_rg_closedloop_parallel(sim_env = sim_env_lowsamp_long, n_sims = n_sims,
                                                     fleet_allocation = fleet_allocation,
                                                     lls_design_type = "all", n_cores = 7)
saveRDS(sim_env_lowsamp_long, here("outputs", "cross_test", "spatial_blocked_scenarios", "spt_rand_OM_lowsamp_long.RDS"))

# Single-region low sample size, low regime (long closed loop year)
sim_env_lowsamp_lowreg_long <- Setup_sim_env(sim_list = sim_list_lowregime)
sim_env_lowsamp_lowreg_long <- add_aggregated_obj_to_simenv(sim_env = sim_env_lowsamp_lowreg_long)
sim_env_lowsamp_lowreg_long <- run_single_rg_closedloop_parallel(sim_env = sim_env_lowsamp_lowreg_long, n_sims = n_sims,
                                                                 fleet_allocation = fleet_allocation,
                                                                 lls_design_type = "all", n_cores = 7)
saveRDS(sim_env_lowsamp_lowreg_long, here("outputs", "cross_test", "spatial_blocked_scenarios", "spt_rand_OM_lowsamp_lowreglong.RDS"))

# Single-region low sample size, low regime (long closed loop year)
sim_env_lowsamp_highreg_long <- Setup_sim_env(sim_list = sim_list_highregime)
sim_env_lowsamp_highreg_long <- add_aggregated_obj_to_simenv(sim_env = sim_env_lowsamp_highreg_long)
sim_env_lowsamp_highreg_long <- run_single_rg_closedloop_parallel(sim_env = sim_env_lowsamp_highreg_long, n_sims = n_sims,
                                                                  fleet_allocation = fleet_allocation,
                                                                  lls_design_type = "all", n_cores = 7)
saveRDS(sim_env_lowsamp_highreg_long, here("outputs", "cross_test", "spatial_blocked_scenarios", "spt_rand_OM_lowsamp_highreglong.RDS"))
