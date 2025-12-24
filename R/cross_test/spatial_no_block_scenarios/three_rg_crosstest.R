# Purpose: To do cross-testing using a 3 region spatial model
# Creator: Matthew LH. Cheng
# 12/23/25


# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)
library(future.apply)
library(progressr)
library(furrr)

# Load in functions
source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "three_rg_em.R"))

# Read in OMs
lowsamp <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

# Run EMs, low sample OM -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 6)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    # get five region data
    asmt_list <- three_rg_em(
      sim_env = lowsamp,
      y = y,
      sim = i,
      srv_idx_se = 0.2,
      age_lag = 1,
      lls_design_type = 'historical',
      srv_wgt = 'numbers',
      fish_wgt = 'numbers',
      UseTagging = 1
    )

    model <- tryCatch(
      {
        # fit model
        model <- fit_model(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          2,
          silent = F
        )

        model$data <- asmt_list$data    # save data
        model$sd_rep <- sdreport(model)
        model
      },
      error = function(e) {
        NA
      }
    )
    p()
    model
  })
})

# Sim, EM list dimension
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_lowsamp.RDS"))


# Process Results ---------------------------------------------------------
model_list_lowsamp <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_lowsamp.RDS"))

# storage containers
ssb_store <- array(NA, dim = c(3, 65, 100))
rec_store <- array(NA, dim = c(3, 65, 100))
dep_store <- array(NA, dim = c(3, 65, 100))
agg_ssb_store <- array(NA, dim = c(65, 100))
agg_rec_store <- array(NA, dim = c(65, 100))
agg_dep_store <- array(NA, dim = c(65, 100))

counter <- 0
for(i in 1:100) {
  if(length(model_list_lowsamp[[i]]) > 1) {
    counter <- counter + 1
    agg_ssb_store[,i] <- (colSums(model_list_lowsamp[[i]]$rep$SSB) - colSums(lowsamp$SSB[,-66,1])) / colSums(lowsamp$SSB[,-66,1])
    agg_rec_store[,i] <- (colSums(model_list_lowsamp[[i]]$rep$Rec) - colSums(lowsamp$Rec[,-66,1])) / colSums(lowsamp$Rec[,-66,1])
  }
}

plot(apply(agg_ssb_store, 1, median, na.rm = T))
plot(apply(agg_rec_store, 1, median, na.rm = T), type = 'l')
