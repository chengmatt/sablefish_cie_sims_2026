# Purpose: To self test the single region EM w/ an extended time series incorporating a change in survey design
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

# Read in OMs
srvchng <- readRDS(here("outputs", "mse_results", "spatial_noblock_scenarios", "five_region_base.RDS"))

# Simulation dimensions
y <- 95 # terminal year
n_sims <- 100 # number of sims

# Current design  -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 10)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    # build single region data (need to iterate to build up data object)
    for(y1 in 65:y) {
      asmt_list <- single_region_em(
        sim_env = srvchng,
        y = y1,
        sim = i,
        srv_idx_se = 0.2,
        age_lag = 1,
        lls_design_type = "current",
        srv_wgt = 'numbers',
        fish_wgt = 'numbers'
      )
    } # end y1 loop

    # fit model
    model <- fit_model(
      asmt_list$data,
      asmt_list$par,
      asmt_list$map,
      NULL,
      2,
      silent = TRUE
    )

    model$data <- asmt_list$data # save data
    model$sd_rep <- sdreport(model)

    p()
    model
  })
})

# Sim, EM list dimension
saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_srvchng_current.RDS"))


# Process Results ---------------------------------------------------------

model_list_srvchng <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios",
                                   "single_region_crosstest_srvchng_current.RDS"))

ssb_em_results <- array(NA, dim = c( length(1:y), n_sims))
rec_em_results <- array(NA, dim = c( length(1:y), n_sims))
ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
rec_om_results <- array(NA, dim = c(length(1:y), n_sims))

for(i in 1:n_sims) {
  # get aggregated OM results
  ssb_om_results[,i] <- colSums(srvchng$SSB[,1:y,i])
  rec_om_results[,i] <- colSums(srvchng$Rec[,1:y,i])

  # get EM results
  ssb_em_results[,i] <- t(model_list_srvchng[[i]]$rep$SSB)
  rec_em_results[,i] <- t(model_list_srvchng[[i]]$rep$Rec)
} # end i


par(mfrow = c(2,2))
# SSB
plot(apply((ssb_em_results - ssb_om_results) / ssb_om_results, 1, median), type = 'l', ylim = c(-0.25, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in SSB')
abline(h = 0, lty = 2, lwd = 3)

# SSB abs
plot(apply(abs((ssb_em_results - ssb_om_results) / ssb_om_results), 1, median), type = 'l', ylim = c(0, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in SSB')

# Rec
plot(apply((rec_em_results - rec_om_results) / rec_om_results, 1, median), type = 'l', ylim = c(-1, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in Rec')
abline(h = 0, lty = 2, lwd = 3)
# Rec abs
plot(apply(abs((rec_em_results - rec_om_results) / rec_om_results), 1, median), type = 'l', ylim = c(0, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in Rec')

# All design -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 10)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    # build single region data (need to iterate to build up data object)
    for(y1 in 65:y) {
      asmt_list <- single_region_em(
        sim_env = srvchng,
        y = y1,
        sim = i,
        srv_idx_se = 0.2,
        age_lag = 1,
        lls_design_type = "all",
        srv_wgt = 'numbers',
        fish_wgt = 'numbers'
      )
    } # end y1 loop

    # fit model
    model <- fit_model(
      asmt_list$data,
      asmt_list$par,
      asmt_list$map,
      NULL,
      2,
      silent = TRUE
    )

    model$data <- asmt_list$data # save data
    model$sd_rep <- sdreport(model)

    p()
    model
  })
})

# Sim, EM list dimension
saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_srvchng_all.RDS"))

# Process Results ---------------------------------------------------------

model_list_srvchng <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios",
                                   "single_region_crosstest_srvchng_all.RDS"))

ssb_em_results <- array(NA, dim = c( length(1:y), n_sims))
rec_em_results <- array(NA, dim = c( length(1:y), n_sims))
ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
rec_om_results <- array(NA, dim = c(length(1:y), n_sims))

for(i in 1:n_sims) {
  # get aggregated OM results
  ssb_om_results[,i] <- colSums(srvchng$SSB[,1:y,i])
  rec_om_results[,i] <- colSums(srvchng$Rec[,1:y,i])

  # get EM results
  ssb_em_results[,i] <- t(model_list_srvchng[[i]]$rep$SSB)
  rec_em_results[,i] <- t(model_list_srvchng[[i]]$rep$Rec)
} # end i


par(mfrow = c(2,2))
# SSB
plot(apply((ssb_em_results - ssb_om_results) / ssb_om_results, 1, median), type = 'l', ylim = c(-0.25, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in SSB')
abline(h = 0, lty = 2, lwd = 3)

# SSB abs
plot(apply(abs((ssb_em_results - ssb_om_results) / ssb_om_results), 1, median), type = 'l', ylim = c(0, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in SSB')

# Rec
plot(apply((rec_em_results - rec_om_results) / rec_om_results, 1, median), type = 'l', ylim = c(-1, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in Rec')
abline(h = 0, lty = 2, lwd = 3)
# Rec abs
plot(apply(abs((rec_em_results - rec_om_results) / rec_om_results), 1, median), type = 'l', ylim = c(0, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in Rec')


# Historical design -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 10)
options(future.globals.maxSize = 15e9)
# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    # build single region data (need to iterate to build up data object)
    for(y1 in 65:y) {
      asmt_list <- single_region_em(
        sim_env = srvchng,
        y = y1,
        sim = i,
        srv_idx_se = 0.2,
        age_lag = 1,
        lls_design_type = "historical",
        srv_wgt = 'numbers',
        fish_wgt = 'numbers'
      )
    } # end y1 loop

    # fit model
    model <- fit_model(
      asmt_list$data,
      asmt_list$par,
      asmt_list$map,
      NULL,
      2,
      silent = TRUE
    )

    model$data <- asmt_list$data # save data
    model$sd_rep <- sdreport(model)

    p()
    model
  })
})

# Sim, EM list dimension
saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_srvchng_historical.RDS"))

# Process Results ---------------------------------------------------------
model_list_srvchng <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios",
                                   "single_region_crosstest_srvchng_historical.RDS"))

ssb_em_results <- array(NA, dim = c( length(1:y), n_sims))
rec_em_results <- array(NA, dim = c( length(1:y), n_sims))
ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
rec_om_results <- array(NA, dim = c(length(1:y), n_sims))

for(i in 1:n_sims) {
  # get aggregated OM results
  ssb_om_results[,i] <- colSums(srvchng$SSB[,1:y,i])
  rec_om_results[,i] <- colSums(srvchng$Rec[,1:y,i])

  # get EM results
  ssb_em_results[,i] <- t(model_list_srvchng[[i]]$rep$SSB)
  rec_em_results[,i] <- t(model_list_srvchng[[i]]$rep$Rec)
} # end i


par(mfrow = c(2,2))
# SSB
plot(apply((ssb_em_results - ssb_om_results) / ssb_om_results, 1, median), type = 'l', ylim = c(-0.25, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in SSB')
abline(h = 0, lty = 2, lwd = 3)

# SSB abs
plot(apply(abs((ssb_em_results - ssb_om_results) / ssb_om_results), 1, median), type = 'l', ylim = c(0, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in SSB')

# Rec
plot(apply((rec_em_results - rec_om_results) / rec_om_results, 1, median), type = 'l', ylim = c(-1, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in Rec')
abline(h = 0, lty = 2, lwd = 3)
# Rec abs
plot(apply(abs((rec_em_results - rec_om_results) / rec_om_results), 1, median), type = 'l', ylim = c(0, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in Rec')

