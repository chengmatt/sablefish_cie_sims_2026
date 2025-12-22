# Purpose: To self test the single region EM prior to the MSE year "R/cross_test/spatial_no_block_scenarios/sglrg_em_crosstest.R"
# across a combination of weighting schemes for the survey and fishery
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
lowsamp <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

# Run EMs, low sample OM (w/o Francis) -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    # get single region data
    asmt_list <- single_region_em(
      sim_env = lowsamp,
      y = y,
      sim = i,
      srv_idx_se = 0.2,
      age_lag = 1,
      lls_design_type = "historical",
      srv_wgt = 'numbers',
      fish_wgt = 'numbers'
    )

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
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp.RDS"))

# Run EMs, low sample OM (w Francis) -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp_francis <- future_lapply(1:n_sims, function(i) {

    # get single region data
    asmt_list <- single_region_em(
      sim_env = lowsamp,
      y = y,
      sim = i,
      srv_idx_se = 0.2,
      age_lag = 1,
      lls_design_type = "historical",
      srv_wgt = 'numbers',
      fish_wgt = 'numbers'
    )

    # fit model
    model_francis <- run_francis(
      asmt_list$data,
      asmt_list$par,
      asmt_list$map,
      NULL,
      n_francis_iter = 10,
      newton_loops = 2
    )

    model <- model_francis$obj # return model object
    model$data <- asmt_list$data # save data
    model$sd_rep <- sdreport(model)

    p()
    model
  })
})

# Sim, EM list dimension
saveRDS(model_list_lowsamp_francis, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp_francis.RDS"))

# Process Results ---------------------------------------------------------
model_list_lowsamp <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp.RDS"))
model_list_lowsamp_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp_francis.RDS"))

ssb_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
rec_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
rec_om_results <- array(NA, dim = c(length(1:y), n_sims))

for(i in 1:n_sims) {
  # get aggregated OM results
  ssb_om_results[,i] <- colSums(lowsamp$SSB[,1:y,i])
  rec_om_results[,i] <- colSums(lowsamp$Rec[,1:y,i])

  # get EM results
  ssb_em_results[1,,i] <- t(model_list_lowsamp[[i]]$rep$SSB)
  rec_em_results[1,,i] <- t(model_list_lowsamp[[i]]$rep$Rec)
  ssb_em_results[2,,i] <- t(model_list_lowsamp_francis[[i]]$rep$SSB)
  rec_em_results[2,,i] <- t(model_list_lowsamp_francis[[i]]$rep$Rec)
} # end i


par(mfrow = c(2,2))
# SSB
plot(apply((ssb_em_results[1,,] - ssb_om_results) / ssb_om_results, 1, median), type = 'l', ylim = c(-0.25, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in SSB')
lines(apply((ssb_em_results[2,,] - ssb_om_results) / ssb_om_results, 1, median), col = 'red', lwd = 3, lty = 2)
abline(h = 0, lty = 2, lwd = 3)
# SSB abs
plot(apply(abs((ssb_em_results[1,,] - ssb_om_results) / ssb_om_results), 1, median), type = 'l', ylim = c(0, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in SSB')
lines(apply(abs((ssb_em_results[2,,] - ssb_om_results) / ssb_om_results), 1, median), col = 'red', lwd = 3, lty = 2)

# Rec
plot(apply((rec_em_results[1,,] - rec_om_results) / rec_om_results, 1, median), type = 'l', ylim = c(-1, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in Rec')
lines(apply((rec_em_results[2,,] - rec_om_results) / rec_om_results, 1, median), col = 'red', lwd = 3, lty = 2)
abline(h = 0, lty = 2, lwd = 3)
# Rec abs
plot(apply(abs((rec_em_results[1,,] - rec_om_results) / rec_om_results), 1, median), type = 'l', ylim = c(0, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in Rec')
lines(apply(abs((rec_em_results[2,,] - rec_om_results) / rec_om_results), 1, median), col = 'red', lwd = 3, lty = 2)

# Single-region models exhibit positive bias, relative to spatial model
# Likely because the spatial heterogeneity in the composition data cannot
# be adqeuately mimiked by a single-region model, even after accounting for
# catch weighitng (composition data looks flatter / smears / confounds the recruitment signal more).
# This then translates into differences in what is "selected/observed"
# by the fishery and survey, with downstream impacts on recruitment. For instance, the
# domestic survey fleet in the spatial model tends to select more younger fish relative to what is
# estimated by the single region model, thus providing a more consistent recruitment signal. So
# during high recruitment events, there are lot more young fish being selected by the survey in the spatial
# model but because the survey is: 1) forced to be asymptotic such as not to select a lot of young fish, and 2)
# the comps don't necessarily fully reflect these recruitment events from the BS when aggregated by weighting spatially,
# the selectivity becomes biased. Because survey selex is connected to the survey index, with a lower selectivity, you
# can't really fit to the index as well given the new influex of recruitments, so you need to then overestimate recruitment
# to fit to index.

# The exception to this is during the more recent yeras, most likely due to
# not having a good enough recruitment signal so it leads to a negative bias
# because you haven't quite seen those recruits just yet - the model
# will eventually retroactively revise those recruitment estimates to be higher.
# The model can't quite fit to those indices in the last couple of years, but is within that
# 95% bar (tho generally always below)

fleet <- 1
sim <- 88
yrs <- 1:65
# plot(apply(model_list_lowsamp[[sim]]$data$ObsSrvAgeComps[,-65,,1,1], 2, mean) /
#        sum(apply(model_list_lowsamp[[sim]]$data$ObsSrvAgeComps[,-65,,1,1], 2, mean)),  ylim = c(0, 0.2))
# lines(apply(lowsamp$SrvIAA[,yrs,,1,fleet,sim], 3, sum) / sum(apply(lowsamp$SrvIAA[,yrs,,1,fleet,sim], 3, sum)))
#
# sum(
#   apply(model_list_lowsamp[[sim]]$data$ObsSrvAgeComps[,-65,,1,1], 2, mean) /
#     sum(apply(model_list_lowsamp[[sim]]$data$ObsSrvAgeComps[,-65,,1,1], 2, mean)) * 1:30
# )
#
# sum(
#   apply(lowsamp$SrvIAA[,yrs,,1,fleet,sim], 3, sum) / sum(apply(lowsamp$SrvIAA[,yrs,,1,fleet,sim], 3, sum)) * 1:30
# )
model_fra <- SPoRC::run_francis(
  asmt_list$data,
  asmt_list$par,
  asmt_list$map, random = NULL, n_francis_iter = 10, newton_loops = 3
)
get_idx_fits_plot(list(asmt_list$data), list(model_fra$obj$rep), 1)
model_fra$recorded_weights %>% view()


