# Purpose: To do cross-testing using a 5 region spatial model
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
source(here("R", "functions", "five_rg_em.R"))

# Read in OMs
lowsamp <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))
highsamp <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_highsamp.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

# Run EMs, low sample OM -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    # get five region data
    asmt_list <- five_rg_em(
      sim_env = lowsamp,
      y = y,
      sim = i,
      srv_idx_se = 0.2,
      age_lag = 1,
      lls_design_type = 'historical',
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
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "five_region_crosstest_lowsamp.RDS"))

# Run EMs, high sample EM -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_highsamp <- future_lapply(1:n_sims, function(i) {

    # get five region data
    asmt_list <- five_rg_em(
      sim_env = highsamp,
      y = y,
      sim = i,
      srv_idx_se = 0.1,
      age_lag = 1,
      lls_design_type = 'all',
      UseTagging = 1
    )

    # Use data every year
    asmt_list$data$UseFishAgeComps[which(asmt_list$data$UseCatch[,,1] == 1)] <- 1
    asmt_list$data$UseFishLenComps[which(asmt_list$data$UseCatch == 1)] <- 1
    asmt_list$data$UseSrvIdx[,,c(1,3)] <- 1
    asmt_list$data$UseSrvAgeComps[,,c(1,3)] <- 1

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
saveRDS(model_list_highsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "five_region_crosstest_highsamp.RDS"))

#
# # # Process Results ---------------------------------------------------------
# model_list_highsamp <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "five_region_crosstest_lowsamp.RDS"))
#
# # storage containers
# agg_ssb_store <- array(NA, dim = c(65, 100))
# agg_rec_store <- array(NA, dim = c(65, 100))
# agg_dep_store <- array(NA, dim = c(65, 100))
#
#
# counter <- 0
# for(i in 1:100) {
#   if(length(model_list_highsamp[[i]]) > 1)agg_ssb_store[,i] <- colSums(model_list_highsamp[[i]]$rep$SSB[,])
#   if(length(model_list_highsamp[[i]]) > 1){
#     agg_rec_store[,i] <- colSums(model_list_highsamp[[i]]$rep$Rec[,])
#     counter <- counter + 1
#   }
#
# }
#
# model_list_highsamp[[1]]$data$ObsSrvIdx[which(model_list_highsamp[[1]]$data$UseSrvIdx == 0)] <- NA
# get_catch_fits_plot(list(model_list_highsamp[[1]]$data), list(model_list_highsamp[[1]]$rep), 1)
# get_idx_fits_plot(list(model_list_highsamp[[1]]$data), list(model_list_highsamp[[1]]$rep), 1)
#
# plot(apply(agg_ssb_store, 1, median))
# lines(colSums(lowsamp$SSB[,-66,1]))
# plot((apply(agg_ssb_store, 1, median) - colSums(lowsamp$SSB[,-66,1])) / colSums(lowsamp$SSB[,-66,1]))
#
# plot(apply(agg_rec_store, 1, median))
# lines(colSums(lowsamp$Rec[,-66,1]))
# plot((apply(agg_rec_store, 1, median) - colSums(lowsamp$Rec[,-66,1])) / colSums(lowsamp$Rec[,-66,1]))
#
#
# ggplot() +
#   geom_line(reshape2::melt(model_list_highsamp[[6]]$rep$Movement),
#             mapping = aes(x = ages, y = value)) +
#   geom_line(reshape2::melt(lowsamp$Movement[,,,,,1]) %>%
#               rename(from = Var1, to = Var2, years = Var3, ages = Var4, sexes = Var5),
#             mapping = aes(x = ages, y = value), color ='red') +
#   ggh4x::facet_grid2(paste("to", to)~paste("from", from), scales = 'free', independent = 'y')
