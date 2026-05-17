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
srvchng <- readRDS(here("outputs", "mse_results", "spatial_noblock_scenarios", "five_region_base.RDS"))

# Simulation dimensions
y <- 95 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

# Run EMs, low sample OM -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    # get three region data (need to iterate to build up data frame)
    for(y1 in 65:y) {
      asmt_list <- three_rg_em(
        sim_env = srvchng,
        y = y1,
        sim = i,
        srv_idx_se = 0.2,
        age_lag = 1,
        lls_design_type = 'current',
        srv_wgt = 'numbers',
        fish_wgt = 'numbers',
        UseTagging = 1
      )
    }

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
saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_srvchng_current.RDS"))

# Process Results ---------------------------------------------------------
# model_list_srvchng <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_srvchng_current.RDS"))
#
# # storage containers
# ssb_store <- array(NA, dim = c(3, 65, 100))
# rec_store <- array(NA, dim = c(3, 65, 100))
# dep_store <- array(NA, dim = c(3, 65, 100))
# agg_ssb_store <- array(NA, dim = c(65, 100))
# agg_rec_store <- array(NA, dim = c(65, 100))
# agg_dep_store <- array(NA, dim = c(65, 100))
#
# counter <- 0
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     counter <- counter + 1
#     agg_ssb_store[,i] <- (colSums(model_list_srvchng[[i]]$rep$SSB) - colSums(srvchng$SSB[,-66,1])) / colSums(srvchng$SSB[,-66,1])
#     agg_rec_store[,i] <- (colSums(model_list_srvchng[[i]]$rep$Rec) - colSums(srvchng$Rec[,-66,1])) / colSums(srvchng$Rec[,-66,1])
#   }
# }
#
# plot(apply(agg_ssb_store, 1, quantile, na.rm = T)[2,], type = 'l', ylim = c(-0.5, 0.5))
# lines(apply(agg_ssb_store, 1, median, na.rm = T), lty = 2)
# lines(apply(agg_ssb_store, 1, quantile, na.rm = T)[4,], type = 'l')
#
#
# plot(apply(agg_rec_store, 1, quantile, na.rm = T)[1,], type = 'l', ylim = c(-2, 2))
# lines(apply(agg_rec_store, 1, median, na.rm = T), lty = 2)
# lines(apply(agg_rec_store, 1, quantile, na.rm = T)[4,], type = 'l')
#
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     if(i == 1) plot(colSums(model_list_srvchng[[i]]$rep$SSB), type = 'l', ylim = c(0, 320))
#     else lines(colSums(model_list_srvchng[[i]]$rep$SSB), type = 'l', ylim = c(0, 320))
#   }
# }
#
# lines(colSums(srvchng$SSB[,-66,1]), lty = 2, col = 'red', lwd = 3)
#
#
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     if(i == 1) plot(colSums(model_list_srvchng[[i]]$rep$Rec[1,,drop = F]), type = 'l', ylim = c(0, 130))
#     else lines(colSums(model_list_srvchng[[i]]$rep$Rec[1,,drop = F]), type = 'l', ylim = c(0, 130))
#   }
# }
#
# lines(colSums(srvchng$Rec[1,,1, drop = F]), lty = 2, col = 'red', lwd = 3)

# All design -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    # get three region data (need to iterate to build up data frame)
    for(y1 in 65:y) {
      asmt_list <- three_rg_em(
        sim_env = srvchng,
        y = y1,
        sim = i,
        srv_idx_se = 0.2,
        age_lag = 1,
        lls_design_type = 'all',
        srv_wgt = 'numbers',
        fish_wgt = 'numbers',
        UseTagging = 1
      )
    }

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
saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_srvchng_all.RDS"))

# Process Results ---------------------------------------------------------
# model_list_srvchng <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_srvchng_all.RDS"))
#
# # storage containers
# ssb_store <- array(NA, dim = c(3, 65, 100))
# rec_store <- array(NA, dim = c(3, 65, 100))
# dep_store <- array(NA, dim = c(3, 65, 100))
# agg_ssb_store <- array(NA, dim = c(65, 100))
# agg_rec_store <- array(NA, dim = c(65, 100))
# agg_dep_store <- array(NA, dim = c(65, 100))
#
# counter <- 0
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     counter <- counter + 1
#     agg_ssb_store[,i] <- (colSums(model_list_srvchng[[i]]$rep$SSB) - colSums(srvchng$SSB[,-66,1])) / colSums(srvchng$SSB[,-66,1])
#     agg_rec_store[,i] <- (colSums(model_list_srvchng[[i]]$rep$Rec) - colSums(srvchng$Rec[,-66,1])) / colSums(srvchng$Rec[,-66,1])
#   }
# }
#
# plot(apply(agg_ssb_store, 1, quantile, na.rm = T)[2,], type = 'l', ylim = c(-0.5, 0.5))
# lines(apply(agg_ssb_store, 1, median, na.rm = T), lty = 2)
# lines(apply(agg_ssb_store, 1, quantile, na.rm = T)[4,], type = 'l')
#
#
# plot(apply(agg_rec_store, 1, quantile, na.rm = T)[1,], type = 'l', ylim = c(-2, 2))
# lines(apply(agg_rec_store, 1, median, na.rm = T), lty = 2)
# lines(apply(agg_rec_store, 1, quantile, na.rm = T)[4,], type = 'l')
#
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     if(i == 1) plot(colSums(model_list_srvchng[[i]]$rep$SSB), type = 'l', ylim = c(0, 320))
#     else lines(colSums(model_list_srvchng[[i]]$rep$SSB), type = 'l', ylim = c(0, 320))
#   }
# }
#
# lines(colSums(srvchng$SSB[,-66,1]), lty = 2, col = 'red', lwd = 3)
#
#
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     if(i == 1) plot(colSums(model_list_srvchng[[i]]$rep$Rec[1,,drop = F]), type = 'l', ylim = c(0, 130))
#     else lines(colSums(model_list_srvchng[[i]]$rep$Rec[1,,drop = F]), type = 'l', ylim = c(0, 130))
#   }
# }
#
# lines(colSums(srvchng$Rec[1,,1, drop = F]), lty = 2, col = 'red', lwd = 3)

# Historical design -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    # get three region data (need to iterate to build up data frame)
    for(y1 in 65:y) {
      asmt_list <- three_rg_em(
        sim_env = srvchng,
        y = y1,
        sim = i,
        srv_idx_se = 0.2,
        age_lag = 1,
        lls_design_type = 'historical',
        srv_wgt = 'numbers',
        fish_wgt = 'numbers',
        UseTagging = 1
      )
    }

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
saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_srvchng_historical.RDS"))

# Process Results ---------------------------------------------------------
# model_list_srvchng <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_srvchng_historical.RDS"))
#
# # storage containers
# ssb_store <- array(NA, dim = c(3, 65, 100))
# rec_store <- array(NA, dim = c(3, 65, 100))
# dep_store <- array(NA, dim = c(3, 65, 100))
# agg_ssb_store <- array(NA, dim = c(65, 100))
# agg_rec_store <- array(NA, dim = c(65, 100))
# agg_dep_store <- array(NA, dim = c(65, 100))
#
# counter <- 0
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     counter <- counter + 1
#     agg_ssb_store[,i] <- (colSums(model_list_srvchng[[i]]$rep$SSB) - colSums(srvchng$SSB[,-66,1])) / colSums(srvchng$SSB[,-66,1])
#     agg_rec_store[,i] <- (colSums(model_list_srvchng[[i]]$rep$Rec) - colSums(srvchng$Rec[,-66,1])) / colSums(srvchng$Rec[,-66,1])
#   }
# }
#
# plot(apply(agg_ssb_store, 1, quantile, na.rm = T)[2,], type = 'l', ylim = c(-0.5, 0.5))
# lines(apply(agg_ssb_store, 1, median, na.rm = T), lty = 2)
# lines(apply(agg_ssb_store, 1, quantile, na.rm = T)[4,], type = 'l')
#
#
# plot(apply(agg_rec_store, 1, quantile, na.rm = T)[1,], type = 'l', ylim = c(-2, 2))
# lines(apply(agg_rec_store, 1, median, na.rm = T), lty = 2)
# lines(apply(agg_rec_store, 1, quantile, na.rm = T)[4,], type = 'l')
#
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     if(i == 1) plot(colSums(model_list_srvchng[[i]]$rep$SSB), type = 'l', ylim = c(0, 320))
#     else lines(colSums(model_list_srvchng[[i]]$rep$SSB), type = 'l', ylim = c(0, 320))
#   }
# }
#
# lines(colSums(srvchng$SSB[,-66,1]), lty = 2, col = 'red', lwd = 3)
#
#
# for(i in 1:100) {
#   if(length(model_list_srvchng[[i]]) > 1) {
#     if(i == 1) plot(colSums(model_list_srvchng[[i]]$rep$Rec[1,,drop = F]), type = 'l', ylim = c(0, 130))
#     else lines(colSums(model_list_srvchng[[i]]$rep$Rec[1,,drop = F]), type = 'l', ylim = c(0, 130))
#   }
# }
#
# lines(colSums(srvchng$Rec[1,,1, drop = F]), lty = 2, col = 'red', lwd = 3)
