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
lowsamp_movesense <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp_tvmove_sens.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

ems_grid <- expand.grid(
  fish_block = c(TRUE, FALSE),
  srv_block = c(TRUE, FALSE)
)

# Run EMs, low sample OM (w/o Francis) -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    # container
    out_i <- vector("list", nrow(ems_grid))

    for(j in 1:nrow(ems_grid)) {
      # get single region data
      asmt_list <- single_region_em(
        sim_env = lowsamp,
        y = y,
        sim = i,
        srv_idx_se = 0.2,
        age_lag = 1,
        lls_design_type = "historical",
        srv_wgt = 'numbers',
        fish_wgt = 'numbers',
        srv_block = ems_grid$srv_block[j],
        fish_block = ems_grid$fish_block[j]
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

      out_i[[j]] <- model
    }

    p()
    out_i

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

# Run EMs, low sample OM (w/o Francis - Sensitivity) -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    # get single region data
    asmt_list <- single_region_em(
      sim_env = lowsamp_movesense,
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
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_movesense.RDS"))

# Run EMs, low sample OM (w Francis - Sensitivity) -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp_francis <- future_lapply(1:n_sims, function(i) {

    # get single region data
    asmt_list <- single_region_em(
      sim_env = lowsamp_movesense,
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
saveRDS(model_list_lowsamp_francis, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_movesense_francis.RDS"))


# Run Single Region Model (w/o Francis Retro) ---------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  retros_wo_francis <- future_lapply(1:n_sims, function(i) {

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

    # fit model via retros
    retros <- do_retrospective(
      8,
      data = asmt_list$data,
      parameters = asmt_list$par,
      mapping = asmt_list$map,
      random = NULL, do_par = F,
      newton_loops = 3,
      do_francis = F
    )
    retros$Sim <- i  # name retros
    p()
    retros

  })
})

retros_wo_francis_df <- data.table::rbindlist(retros_wo_francis)
write.csv(retros_wo_francis_df, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_cross_test_retro.csv"))

# Run Single Region Model (w Francis Retro) ---------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  retros_w_francis <- future_lapply(1:n_sims, function(i) {

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

    # fit model via retros
    retros <- do_retrospective(
      8,
      data = asmt_list$data,
      parameters = asmt_list$par,
      mapping = asmt_list$map,
      random = NULL, do_par = F,
      newton_loops = 3,
      do_francis = T,
      n_francis_iter = 10
    )
    retros$Sim <- i  # name retros
    p()
    retros

  })
})

retros_w_francis_df <- data.table::rbindlist(retros_w_francis)
write.csv(retros_w_francis_df, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_cross_test_retro_francis.csv"))


# Process Results ---------------------------------------------------------
model_list_lowsamp <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp.RDS"))
model_list_lowsamp_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp_francis.RDS"))

ssb_em_results <- array(NA, dim = c(length(1:y), n_sims, nrow(ems_grid)))
rec_em_results <- array(NA, dim = c(length(1:y), n_sims, nrow(ems_grid)))
srv_em_results <- array(NA, dim = c(length(1:y), n_sims, nrow(ems_grid)))
ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
rec_om_results <- array(NA, dim = c(length(1:y), n_sims))
srv_om_results <- array(NA, dim = c(length(1:y), n_sims))

for (i in 1:n_sims) {
  # OM (aggregated across regions)
  ssb_om_results[, i] <- colSums(lowsamp$SSB[, 1:y, i])
  rec_om_results[, i] <- colSums(lowsamp$Rec[, 1:y, i])
  srv_om_results[, i] <- lowsamp$Agg_ObsSrvIdx[,1:y,1,i]

  # EM — all 4 models
  for (m in 1:nrow(ems_grid)) {
    ssb_em_results[ , i, m] <- t(model_list_lowsamp[[i]][[m]]$rep$SSB)
    rec_em_results[ , i, m] <- t(model_list_lowsamp[[i]][[m]]$rep$Rec)
    srv_em_results[, i, m] <- t(model_list_lowsamp[[i]][[m]]$rep$PredSrvIdx[,,1])
  }

  # Model 5 (Francis)
  # ssb_em_results[5, , i] <- t(model_list_5[[i]]$rep$SSB)
  # rec_em_results[5, , i] <- t(model_list_5[[i]]$rep$Rec)
}

# Summarize biases
model_labels <- c("none", "srv_block", "fish_block", "fish_srv_block")
ssb_rel_bias <- sweep(
  sweep(ssb_em_results, c(1, 2), ssb_om_results, FUN = "-"),
  c(1, 2), ssb_om_results, FUN = "/"
)

rec_rel_bias <- sweep(
  sweep(rec_em_results, c(1, 2), rec_om_results, FUN = "-"),
  c(1, 2), rec_om_results, FUN = "/"
)

srv_rel_bias <- sweep(
  sweep(srv_em_results, c(1, 2), srv_om_results, FUN = "-"),
  c(1, 2), srv_om_results, FUN = "/"
)

summarise_bias <- function(bias_array, quantity_name) {
  expand.grid(year = 1:y, sim = 1:n_sims, model = 1:nrow(ems_grid)) |>
    mutate(
      rel_bias = as.vector(bias_array),
      model    = factor(model, labels = model_labels),
      quantity = quantity_name
    )
}

bias_df <- bind_rows(
  summarise_bias(ssb_rel_bias, "SSB"),
  summarise_bias(rec_rel_bias, "Recruitment"),
  summarise_bias(srv_rel_bias, "Survey"),

)

# Quantile summary across sims
bias_summary <- bias_df |>
  group_by(model, year, quantity) |>
  summarise(
    med   = median(rel_bias),
    lo50  = quantile(rel_bias, 0.25),
    hi50  = quantile(rel_bias, 0.75),
    lo90  = quantile(rel_bias, 0.05),
    hi90  = quantile(rel_bias, 0.95),
    .groups = "drop"
  )

ggplot(bias_summary, aes(x = year, colour = model, fill = model)) +
  # median
  geom_line(aes(y = med), linewidth = 0.8) +
  facet_wrap(~quantity, ncol = 1, scales = "free_y",
             labeller = labeller(quantity = c(SSB = "SSB", Recruitment = "Recruitment"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 15)

ggplot(bias_summary %>% filter(year >= 50), aes(x = year, colour = model, fill = model)) +
  # median
  geom_line(aes(y = med), linewidth = 0.8) +
  facet_wrap(~quantity, ncol = 1, scales = "free_y",
             labeller = labeller(quantity = c(SSB = "SSB", Recruitment = "Recruitment"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 15)

# Single region model blocking doesn't really decrease bias
# Think a bit part of the bias from the single region models is aggregating data and trying to
# fit data that are not really data (i.e., carry over values, etc)
# Trying to fit carry over values with time blocks leads to underestiamting recruitment.

# Look at estimated survey q
srv_q_results <- array(NA, dim = c(n_sims, nrow(ems_grid)))
for (i in 1:n_sims) {
  for (m in 1:nrow(ems_grid)) {
    srv_q_results[i, m] <- unique(as.vector(model_list_lowsamp[[i]][[m]]$rep$srv_q))[1]
  }
}
apply(srv_q_results, 2, mean)

# Look at survey age compositions being fitted to (no blocking)
par(mfrow = c(1,2))
plot(as.vector(apply(model_list_lowsamp[[i]][[1]]$data$ObsSrvAgeComps[1,1:64,,,1], 2:3, mean))
     / sum(as.vector(apply(model_list_lowsamp[[i]][[1]]$data$ObsSrvAgeComps[1,1:64,,,1], 2:3, mean))), ylim = c(0,0.15),
     ylab = 'Proportion (1960:2025)', xlab = 'Age-Sex')
lines(as.vector(apply(apply(lowsamp$SrvIAA[,1:64,,,1,1], 2:4, sum), 2:3, mean)) / sum(as.vector(apply(apply(lowsamp$SrvIAA[,1:64,,,1,1], 2:4, sum), 2:3, mean))))
legend("topright", legend = c("Obs", "True"),
       pch = c(1, NA), lty = c(NA, 1))

# Look at survey age compositions being fitted to (blocking)
plot(as.vector(apply(model_list_lowsamp[[i]][[1]]$data$ObsSrvAgeComps[1,55:64,,,1], 2:3, mean))
     / sum(as.vector(apply(model_list_lowsamp[[i]][[1]]$data$ObsSrvAgeComps[1,55:64,,,1], 2:3, mean))),
     ylim = c(0,0.15), ylab = 'Proportion (2015:2025)', xlab = 'Age-Sex')
lines(as.vector(apply(apply(lowsamp$SrvIAA[,55:64,,,1,1], 2:4, sum), 2:3, mean)) / sum(as.vector(apply(apply(lowsamp$SrvIAA[,55:64,,,1,1], 2:4, sum), 2:3, mean))))
legend("topright", legend = c("Obs", "True"),
       pch = c(1, NA), lty = c(NA, 1))


# Sensitivity Run Results -------------------------------------------------
model_list_lowsamp <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_movesense.RDS"))
model_list_lowsamp_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_movesense_francis.RDS"))

ssb_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
rec_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
rec_om_results <- array(NA, dim = c(length(1:y), n_sims))

for(i in 1:n_sims) {
  # get aggregated OM results
  ssb_om_results[,i] <- colSums(lowsamp_movesense$SSB[,1:y,i])
  rec_om_results[,i] <- colSums(lowsamp_movesense$Rec[,1:y,i])

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




# Retrospective Results ---------------------------------------------------
retro_nofrancis <- read_csv(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_cross_test_retro.csv"))
retro_wfrancis <- read_csv(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_cross_test_retro_francis.csv"))

# Summarize info
retro_nofrancis_sum <- retro_nofrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'No Francis') %>%
  ungroup()


# Summarize info
retro_wfrancis_sum <- retro_wfrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'W Francis') %>%
  ungroup()

# get relative difference
retro_nofrancis_rd <- get_retrospective_relative_difference(retro_nofrancis_sum) %>% mutate(type = 'No Francis')
retro_francis_rd <- get_retrospective_relative_difference(retro_wfrancis_sum)%>% mutate(type = 'W Francis')

# get mohn's rho
mohns_rho_nofrancis <- retro_nofrancis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'No Francis')
mohns_rho_wfrancis <- retro_francis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'W Francis')

# bind
retro_sgl_region <- rbind(retro_wfrancis_sum, retro_nofrancis_sum)
mohns_rho <- rbind(mohns_rho_nofrancis, mohns_rho_wfrancis)
ret_df <- rbind(retro_nofrancis_rd, retro_francis_rd)

# Generally look the same
ggplot() +
  geom_hline(yintercept = 0,
                      lty = 2,
                      lwd = 1.3) +
  geom_line(
    ret_df,
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      color = as.numeric(peel)
    ),
    lwd = 1.3
  ) +
  geom_point(
    ret_df %>% dplyr::filter(peel ==
                               max(Year) - Year),
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      fill = as.numeric(peel)
    ),
    pch = 21,
    size = 6
  ) +
  geom_text(
    mohns_rho,
    mapping = aes(
      x = -Inf,
      y = Inf,
      label = paste("Mohns Rho:", round(rho, 4))
    ),
    hjust = -0.3,
    vjust = 3,
    size = 5
  )+
  guides(color = guide_colourbar(barwidth = 15, barheight = 1.3)) +
  labs(x = "Year",
                y = "Relative Difference from Terminal Year",
                color = "Retrospective Year",
                fill = "Retrospective Year") +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  facet_grid(type ~ Type, scales = "free") +
  theme_sablefish() +
  theme(legend.position = "top")

ggplot() +
  geom_line(
    retro_sgl_region %>%
      filter(peel != 0),
    mapping = aes(
      x = Year,
      y = value,
      group = peel,
      color = peel
    ),
    lwd = 1
  ) +
  geom_line(
    retro_sgl_region %>%
      filter(peel == 0),
    mapping = aes(x = Year, y = value),
    lty = 2,
    lwd = 1
  ) +
  # geom_ribbon(
  #   retro_sgl_region %>%
  #     filter(peel != 0),
  #   mapping = aes(
  #     x = Year,
  #     y = value,
  #     group = peel,
  #     ymin = lwr_95,
  #     ymax = upr_95,
  #     fill = peel
  #   ),
  #   color = NA, alpha = 0.2
  # ) +
  geom_text(mohns_rho,
            mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
            hjust = -0.3, vjust = 3, size = 5) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +

  coord_cartesian(ylim = c(0, NA)) +
  facet_grid(Type ~ type, scales = "free_y") +
  theme_sablefish() +
  labs(x = "Year", y = "Value", color = "Peel")


# Compare Retro to True Values --------------------------------------------
# Input true values into peel == 0 for evaluation
true_ssb <- colSums(lowsamp$SSB[,-66,1])
true_rec <- colSums(lowsamp$Rec[,-66,1])
true_values <- c(true_ssb, true_rec)

retro_nofrancis[retro_nofrancis$peel == 0, ]$value <- rep(true_values, length(unique(retro_nofrancis$Sim)))
retro_wfrancis[retro_wfrancis$peel == 0, ]$value <- rep(true_values, length(unique(retro_wfrancis$Sim)))

# Summarize info
retro_nofrancis_sum <- retro_nofrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'No Francis') %>%
  ungroup()


# Summarize info
retro_wfrancis_sum <- retro_wfrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'W Francis') %>%
  ungroup()

# get relative difference
retro_nofrancis_rd <- get_retrospective_relative_difference(retro_nofrancis_sum) %>% mutate(type = 'No Francis')
retro_francis_rd <- get_retrospective_relative_difference(retro_wfrancis_sum)%>% mutate(type = 'W Francis')

# get mohn's rho
mohns_rho_nofrancis <- retro_nofrancis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'No Francis')
mohns_rho_wfrancis <- retro_francis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'W Francis')

# bind
retro_sgl_region <- rbind(retro_wfrancis_sum, retro_nofrancis_sum)
mohns_rho <- rbind(mohns_rho_nofrancis, mohns_rho_wfrancis)
ret_df <- rbind(retro_nofrancis_rd, retro_francis_rd)
# Generally look the same
ggplot() +
  geom_hline(yintercept = 0,
             lty = 2,
             lwd = 1.3) +
  geom_line(
    ret_df,
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      color = as.numeric(peel)
    ),
    lwd = 1.3
  ) +
  geom_point(
    ret_df %>% dplyr::filter(peel ==
                               max(Year) - Year),
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      fill = as.numeric(peel)
    ),
    pch = 21,
    size = 6
  ) +
  geom_text(
    mohns_rho,
    mapping = aes(
      x = -Inf,
      y = Inf,
      label = paste("Mohns Rho:", round(rho, 4))
    ),
    hjust = -0.3,
    vjust = 3,
    size = 5
  )+
  guides(color = guide_colourbar(barwidth = 15, barheight = 1.3)) +
  labs(x = "Year",
       y = "Relative Difference from Terminal Year",
       color = "Retrospective Year",
       fill = "Retrospective Year") +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  facet_grid(type ~ Type, scales = "free") +
  theme_sablefish() +
  theme(legend.position = "top")

ggplot() +
  geom_line(
    retro_sgl_region %>%
      filter(peel != 0),
    mapping = aes(
      x = Year,
      y = value,
      group = peel,
      color = peel
    ),
    lwd = 1
  ) +
  geom_line(
    retro_sgl_region %>%
      filter(peel == 0),
    mapping = aes(x = Year, y = value),
    lty = 2,
    lwd = 1
  ) +
  geom_ribbon(
    retro_sgl_region %>%
      filter(peel != 0),
    mapping = aes(
      x = Year,
      y = value,
      group = peel,
      ymin = lwr_95,
      ymax = upr_95,
      fill = peel
    ),
    color = NA, alpha = 0.2
  ) +
  geom_text(mohns_rho,
            mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
            hjust = -0.3, vjust = 3, size = 5) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_grid(Type ~ type, scales = "free_y") +
  theme_sablefish() +
  labs(x = "Year", y = "Value", color = "Peel")
