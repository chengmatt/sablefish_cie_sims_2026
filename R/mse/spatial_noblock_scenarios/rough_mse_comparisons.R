# Purpose: To roughly compare MSE performance across models
# Creator: Matthew LH. Cheng
# 1/22/25


# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)

source(here("R", "functions", "mse_functions.R"))

# Read in MSEs (base recruitment)
file_dir <- here("outputs", "mse_results", "spatial_noblock_scenarios")
sgl_rg <- readRDS(here(file_dir, "single_region_base.RDS"))
faa <- readRDS(here(file_dir, "faa_base.RDS"))
five_rg <- readRDS(here(file_dir, 'five_region_base.RDS'))

# Dimensions
n_sims <- 150
st_yr <- 65
end_yr <- 95

# combine to a list
mse_list <- list(sgl_rg, faa, five_rg)

# Storage containers
ssb_results <- array(NA, dim = c(length(mse_list), dim(sgl_rg$SSB)))
rec_results <- array(NA, dim = c(length(mse_list), dim(sgl_rg$SSB)))
catch_results <- array(NA, dim = c(length(mse_list), dim(sgl_rg$TrueCatch)))

# Loop through to get results
for(i in 1:length(mse_list)) {
  ssb_results[i,,,] <- mse_list[[i]]$SSB
  rec_results[i,,,] <- mse_list[[i]]$Rec
  catch_results[i,,,,] <- mse_list[[i]]$TrueCatch
}

# Munge results
# regional ssb
ssb_rg_df <- reshape2::melt(ssb_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# global ssb
ssb_global_df <- reshape2::melt(ssb_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# regional recruitment
rec_rg_df <- reshape2::melt(rec_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# global recruitment
rec_global_df <- reshape2::melt(rec_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# regional catch
catch_rg_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, region, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# global catch
catch_global_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# # Compute reference points from five region model for comparison
# b40_mat <- array(NA, dim = c(length(mse_list), length(st_yr:end_yr), n_sims))
# for(i in 1:n_sims) {
#   for(y in st_yr:end_yr) {
#     for(j in 1:length(mse_list)) {
#       if(j %in% c(1:2) && !is.null(mse_list[[j]]$models[[i]])) {
#         # get true b40 for every sim
#         ref_pts <- get_closed_loop_reference_points(
#           use_true_values = TRUE,
#           sim_env = mse_list[[j]],
#           asmt_data = NULL,
#           asmt_rep = NULL,
#           y = y,
#           sim = i,
#
#           # single region reference points
#           reference_points_opt = list(
#             n_avg_yrs = 1,
#             SPR_x = 0.4,
#             calc_rec_st_yr = 20,
#             rec_age = 2,
#             type = 'multi_region',
#             what = "global_SPR",
#             B_x = 0.4
#           ),
#           n_proj_yrs = 2
#         )
#       }
#
#       # sum and save results
#       b40_mat[j,y-st_yr + 1,i] <- sum(ref_pts$b_ref_pt[,2])
#     }
#   }
# }

# Terminal Year Bias
ssb_bias_results <- array(NA, dim = c(2, 95, 150))

# Loop through to get results
for(i in 1:2) { # sgl, faa
  for(j in 1:150) { # 150 sims
    if((sum(mse_list[[i]]$models[[j]]$rep$SSB)) != 0) ssb_bias_results[i,,j] <- (mse_list[[i]]$models[[j]]$rep$SSB - colSums(mse_list[[i]]$SSB[,,j])) / colSums(mse_list[[i]]$SSB[,,j])
  } # end j loop
} # end i loop

reshape2::melt(ssb_bias_results) %>%
  group_by(Var2, Var1) %>%
  mutate(median = median(value, na.rm = T)) %>%
  ggplot() +
  geom_line(aes(x = Var2 + 1959, y = value, group = Var3)) +
  geom_line(aes(x = Var2 + 1959, y = median, group = Var1), col = 'green', lwd = 3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 2, col = 'red') +
  facet_wrap(~Var1, labeller = labeller(Var1 = c("1" = "sgl", "2" = "faa"))) +
  labs(x = 'Year', y = 'SSB Bias')


# Compare Time Series -------------------------------------------------------------------

# Regional SSB
ggplot() +
  geom_line(ssb_rg_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(ssb_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1.3) +
  geom_ribbon(ssb_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.25, color = NA) +
  facet_wrap(~region, nrow = 1, scales = 'free') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  ggthemes::scale_color_solarized() +
  ggthemes::scale_fill_solarized() +
  labs(x = 'Year', y = 'SSB', fill = 'Model', color = 'Model') +
  theme_sablefish()

# Global SSB
ggplot() +
  geom_line(ssb_global_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(ssb_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  geom_ribbon(ssb_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)),
              alpha = 0.25, color = NA) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  theme_bw(base_size = 18) +
  ggthemes::scale_color_solarized() +
  ggthemes::scale_fill_solarized() +
  labs(x = 'Year', y = 'SSB', fill = 'Model', color = 'Model')

# Regional Catch
ggplot() +
  geom_line(catch_rg_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(catch_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1.3) +
  geom_ribbon(catch_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.25, color = NA) +
  facet_wrap(~region, nrow = 1, scales = 'free') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  ggthemes::scale_color_solarized() +
  ggthemes::scale_fill_solarized() +
  labs(x = 'Year', y = 'Catch', fill = 'Model', color = 'Model') +
  theme_sablefish()

# Global Catch
ggplot() +
  geom_line(catch_global_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(catch_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  geom_ribbon(catch_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)),
              alpha = 0.25, color = NA) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  theme_bw(base_size = 18) +
  ggthemes::scale_color_solarized() +
  ggthemes::scale_fill_solarized() +
  labs(x = 'Year', y = 'Catch', fill = 'Model', color = 'Model')

# Regional Recruitment
ggplot() +
  geom_line(rec_rg_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(rec_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1.3) +
  geom_ribbon(rec_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.25, color = NA) +
  facet_wrap(~region, nrow = 1, scales = 'free') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  ggthemes::scale_color_solarized() +
  ggthemes::scale_fill_solarized() +
  labs(x = 'Year', y = 'Recruitment', fill = 'Model', color = 'Model') +
  theme_sablefish()


# Compare Apportionment ---------------------------------------------------

# Figure out apportionment for each model
get_apportionment <- function(sim_env, lls_design_type, y, sim, rolling_avg_yrs = 5) {
  # get survey apportionment (survey biomass)
  yr_lag <- rolling_avg_yrs - 1 # how many years to rollying avg across
  tmp_Srv_BiomIA <- sim_env$SrvIAA[,((y - yr_lag): y),,,1,sim] * sim_env$WAA_srv[,((y - yr_lag): y),,,1,sim]
  tmp_srvidx_biom <- apply(tmp_Srv_BiomIA, 1:2, sum) # get biomass
  apportionment <- get_single_region_survey_apportionment(feedback_start_yr = sim_env$feedback_start_yr,
                                                          n_yrs = sim_env$n_yrs,
                                                          y = y,
                                                          srv_idx = tmp_srvidx_biom, # using true survey biomass values to do apportionment
                                                          rolling_avg_yrs = rolling_avg_yrs,
                                                          lls_design_type = lls_design_type)
  return(apportionment)
}

# loop through to get apportionment
apportionment <- array(NA, dim = c(3, 5, length(st_yr:end_yr), n_sims))
for(y in st_yr:end_yr) {
  for(sim in 1:n_sims) {
    apportionment[1,,y - st_yr + 1,sim] <- get_apportionment(sim_env = sgl_rg, 'current', y = y, sim = sim)
    apportionment[2,,y - st_yr + 1,sim] <- get_apportionment(sim_env = faa, 'current', y = y, sim = sim)
    apportionment[3,,y - st_yr + 1,sim] <- rowSums(five_rg$TrueCatch[,y,,sim]) / sum(five_rg$TrueCatch[,y,,sim])
  }
}

# Plot apportionment
apportionment_df <- reshape2::melt(apportionment) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'five_rg'
  )) %>%
  filter(year != 1)

# summarize apportionment
approtionment_sum <- apportionment_df %>%
  group_by(region, year, model) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T))

ggplot(approtionment_sum, aes(x = year + 65, y = median,
                              ymin = lwr, ymax = upr,
                              color = factor(model), fill = factor(model))) +
  geom_line(lwd = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  facet_wrap(~region, nrow = 1) +
  theme_bw(base_size = 15) +
  ggthemes::scale_color_solarized() +
  ggthemes::scale_fill_solarized() +
  labs(x = 'Year', y = 'Apportionment', fill = 'Model', color = 'Model')


