# Purpose: To roughly compare MSE performance across models
# Creator: Matthew LH. Cheng
# 1/22/25

# Helper functions --------------------------------------------------------
combine_results <- function(all_results) {
  setNames(lapply(names(all_results[[1]]), function(nm) {
    slices <- lapply(all_results, `[[`, nm)
    if (nm == "models") {
      n_cl_yrs <- length(slices[[1]])
      lapply(1:n_cl_yrs, function(cl) {
        do.call(c, lapply(slices, `[[`, cl))
      })
    } else if (is.array(slices[[1]])) {
      abind::abind(slices, along = length(dim(slices[[1]])))
    } else {
      do.call(c, slices)
    }
  }), names(all_results[[1]]))
}

summarize_bias_ts <- function(arr) {
  reshape2::melt(arr) %>%
    rename(model = Var1, cl_yr = Var2, year = Var3, sim = Var4) %>%
    group_by(model, cl_yr, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr = quantile(value, 0.025, na.rm = TRUE),
              upr = quantile(value, 0.975, na.rm = TRUE)) %>%
    mutate(model = case_when(
      model == 1 ~ 'sgl',
      model == 2 ~ 'faa',
      model == 3 ~ 'three_rg'
    ),
    model = factor(model, levels = c('sgl', 'faa', 'three_rg')),
    cl_yr = cl_yr + st_yr + 1959,
    year  = year + 1960)
}

summarize_bias_ts_spt <- function(arr) {
  reshape2::melt(arr) %>%
    rename(Region = Var1, cl_yr = Var2, year = Var3, sim = Var4) %>%
    group_by(Region, cl_yr, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr = quantile(value, 0.025, na.rm = TRUE),
              upr = quantile(value, 0.975, na.rm = TRUE)) %>%
    mutate(Region = case_when(
      Region == 1 ~ 'BS',
      Region == 2 ~ 'AI',
      Region == 3 ~ 'GOA'
    ),
    Region = factor(Region, levels = c('BS', 'AI', 'GOA')),
    cl_yr = cl_yr + st_yr + 1959,
    year  = year + 1960)
}

bias_ts_plot <- function(df, title) {
  ggplot(df, aes(x = year, y = median, color = cl_yr, fill = cl_yr, group = cl_yr)) +
    geom_line(lwd = 1) +
    facet_wrap(~model) +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(x = 'Year', y = 'Relative Error', title = title) +
    theme_bw(base_size = 13)
}


bias_spt_ts_plot <- function(df, title) {
  ggplot(df, aes(x = year, y = median, color = cl_yr, fill = cl_yr, group = cl_yr)) +
    geom_line(lwd = 1) +
    facet_wrap(~Region) +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(x = 'Year', y = 'Relative Error', title = title) +
    theme_bw(base_size = 13)
}

# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)
library(gganimate)
library(patchwork)

source(here("R", "functions", "mse_functions.R"))
file_dir <- here("outputs", "mse_results", "spatial_noblock_scenarios")
fig_dir  <- here("outputs", "mse_results", "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Read in results ---------------------------------------------------------
files <- list.files(file_dir, full.names = TRUE)
sgl_rg_files <- files[stringr::str_detect(files, "single_region_base")]
all_sgl_rg_results <- lapply(sgl_rg_files, readRDS)
sgl_rg <- combine_results(all_sgl_rg_results)

files <- list.files(file_dir, full.names = TRUE)
faa_files <- files[stringr::str_detect(files, "faa_base")]
all_faa_results <- lapply(faa_files, readRDS)
faa <- combine_results(all_faa_results)

files <- list.files(file_dir, full.names = TRUE)
three_rg_files <- files[stringr::str_detect(files, "three_region_base")]
all_three_rg_results <- lapply(three_rg_files, readRDS)
three_rg <- combine_results(all_three_rg_results)

five_rg <- readRDS(here(file_dir, 'five_region_base.RDS'))

## Process MSE results ---------------------------------------------------------
n_sims <- 100
st_yr <- 65
end_yr <- 95
n_regions <- 5
n_fish_fleets <- 2

mse_list <- list(sgl_rg, faa, three_rg, five_rg)

ssb_results   <- array(NA, dim = c(length(mse_list), n_regions, end_yr, n_sims))
rec_results   <- array(NA, dim = c(length(mse_list), n_regions, end_yr, n_sims))
catch_results <- array(NA, dim = c(length(mse_list), n_regions, end_yr, n_fish_fleets, n_sims))
keep_idx      <- vector("list", length(mse_list))

for (i in 1:length(mse_list)) {
  ssb_all       <- mse_list[[i]]$SSB
  rec_all       <- mse_list[[i]]$Rec
  truecatch_all <- mse_list[[i]]$TrueCatch
  total_sims    <- dim(ssb_all)[3]
  converged     <- which(apply(ssb_all, 3, sum) != 0)
  keep          <- converged[1:n_sims]
  keep_idx[[i]] <- converged[1:n_sims]
  ssb_results[i,,,]    <- ssb_all[,,keep]
  rec_results[i,,,]    <- rec_all[,,keep]
  catch_results[i,,,,] <- truecatch_all[,,,keep]
}

### Munge MSE Results -------------------------------------------------------

ssb_rg_df <- reshape2::melt(ssb_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
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
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg', 'five_rg')),
  year = year + 1960)

ssb_global_df <- reshape2::melt(ssb_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg', 'five_rg')),
  year = year + 1960)

rec_rg_df <- reshape2::melt(rec_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
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
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg', 'five_rg')),
  year = year + 1960)

rec_global_df <- reshape2::melt(rec_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg', 'five_rg')),
  year = year + 1960)

catch_rg_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, region, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
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
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg', 'five_rg')),
  year = year + 1960)

catch_global_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg', 'five_rg')),
  year = year + 1960)

apport_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, region, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year, sim) %>%
  mutate(value = value / sum(value)) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
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
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg', 'five_rg')),
  year = year + 1960)

#### Compare MSE Time Series -------------------------------------------------------------------

##### SSB ---------------------------------------------------------------------
ggplot() +
  geom_ribbon(ssb_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(ssb_rg_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(ssb_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  facet_wrap(~region, nrow = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'SSB', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, 'ssb_regional_base.png'), width = 14, height = 7, dpi = 300)

ggplot() +
  geom_ribbon(ssb_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(ssb_global_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(ssb_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  theme_bw(base_size = 18) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'SSB', fill = 'Model', color = 'Model')
ggsave(file.path(fig_dir, 'ssb_global_base.png'), width = 8, height = 5, dpi = 300)

##### Catch -------------------------------------------------------------------
ggplot() +
  geom_ribbon(catch_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(catch_rg_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(catch_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  facet_wrap(~region, nrow = 1, scales = 'free') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Catch', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, 'catch_regional_base.png'), width = 14, height = 7, dpi = 300)

ggplot() +
  geom_ribbon(catch_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(catch_global_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(catch_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  theme_bw(base_size = 18) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Catch', fill = 'Model', color = 'Model')
ggsave(file.path(fig_dir, 'catch_global_base.png'), width = 8, height = 5, dpi = 300)

##### Recruitment -------------------------------------------------------------
ggplot() +
  geom_ribbon(rec_rg_df %>% filter(year >= st_yr + 1960, model == 'sgl'), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(rec_rg_df %>% filter(year <= st_yr + 1960, model == 'sgl'), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(rec_rg_df %>% filter(year >= st_yr + 1960, model == 'sgl'), mapping = aes(x = year, y = median), lwd = 1) +
  facet_wrap(~region, nrow = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Recruitment', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, 'rec_regional_base.png'), width = 18, height = 7, dpi = 300)

ggplot() +
  geom_ribbon(rec_global_df %>% filter(year >= st_yr + 1960, model == 'sgl'), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(rec_global_df %>% filter(year <= st_yr + 1960, model == 'sgl'), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(rec_global_df %>% filter(year >= st_yr + 1960, model == 'sgl'), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  theme_bw(base_size = 18) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = 'none') +
  labs(x = 'Year', y = 'Recruitment', fill = 'Model', color = 'Model')
ggsave(file.path(fig_dir, 'rec_global_base.png'), width = 8, height = 5, dpi = 300)

##### Apportionment -----------------------------------------------------------
ggplot() +
  geom_ribbon(apport_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(apport_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(apport_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  facet_wrap(~region, nrow = 1) +
  geom_vline(xintercept = 2025, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Apportionment', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, 'apportionment_base.png'), width = 14, height = 7, dpi = 300)

## Process Estimation Bias Results ------------------------------------------------
n_cl_yrs <- end_yr - st_yr + 1

dep_est_results      <- ssb_est_results      <- totbiom_est_results <- rec_est_results <- array(NA, dim = c(3, n_cl_yrs, end_yr, n_sims))
dep_true_results      <- ssb_true_results     <- totbiom_true_results <- rec_true_results <- array(NA, dim = c(3, end_yr, n_sims))
q_est_results        <- r0_est_results <- array(NA, dim = c(3, n_cl_yrs, n_sims))
q_true_results       <- r0_true_results <- array(NA, dim = c(3, n_sims))
dep_bias_results      <- ssb_bias_results     <- totbiom_bias_results <- rec_bias_results <- array(NA, dim = c(3, n_cl_yrs, end_yr, n_sims))
q_bias_results       <- r0_bias_results <- array(NA, dim = c(3, n_cl_yrs, n_sims))
spt_rec_bias_results     <- spt_rec_est_results <- spt_ssb_bias_results     <- spt_ssb_est_results      <- array(NA, dim = c(3, n_cl_yrs, end_yr, n_sims))
spt_rec_true_results <- spt_ssb_true_results     <- array(NA, dim = c(3, end_yr, n_sims))


for (i in 1:3) {

  # extract quants
  keep          <- keep_idx[[i]]
  ssb_true      <- apply(mse_list[[i]]$SSB[,,keep], 3, colSums)
  dep_true      <- ssb_true / ssb_true[1,]
  totbiom_true  <- apply(mse_list[[i]]$Total_Biom[,,keep], 3, colSums)
  rec_true      <- apply(mse_list[[i]]$Rec[,,keep], 3, colSums)
  q_true        <- mse_list[[i]]$srv_q[1,1,1,keep]
  r0_true       <- colSums(mse_list[[i]]$R0[,1,keep])

  # store true quants
  ssb_true_results[i,,]     <- ssb_true
  totbiom_true_results[i,,] <- totbiom_true
  rec_true_results[i,,]     <- rec_true
  q_true_results[i,]        <- q_true
  r0_true_results[i,]       <- r0_true
  dep_true_results[i,,]       <- dep_true

  # if spatial model
  if(i == 3) {
    # ssb
    spt_ssb_true <- mse_list[[i]]$SSB[,,keep]
    spt_ssb_true_results[1:2,,] <- mse_list[[i]]$SSB[1:2,,keep] # bs and ai
    spt_ssb_true_results[3,,] <- apply(mse_list[[i]]$SSB[3:5,,keep], c(2,3), sum) # goa

    # rec
    spt_rec_true <- mse_list[[i]]$Rec[,,keep]
    spt_rec_true_results[1:2,,] <- mse_list[[i]]$Rec[1:2,,keep] # bs and ai
    spt_rec_true_results[3,,] <- apply(mse_list[[i]]$Rec[3:5,,keep], c(2,3), sum) # goa
  }

  for (cl in 1:n_cl_yrs) {

    # get estimates
    yr_idx      <- 1:(st_yr + cl - 1)
    ssb_est     <- sapply(1:length(keep), \(s) colSums(mse_list[[i]]$models[[cl]][[keep[s]]]$rep$SSB))
    dep_est     <- ssb_est / ssb_est[1,]
    totbiom_est <- sapply(1:length(keep), \(s) colSums(mse_list[[i]]$models[[cl]][[keep[s]]]$rep$Total_Biom))
    rec_est     <- sapply(1:length(keep), \(s) colSums(mse_list[[i]]$models[[cl]][[keep[s]]]$rep$Rec))
    if (i == 2) q_est <- NA # faa model
    else q_est  <- sapply(1:length(keep), \(s) mse_list[[i]]$models[[cl]][[keep[s]]]$rep$srv_q[1,1,1])
    r0_est      <- sapply(1:length(keep), \(s) mse_list[[i]]$models[[cl]][[keep[s]]]$rep$R0)

    # store and summarize results
    ssb_est_results[i, cl, yr_idx,]      <- ssb_est
    dep_est_results[i, cl, yr_idx,]      <- dep_est
    totbiom_est_results[i, cl, yr_idx,]  <- totbiom_est
    rec_est_results[i, cl, yr_idx,]      <- rec_est
    q_est_results[i, cl,]                <- q_est
    r0_est_results[i, cl,]               <- r0_est
    ssb_bias_results[i, cl, yr_idx,]     <- (ssb_est - ssb_true[yr_idx,]) / ssb_true[yr_idx,]
    dep_bias_results[i, cl, yr_idx,]     <- (dep_est - dep_true[yr_idx,]) / dep_true[yr_idx,]
    totbiom_bias_results[i, cl, yr_idx,] <- (totbiom_est - totbiom_true[yr_idx,]) / totbiom_true[yr_idx,]
    rec_bias_results[i, cl, yr_idx,]     <- (rec_est - rec_true[yr_idx,]) / rec_true[yr_idx,]
    q_bias_results[i, cl,]               <- (q_est - q_true) / q_true
    r0_bias_results[i, cl,]              <- (r0_est - r0_true) / r0_true

    # spatial model regional ssb results
    if(i == 3) {

      # get spatial ssb estimates
      spt_ssb_est <- abind::abind(
        lapply(1:length(keep), \(s) mse_list[[i]]$models[[cl]][[keep[s]]]$rep$SSB),
        along = length(dim(mse_list[[i]]$models[[cl]][[keep[1]]]$rep$SSB)) + 1
      )
      spt_ssb_est_results[, cl, yr_idx,]      <- spt_ssb_est
      spt_ssb_bias_results[, cl, yr_idx,]     <- (spt_ssb_est - spt_ssb_true_results[,yr_idx,]) / spt_ssb_true_results[,yr_idx,]

      # get spatial rec estimates
      spt_rec_est <- abind::abind(
        lapply(1:length(keep), \(s) mse_list[[i]]$models[[cl]][[keep[s]]]$rep$Rec),
        along = length(dim(mse_list[[i]]$models[[cl]][[keep[1]]]$rep$Rec)) + 1
      )
      spt_rec_est_results[, cl, yr_idx,]      <- spt_rec_est
      spt_rec_bias_results[, cl, yr_idx,]     <- (spt_rec_est - spt_rec_true_results[,yr_idx,]) / spt_rec_true_results[,yr_idx,]

    }
  }
}

ssb_bias_df     <- summarize_bias_ts(ssb_bias_results)
dep_bias_df     <- summarize_bias_ts(dep_bias_results)
totbiom_bias_df <- summarize_bias_ts(totbiom_bias_results)
rec_bias_df     <- summarize_bias_ts(rec_bias_results)
spt_ssb_bias_df     <- summarize_bias_ts_spt(spt_ssb_bias_results)
spt_rec_bias_df     <- summarize_bias_ts_spt(spt_rec_bias_results)

scalar_bias_df <- rbind(
  reshape2::melt(q_bias_results)  %>% mutate(metric = 'q'),
  reshape2::melt(r0_bias_results) %>% mutate(metric = 'R0')
) %>%
  rename(model = Var1, cl_yr = Var2, sim = Var3) %>%
  group_by(model, cl_yr, metric) %>%
  summarize(median = median(value, na.rm = TRUE),
            lwr = quantile(value, 0.025, na.rm = TRUE),
            upr = quantile(value, 0.975, na.rm = TRUE)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg'
  ),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg')),
  cl_yr = cl_yr + st_yr + 1959)

### EM Time Series Bias -----------------------------------------------------
bias_ts_plot(ssb_bias_df, 'SSB Bias')
ggsave(file.path(fig_dir, 'ssb_bias_base.png'), width = 10, height = 4, dpi = 300)

bias_ts_plot(dep_bias_df, 'Depletion Bias')
ggsave(file.path(fig_dir, 'dep_bias_base.png'), width = 10, height = 4, dpi = 300)

bias_ts_plot(totbiom_bias_df, 'Total Biomass Bias')
ggsave(file.path(fig_dir, 'totbiom_bias_base.png'), width = 10, height = 4, dpi = 300)

bias_ts_plot(rec_bias_df, 'Recruitment Bias')
ggsave(file.path(fig_dir, 'rec_bias_base.png'), width = 10, height = 4, dpi = 300)

bias_spt_ts_plot(spt_ssb_bias_df, "Spatial SSB Bias")
ggsave(file.path(fig_dir, 'spt_ssb_bias_base.png'), width = 10, height = 4, dpi = 300)

bias_spt_ts_plot(spt_rec_bias_df, "Spatial Recruitment Bias")
ggsave(file.path(fig_dir, 'spt_rec_bias_base.png'), width = 10, height = 4, dpi = 300)


### EM Scalar Bias ----------------------------------------------------------
ggplot(scalar_bias_df, aes(x = cl_yr, y = median)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line(lwd = 1) +
  facet_grid(metric ~ model) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Year', y = 'Relative Error', title = 'R0 and Survey q Bias') +
  theme_bw(base_size = 13)

# Q and R0 bias from sgl model, likely due to carry over values, manifesting as time-varying q

ggsave(file.path(fig_dir, 'scalar_bias_base.png'), width = 10, height = 6, dpi = 300)

### Animate EM Time Series Bias (Aggregated) ---------------------------------------------
build_tracking_df <- function(est_arr, true_arr, bias_arr, metric_label) {
  est_df <- reshape2::melt(est_arr) %>%
    rename(model = Var1, cl_yr = Var2, year = Var3, sim = Var4) %>%
    group_by(model, cl_yr, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr = quantile(value, 0.025, na.rm = TRUE),
              upr = quantile(value, 0.975, na.rm = TRUE)) %>%
    mutate(type = 'Estimate')
  true_df <- reshape2::melt(true_arr) %>%
    rename(model = Var1, year = Var2, sim = Var3) %>%
    group_by(model, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr = quantile(value, 0.025, na.rm = TRUE),
              upr = quantile(value, 0.975, na.rm = TRUE)) %>%
    crossing(cl_yr = 1:n_cl_yrs) %>%
    mutate(type = 'Truth')
  bias_df <- reshape2::melt(bias_arr) %>%
    rename(model = Var1, cl_yr = Var2, year = Var3, sim = Var4) %>%
    group_by(model, cl_yr, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr = quantile(value, 0.025, na.rm = TRUE),
              upr = quantile(value, 0.975, na.rm = TRUE)) %>%
    mutate(type = 'Estimate')
  bind_rows(
    rbind(est_df, true_df) %>% mutate(panel = metric_label),
    bias_df %>% mutate(panel = paste0(metric_label, ' RE'))
  )
}

panel_levels <- c('SSB', 'SSB RE', 'Recruitment', 'Recruitment RE', 'Total Biomass', 'Total Biomass RE')

combined_anim_df <- bind_rows(
  build_tracking_df(ssb_est_results,     ssb_true_results,     ssb_bias_results,     'SSB'),
  build_tracking_df(rec_est_results,     rec_true_results,     rec_bias_results,     'Recruitment'),
  build_tracking_df(totbiom_est_results, totbiom_true_results, totbiom_bias_results, 'Total Biomass')
) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg'
  ),
  model = factor(model, levels = c('sgl', 'faa', 'three_rg')),
  year  = year + 1960,
  cl_yr = cl_yr + st_yr + 1959,
  panel = factor(panel, levels = panel_levels)) %>%
  filter(type == 'Truth' | year <= cl_yr)

hline_df <- combined_anim_df %>%
  filter(grepl('RE', panel)) %>%
  distinct(panel, model, cl_yr)

p <- ggplot(combined_anim_df, aes(x = year, y = median, color = type, fill = type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(lwd = 0.8) +
  geom_vline(aes(xintercept = cl_yr), linetype = 'dashed', color = 'grey40') +
  geom_hline(data = hline_df, aes(yintercept = 0), linetype = 'dashed', color = 'black') +
  ggh4x::facet_grid2(panel ~ model, scales = 'free', independent = 'y') +
  scale_color_manual(values = c('Truth' = 'black', 'Estimate' = 'steelblue')) +
  scale_fill_manual(values  = c('Truth' = 'black', 'Estimate' = 'steelblue')) +
  labs(x = 'Year', y = NULL, title = 'Assessment Year: {frame_time}', color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  transition_time(cl_yr) +
  ease_aes('linear')

animate(p, nframes = n_cl_yrs, fps = 4, width = 1000, height = 1000)
anim_save(file.path(fig_dir, 'all_tracking_base.gif'))

### Animate Spatial Model Regional SSB/Rec Bias --------------------------------

build_spatial_tracking_df <- function(est_arr, true_arr, bias_arr, metric_label) {
  est_df <- reshape2::melt(est_arr) %>%
    rename(region = Var1, cl_yr = Var2, year = Var3, sim = Var4) %>%
    group_by(region, cl_yr, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr    = quantile(value, 0.025, na.rm = TRUE),
              upr    = quantile(value, 0.975, na.rm = TRUE),
              .groups = 'drop') %>%
    mutate(type = 'Estimate')

  true_df <- reshape2::melt(true_arr) %>%
    rename(region = Var1, year = Var2, sim = Var3) %>%
    group_by(region, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr    = quantile(value, 0.025, na.rm = TRUE),
              upr    = quantile(value, 0.975, na.rm = TRUE),
              .groups = 'drop') %>%
    crossing(cl_yr = 1:n_cl_yrs) %>%
    mutate(type = 'Truth')

  bias_df <- reshape2::melt(bias_arr) %>%
    rename(region = Var1, cl_yr = Var2, year = Var3, sim = Var4) %>%
    group_by(region, cl_yr, year) %>%
    summarize(median = median(value, na.rm = TRUE),
              lwr    = quantile(value, 0.025, na.rm = TRUE),
              upr    = quantile(value, 0.975, na.rm = TRUE),
              .groups = 'drop') %>%
    mutate(type = 'Estimate')

  bind_rows(
    rbind(est_df, true_df) %>% mutate(panel = metric_label),
    bias_df               %>% mutate(panel = paste0(metric_label, ' RE'))
  )
}

spt_panel_levels <- c('SSB', 'SSB RE', 'Recruitment', 'Recruitment RE')

spt_combined_anim_df <- bind_rows(
  build_spatial_tracking_df(spt_ssb_est_results, spt_ssb_true_results, spt_ssb_bias_results, 'SSB'),
  build_spatial_tracking_df(spt_rec_est_results, spt_rec_true_results, spt_rec_bias_results, 'Recruitment')
) %>%
  mutate(
    region = case_when(
      region == 1 ~ 'BS',
      region == 2 ~ 'AI',
      region == 3 ~ 'GOA'
    ),
    region = factor(region, levels = c('BS', 'AI', 'GOA')),
    year   = year  + 1960,
    cl_yr  = cl_yr + st_yr + 1959,
    panel  = factor(panel, levels = spt_panel_levels)
  ) %>%
  filter(type == 'Truth' | year <= cl_yr)

spt_hline_df <- spt_combined_anim_df %>%
  filter(grepl('RE', panel)) %>%
  distinct(panel, region, cl_yr)

p_spt <- ggplot(spt_combined_anim_df,
                aes(x = year, y = median, color = type, fill = type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(lwd = 0.8) +
  geom_vline(aes(xintercept = cl_yr), linetype = 'dashed', color = 'grey40') +
  geom_hline(data = spt_hline_df, aes(yintercept = 0),
             linetype = 'dashed', color = 'black') +
  ggh4x::facet_grid2(panel ~ region, scales = 'free', independent = 'y') +
  scale_color_manual(values = c('Truth' = 'black', 'Estimate' = 'steelblue')) +
  scale_fill_manual(values  = c('Truth' = 'black', 'Estimate' = 'steelblue')) +
  labs(x = 'Year', y = NULL,
       title = 'Spatial Model (three_rg) — Assessment Year: {frame_time}',
       color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  transition_time(cl_yr) +
  ease_aes('linear')

animate(p_spt, nframes = n_cl_yrs, fps = 4, width = 900, height = 700)
anim_save(file.path(fig_dir, 'spt_regional_tracking_base.gif'))


# Publication Plots -------------------------------------------------------
ssb_global_df <- ssb_global_df %>% mutate(model = factor(model, levels = c("five_rg", 'sgl', 'faa', 'three_rg')))
catch_global_df <- catch_global_df %>% mutate(model = factor(model, levels = c("five_rg", 'sgl', 'faa', 'three_rg')))

# SSB
ssb_base_plot <- ggplot() +
  geom_ribbon(ssb_global_df %>% filter(year >= st_yr + 1960),
              mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model),
                            fill = factor(model)), alpha = 0.1, color = NA) +
  geom_line(ssb_global_df %>% filter(year >= 60 + 1960, year <= st_yr + 1960), mapping = aes(x = year, y = median,
                                                                                             lty = factor(model)), color = 'black', lwd = 1.3) +
  geom_line(ssb_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median,
                                                                          color = factor(model), lty = factor(model)), lwd = 1.3) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.075, 0.85), legend.background = element_blank()) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'SSB (kt)', fill = 'Model', color = 'Model', lty = 'Model')


# Catch
catch_base_plot <- ggplot() +
  geom_ribbon(catch_global_df %>% filter(year >= st_yr + 1960),
              mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model),
                            fill = factor(model)), alpha = 0.1, color = NA) +
  geom_line(catch_global_df %>% filter(year >= 60 + 1960, year <= st_yr + 1960), mapping = aes(x = year, y = median,
                                                                                               lty = factor(model)), color = 'black', lwd = 1.3) +
  geom_line(catch_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median,
                                                                            color = factor(model), lty = factor(model)), lwd = 1.3) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.075, 0.85), legend.background = element_blank()) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Catch (kt)', fill = 'Model', color = 'Model', lty = 'Model')

# SSB Bias
cl_bias_plot <- bias_ts_plot(df = ssb_bias_df, title = 'SSB Bias') +
  geom_line(lwd = 1, alpha = 0.85) +
  geom_point(ssb_bias_df %>% filter(cl_yr == year),
             mapping = aes(x = year, y = median, fill = cl_yr),
             pch = 21, size = 4, color = 'black', alpha = 0.85) +
  labs(x = 'Year', y = 'Relative Error in SSB', title = '', color = '', fill = '') +
  scale_color_viridis_c(option = 'magma', direction = -1) +
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  facet_grid(model~.) +
  theme_bw(base_size = 18) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(legend.position = c(0.875, 0.08),
        legend.background = element_blank())

# Align all 3 individual plots first
aligned <- cowplot::align_plots(
  ssb_base_plot,
  catch_base_plot,
  cl_bias_plot,
  align = 'hv',
  axis = 'tblr'
)

# Plot
left_col <- cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 1,
                               labels = c('A', 'B'), label_size = 30, hjust = -0.55, vjust =  2.5)
cowplot::plot_grid(left_col, aligned[[3]],
                   rel_widths = c(0.7, 0.3), labels = c('', 'C'),
                   label_size = 30, hjust = -0.55, vjust =  2.5)

ggsave(file.path(fig_dir, 'mse_base.png'), width = 18, height = 15, dpi = 300)

