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
      model == 3 ~ 'spat'
    ),
    model = factor(model, levels = c('sgl', 'faa', 'spat')),
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
    theme_bw(base_size = 18)
}


bias_spt_ts_plot <- function(df, title) {
  ggplot(df, aes(x = year, y = median, color = cl_yr, fill = cl_yr, group = cl_yr)) +
    geom_line(lwd = 1) +
    facet_wrap(~Region) +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(x = 'Year', y = 'Relative Error', title = title) +
    theme_bw(base_size = 18)
}

# Rotated densities of the terminal year distribution, returned as a narrow companion
# panel to sit immediately right of a time series plot. df is simulation level
# (model, year, sim, value) and y limits are lifted off main_plot so the strip shares
# its scale exactly. Densities are trimmed to trim_q (matching the ribbon interval, so
# tail sims do not stretch the axis) and share a common scaling, so relative peak
# heights are comparable across models.
terminal_density_panel <- function(df, term_yr, main_plot, alpha = 0.35,
                                   trim_q = c(0.025, 0.975), base_size = 18) {

  dens_df <- df %>%
    filter(year == term_yr) %>%
    group_by(model) %>%
    reframe({
      d <- density(value, na.rm = TRUE,
                   from = quantile(value, trim_q[1], na.rm = TRUE),
                   to   = quantile(value, trim_q[2], na.rm = TRUE))
      tibble(value = d$x, dens = d$y)
    }) %>%
    mutate(dens = dens / max(dens))

  ggplot(dens_df) +
    geom_ribbon(mapping = aes(y = value, xmin = 0, xmax = dens, fill = factor(model)),
                orientation = 'y', alpha = alpha, color = NA) +
    geom_path(mapping = aes(y = value, x = dens, color = factor(model)), lwd = 0.85) +
    coord_cartesian(ylim = ggplot_build(main_plot)$layout$panel_params[[1]]$y.range,
                    xlim = c(0, 1.05), expand = FALSE) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    theme_bw(base_size = base_size) +
    theme(legend.position = 'none',
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0)) # sit tight against the time series
}

# Setup -------------------------------------------------------------------
library(here)
library(tidyverse)
library(SPoRC)
library(gganimate)
library(patchwork)

source(here("R", "functions", "mse_functions.R"))
file_dir <- here("outputs", "mse_results", "apportionment")
fig_dir  <- here("outputs", "mse_results", "apportionment", "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Read in results ---------------------------------------------------------
carry <- readRDS(here(file_dir, 'carryover.RDS'))
fixed <- readRDS(here(file_dir, '2024_values.RDS'))
eq <- readRDS(here(file_dir, 'equal.RDS'))

## Process MSE results ---------------------------------------------------------
n_sims <- 100
st_yr <- 65
end_yr <- 95
n_regions <- 5
n_fish_fleets <- 2

mse_list <- list(carry, fixed, eq)

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
    model == 1 ~ 'carryover',
    model == 2 ~ 'fixed',
    model == 3 ~ 'equal'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('carryover', 'fixed', 'equal')),
  year = year + 1960)

ssb_global_sim_df <- reshape2::melt(ssb_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value), .groups = 'drop') %>%
  mutate(model = case_when(
    model == 1 ~ 'carryover',
    model == 2 ~ 'fixed',
    model == 3 ~ 'equal'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('carryover', 'fixed', 'equal')),
  year = year + 1960)

ssb_global_df <- ssb_global_sim_df %>%
  group_by(model, period, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T), .groups = 'drop')

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
    model == 1 ~ 'carryover',
    model == 2 ~ 'fixed',
    model == 3 ~ 'equal'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('carryover', 'fixed', 'equal')),
  year = year + 1960)

catch_global_sim_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value), .groups = 'drop') %>%
  mutate(model = case_when(
    model == 1 ~ 'carryover',
    model == 2 ~ 'fixed',
    model == 3 ~ 'equal'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('carryover', 'fixed', 'equal')),
  year = year + 1960)

catch_global_df <- catch_global_sim_df %>%
  group_by(model, period, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T), .groups = 'drop')

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
    model == 1 ~ 'carryover',
    model == 2 ~ 'fixed',
    model == 3 ~ 'equal'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'),
  model = factor(model, levels = c('carryover', 'fixed', 'equal')),
  year = year + 1960)


##### SSB ---------------------------------------------------------------------
ssb_local <- ggplot() +
  geom_ribbon(ssb_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(ssb_rg_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(ssb_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  facet_wrap(~region, ncol = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'SSB (kt)', fill = 'Apportionment Type', color = 'Apportionment Type') +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none')
ggsave(file.path(fig_dir, 'ssb_regional_base.png'), ssb_local, width = 14, height = 14, dpi = 300)

ssb_global <- ggplot() +
  geom_ribbon(ssb_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(ssb_global_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(ssb_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.19, 0.9),
        legend.background = element_blank()) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'SSB (kt)', fill = 'Apportionment Type', color = 'Apportionment Type')
ggsave(file.path(fig_dir, 'ssb_global_base.png'), ssb_global, width = 8, height = 8, dpi = 300)

# recruitment plot combined
ssb_combined <- cowplot::plot_grid(ssb_global, ssb_local, rel_widths = c(0.75, 0.25), labels = c('A', 'B'), label_size = 20)
ggsave(file.path(fig_dir, 'ssb_combined_base.png'), ssb_combined, width = 15, height = 10, dpi = 300)


##### Catch -------------------------------------------------------------------
catch_local <- ggplot() +
  geom_ribbon(catch_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(catch_rg_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(catch_rg_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  facet_wrap(~region, ncol = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Catch (kt)', fill = 'Apportionment Type', color = 'Apportionment Type') +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none')
ggsave(file.path(fig_dir, 'catch_regional_base.png'), catch_local, width = 14, height = 14, dpi = 300)

catch_global <- ggplot() +
  geom_ribbon(catch_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(catch_global_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(catch_global_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 2025, lty = 2) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.19, 0.9),
        legend.background = element_blank()) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Catch (kt)', fill = 'Apportionment Type', color = 'Apportionment Type')
ggsave(file.path(fig_dir, 'catch_global_base.png'), catch_global, width = 8, height = 8, dpi = 300)

# recruitment plot combined
catch_combined <- cowplot::plot_grid(catch_global, catch_local, rel_widths = c(0.75, 0.25), labels = c('A', 'B'), label_size = 20)
ggsave(file.path(fig_dir, 'catch_combined_base.png'), catch_combined, width = 15, height = 10, dpi = 300)

##### Apportionment -----------------------------------------------------------
ggplot() +
  geom_ribbon(apport_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.15, color = NA) +
  geom_line(apport_df %>% filter(year <= st_yr + 1960), mapping = aes(x = year, y = median), color = 'black', lwd = 1) +
  geom_line(apport_df %>% filter(year >= st_yr + 1960), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  facet_wrap(~region, nrow = 1) +
  geom_vline(xintercept = 2025, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Apportionment', fill = 'Apportionment Type', color = 'Apportionment Type') +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top')
ggsave(file.path(fig_dir, 'apportionment_base.png'), width = 14, height = 7, dpi = 300)
