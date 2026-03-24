# Purpose: To compare the simulation results of a single-region, FAA, and three-region model under historical survey design.
# Creator: Matthew LH. Cheng
# Date 12/24/25

# Helper functions --------------------------------------------------------
conv_check <- function(sd_rep) {
  if(sd_rep$pdHess == T &  max(abs(sd_rep$gradient.fixed)) <= 0.01) TRUE
  else FALSE
}

melt_om <- function(arr, is_spatial = FALSE) {
  df <- reshape2::melt(arr)
  if (is_spatial)
    rename(df, Region = Var1, Year = Var2, Sim = Var3, OM = value)
  else
    rename(df, Year = Var1, Sim = Var2, OM = value)
}

melt_em <- function(store, om_df, is_spatial = FALSE, is_spt_store = FALSE) {
  df <- reshape2::melt(store)
  if (is_spt_store) {
    df <- rename(df, Region = Var1, Year = Var2, Sim = Var3, Data = Var4)
    join_keys <- c("Region", "Year", "Sim")
    df <- df %>%
      mutate(Data = ifelse(Data == 1, "low", "high"))
  } else {
    df <- df %>%
      left_join(model_names, by = c("Var1" = "model")) %>%
      rename(Model = Var1, Year = Var2, Sim = Var3)
    join_keys <- c("Year", "Sim")
  }
  df %>%
    left_join(om_df, by = join_keys) %>%
    mutate(RE = (value - OM) / OM, Abs_RE = abs(RE))
}

converged <- function(mod) length(mod) > 1
get_agg  <- function(mod, qty, multi) {
  x <- mod$rep[[qty]]
  if (multi) colSums(x) else x
}

HCR_threshold <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp

  if (stock_status >= 1) {
    f <- frp
  } else if (stock_status > alpha) {
    f <- frp * (stock_status - alpha) / (1 - alpha)
  } else {
    f <- 0
  }

  return(f)
}

melt_ref_pt <- function(arr) {
  reshape2::melt(arr) %>%
    left_join(model_names, by = c("Var1" = "model")) %>%
    rename(Model = Var1, Sim = Var2) %>%
    mutate(model_name = factor(model_name, levels = model_levels))
}


# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)

# Historical Survey Design (1960 - 2024) ----------------------------------

# Read in estimation models
sgl_rg <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp.RDS"))
faa <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_firsthalf.RDS")) # Model 19
three_rg_low <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_lowsamp.RDS"))
three_rg_high <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_highsamp.RDS"))

# Read in operating model
oms <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Reverse order of FAA models
faa_model <- vector("list", 100)
for(i in 1:100) faa_model[[i]] <- faa[[i]][[19]]

# Dimensions
n_yrs <- 65
n_sims <- 100
n_ages <- 30
n_sexes <- 2
n_rg <- 5

### Process Results ---------------------------------------------------------

#### Time Series -------------------------------------------------------------
multi_rg <- c(FALSE, FALSE, TRUE, TRUE)

model_names <- data.frame(
  model      = 1:4,
  model_name = c("sgl", "faa", "three_rg_low", "three_rg_high")
)

n_models <- nrow(model_names)
n_rg     <- 3   # BS, AI, GOA combined

# Storage
agg_ssb_store <- array(NA, dim = c(n_models, n_yrs, n_sims))
agg_rec_store <- array(NA, dim = c(n_models, n_yrs, n_sims))
agg_dep_store <- array(NA, dim = c(n_models, n_yrs, n_sims))
agg_ssb_om    <- array(NA, dim = c(n_yrs, n_sims))
agg_rec_om    <- array(NA, dim = c(n_yrs, n_sims))
agg_dep_om    <- array(NA, dim = c(n_yrs, n_sims))
conv_store    <- array(NA, dim = c(n_models, n_sims))
spt_ssb_store <- array(NA, dim = c(n_rg, n_yrs, n_sims, 2))
spt_rec_store <- array(NA, dim = c(n_rg, n_yrs, n_sims, 2))
spt_dep_store <- array(NA, dim = c(n_rg, n_yrs, n_sims, 2))
spt_ssb_om    <- array(NA, dim = c(n_rg, n_yrs, n_sims))
spt_rec_om    <- array(NA, dim = c(n_rg, n_yrs, n_sims))
spt_dep_om    <- array(NA, dim = c(n_rg, n_yrs, n_sims))

for (i in 1:n_sims) {

  # OM aggregates
  agg_ssb_om[, i] <- colSums(oms$SSB[, 1:n_yrs, i])
  agg_rec_om[, i] <- colSums(oms$Rec[, 1:n_yrs, i])
  agg_dep_om[, i] <- agg_ssb_om[, i] / agg_ssb_om[1, i]

  # EM aggregates
  mods <- list(sgl_rg[[i]], faa_model[[i]], three_rg_low[[i]], three_rg_high[[i]])
  for (m in seq_along(mods)) {
    if (!converged(mods[[m]])) next
    ssb <- get_agg(mods[[m]], "SSB", multi_rg[m])
    rec <- get_agg(mods[[m]], "Rec", multi_rg[m])
    agg_ssb_store[m, , i] <- ssb
    agg_rec_store[m, , i] <- rec
    agg_dep_store[m, , i] <- ssb / ssb[1]
  }

  # OM spatial (BS / AI / GOA combined)
  spt_ssb_om[1, , i] <- oms$SSB[1, 1:n_yrs, i]
  spt_ssb_om[2, , i] <- oms$SSB[2, 1:n_yrs, i]
  spt_ssb_om[3, , i] <- colSums(oms$SSB[3:5, 1:n_yrs, i])
  spt_rec_om[1, , i] <- oms$Rec[1, 1:n_yrs, i]
  spt_rec_om[2, , i] <- oms$Rec[2, 1:n_yrs, i]
  spt_rec_om[3, , i] <- colSums(oms$Rec[3:5, 1:n_yrs, i])
  for (r in 1:n_rg) spt_dep_om[r, , i] <- spt_ssb_om[r, , i] / spt_ssb_om[r, 1, i]

  # EM spatial (three-region low=1, high=2)
  three_rg_mods <- list(three_rg_low[[i]], three_rg_high[[i]])
  for (k in 1:2) {
    if (!converged(three_rg_mods[[k]])) next
    ssb <- three_rg_mods[[k]]$rep$SSB[, 1:n_yrs]
    spt_ssb_store[, , i, k] <- ssb
    spt_rec_store[, , i, k] <- three_rg_mods[[k]]$rep$Rec[, 1:n_yrs]
    spt_dep_store[, , i, k] <- ssb / three_rg_mods[[k]]$rep$SSB[, 1]
  }

} # end i loop

# OM dataframes
agg_ssb_om_df <- melt_om(agg_ssb_om)
agg_rec_om_df <- melt_om(agg_rec_om)
agg_dep_om_df <- melt_om(agg_dep_om)
spt_ssb_om_df <- melt_om(spt_ssb_om, is_spatial = TRUE)
spt_rec_om_df <- melt_om(spt_rec_om, is_spatial = TRUE)
spt_dep_om_df <- melt_om(spt_dep_om, is_spatial = TRUE)

# EM dataframes
agg_ssb_df <- melt_em(agg_ssb_store, agg_ssb_om_df)
agg_rec_df <- melt_em(agg_rec_store, agg_rec_om_df)
agg_dep_df <- melt_em(agg_dep_store, agg_dep_om_df)
spt_ssb_df <- melt_em(spt_ssb_store, spt_ssb_om_df, is_spt_store = TRUE)
spt_rec_df <- melt_em(spt_rec_store, spt_rec_om_df, is_spt_store = TRUE)
spt_dep_df <- melt_em(spt_dep_store, spt_dep_om_df, is_spt_store = TRUE)

#### Reference Points and Catch Advice --------------------------------------------------------
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

true_five_rg_f40 <- get_key_quants(
  data = list(om_values$data), rep = list(om_values$rep),
  reference_points_opt = list(SPR_x = 0.4, t_spwn = 0, sex_ratio_f = 0.5,
                              calc_rec_st_yr = 20, rec_age = 2,
                              type = "multi_region", what = "global_SPR"),
  proj_model_opt = list(n_proj_yrs = 2, HCR_function = HCR_threshold,
                        recruitment_opt = "mean_rec", fmort_opt = "HCR_global",
                        n_avg_yrs = 1), 1)[[1]]

# Per-model configuration
mod_cfg <- list(
  sgl        = list(obj = sgl_rg,       multi = FALSE),
  faa        = list(obj = faa_model,    multi = FALSE),
  three_rg_low  = list(obj = three_rg_low,  multi = TRUE),
  three_rg_high = list(obj = three_rg_high, multi = TRUE)
)

ref_pt_opt_sgl   <- list(SPR_x = 0.4, t_spwn = 0, sex_ratio_f = 0.5,
                         calc_rec_st_yr = 20, rec_age = 2,
                         type = "single_region", what = "SPR")
ref_pt_opt_multi <- list(SPR_x = 0.4, t_spwn = 0, sex_ratio_f = 0.5,
                         calc_rec_st_yr = 20, rec_age = 2,
                         type = "multi_region", what = "global_SPR")
proj_opt_sgl   <- list(n_proj_yrs = 2, HCR_function = HCR_threshold,
                       recruitment_opt = "mean_rec", fmort_opt = "HCR",        n_avg_yrs = 1)
proj_opt_multi <- list(n_proj_yrs = 2, HCR_function = HCR_threshold,
                       recruitment_opt = "mean_rec", fmort_opt = "HCR_global", n_avg_yrs = 1)

f40_values <- array(NA, dim = c(n_models, n_sims))
b40_values <- array(NA, dim = c(n_models, n_sims))
abc_values <- array(NA, dim = c(n_models, n_sims))

for (sim in 1:n_sims) {
  for (m in seq_along(mod_cfg)) {

    # get model object
    mod <- mod_cfg[[m]]$obj[[sim]]
    if (!converged(mod)) next
    multi <- mod_cfg[[m]]$multi

    # get key quants
    kq <- get_key_quants(
      data = list(mod$data), rep = list(mod$rep),
      reference_points_opt = if (multi) ref_pt_opt_multi else ref_pt_opt_sgl,
      proj_model_opt       = if (multi) proj_opt_multi   else proj_opt_sgl,
      1)[[1]]

    # get reference points
    f40_values[m, sim] <- if (multi) unique(kq$F_Ref_Pt)  else kq$F_Ref_Pt
    b40_values[m, sim] <- if (multi) sum(kq$B_Ref_Pt)     else kq$B_Ref_Pt
    abc_values[m, sim] <- if (multi) sum(kq$Catch_Advice) else kq$Catch_Advice
  }
}

# munge reference points
model_levels <- c("sgl", "faa", "three_rg_low", "three_rg_high")
f40_df <- melt_ref_pt(f40_values)
b40_df <- melt_ref_pt(b40_values)
abc_df <- melt_ref_pt(abc_values)


#### Movement ----------------------------------------------------------------
# Read in five region model
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

# Extract movement estimates
move_df <- data.frame()
for(i in 1:100) {
  if(length(three_rg_high[[i]]) > 1) high_move <- reshape2::melt(three_rg_high[[i]]$rep$Movement) %>% mutate(Sim = i, Type = 'high')
  if(length(three_rg_low[[i]]) > 1) low_move <- reshape2::melt(three_rg_low[[i]]$rep$Movement) %>% mutate(Sim = i, Type = 'low')
  move_df <- rbind(move_df, high_move, low_move)
}

# Melt the OM movement array
mov_long <- reshape2::melt(om_values$rep$Movement)

# Collapse destinations 3:5 → 3 for from = 3
mov_3 <- mov_long %>%
  filter(from == 3) %>%
  mutate(to_collapsed = ifelse(to %in% 3:5, 3, to)) %>%
  group_by(from, to_collapsed, years, ages, sexes) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  rename(to = to_collapsed)

# Keep origins 1 and 2
mov_12 <- mov_long %>% filter(from %in% 1:2, to %in% 1:3)

# Back out destination 3 for from = 1 and 2
mov_12_backed_out <- mov_12 %>%
  filter(to != 3) %>%                  # keep only known destinations
  group_by(from, years, ages, sexes) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(to = 3, value = pmax(0, 1 - value))   # remainder goes to to = 3

# Combine everything
mov_final <- bind_rows(
  mov_12 %>% filter(to != 3),  # from 1/2 → 1/2
  mov_12_backed_out,            # from 1/2 → 3
  mov_3                         # from 3 → 1/2/3 (collapsed)
)

move_joined_df <- move_df %>%
  left_join(mov_final %>% rename(true = value),
            by = c("from", 'to', 'years', 'ages', 'sexes'))

# Summarise to quantile ribbons
move_ribbon <- move_joined_df %>%
  filter(years == 1, sexes == 1) %>%
  group_by(from, to, ages, Type) %>%
  summarise(
    lo     = quantile(value, 0.025, na.rm = TRUE),
    med    = quantile(value, 0.500, na.rm = TRUE),
    hi     = quantile(value, 0.975, na.rm = TRUE),
    true   = first(true),
    .groups = "drop"
  ) %>%
  mutate(from =
           case_when(from == 1 ~ "BS+AI+WGOA",
                     from == 2 ~ 'CGOA',
                     from == 3 ~ 'EGOA'),
         to =
           case_when(to == 1 ~ "BS+AI+WGOA",
                     to == 2 ~ 'CGOA',
                     to == 3 ~ 'EGOA'))


### Model Convergence -------------------------------------------------------

# Check convergence
mod_list <- list(sgl_rg, faa_model, three_rg_low, three_rg_high)

for (i in 1:n_sims) {
  for (m in seq_along(mod_list)) {
    if (!converged(mod_list[[m]][[i]])) next
    conv_store[m, i] <- conv_check(sd_rep = mod_list[[m]][[i]]$sd_rep)
  }
}

# figure out convergence rates
apply(conv_store, 1, sum, na.rm = T)

fig_dir_ct <- here("outputs", "cross_test", "spatial_no_block_scenarios", "figs")
dir.create(fig_dir_ct, showWarnings = FALSE, recursive = TRUE)

### SSB ---------------------------------------------------------------------
#### Aggregate ---------------------------------------------------------------

# Absolute
ggplot() +
  geom_line(agg_ssb_df, mapping = aes(x = Year, y = value, group = Sim)) +
  geom_line(agg_ssb_df, mapping = aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~model_name) +
  ylim(0, NA) +
  theme_bw(base_size = 15) +
  labs(y = 'SSB')
ggsave(file.path(fig_dir_ct, 'agg_ssb_abs_historical.png'), width = 10, height = 4, dpi = 300)

# Relative Error
ggplot() +
  geom_line(agg_ssb_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping = aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'agg_ssb_re_historical.png'), width = 10, height = 4, dpi = 300)

#### Spatial Estimates -------------------------------------------------------

# Absolute
ggplot() +
  geom_line(spt_ssb_df, mapping = aes(x = Year, y = value, group = Sim)) +
  geom_line(spt_ssb_df %>% filter(Sim == 1), mapping = aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_grid(Region~Data, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15) +
  labs(y = 'SSB')
ggsave(file.path(fig_dir_ct, 'spt_ssb_abs_historical.png'), width = 8, height = 6, dpi = 300)

# Relative Error
ggplot() +
  geom_line(spt_ssb_df %>%
              group_by(Year, Region, Data) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping = aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15) +
  facet_grid(Region~Data, scales = 'free')
ggsave(file.path(fig_dir_ct, 'spt_ssb_re_historical.png'), width = 8, height = 6, dpi = 300)

### Recruitment -------------------------------------------------------------
#### Aggregate ---------------------------------------------------------------

# Absolute
ggplot() +
  geom_line(agg_rec_df, mapping = aes(x = Year, y = value, group = Sim)) +
  geom_line(agg_rec_df %>% filter(Sim == 1), mapping = aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~model_name) +
  ylim(0, NA) +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'agg_rec_abs_historical.png'), width = 10, height = 4, dpi = 300)

# Relative Error
ggplot() +
  geom_line(agg_rec_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping = aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'agg_rec_re_historical.png'), width = 10, height = 4, dpi = 300)

#### Spatial Estimates -------------------------------------------------------

# Absolute
ggplot() +
  geom_line(spt_rec_df, mapping = aes(x = Year, y = value, group = Sim)) +
  geom_line(spt_rec_df %>% filter(Sim == 1), mapping = aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_grid(Region~Data, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'spt_rec_abs_historical.png'), width = 8, height = 6, dpi = 300)

# Relative Error
ggplot() +
  geom_line(spt_rec_df %>%
              group_by(Year, Region, Data) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping = aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  theme_bw(base_size = 15) +
  facet_grid(Region~Data, scales = 'free')
ggsave(file.path(fig_dir_ct, 'spt_rec_re_historical.png'), width = 8, height = 6, dpi = 300)

### Depletion ---------------------------------------------------------------
#### Aggregate ---------------------------------------------------------------

# Absolute
ggplot() +
  geom_line(agg_dep_df, mapping = aes(x = Year, y = value, group = Sim)) +
  geom_line(agg_dep_df %>% filter(!str_detect(model_name, "_f"), Sim == 1), mapping = aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~model_name) +
  ylim(0, NA) +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'agg_dep_abs_historical.png'), width = 10, height = 4, dpi = 300)

# Relative Error
ggplot() +
  geom_line(agg_dep_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping = aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'agg_dep_re_historical.png'), width = 10, height = 4, dpi = 300)

#### Spatial Estimates -------------------------------------------------------

# Absolute
ggplot() +
  geom_line(spt_dep_df, mapping = aes(x = Year, y = value, group = Sim)) +
  geom_line(spt_dep_df %>% filter(Sim == 1), mapping = aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_grid(Region~Data, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'spt_dep_abs_historical.png'), width = 8, height = 6, dpi = 300)

# Relative Error
ggplot() +
  geom_line(spt_dep_df %>%
              group_by(Year, Region, Data) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping = aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  facet_grid(Region~Data, scales = 'free')
ggsave(file.path(fig_dir_ct, 'spt_dep_re_historical.png'), width = 8, height = 6, dpi = 300)

### Reference Points and Catch Advice ---------------------------------------

# Absolute F40
ggplot(f40_df %>%
         group_by(model_name) %>%
         summarize(median = median(value, na.rm = T),
                   lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = unique(true_five_rg_f40$F_Ref_Pt), lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'F40') +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'f40_abs_historical.png'), width = 6, height = 4, dpi = 300)

# Relative Error F40
ggplot(f40_df %>%
         mutate(RE = (value - unique(true_five_rg_f40$F_Ref_Pt)) / unique(true_five_rg_f40$F_Ref_Pt)) %>%
         group_by(model_name) %>%
         summarize(median = median(RE, na.rm = T),
                   lwr = quantile(RE, 0.025, na.rm = T),
                   upr = quantile(RE, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'F40 Relative Error') +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'f40_re_historical.png'), width = 6, height = 4, dpi = 300)

# Absolute B40
ggplot(b40_df %>%
         group_by(model_name) %>%
         summarize(median = median(value, na.rm = T),
                   lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = sum(true_five_rg_f40$B_Ref_Pt), lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'B40') +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'b40_abs_historical.png'), width = 6, height = 4, dpi = 300)

# Relative Error B40
ggplot(b40_df %>%
         mutate(RE = (value - sum(true_five_rg_f40$B_Ref_Pt)) / sum(true_five_rg_f40$B_Ref_Pt)) %>%
         group_by(model_name) %>%
         summarize(median = median(RE, na.rm = T),
                   lwr = quantile(RE, 0.025, na.rm = T),
                   upr = quantile(RE, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'B40 Relative Error') +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'b40_re_historical.png'), width = 6, height = 4, dpi = 300)

# Absolute ABC
ggplot(abc_df %>%
         group_by(model_name) %>%
         summarize(median = median(value, na.rm = T),
                   lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = sum(true_five_rg_f40$Catch_Advice), lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'ABC') +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'abc_abs_historical.png'), width = 6, height = 4, dpi = 300)

# Relative Error ABC
ggplot(abc_df %>%
         mutate(RE = (value - sum(true_five_rg_f40$Catch_Advice)) / sum(true_five_rg_f40$Catch_Advice)) %>%
         group_by(model_name) %>%
         summarize(median = median(RE, na.rm = T),
                   lwr = quantile(RE, 0.025, na.rm = T),
                   upr = quantile(RE, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'ABC Relative Error') +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'abc_re_historical.png'), width = 6, height = 4, dpi = 300)

### Spatial Model Movement Bias ---------------------------------------------
ggplot(move_ribbon, aes(x = ages, fill = Type, color = Type)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, color = NA) +
  geom_line(aes(y = med)) +
  geom_line(aes(y = true), lty = 2, color = "black", lwd = 1.3) +
  ggh4x::facet_grid2(paste("to", to) ~ paste("from", from),
                     scales = "free", independent = "all") +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'movement_bias_historical.png'), width = 10, height = 8, dpi = 300)
