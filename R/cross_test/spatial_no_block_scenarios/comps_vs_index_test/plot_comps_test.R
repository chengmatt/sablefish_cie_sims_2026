# Purpose: Summarize and plot the comps vs. index test (single region variants + 26.2_FAA-like FAA EM,
#          fit to LLS data as run and to complete LLS data; indices summed from regional observations).
#          Error is relative to the five-region OM truth, averaged over 2010-2024, per simulation.
# Creator: Matthew LH. Cheng
# Date: 11/10/26

# Setup -------------------------------------------------------------------
library(here)
library(tidyverse)

out_dir <- here("outputs", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test")
fig_dir <- file.path(out_dir, "figs"); dir.create(fig_dir, showWarnings = FALSE)

om <- readRDS(here("outputs", "mse_results", "spatial_no_block_scenarios", "five_region_base.RDS"))
sgl <- readRDS(file.path(out_dir, "single_region_comps_test.RDS"))
faa <- readRDS(file.path(out_dir, "faa_comps_test.RDS"))
yrs <- 51:65           # 2010-2024
old_ages <- 9:30       # model ages 10+

# Per-simulation errors ---------------------------------------------------
err <- function(x, i) {
  tr_naa <- apply(om$NAA[, yrs, , , i], 2:3, sum)
  c(SSB = mean(x$SSB[yrs] / colSums(om$SSB[, yrs, i]) - 1),
    Old = mean(rowSums(x$NAA[yrs, old_ages]) / rowSums(tr_naa[, old_ages]) - 1))
}
collect <- function(lst) bind_rows(lapply(names(lst), function(i) bind_rows(lapply(names(lst[[i]]), function(k) {
  x <- lst[[i]][[k]]; if (is.null(x) || x$maxgrad > 1e-3) return(NULL)   # drop non-converged fits
  e <- err(x, as.integer(i)); data.frame(sim = as.integer(i), variant = k, SSB = e[["SSB"]], Old = e[["Old"]])
}))))
d <- bind_rows(collect(sgl), collect(faa))

# scenario key: model, LLS index, LLS age comps; 'scen' is the plotting label for the main figures
key <- tribble(~variant, ~model, ~comps, ~index, ~scen,
  "hist",                   "Single region", "missing",  "Carryover", "sr_asrun",
  "full_idx_partial_comps", "Single region", "missing",  "Complete",  NA,
  "goa_comps_carry_idx",    "Single region", "goa",      "Carryover", "sr_goa",
  "goa_comps_full_idx",     "Single region", "goa",      "Complete",  NA,
  "bs_comps_carry_idx",     "Single region", "bs",       "Carryover", "sr_bs",
  "bs_comps_full_idx",      "Single region", "bs",       "Complete",  NA,
  "carry_idx_full_comps",   "Single region", "complete", "Carryover", "sr_complete",
  "all",                    "Single region", "complete", "Complete",  NA,
  "faa_1fish",              "26.2_FAA",      "asrun",    "As run",    "faa_asrun",
  "faa_1fish_all",          "26.2_FAA",      "complete", "Complete",  "faa_complete")
d <- left_join(d, key, by = "variant")
write_csv(d, file.path(out_dir, "comps_test_errors_by_sim.csv"))
med <- d %>% group_by(model, comps, index, variant) %>% summarise(n = n(), SSB_median = median(SSB), Old_median = median(Old), .groups = "drop")
write_csv(med, file.path(out_dir, "comps_test_medians.csv"))
print(med, n = Inf)

# Plot setup ----------------------------------------------------------------
# main scenarios: single region on the carryover index (as in 26.1_1Reg) with LLS age comps from different regions,
# and 26.2_FAA fit to the regional LLS data as surveyed or with every region surveyed every year
scen_lab <- c(sr_asrun = "Comps: GOA +\nBS or AI\n(as run)", sr_goa = "Comps:\nGOA only", sr_bs = "Comps:\nBS only",
              sr_complete = "Comps:\nall regions", faa_asrun = "As surveyed\n(BS / AI every\nother year)", faa_complete = "All regions\nevery year")
scen_long <- c(sr_asrun = "Single region, carryover index; comps from GOA + BS or AI (as run)",
               sr_goa = "Single region, carryover index; comps from GOA only",
               sr_bs = "Single region, carryover index; comps from BS only",
               sr_complete = "Single region, carryover index; comps from all regions",
               faa_asrun = "26.2_FAA, LLS data as surveyed (BS / AI every other year)",
               faa_complete = "26.2_FAA, LLS data from all regions every year")
scen_col <- c(sr_asrun = "#E69F00", sr_goa = "#D55E00", sr_bs = "#CC79A7", sr_complete = "#009E73",
              faa_asrun = "#0072B2", faa_complete = "#56B4E9")   # ggthemes colorblind
scen_lvl <- names(scen_lab)
group_lab <- c("Single Region: Carryover LLS Index (as in 26.1_1Reg)", "26.2_FAA")

# Main figure: single region on the carryover index with comps from different regions, and 26.2_FAA as surveyed ------
dm <- d %>% filter(scen %in% c("sr_asrun", "sr_goa", "sr_bs", "sr_complete", "faa_asrun")) %>%
  mutate(scen = factor(scen, scen_lvl), group = factor(ifelse(model == "26.2_FAA", group_lab[2], group_lab[1]), group_lab))
p <- ggplot(dm, aes(x = scen, y = SSB, fill = scen)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  geom_boxplot(outliers = FALSE) +
  facet_grid(~group, scales = 'free_x', space = 'free_x') +
  scale_x_discrete(labels = scen_lab) +
  scale_fill_manual(values = scen_col) +
  labs(x = NULL, y = 'SSB Relative Error\n(Estimated / True - 1, 2010-2024)') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none')
ggsave(file.path(fig_dir, "comps_test_main.png"), p, width = 13, height = 6.5, dpi = 300)

# Backup figure: each model fit to the LLS data as run and with all regions surveyed every year, SSB and fish 10+ ----
bx_lab <- c(missing = "Comps: GOA +\nBS or AI (as run)", goa = "Comps:\nGOA only", bs = "Comps:\nBS only", complete = "Comps:\nall regions", faa = "Regional LLS\nindices and comps")
db <- d %>%
  mutate(x = factor(ifelse(model == "26.2_FAA", "faa", comps), names(bx_lab)),
         data = factor(ifelse(index %in% c("Carryover", "As run"), "asrun", "all"), c("asrun", "all")),
         group = factor(ifelse(model == "26.2_FAA", "26.2_FAA", "Single Region"), c("Single Region", "26.2_FAA"))) %>%
  pivot_longer(c(SSB, Old), names_to = "metric") %>%
  mutate(metric = factor(metric, c("SSB", "Old"), c("SSB", "Numbers of Fish Age 10+")))
p <- ggplot(db, aes(x = x, y = value, fill = data)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  geom_boxplot(outliers = FALSE) +
  facet_grid(metric ~ group, scales = 'free', space = 'free_x') +
  scale_x_discrete(labels = bx_lab) +
  scale_fill_manual(values = c(asrun = "#E69F00", all = "#56B4E9"),
                    labels = c(asrun = "As run (single region: carryover index; 26.2_FAA: as surveyed)",
                               all = "All regions every year (single region: index only; 26.2_FAA: index and comps)")) +
  labs(x = NULL, y = 'Relative Error (Estimated / True - 1, 2010-2024)', fill = 'LLS Data') +
  theme_bw(base_size = 15) +
  theme(legend.position = 'top') +
  guides(fill = guide_legend(ncol = 1))
ggsave(file.path(fig_dir, "comps_test_all_scenarios.png"), p, width = 13, height = 8.5, dpi = 300)

# Time series figures: SSB by year, main scenarios ------------------------------------------
yrs_ts <- 1:65                                   # 1960-2024
tr_ssb <- colSums(om$SSB[, yrs_ts, 1])          # OM truth is identical across sims through 2024
ts_collect <- function(lst) bind_rows(lapply(names(lst), function(i) bind_rows(lapply(names(lst[[i]]), function(k) {
  x <- lst[[i]][[k]]; if (is.null(x) || x$maxgrad > 1e-3) return(NULL)
  data.frame(sim = as.integer(i), variant = k, year = 1959 + yrs_ts, SSB = x$SSB[yrs_ts])
}))))
dts <- bind_rows(ts_collect(sgl), ts_collect(faa)) %>% left_join(key, by = "variant") %>% filter(!is.na(scen)) %>%
  mutate(truth = tr_ssb[year - 1959], err = SSB / truth - 1, scen = factor(scen, scen_lvl))
tsm <- dts %>% group_by(scen, year, truth) %>%
  summarise(med = median(SSB), lo = quantile(SSB, 0.1), hi = quantile(SSB, 0.9),
            emed = median(err), elo = quantile(err, 0.1), ehi = quantile(err, 0.9), .groups = "drop")
write_csv(tsm, file.path(out_dir, "comps_test_ssb_timeseries.csv"))

# SSB vs. truth, one panel per scenario (lines = simulations, red = OM)
p <- ggplot() +
  geom_line(dts, mapping = aes(x = year, y = SSB, group = sim), alpha = 0.2) +
  geom_line(tsm, mapping = aes(x = year, y = truth), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~scen, nrow = 2, labeller = as_labeller(gsub("; ", "\n", scen_long))) +
  ylim(0, NA) +
  labs(x = 'Year', y = 'SSB (kt)') +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir, "comps_test_ssb_timeseries.png"), p, width = 13, height = 8, dpi = 300)

# Relative SSB error by year (median over simulations)
p <- ggplot(tsm, aes(x = year, y = emed, color = scen)) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  scale_color_manual(values = scen_col, labels = scen_long) +
  labs(x = 'Year', y = 'SSB Relative Error\n(Median of 100 Simulations)', color = NULL) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'top') +
  guides(color = guide_legend(ncol = 2))
ggsave(file.path(fig_dir, "comps_test_error_timeseries.png"), p, width = 15, height = 7, dpi = 300)

# Single region under the current alternating LLS design applied historically (run_single_region_alt_design.R) ------
alt_file <- file.path(out_dir, "single_region_alt_design.RDS")
if (file.exists(alt_file)) {
  alt <- readRDS(alt_file)
  da <- bind_rows(lapply(names(sgl), function(i) {
    L <- list(hist = sgl[[i]]$hist, alt = alt[[i]]$alt_both)
    bind_rows(lapply(names(L), function(k) { x <- L[[k]]; if (is.null(x) || x$maxgrad > 1e-3) return(NULL)
      data.frame(sim = as.integer(i), design = k, year = 1959 + yrs_ts, err = x$SSB[yrs_ts] / tr_ssb - 1) }))
  }))
  print(da %>% group_by(design, sim) %>% summarise(avg = mean(err[year %in% (1959 + yrs)]), y2024 = err[year == 2024], .groups = "drop") %>%
          group_by(design) %>% summarise(n = n(), avg = round(100 * median(avg)), y2024 = round(100 * median(y2024))))
  am <- da %>% group_by(design, year) %>% summarise(err = median(err), .groups = "drop")
  p <- ggplot(am %>% filter(year >= 1990), aes(x = year, y = err, color = design)) +
    geom_line(lwd = 1.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
    scale_color_manual(values = c(hist = "#E69F00", alt = "#000000"), breaks = c("hist", "alt"),
                       labels = c(hist = "Historical design (as run): GOA every year, BS / AI in alternate years",
                                  alt = "Current design from 1997: GOA and BS + AI in alternate years")) +
    labs(x = 'Year', y = 'SSB Relative Error\n(Median of 100 Simulations)', color = 'Single Region, LLS Design') +
    theme_bw(base_size = 15) +
    theme(legend.position = 'top') +
    guides(color = guide_legend(ncol = 1))
  ggsave(file.path(fig_dir, "sr_alternating_design.png"), p, width = 12, height = 6.5, dpi = 300)
}

# Carryover lag over time: single region with the as-run comps, complete vs. carryover index --------------------------
# (historical design = `hist`; current design = `alt_index` from run_single_region_alt_design.R; complete = `full_idx_partial_comps`)
if (file.exists(alt_file) && !is.null(alt[[1]]$alt_index)) {
  source(here("R", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test", "index_helpers.R"))   # carry_hist(), carry_alt()
  des_col <- c(complete = "#56B4E9", hist = "#E69F00", alt = "#000000")
  des_lab <- c(complete = "Complete index (all regions every year)", hist = "Carryover index, historical design (as run)",
               alt = "Carryover index, current design from 1997")
  dc <- bind_rows(lapply(names(sgl), function(i) {
    L <- list(complete = sgl[[i]]$full_idx_partial_comps, hist = sgl[[i]]$hist, alt = alt[[i]]$alt_index)
    if (any(sapply(L, function(x) is.null(x) || x$maxgrad > 1e-3))) return(NULL)
    bind_rows(lapply(names(L), function(k) data.frame(sim = as.integer(i), index = k, year = 1959 + yrs_ts, err = L[[k]]$SSB[yrs_ts] / tr_ssb - 1)))
  }))
  # within-sim carryover effect on SSB error: carryover index minus complete index (same comps, same data otherwise)
  dd <- dc %>% pivot_wider(names_from = index, values_from = err) %>%
    transmute(sim, year, hist = hist - complete, alt = alt - complete) %>% pivot_longer(c(hist, alt), names_to = "index", values_to = "diff") %>%
    group_by(index, year) %>% summarise(med = median(diff), lo = quantile(diff, 0.1), hi = quantile(diff, 0.9), .groups = "drop")
  # carryover index error vs. the true index (OM truth, identical across sims through 2024), 1990-2023
  tri <- om$TrueSrvIdx[, 1:65, , 1]; iy <- 31:64
  di <- tibble(year = 1959 + iy, true = colSums(tri[, iy, 1]), hist = colSums(carry_hist(tri)[, iy, 1]) / true - 1,
               alt = colSums(carry_alt(tri)[, iy, 1]) / true - 1) %>% pivot_longer(c(hist, alt), names_to = "index", values_to = "err")
  ssb_med <- dc %>% group_by(index, year) %>% summarise(err = median(err), .groups = "drop") %>% filter(year >= 1990)

  print(ssb_med %>% filter(year %in% c(2017, 2020, 2022, 2024)) %>% mutate(err = round(100 * err)) %>% pivot_wider(names_from = index, values_from = err))
  ssb_plot <- function(keep) {
    ggplot(ssb_med %>% filter(index %in% keep) %>% mutate(index = factor(index, keep)), aes(x = year, y = err, color = index)) +
      geom_line(lwd = 1.3) +
      geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
      scale_color_manual(values = des_col, labels = des_lab) +
      labs(x = 'Year', y = 'SSB Relative Error\n(Median of 100 Simulations)', color = NULL) +
      theme_bw(base_size = 15) + theme(legend.position = 'top') + guides(color = guide_legend(ncol = 1))
  }
  # (1) historical design: SSB error over time, complete vs. carryover index (terminal years biased low as the index lags)
  ggsave(file.path(fig_dir, "carryover_ssb_terminal_hist.png"), ssb_plot(c("complete", "hist")), width = 10, height = 6.5, dpi = 300)
  # (2) current design added: the terminal-year low bias is exacerbated
  ggsave(file.path(fig_dir, "carryover_ssb_terminal_current.png"), ssb_plot(c("complete", "hist", "alt")), width = 10, height = 6.5, dpi = 300)
}
