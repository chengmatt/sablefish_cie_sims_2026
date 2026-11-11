# Purpose: Carryover LLS index vs. the true Alaska-wide LLS index over 1990-2023 (the survey years fit in the
#          simulation test), under the historical design (as in 26.1_1Reg; one of BS / AI carried each year) and the
#          current alternating design applied from 1997 (BS + AI or GOA carried each year). Same carryover rules as the
#          fitted runs (index_helpers.R). OM truth only (no observation error, no model fits); the OM's historical
#          period is identical across simulations.
# Creator: Matthew LH. Cheng
# Date: 11/11/26

# Setup -------------------------------------------------------------------
library(here)
library(tidyverse)
library(patchwork)

source(here("R", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test", "index_helpers.R"))   # carry_hist(), carry_alt()

out_dir <- here("outputs", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test")
fig_dir <- file.path(out_dir, "figs"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
om <- readRDS(here("outputs", "mse_results", "spatial_no_block_scenarios", "five_region_base.RDS"))

tr <- om$TrueSrvIdx[, 1:65, , 1]   # region x year x fleet, sim 1 (identical across sims through 2024)
yrs <- 31:64                       # 1990-2023 (no survey in 2024)

# Carryover index by design -------------------------------------------------
idx <- function(x) colSums(x[, yrs, 1])
d <- tibble(year = 1959 + yrs, true = idx(tr), hist = idx(carry_hist(tr)), alt = idx(carry_alt(tr))) %>%
  mutate(err_hist = hist / true - 1, err_alt = alt / true - 1)
# share of the true index that is carried over each year (regions whose value was not observed that year)
carried_share <- function(cx) colSums(tr[, yrs, 1] * (cx[, yrs, 1] != tr[, yrs, 1])) / colSums(tr[, yrs, 1])
d <- d %>% mutate(share_hist = carried_share(carry_hist(tr)), share_alt = carried_share(carry_alt(tr)))
write_csv(d, file.path(out_dir, "carryover_index_lag.csv"))

# Summaries ---------------------------------------------------------------
cat("mean share of the index carried, 1997-2023: historical", round(100 * mean(d$share_hist[d$year >= 1997])),
    "%; current design", round(100 * mean(d$share_alt[d$year >= 1997])), "%\n")
print(d %>% mutate(period = case_when(year %in% 1997:2015 ~ "1997-2015 (stable / slow decline)", year %in% 2017:2023 ~ "2017-2023 (rapid rise)")) %>%
        filter(!is.na(period)) %>% group_by(period) %>%
        summarise(hist_mean_err = round(100 * mean(err_hist), 1), alt_mean_err = round(100 * mean(err_alt), 1),
                  hist_min = round(100 * min(err_hist)), alt_min = round(100 * min(err_alt))))

# Figure ------------------------------------------------------------------
# top: what the population is doing (true Alaska-wide LLS index); bottom: carryover index error under each design
de <- d %>% select(year, err_hist, err_alt) %>% pivot_longer(-year, names_to = "design", values_to = "err") %>%
  mutate(design = factor(design, c("err_hist", "err_alt")))
pA <- ggplot(d, aes(x = year, y = true)) +
  geom_line(lwd = 1.3) +
  labs(x = NULL, y = 'True LLS Index') +
  theme_bw(base_size = 15)
pB <- ggplot(de, aes(x = year, y = err, color = design)) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  scale_color_manual(values = c(err_hist = "#E69F00", err_alt = "#000000"),
                     labels = c(err_hist = "Historical design (as run): one of BS / AI carried",
                                err_alt = "Current design from 1997: GOA or BS + AI carried")) +
  labs(x = 'Year', y = 'Relative Error\n(Carryover vs. True Index)', color = NULL) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(ncol = 1))
p <- pA / pB
ggsave(file.path(fig_dir, "carryover_index_lag.png"), p, width = 10, height = 7.5, dpi = 300)

# Historical design only (for the first carryover slide) ---------------------------------------------------------
pB_hist <- ggplot(de %>% filter(design == "err_hist"), aes(x = year, y = err)) +
  geom_line(lwd = 1.3, color = "#E69F00") +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  labs(x = 'Year', y = 'Relative Error\n(Carryover vs. True Index)') +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir, "carryover_index_lag_hist.png"), pA / pB_hist, width = 10, height = 7, dpi = 300)

# True LLS index by region relative to 2010: where the large 2014-2020 year classes showed up ------------------
ti <- om$TrueSrvIdx[, 1:65, 1, 1]   # identical across sims through 2024
ri <- tibble(year = 1959 + 1:65, BS = ti[1, ], AI = ti[2, ], GOA = colSums(ti[3:5, ]), `Alaska-wide` = colSums(ti)) %>%
  pivot_longer(-year, names_to = "region") %>% group_by(region) %>% mutate(value = value / value[year == 2010]) %>% ungroup() %>%
  mutate(region = factor(region, c("BS", "AI", "GOA", "Alaska-wide"))) %>% filter(year >= 1990)
p <- ggplot(ri, aes(x = year, y = value, color = region, lwd = region)) +
  geom_line() +
  scale_color_manual(values = c(BS = "#56B4E9", AI = "#009E73", GOA = "#E69F00", `Alaska-wide` = "black")) +
  scale_linewidth_manual(values = c(BS = 1.1, AI = 1.1, GOA = 1.1, `Alaska-wide` = 1.8)) +
  labs(x = 'Year', y = 'True LLS Index (Relative to 2010)', color = 'Region', lwd = 'Region') +
  theme_bw(base_size = 15) +
  theme(legend.position = 'top')
ggsave(file.path(fig_dir, "true_regional_lls_index.png"), p, width = 10, height = 6, dpi = 300)
