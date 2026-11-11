# Purpose: FAA EM structured like operational Model 26.2_FAA (regional BS / AI / GOA LLS fleets, a single
#          fixed gear fishery fleet + trawl, logistic selectivity, no interpolation) fit to the same OM, for
#          comparison with the single region comps test. Fit twice per sim: LLS data as run (historical design:
#          GOA every year, BS / AI alternating, unsurveyed years left missing) and complete LLS data (all regions
#          surveyed every year: all indices and all age comps). Survey indices are summed from the OM's regional
#          observations (index_helpers.R), so the GOA index (3 OM regions) is more precise than BS or AI.
# Creator: Matthew LH. Cheng
# Date: 11/10/26
# Notes: Requires SPoRC at commit 72b8426a (pre-seasonality); see README. Run in its own R session:
#        faa_em() replaces the aggregated data arrays in the OM environment.

# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)
library(parallel)

source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "faa_em.R"))
source(here("R", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test", "index_helpers.R"))

sims <- 1:100
ncores <- 3
Y <- 65

out_dir <- here("outputs", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

om <- readRDS(here("outputs", "mse_results", "spatial_no_block_scenarios", "five_region_base.RDS"))
n_sims <- om$n_sims   # needed as a global by add_aggregated_obj_to_simenv()
sp <- expand.grid(sex = 1:2, par = 1:2)

# Run ---------------------------------------------------------------------
variants <- data.frame(design = c("historical", "all"), name = c("faa_1fish", "faa_1fish_all"))

run_one <- function(i) {
  out <- list()
  for (k in seq_len(nrow(variants))) {
    des <- variants$design[k]
    set.seed(3000 + i)                     # same seed (same data) for every variant
    a <- faa_em(sim_env = om, y = Y, sim = i, srv_idx_se = 0.2, age_lag = 1, lls_design_type = des,
                faa_n_fish_fleets = 2, faa_n_srv_fleets = 4, srv_wgt = "numbers", fish_wgt = "numbers",
                fish_sel_blocks = paste0("none_Fleet_", 1:2),
                fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),
                fish_fixed_sel_pars_spec = rep("est_all", 2),
                fish_selex_prior = build_selex_prior(fleets_blocks = list("1" = 1, "2" = 1), sex_par = sp, region = 1, mu = 3.5, sd = 2),
                srv_sel_blocks = paste0("none_Fleet_", 1:4),
                srv_sel_model = paste0("logist1_Fleet_", 1:4),
                srv_fixed_sel_pars_spec = rep("est_all", 4),
                srv_selex_prior = build_selex_prior(fleets_blocks = list("1" = 1, "2" = 1, "3" = 1, "4" = 1), sex_par = sp, region = 1, mu = 3.5, sd = 2))
    a <- rebuild_faa_index(a, om, i, Y)   # indices summed from the regional observations (see index_helpers.R)
    fit <- tryCatch(fit_model(a$data, a$par, a$map, NULL, 2, silent = TRUE), error = function(e) NULL)
    out[[variants$name[k]]] <- if (is.null(fit)) NULL else list(SSB = as.vector(fit$rep$SSB[1, 1:Y]),
                                                               NAA = apply(fit$rep$NAA[1, 1:Y, , ], 1:2, sum),
                                                               maxgrad = max(abs(fit$gr())))
  }
  message("sim ", i, " done")
  out
}

res <- mclapply(sims, run_one, mc.cores = ncores)
names(res) <- sims
saveRDS(res, file.path(out_dir, "faa_comps_test.RDS"))
