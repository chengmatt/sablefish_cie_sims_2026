# Purpose: Shared LLS index helpers for the comps vs. index test.
#          (1) Carryover rules for the single region index: the historical design (as in agg_data_to_single_rg()) and
#              the current alternating design applied from 1997 (BS + AI in odd years, GOA in even years).
#          (2) FAA survey indices built from the OM's regional observed indices (each region observed with SE 0.2), so
#              the GOA index (3 OM regions) is more precise (SE ~0.12), as a design-based index would be. The SE given
#              to the EM is the delta-method CV of the sum. Replaces the default in agg_data_to_faa(), which adds a
#              fresh SE 0.2 error to the true summed index. (The single region EM keeps the default index.)
# Creator: Matthew LH. Cheng
# Date: 11/11/26

# Carryover rules (x: region x year x fleet; regions 1 = BS, 2 = AI, 3:5 = GOA; model year 31 = 1990) ----------------
bsai_year <- function(y) y >= 38 && y %% 2 == 0   # current design: odd calendar years from 1997, BS + AI surveyed
goa_year  <- function(y) y >= 38 && y %% 2 == 1   # current design: even calendar years from 1998, GOA surveyed

# historical design, as in agg_data_to_single_rg(): GOA every year, BS in odd / AI in even years from 1997 / 1996; an
# unsurveyed BS / AI year takes the last observed value (earlier years take the first observed value)
carry_hist <- function(x) {
  x[1, 31:37, 1] <- x[1, 38, 1]; x[1, seq(39, 63, 2), 1] <- x[1, seq(39, 63, 2) - 1, 1]
  x[2, 31:36, 1] <- x[2, 37, 1]; x[2, seq(38, 64, 2), 1] <- x[2, seq(38, 64, 2) - 1, 1]
  x
}
# current alternating design from 1997 (1990-1996 as historical): the unsurveyed half takes the previous year's value
carry_alt <- function(x) {
  x[1, 31:37, 1] <- x[1, 38, 1]; x[2, 31:36, 1] <- x[2, 37, 1]
  for (y in 38:dim(x)[2]) {
    if (bsai_year(y)) x[3:5, y, 1] <- x[3:5, y - 1, 1]
    if (goa_year(y))  x[1:2, y, 1] <- x[1:2, y - 1, 1]
  }
  x
}

# delta-method CV of a sum of independent regional indices (obs, se: region x year)
sum_cv <- function(obs, se) sqrt(colSums((obs * se)^2)) / colSums(obs)

# FAA EM with 4 survey fleets: BS LLS, AI LLS, GOA LLS (OM regions 3-5), and OM survey fleet 3 summed over all regions
rebuild_faa_index <- function(a, om, i, y) {
  regs <- list(1, 2, 3:5, 1:5); fl <- c(1, 1, 1, 3)
  for (k in 1:4) {
    O  <- matrix(om$ObsSrvIdx[regs[[k]], 1:y, fl[k], i], nrow = length(regs[[k]]))
    SE <- matrix(om$ObsSrvIdx_SE[regs[[k]], 1:y, fl[k]], nrow = length(regs[[k]]))
    a$data$ObsSrvIdx[1, 1:y, k]    <- colSums(O)
    a$data$ObsSrvIdx_SE[1, 1:y, k] <- sum_cv(O, SE)
  }
  a
}
