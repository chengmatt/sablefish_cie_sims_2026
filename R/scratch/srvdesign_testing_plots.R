# Process Results ---------------------------------------------------------
library(here)
library(tidyverse)
srvchng <- readRDS(here("outputs", "mse_results", "spatial_noblock_scenarios", "five_region_base.RDS"))

model_list_srvchng_current <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios",
                                   "single_region_crosstest_srvchng_current.RDS"))
model_list_srvchng_hist <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios",
                                   "single_region_crosstest_srvchng_historical.RDS"))
model_list_srvchng_all <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios",
                                        "single_region_crosstest_srvchng_all.RDS"))

n_sims = 100
sim = 1
i = 1
q_current <- do.call(rbind, lapply(1:n_sims, function(i) model_list_srvchng_current[[i]]$rep$srv_q[1,1,1]))
q_hist <- do.call(rbind, lapply(1:n_sims, function(i) model_list_srvchng_all[[i]]$rep$srv_q[1,1,1]))

# catchability difference ...
(median(q_current) - srvchng$srv_q[1,1,1,1]) / srvchng$srv_q[1,1,1,1]
(median(q_hist) - srvchng$srv_q[1,1,1,1]) / srvchng$srv_q[1,1,1,1]


sim = 5
plot(t(model_list_srvchng_current[[sim]]$rep$Rec))
lines(t(model_list_srvchng_hist[[sim]]$rep$Rec))

plot(model_list_srvchng_current[[sim]]$rep$NAA[1,65,,1])
lines(model_list_srvchng_hist[[sim]]$rep$NAA[1,65,,1])



sim = 12
plot(apply(model_list_srvchng_current[[sim]]$rep$Fmort, 2, sum))
lines(apply(model_list_srvchng_hist[[sim]]$rep$Fmort, 2, sum))

plot(model_list_srvchng_current[[sim]]$rep$fish_sel[1,80,,1,1])
lines(model_list_srvchng_hist[[sim]]$rep$fish_sel[1,80,,1,1])

# median fishery selectivity across sims for current vs historical
sel_current <- do.call(rbind, lapply(1:n_sims, function(i) model_list_srvchng_current[[i]]$rep$srv_sel[1,66,,1,1]))
sel_hist    <- do.call(rbind, lapply(1:n_sims, function(i) model_list_srvchng_hist[[i]]$rep$srv_sel[1,66,,1,1]))

naa_current <- do.call(rbind, lapply(1:n_sims, function(i) model_list_srvchng_current[[i]]$rep$NAA[1,80,,1]))
naa_hist    <- do.call(rbind, lapply(1:n_sims, function(i) model_list_srvchng_hist[[i]]$rep$NAA[1,80,,1]))

plot(apply(naa_current, 2, median), type = 'l', col = 'red', xlab = 'Age', ylab = 'NAA')
lines(apply(naa_hist, 2, median), col = 'blue')
legend('topleft', legend = c('current', 'historical'), col = c('red', 'blue'), lty = 1)

model_list <- list( model_list_srvchng_current, model_list_srvchng_hist, model_list_srvchng_all)
y <- 95
n_sims = 100
ssb_em_results <- array(NA, dim = c(3, length(1:y), n_sims))
rec_em_results <- array(NA, dim = c(3, length(1:y), n_sims))
ssb_om_results <- matrix(NA, nrow = length(1:y), ncol = n_sims)
rec_om_results <- matrix(NA, nrow = length(1:y), ncol = n_sims)

for(i in 1:n_sims) {
  ssb_om_results[,i] <- colSums(srvchng$SSB[,1:y,i])
  rec_om_results[,i] <- colSums(srvchng$Rec[,1:y,i])
}

for(j in 1:3) {
  for(i in 1:n_sims) {
    if(sum(length(model_list[[j]][[i]])) > 1) {
      ssb_em_results[j,,i] <- (as.vector(model_list[[j]][[i]]$rep$SSB) - ssb_om_results[,i]) / ssb_om_results[,i]
      rec_em_results[j,,i] <- (as.vector(model_list[[j]][[i]]$rep$Rec) - rec_om_results[,i]) / rec_om_results[,i]
    }
  }
}

ggplot(reshape2::melt(  ssb_em_results) %>%
         group_by(Var1, Var2) %>%
         summarize(lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T),
                   median = median(value, na.rm = T), na.rm = T) %>%
         mutate(Var1 = case_when(Var1 == 1 ~ 'curr',
                                 Var1 == 2 ~ 'hist',
                                 Var1 == 3 ~ 'all')),
       mapping = aes(x = Var2, y = median, ymin = lwr, ymax = upr,
                     fill = factor(Var1), color = factor(Var1))) +
  geom_line() +
  geom_ribbon(alpha = 0.3, color = NA)

ggplot(reshape2::melt(  rec_em_results) %>%
         group_by(Var1, Var2) %>%
         summarize(lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T),
                   median = median(value, na.rm = T)) %>%
         mutate(Var1 = case_when(Var1 == 1 ~ 'hist',
                                 Var1 == 2 ~ 'all')),
       mapping = aes(x = Var2, y = median,
                     # ymin = lwr, ymax = upr,
                     fill = factor(Var1), color = factor(Var1))) +
  geom_line()
  # geom_ribbon(alpha = 0.3, color = NA)

x <- (srvchng$feedback_start_yr + 1):srvchng$n_yrs
odd_yrs <- x[x %% 2 == 1]
even_yrs <- x[x %% 2 == 0]

agg_true_srvidx_store <- array(NA, dim = c(2, 30, 150))  # design x years 66:95 x sims
designs <- c('current', 'historical')

for(d in 1:2) {
  lls_design_type <- designs[d]
  for(i in 1:150) {
    for(y in 66:95) {
      if(lls_design_type == 'current') {
        if(y == srvchng$feedback_start_yr + 1) {
          bs_true_srvidx <- srvchng$TrueSrvIdx[1,y-2,1,i]
          ai_true_srvidx <- srvchng$TrueSrvIdx[2,y-3,1,i]
          agg_true_srvidx <- srvchng$TrueSrvIdx[,y,,i]
          agg_true_srvidx[1:2,1] <- c(bs_true_srvidx, ai_true_srvidx)
        } else {
          prev_true_srvidx <- srvchng$TrueSrvIdx[,y-1,,i]
          agg_true_srvidx <- srvchng$TrueSrvIdx[,y,,i]
          if(y %in% even_yrs) agg_true_srvidx[1:2,1] <- prev_true_srvidx[1:2,1]
          if(y %in% odd_yrs) agg_true_srvidx[3:5,1] <- prev_true_srvidx[3:5,1]
        }
      }

      if(lls_design_type == 'historical') {
        if(y == srvchng$feedback_start_yr + 1) {
          ai_true_srvidx <- srvchng$TrueSrvIdx[2,y-3,1,i]
          agg_true_srvidx <- srvchng$TrueSrvIdx[,y,,i]
          agg_true_srvidx[2,1] <- ai_true_srvidx
        } else {
          prev_true_srvidx <- srvchng$TrueSrvIdx[,y-1,,i]
          agg_true_srvidx <- srvchng$TrueSrvIdx[,y,,i]
          if(y %in% even_yrs) agg_true_srvidx[2,1] <- prev_true_srvidx[2,1]
          if(y %in% odd_yrs) agg_true_srvidx[1,1] <- prev_true_srvidx[1,1]
        }
      }

      agg_true_srvidx_store[d, y-65, i] <- sum(agg_true_srvidx[,1])

    }
  }
}

colSums(agg_true_srvidx_store)



# compare sd across years for each design
sd_current <- apply(agg_true_srvidx_store[1,,], 1, sd)
sd_hist <- apply(agg_true_srvidx_store[2,,], 1, sd)

plot(66:95, sd_current, type = 'l', col = 'red', ylim = range(c(sd_current, sd_hist)),
     xlab = 'Year', ylab = 'SD of aggregated true survey index')
lines(66:95, sd_hist, col = 'blue')
legend('topleft', legend = c('current', 'historical'), col = c('red', 'blue'), lty = 1)

# true index = sum across ALL regions every year (no imputation)
true_agg_srvidx <- matrix(NA, nrow = 30, ncol = 150)
for(i in 1:150) {
  for(y in 66:95) {
    true_agg_srvidx[y-65, i] <- sum(srvchng$TrueSrvIdx[,y,1,i])
  }
}

# relative bias for each design
bias_current  <- (agg_true_srvidx_store[1,,] - true_agg_srvidx) / true_agg_srvidx
bias_hist     <- (agg_true_srvidx_store[2,,] - true_agg_srvidx) / true_agg_srvidx

# plot median + quantiles
plot(66:95, apply(bias_current, 1, median), type = 'l', col = 'red',
     ylim = range(c(bias_current, bias_hist)),
     xlab = 'Year', ylab = 'Relative bias in aggregated survey index')
lines(66:95, apply(bias_hist, 1, median), col = 'blue')
polygon(c(66:95, 95:66), c(apply(bias_current, 1, quantile, 0.25), rev(apply(bias_current, 1, quantile, 0.75))),
        col = adjustcolor('red', 0.2), border = NA)
polygon(c(66:95, 95:66), c(apply(bias_hist, 1, quantile, 0.25), rev(apply(bias_hist, 1, quantile, 0.75))),
        col = adjustcolor('blue', 0.2), border = NA)
abline(h = 0, lty = 2)
legend('topleft', legend = c('current', 'historical'), col = c('red', 'blue'), lty = 1)

# store srv_prop for GOA regions (3:5) across projection years for both designs
srv_prop_store <- array(NA, dim = c(2, 30, 5, 150))  # design x years 66:95 x regions x sims

designs <- c('current', 'historical')

for(d in 1:2) {
  lls_design_type <- designs[d]
  for(i in 1:150) {
    # compute srv_prop using numbers weighting (same as in agg_data_to_single_rg)
    true_srvidx <- srvchng$TrueSrvIdx[,1:95,1,i]  # fleet 1
    total_srv <- apply(true_srvidx, 2, sum)
    srv_prop <- sweep(true_srvidx, 2, total_srv, "/")

    # zero out historical pre-feedback unsampled years
    srv_prop[1, c(31:37, seq(39,65,2))] <- 0
    srv_prop[2, c(31:36, seq(38,64,2))] <- 0
    srv_prop[, 65] <- 0  # no survey in 2024

    # zero out projection years based on design
    for(yr in 66:95) {
      if(lls_design_type == 'current') {
        if(yr %in% even_yrs) srv_prop[1:2, yr] <- 0
        if(yr %in% odd_yrs)  srv_prop[3:5, yr] <- 0
      }
      if(lls_design_type == 'historical') {
        if(yr %in% even_yrs) srv_prop[2, yr] <- 0
        if(yr %in% odd_yrs)  srv_prop[1, yr] <- 0
      }
      srv_prop[, yr] <- srv_prop[, yr] / sum(srv_prop[, yr])  # renormalize
      srv_prop_store[d, yr-65, , i] <- srv_prop[, yr]
    }
  }
}

# median GOA weight (regions 3:5) per year per design
goa_wgt_current <- apply(srv_prop_store[1,, 3:5,], c(1), function(x) median(rowSums(matrix(x, nrow=30))))
goa_wgt_hist    <- apply(srv_prop_store[2,, 3:5,], c(1), function(x) median(rowSums(matrix(x, nrow=30))))

# cleaner: sum GOA regions first then take median across sims
goa_wgt_current <- apply(apply(srv_prop_store[1,,3:5,], c(1,3), sum), 1, median)
goa_wgt_hist    <- apply(apply(srv_prop_store[2,,3:5,], c(1,3), sum), 1, median)

plot(66:95, goa_wgt_current, type = 'l', col = 'red',
     ylim = c(0,1), xlab = 'Year', ylab = 'Median GOA weight in survey comps')
lines(66:95, goa_wgt_hist, col = 'blue')
abline(h = 0.5, lty = 2)
legend('topleft', legend = c('current', 'historical'), col = c('red', 'blue'), lty = 1)
