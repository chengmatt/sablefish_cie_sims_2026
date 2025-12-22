# Figure out what the FAA selex curves look like
nll_caa <- function(data, pars) {
  RTMB::getAll(data, pars)
  f_ar = exp(ln_f_ar)
  aa_est = naa_total * (f_ar / (f_ar + m)) * (1 - exp(- (f_ar + m)))
  nll = sum((log(aa_r) - log(aa_est))^2)
  REPORT(f_ar)
  REPORT(aa_est)
  return(nll)
}

nll_srv_aa <- function(data, pars) {
  RTMB::getAll(data, pars)
  f_ar = exp(ln_f_ar)
  aa_est = naa_total * f_ar
  nll = sum((log(aa_r) - log(aa_est))^2)
  REPORT(f_ar)
  REPORT(aa_est)
  return(nll)
}


# get faa from fishery
get_faa_sel <- function(fleet, sex, yrs = 1:65) {
  # BS
  aa_r <- colMeans(om_values$CAA[1,yrs,,sex,fleet,1])
  # aa_r <- apply(om_values$CAA[,yrs,,sex,fleet,1], 3, mean)
  naa_total <- apply(om_values$NAA[,c(yrs),,sex,1], 3, mean)

  data <- list(
    aa_r = aa_r,
    naa_total = naa_total,
    m = m
  )

  pars <- list(
    ln_f_ar = rep(log(0.1), length(aa_r))
  )

  adfun <- MakeADFun(SPoRC:::cmb(nll_caa, data), pars, random = NULL)
  optim <- stats::nlminb(adfun$par, adfun$fn, adfun$gr,
                         control =  list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15))
  sd <- sdreport(adfun)
  rep <- adfun$report(adfun$env$last.par.best)

  par(mfrow = c(3,2), cex = 0.8)
  # look at selex
  plot(rep$f_ar / max(rep$f_ar), type = 'l', ylim = c(0,1), xlab = 'Age', ylab = 'Selex', main = 'BS')

  # compare CAA
  plot(rep$aa_est, type = 'l', lty = 2, ylab = 'CAA', xlab = 'Age', main = 'BS')
  lines(data$aa_r, col = 'red')

  # AI
  aa_r <- colSums(om_values$CAA[2,yrs,,sex,fleet,1])
  naa_total <- apply(om_values$NAA[,yrs,,sex,1], 3, sum)

  data <- list(
    aa_r = aa_r,
    naa_total = naa_total,
    m = m
  )

  pars <- list(
    ln_f_ar = rep(log(0.1), length(aa_r))
  )

  adfun <- MakeADFun(SPoRC:::cmb(nll_caa, data), pars, random = NULL)
  optim <- stats::nlminb(adfun$par, adfun$fn, adfun$gr,
                         control =  list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15))
  sd <- sdreport(adfun)
  rep <- adfun$report(adfun$env$last.par.best)

  # look at selex
  plot(rep$f_ar / max(rep$f_ar), type = 'l', ylim = c(0,1), xlab = 'Age', ylab = 'Selex', main = 'AI')

  # compare CAA
  plot(rep$aa_est, type = 'l', lty = 2, ylab = 'CAA', xlab = 'Age', main = 'AI')
  lines(data$aa_r, col = 'red')

  # GOA
  aa_r <- apply(om_values$CAA[3:5,,,sex,fleet,1], 3, mean)
  naa_total <- apply(om_values$NAA[,-66,,sex,1], 3, mean)

  data <- list(
    aa_r = aa_r,
    naa_total = naa_total,
    m = m
  )

  pars <- list(
    ln_f_ar = rep(log(0.1), length(aa_r))
  )

  adfun <- MakeADFun(SPoRC:::cmb(nll_caa, data), pars, random = NULL)
  optim <- stats::nlminb(adfun$par, adfun$fn, adfun$gr,
                         control =  list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15))
  sd <- sdreport(adfun)
  rep <- adfun$report(adfun$env$last.par.best)

  # look at selex
  plot(rep$f_ar / max(rep$f_ar), type = 'l', ylim = c(0,1), xlab = 'Age', ylab = 'Selex', main = 'GOA')

  # compare CAA
  plot(rep$aa_est, type = 'l', lty = 2, ylab = 'CAA', xlab = 'Age', main = 'GOA')
  lines(data$aa_r, col = 'red')
}


# get srvaa from survey
get_srvaa_sel <- function(fleet, sex, yrs = 1:65) {
  # BS
  aa_r <- colMeans(om_values$SrvIAA[1,yrs,,sex,fleet,1])
  naa_total <- apply(om_values$NAA[,yrs,,sex,1], 3, mean)

  data <- list(
    aa_r = aa_r,
    naa_total = naa_total,
    m = m
  )

  pars <- list(
    ln_f_ar = rep(log(0.1), length(aa_r))
  )

  adfun <- MakeADFun(SPoRC:::cmb(nll_srv_aa, data), pars, random = NULL)
  optim <- stats::nlminb(adfun$par, adfun$fn, adfun$gr,
                         control =  list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15))
  sd <- sdreport(adfun)
  rep <- adfun$report(adfun$env$last.par.best)

  par(mfrow = c(3,2), cex = 0.8)
  # look at selex
  plot(rep$f_ar / max(rep$f_ar), type = 'l', ylim = c(0,1), xlab = 'Age', ylab = 'Selex', main = 'BS')

  # compare SrvIAA
  plot(rep$aa_est, type = 'l', lty = 2, ylab = 'SrvIAA', xlab = 'Age', main = 'BS')
  lines(data$aa_r, col = 'red')

  # AI
  aa_r <- colSums(om_values$SrvIAA[2,yrs,,sex,fleet,1])
  naa_total <- apply(om_values$NAA[,yrs,,sex,1], 3, sum)

  data <- list(
    aa_r = aa_r,
    naa_total = naa_total,
    m = m
  )

  pars <- list(
    ln_f_ar = rep(log(0.1), length(aa_r))
  )

  adfun <- MakeADFun(SPoRC:::cmb(nll_srv_aa, data), pars, random = NULL)
  optim <- stats::nlminb(adfun$par, adfun$fn, adfun$gr,
                         control =  list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15))
  sd <- sdreport(adfun)
  rep <- adfun$report(adfun$env$last.par.best)

  # look at selex
  plot(rep$f_ar / max(rep$f_ar), type = 'l', ylim = c(0,1), xlab = 'Age', ylab = 'Selex', main = 'AI')

  # compare SrvIAA
  plot(rep$aa_est, type = 'l', lty = 2, ylab = 'SrvIAA', xlab = 'Age', main = 'AI')
  lines(data$aa_r, col = 'red')

  # GOA
  aa_r <- apply(om_values$SrvIAA[3:5,,,sex,fleet,1], 3, mean)
  naa_total <- apply(om_values$NAA[,-66,,sex,1], 3, mean)

  data <- list(
    aa_r = aa_r,
    naa_total = naa_total,
    m = m
  )

  pars <- list(
    ln_f_ar = rep(log(0.1), length(aa_r))
  )

  adfun <- MakeADFun(SPoRC:::cmb(nll_srv_aa, data), pars, random = NULL)
  optim <- stats::nlminb(adfun$par, adfun$fn, adfun$gr,
                         control =  list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15))
  sd <- sdreport(adfun)
  rep <- adfun$report(adfun$env$last.par.best)

  # look at selex
  plot(rep$f_ar / max(rep$f_ar), type = 'l', ylim = c(0,1), xlab = 'Age', ylab = 'Selex', main = 'GOA')

  # compare SrvIAA
  plot(rep$aa_est, type = 'l', lty = 2, ylab = 'SrvIAA', xlab = 'Age', main = 'GOA')
  lines(data$aa_r, col = 'red')
}

# Load packages
library(RTMB)
library(here)

# Read in model output
om_values <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))
m <- unique(as.vector(om_values$natmort))

# Fishery
get_faa_sel(1, 1, 1:19) # fleet 1, sex 1
get_faa_sel(1, 2, 1:19) # fleet 1, sex 2
get_faa_sel(2, 1, 1:19) # fleet 2, sex 1
get_faa_sel(2, 2, 1:19) # fleet 2, sex 2

get_faa_sel(1, 1, 20:39) # fleet 1, sex 1
get_faa_sel(1, 2, 20:39) # fleet 1, sex 2
get_faa_sel(2, 1, 20:39) # fleet 2, sex 1
get_faa_sel(2, 2, 20:39) # fleet 2, sex 2

get_faa_sel(1, 1, 40:65) # fleet 1, sex 1
get_faa_sel(1, 2, 40:65) # fleet 1, sex 2
get_faa_sel(2, 1, 40:65) # fleet 2, sex 1
get_faa_sel(2, 2, 40:65) # fleet 2, sex 2

# Survey
get_srvaa_sel(1, 1, 1:19) # fleet 1, sex 1
get_srvaa_sel(1, 2, 1:19) # fleet 1, sex 2
get_srvaa_sel(3, 1, 1:19) # fleet 2, sex 1
get_srvaa_sel(3, 2, 1:19) # fleet 2, sex 2

get_srvaa_sel(1, 1, 20:39) # fleet 1, sex 1
get_srvaa_sel(1, 2, 20:39) # fleet 1, sex 2
get_srvaa_sel(3, 1, 20:39) # fleet 2, sex 1
get_srvaa_sel(3, 2, 20:39) # fleet 2, sex 2

get_srvaa_sel(1, 1, 40:65) # fleet 1, sex 1
get_srvaa_sel(1, 2, 40:65) # fleet 1, sex 2
get_srvaa_sel(3, 1, 40:65) # fleet 2, sex 1
get_srvaa_sel(3, 2, 40:65) # fleet 2, sex 2
