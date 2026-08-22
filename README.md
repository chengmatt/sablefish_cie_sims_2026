# sablefish_cie_sims_2026

Simulation testing (cross tests and MSE closed loops) for the 2026 sablefish CIE review, built on the SPoRC package.

## Layout

Each major folder has its own README describing what every script and result file is: `R/functions`, `R/cross_test`, `R/mse`, `R/scratch`, `R/old`, `data`, and `outputs`.

Scripts under `R/` and results under `outputs/` mirror each other. Within each stage, work is grouped by scenario family, and within a family by estimation model (EM). The current family is `spatial_no_block_scenarios` (conditioned OM without selectivity blocks).

```
R/
  functions/       shared EM builders and MSE utilities, sourced by every run script
    mse_functions.R, single_region_em.R, three_rg_em.R, five_rg_em.R, faa_em.R
  cross_test/
    spatial_no_block_scenarios/
      crosstest_sims.R                 simulates the OM datasets (spt_rand_OM_*.RDS) used by all EMs
      single_region_em/ three_rg_em/ five_rg_em/ faa_em/
                                       per-EM cross test runs (baseline and srvchng variants)
      cross_testing_comparisons_*.R    comparison figures across EMs
      cross_test_aggregation_test.R    multi-EM aggregation experiment
      em_cross_test_goa_crash.R        multi-EM GOA crash experiment
  mse/
    spatial_no_block_scenarios/
      single_region_em/ three_rg_em/ five_rg_em/ faa_em/
                                       per-EM closed loops (base, highregime, lowregime recruitment)
      mse_comparisons_*.R              comparison figures across EMs, one per recruitment scenario
  scratch/         active exploratory scripts only
  old/             superseded blocked-era scenario scripts (spatial_blocked_scenarios,
                   spatial_blocked_selex_scenarios); outputs from that era live in
                   outputs/*/old/. Frozen: paths inside are not kept up to date.

data/
  spatial_outputs/ conditioned OM model results read by crosstest_sims.R and the MSE scripts

outputs/
  cross_test/
    spatial_no_block_scenarios/        simulated OM datasets and EM cross test results
    figs/                              cross test comparison figures (_historical, _srvchng)
    old/                               blocked-era results
  mse_results/                         (gitignored)
    spatial_no_block_scenarios/        closed loop results per EM and recruitment scenario
    figs/                              MSE comparison figures (_base, _highregimerec, _lowregimerec)
    old/                               blocked-era results
```

## Conventions

- Every script resolves paths with `here()` from the project root, so scripts can be moved within `R/` without breaking.
- Run scripts source what they need from `R/functions/`; nothing sources a run script.
- Result files are named `<em>_<scenario>[_<part>].RDS`; figures carry the scenario as a suffix.
- Superseded work is parked in `old/` next to where its replacement lives, never deleted.
