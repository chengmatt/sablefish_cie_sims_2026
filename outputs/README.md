# outputs

Results, mirroring `R/`: `cross_test/` holds fits without management feedback, `mse_results/` holds closed-loop results. Everything is regenerable from the scripts in `R/`. All `.RDS` here are gitignored (and `mse_results/` entirely), so only figures and csvs under `cross_test/` are tracked.

## cross_test/

- `spatial_no_block_scenarios/`: the working results.
  - `spt_rand_OM_*.RDS`: simulated OM datasets from `crosstest_sims.R` (`lowsamp`/`highsamp` = composition sample size; `tvmove_sens` = time-varying movement).
  - `<em>_crosstest_*.RDS`: fitted EMs. Suffixes: `_francis` = Francis-reweighted; `_srvchng_{historical,current,all}` = survey design eras; `_movesense` = movement sensitivity; `_firstpart`..`_sixthpart` = chunked FAA batches; `_sens_dt` = the GOA crash experiment; `*_aggregation_test` = aggregation scheme experiments.
  - `*_retro*.csv`: retrospective (Mohn's rho) summaries.
- `figs/`: comparison figures from the `cross_testing_comparisons_*` scripts. Suffix `_historical` or `_srvchng` says which comparison script made it; prefixes `agg_` vs `spt_` distinguish aggregate and spatially resolved quantities.
- `old/`: blocked-era results (scripts in `R/old/`).

## mse_results/

- `spatial_no_block_scenarios/`: closed-loop output as `<em>_<scenario>_<i>.RDS`, where scenario is `base`, `highregimerec`, or `lowregimerec` and `i` is a simulation batch index. `faa_*_tb_*` = FAA EM with time-blocked survey selectivity. The five-region EM saves one file per scenario (`five_region_<scenario>.RDS`).
- `figs/`: comparison figures from the `mse_comparisons_*` scripts, suffixed by scenario. Gifs are the animated OM vs EM tracking plots.
- `old/`: blocked-era results.
