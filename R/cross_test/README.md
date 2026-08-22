# R/cross_test

Self tests and cross tests: fit each EM to data simulated from the conditioned five-region OM, without management feedback, to check for bias before trusting the EMs inside the MSE.

Run order within `spatial_no_block_scenarios/`:

1. `crosstest_sims.R` simulates the OM datasets every EM consumes. Writes the `spt_rand_OM_*.RDS` files in `outputs/cross_test/spatial_no_block_scenarios/`: `lowsamp` and `highsamp` (low/high composition sample sizes), `lowsamp_tvmove_sens` (time-varying movement sensitivity), `lowsamp_dt`, and `highsamp_aggtest`.
2. Per-EM fit scripts, one folder per EM:
   - `single_region_em/sglrg_em_crosstest.R`: single-region EM self test across survey/fishery weighting schemes, plus movement-sensitivity refits. Writes `single_region_crosstest_lowsamp*.RDS` and `single_region_crosstest_movesense*.RDS`.
   - `single_region_em/sgl_rg_em_crosstest_srvchg.R`: single-region EM on an extended time series with a change in survey design. Writes the `single_region_crosstest_srvchng_{historical,current,all}.RDS` set (historical design, current/future design, all years).
   - `three_rg_em/three_rg_crosstest.R` and `three_rg_crosstest_srvchng.R`: same two ideas for the three-region EM.
   - `five_rg_em/five_rg_crosstest.R`: five-region EM cross test (`five_region_crosstest_{low,high}samp.RDS`).
   - `faa_em/faa_em_cross_test.R`: FAA EM self test. Large runs are saved in chunks (`faa_crosstest_lowsamp_firstpart.RDS` through `_sixthpart`), with `finalFAA` marking the final FAA configuration. `faa_crosstest_srvchng.R` is the survey-design-change version.
3. Cross-EM experiments at this level:
   - `cross_test_aggregation_test.R`: how aggregation scheme choices (which error sources, age lag on/off, simulate-then-aggregate vs aggregate-then-simulate) affect the single-region and FAA EMs. Writes the `sgl_error*_aggregation_test.RDS` family.
   - `em_cross_test_goa_crash.R`: crashes the GOA sub-stock in the OM and checks how the single-region and FAA EMs respond. Writes the `*_sens_dt*.RDS` files.
4. Comparison figures: `cross_testing_comparisons_historical.R` and `cross_testing_comparisons_srvchng.R` read the fitted results and write the pngs in `outputs/cross_test/figs/` (suffix `_historical` or `_srvchng` matches the script).

A `_francis` suffix anywhere means the fit was Francis-reweighted before saving.
