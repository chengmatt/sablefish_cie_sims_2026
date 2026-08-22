# R/functions

Shared machinery sourced by every run script. Nothing in here is run directly.

## EM builders

Each of these constructs a complete SPoRC input list (and fit) for one estimation model, from data simulated by the five-region OM:

- `single_region_em.R`: single-region EM. Data are aggregated across regions before fitting.
- `faa_em.R`: single-region, fleet-as-area (FAA) EM. Regions enter as separate fleets rather than being pooled.
- `three_rg_em.R`: three-region EM (Bering Sea, Aleutian Islands, Gulf of Alaska), aggregated down from the five-region OM.
- `five_rg_em.R`: five-region EM matching the OM structure ("true" spatial configuration).

## mse_functions.R

Everything else, grouped roughly as:

- Closed-loop runners: `run_single_rg_closedloop_parallel()`, `run_faa_closedloop_parallel()`, `run_three_rg_closedloop_parallel()`, `run_five_rg_closedloop_parallel()` and their per-simulation `_i` workers. These drive the annual assess/project/harvest feedback loop and are what the `R/mse` scripts call.
- Data aggregation: `agg_data_to_single_rg()`, `agg_data_to_faa()`, `agg_data_to_three_rg()` collapse five-region simulated data to each EM's resolution; `add_aggregated_obj_to_simenv()` attaches the result to a simulation environment.
- Data availability switches: `single_region_use_indicators()`, `faa_use_indicators()`, `threerg_use_indicators()`, `fiverg_use_indicators()` control which data sources an EM sees in a given year.
- Harvest policy: `HCR_threshold()` (NPFMC-style threshold rule), `solve_multifleet_F()` (F that achieves target catch across fleets), `get_population_projection()` (2-year projection for ABC), `get_single_region_survey_apportionment()` and `get_three_rg_survey_apportionment()` (survey-based regional apportionment of the ABC).
- `build_selex_prior()`: selectivity priors for blocked/fleeted configurations.
