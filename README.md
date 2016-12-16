# Master's Thesis in Finance

![moleman](moleman.jpg)

The source data is in `data/source` and downloaded from Ken French's website

## Generate Data and Estimate Models

1. `do_load.R` is used to translate source data into weekly returns
2. `do_garch_fitting.R` is used to fit the battery of GARCH models
3. `do_garch_choose.R` picks the best GARCH model according to our criteria
   and spits out uniforms from each marginal distribution used to fit the
   copula
4. `do_fit_copula2.R` fits copula models to uniform residuals.

This will generate model estimates and filtered data.

## Simulating

To simulate distributions, run `do_simulate_copula_full.R`

## MV and CDB optimization

See `do_mv.R` and `do_cdb2.R`

Use `do_cdb_eval.R` to get results of the cool stuff.

## Other stuff

The rest is mostly robustness checks and outputting graphs.
