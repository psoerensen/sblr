# BLR output and raw schemas

Every formatted fit has `family`, `model`, `operator`, `input`, `data`,
`model_parameters`, `diagnostics`, `convergence`, `convergence_traces`,
`chains`, and `memory_estimate`. Required structural fields remain present and
may be `NULL`; model-specific scientific fields are present only where defined.

Common scientific fields are marker × trait `bm`, `dm`, `b_final`, and
`d_final`; iteration × trait `vbs`, `vgs`, `ves`, `vle`, and `vld`; and
model-specific `beta_final`, `component_final`, and
`component_probabilities`. MT covariance fields are `cov_b_mean`,
`cov_g_mean`, `cov_e_mean`, `cov_b_final`, `cov_g_final`, and `cov_e_final`.
Probability fields are `pi_final`, `pi_mean`, and `pi_trace`; BayesRC pattern
and annotation parameters live under `model_parameters`.

`b_final` and `bm` are effective effects. `d_final` is binary activity for
BayesC and a mixture state only where the documented model output explicitly
defines it; `dm` is posterior non-null probability. No legacy field aliases
are part of the formatted contract.

Both `stblr_raw` and `mtblr_raw` retain raw schema version 1. Model semantics
version 2 is mandatory. Raw objects are named, validated, and converted by one
family-specific formatter; positional fallback is unsupported.
