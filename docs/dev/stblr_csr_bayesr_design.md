# Exact CSR BayesR Backend Design

## Executive Summary

Plain summary-statistics CSR BayesR does not exist in the package today. The
recommended implementation is a new exact, non-scheduled backend:

```text
src/st_cpg_omp_csr_bayesr.cpp
```

exporting:

```cpp
stblr_cpg_omp_csr_bayesr(...)
```

The implementation should use `src/st_cpg_omp_csr.cpp` as the primary template
for exact CSR likelihood updates, residual handling, chains, seed handling,
`keep_chains`, and standard output layout. Existing SBayesRC-style CSR files
should be used only as references for categorical mixture math and component
probability accumulation. The BED BayesR backend should define the BayesR output
contract: standard `dm = P(component > 0)`, marker-by-component `comp_prob`, and
`dm_component_mean` for posterior mean component index.

Do not include LD-swap in the first CSR BayesR version. BayesR swap moves need a
separate Metropolis-Hastings design because the state includes both marker
effect and mixture component.

## Code Inventory

### Existing Exact CSR BayesC

`src/st_cpg_omp_csr.cpp` is the best architectural template. It provides:

- exact full-sweep summary-statistics CSR marker updates;
- shared LD reading through `read_and_build_st_ld_csr()`;
- residual rebuild/update logic through `st_csr_common.h`;
- multi-chain task mapping through `st_chain_utils.h`;
- `nchains`, `keep_chains`, and `chain_seeds`;
- standard chain summaries:
  `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`, `dm_max`;
- compact per-chain output when `keep_chains = TRUE`;
- LD-swap, which should not be copied into the first BayesR backend.

`src/st_cpg_omp_csr_scheduled.cpp` is useful only as a later reference for a
scheduled BayesR backend. It supports `nchains` and standard summaries but
rejects `keep_chains = TRUE` and has scheduler-specific logic that should not
be part of the exact BayesR implementation.

### Existing BED BayesR

`src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp` defines the current
BayesR convention:

- component `0` is null;
- components `1..K-1` are non-null;
- `pi` is a K-vector updated by `Dirichlet(alpha + counts)`;
- standard `dm` is posterior `P(component > 0)`;
- `pip_k` is accumulated per component;
- `component_mean` preserves posterior mean component index;
- chain summaries are standard `bm`/`dm` summaries.

This file is useful for output and Dirichlet-update conventions, but its marker
likelihood is BED/individual-data specific and scheduled. It should not be the
primary CSR implementation template.

### Existing SBayesRC-Style CSR

The BayesR-like CSR files are:

- `src/st_sbayesrc_omp_csr.cpp`
- `src/st_sbayesrc_omp_csr_annot.cpp`
- `src/st_sbayresrc_omp_csr.cpp`

They use the same summary-statistics CSR residual representation as
`st_cpg_omp_csr.cpp`: `wy`, `ww`, `yy`, `ld_prefix`, `b`, `r`, and sparse LD
through `st_csr_common.h`. They also contain useful BayesR-like pieces:

- categorical component sampling by log posterior probability;
- validation that `gamma[0]` or `mixture_var[0]` is zero and non-null scales are
  positive;
- `sampleB_*` logic using `b_i^2 / gamma_k` or `b_i^2 / mixture_var_k`;
- `dm = component > 0`;
- component probability accumulation into marker-by-component matrices;
- component-count summaries;
- Dirichlet updates for class-specific mixture weights in annotation-class
  variants.

They are not plain BayesR backends. They include annotation-dependent
probabilities, annotation-class priors, learned alpha coefficients, or
class-specific mixture weights. They also do not follow the recently harmonized
multi-chain CSR architecture.

### Shared Headers and R Contracts

`src/st_chain_utils.h` should be reused for task and seed helpers:

- `stblr_task_trait(task, nchains)`
- `stblr_task_chain(task, nchains)`
- `stblr_num_chain_tasks(nt, nchains)`
- `stblr_trait_seed(seed, trait)`
- `stblr_chain_seed(seed, trait, chain)`
- `stblr_seed_with_chain_base(chain_seed, trait)`

`src/st_csr_common.h` should be reused for:

- `STLDCSR`
- `read_and_build_st_ld_csr()`
- `rebuild_residual_st_csr()`
- `sampleE_ST_csr()`
- `computeG_ST_csr()`

`src/cpg_samplers.h` does not currently expose reusable BayesR kernels.

R formatting and downstream compatibility depend on:

- `.format_stblr_fit()` for standard ST-BLR fields and chain summaries;
- `.format_stblr_bayesr_fit()` for BayesR component probabilities in the BED
  layout;
- `format_sbayesrc_csr_fit()` for existing SBayesRC marker-by-component
  formatting;
- `check_stblr_backend_consistency()` for standard `dm`/`bm` and chain summary
  validation;
- `extract_stblr_finemap_loci()` using `fit$dm`, `fit$bm`, and optional
  `dm_*`/`bm_*` summaries.

## Reusable Pieces

The best source for each implementation part is:

- CSR residual updates: `src/st_cpg_omp_csr.cpp` and `src/st_csr_common.h`.
- BayesR mixture sampling: `sampleBeta_SBayesRC_ST_csr()` from SBayesRC files,
  adapted to a global K-vector `pi` rather than marker/class-specific rows.
- Effect variance update: SBayesRC `sampleB_*` logic, adapted to the exact CSR
  chain structure and plain `mixture_var`/`gamma`.
- Chain aggregation: `src/st_cpg_omp_csr.cpp`.
- Seed/task mapping: `src/st_chain_utils.h`.
- Output convention: BED BayesR plus BayesC CSR chain slot conventions.
- Component probability tracking: SBayesRC `comp_prob_accum` and BED BayesR
  `pip_k` convention.

BayesR update kernels are not currently separated into reusable headers. The
first implementation should copy/adapt the smallest required helpers into the
new file. After the plain CSR backend is tested, common BayesR helpers can be
extracted if duplication becomes a real maintenance problem.

## Target Plain CSR BayesR Model

For trait `t` and marker `i`, let component `z_it` be in `0..K-1`.

- Component `0` is the exact null: `b_it = 0`.
- Components `1..K-1` are non-null normal components.
- Conditional on component `k > 0`:

```text
b_it | z_it = k ~ N(0, vb_t * mixture_var[k])
```

- `mixture_var[0] = 0`.
- `mixture_var[k] > 0` for `k > 0`.
- Recommended default scale grid should follow existing package convention:
  `c(0, 0.01, 0.1, 1)` for BayesR-like CSR work. If a smaller GCTB-style grid
  is desired, expose it as a user-supplied option rather than hard-coding it.
- Each trait and chain has a global K-vector `pi_t`.
- If `updatePi = TRUE`, update `pi_t` using:

```text
pi_t ~ Dirichlet(alpha + component_counts_t)
```

- If `updatePi = FALSE`, keep the supplied `pi` fixed.
- Residual variance `ve_t` should use the same `sampleE_ST_csr()` update as
  exact CSR BayesC.
- Marker effect variance `vb_t` should be global per trait, with active
  component scales applied by `mixture_var[k]`.
- Standard `dm_i,t` is:

```text
P(z_it > 0)
```

- Standard `bm_i,t` is posterior mean marker effect.
- `comp_prob[[trait]][i, k]` is posterior `P(z_it = k)`.
- `dm_component_mean_i,t` is posterior mean component index and should be
  returned as a BayesR-specific diagnostic.

## C++ Backend Design

### File and Export

Recommended file:

```text
src/st_cpg_omp_csr_bayesr.cpp
```

Recommended export:

```cpp
// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_bayesr(...)
```

### Rcpp Signature

The signature should mirror exact CSR BayesC where possible:

```cpp
stblr_cpg_omp_csr_bayesr(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> comp_init,
  bool use_comp_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<double> pi,
  std::vector<double> mixture_var,
  std::vector<double> alpha,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds,
  bool updateLDswap = false
)
```

`updateLDswap` should be accepted only if useful for wrapper symmetry, but the
first implementation should reject `TRUE` with a clear error.

### Argument Validation

Validate:

- dimensions of `wy`, `ww`, `yy`, `b_init`, priors, `B`, `E`, and `n` exactly
  as in `st_cpg_omp_csr.cpp`;
- `mixture_var` length at least 2;
- `mixture_var[0] == 0`;
- all non-null `mixture_var` values are finite and positive;
- `pi` and `alpha` have length `K`;
- `pi` is finite, non-negative, and has positive sum; normalize internally if
  the sum is not exactly one;
- `alpha` is finite and positive;
- `comp_init` has length `nt` and each entry length `m` when enabled;
- component states are integers in `0..K-1`;
- `nchains >= 1`;
- `chain_seeds` is empty or length `nchains`;
- `keep_chains` is supported from the start for exact CSR BayesR;
- `updateLDswap = TRUE` is rejected.

### Marker Update

Adapt SBayesRC categorical update logic to a global `pi_t`:

1. Remove the old effect contribution from residual:

```text
r_i += ww_i * b_old
r_j += x_ij * b_old for LD neighbors
```

2. Compute log posterior mass for all components.
3. Component `0` has effect zero.
4. For component `k > 0`, use prior variance `vb_t * mixture_var[k]`.
5. Sample the component using a normalized softmax / inverse CDF.
6. If `k > 0`, sample the marker effect from the conditional normal.
7. Add the new effect contribution back to the residual.

The exact formulas should be taken from `sampleBeta_SBayesRC_ST_csr()` and
checked against `sampleBetaC_ST_csr()` for residual update sign and `ww`/`r`
interpretation.

### Variance and Pi Updates

- Use `sampleE_ST_csr()` unchanged for `ve_t`.
- Adapt SBayesRC `sampleB_*`:

```text
ssb_scaled = sum_{i:z_i>0} b_i^2 / mixture_var[z_i]
vb_t ~ scaled inverse-chi-square using ssb_scaled + nub * ssb_prior_t
```

- Use BED BayesR-style global Dirichlet update:

```text
pi_t[k] = Gamma(alpha[k] + count_k, 1) / sum_k Gamma(...)
```

### Chain Behavior

Use the exact CSR BayesC mapping:

```text
task = trait * nchains + chain
trait = stblr_task_trait(task, nchains)
chain = stblr_task_chain(task, nchains)
```

Seed behavior should match BayesC CSR:

- explicit `chain_seeds`: `chain_seeds[chain] + 1000003 * (trait + 1)`;
- single chain with no explicit seeds: `seed + 1000003 * (trait + 1)`;
- multiple chains with no explicit seeds:
  `seed + 1000003 * (trait + 1) + 9176 * (chain + 1)`.

OpenMP should run over `nt * nchains` tasks with
`stblr_num_threads_for_tasks(ncores, ntasks)`.

### Output Layout

Use the standard BayesC CSR base layout where possible, but do not reuse slot 22
for LD-swap diagnostics because LD-swap is unsupported initially. A safe
BayesR-specific layout is:

- `0`: `bm`
- `1`: `dm`, standard `P(component > 0)`
- `2`: `wy`
- `3`: `r`
- `4`: `b`
- `5`: `component`, final zero-based component index as double
- `6`: marker index
- `7`: `vbs`
- `8`: `vgs`
- `9`: `ves`
- `10`: `covb`
- `11`: `covg`
- `12`: `cove`
- `13`: `vb`
- `14`: `vg`
- `15`: `ve`
- `16`: `pi`, final mixture weights, length `K` per trait
- `17`: `pim`, posterior mean mixture weights, length `K` per trait
- `18`: diagnostics/CPO-compatible slot if retained
- `19`: marker diagnostics if retained
- `20`: `vle`
- `21`: `vld`
- `22`: `comp_prob`, flattened marker-major `m x K` per trait
- `23`: `bm_sd`
- `24`: `bm_min`
- `25`: `bm_max`
- `26`: `dm_sd`
- `27`: `dm_min`
- `28`: `dm_max`
- `29`: `dm_component_mean`
- `30`: `ncomp`, posterior mean component counts, length `K` per trait

If `keep_chains = TRUE`, append compact chain details after these fixed slots:

- `31`: chain-level `dm`, chain-major `nchains * m`
- `32`: chain-level `bm`, chain-major `nchains * m`
- optionally `33`: chain-level `component_mean`, chain-major `nchains * m`

Do not overload slots `23:28`; these should remain the standard chain summaries
recognized by `.format_stblr_fit()`.

### Aggregation Rules

For each trait and marker:

- `bm` = arithmetic mean of per-chain posterior mean effects;
- `dm` = arithmetic mean of per-chain `P(component > 0)`;
- `comp_prob` = arithmetic mean of per-chain component probabilities;
- `dm_component_mean` = arithmetic mean of per-chain posterior mean component
  index;
- `bm_sd`/`dm_sd` = sample SD across chain-level summaries with denominator
  `nchains - 1` when `nchains > 1`, otherwise zero;
- min/max are per-marker minima/maxima of chain-level summaries;
- `pi` and `pim` = arithmetic means across chains;
- traces can follow BayesC CSR convention by averaging chain traces.

## LD-Swap Decision

LD-swap should be deferred for the first exact CSR BayesR backend.

Reasons:

- BayesR state is `(effect, component)`, not only `(effect, included)`.
- Swapping only included status/effect is not a full state exchange.
- Swapping effect and component together changes component prior and variance
  terms.
- A correct Metropolis-Hastings ratio must account for proposal probabilities,
  component probabilities, variance multipliers, and any future marker-specific
  priors.
- Plain exact BayesR CSR can be validated against SBayesRC/BED conventions
  faster without this method change.

Initial behavior should be:

```text
if (updateLDswap) stop("BayesR CSR LD-swap is not yet supported.")
```

## R Wrapper Design

Add an internal experimental helper first:

```r
.stblr_csr_bayesr_experimental(...)
```

Do not extend `stblr_csr()` with `model = "bayesr"` in the first implementation.
The current public `stblr_csr()` API is BayesC-oriented and has no general
`model`/`prior` convention. A public `stblr_csr_bayesr()` wrapper can be added
after native behavior, formatting, and tests are stable.

The internal helper should:

- accept `stats`, `ld_prefix`, priors, `mixture_var`, `pi`, `alpha`,
  `nchains`, `keep_chains`, `chain_seeds`, and MCMC controls;
- create `b_init`, `comp_init`, prior matrices, and marker/trait names using
  the same conventions as `stblr_csr()`;
- reject LD-swap arguments;
- call `stblr_cpg_omp_csr_bayesr()`;
- call a BayesR CSR formatter;
- set metadata:

```r
fit$input$model = "bayesr"
fit$input$backend = "csr_bayesr"
fit$input$scheduled = FALSE
fit$input$nchains = nchains
fit$input$keep_chains = keep_chains
fit$input$mixture_var = mixture_var
fit$input$alpha = alpha
```

## Formatter Design

Add a CSR-capable BayesR formatter rather than forcing the BED formatter to
handle every layout immediately:

```r
.format_stblr_csr_bayesr_fit(...)
```

Recommended behavior:

1. Call `.format_stblr_fit()` for slots `0:21` and standard chain summaries
   `23:28`.
2. Parse raw slot `23` in R one-based indexing only if the C++ layout places
   `comp_prob` at zero-based slot `22`.
3. Format `comp_prob` as a named list by trait, each `m x K` with columns
   `component_0`, `component_1`, ...
4. Set `fit$dm <- 1 - comp_prob[[trait]][, "component_0"]`.
5. Parse `dm_component_mean` from the chosen slot.
6. Parse `ncomp` if returned.
7. Validate `comp_prob` values are finite, within `[0, 1]`, and row sums are
   approximately one.

Longer term, `.format_stblr_bayesr_fit()` and `.format_stblr_csr_bayesr_fit()`
can share a small internal helper for component-probability formatting. Keep
BayesC `.format_stblr_fit()` stable.

Formatted CSR BayesR fits should pass:

```r
check_stblr_backend_consistency(fit, require_chain_summaries = TRUE)
```

and should work with:

```r
extract_stblr_finemap_loci(fit = fit, Glist = Glist, credible_sets = FALSE)
```

because the extractor relies on standard `dm`, `bm`, and optional chain summary
fields.

## Testing Plan

Formatter-level tests:

- mocked raw C++ output with one trait and three markers;
- `K = 3` or `K = 4`;
- `dm = 1 - component_0`;
- `comp_prob` dimensions are `m x K`;
- component probabilities are finite, in `[0, 1]`, and row-sum to one;
- `dm_component_mean` is exposed;
- standard chain summaries have matching dimensions/dimnames.

Tiny CSR fixture tests:

- `nchains = 1`;
- `nchains = 2`;
- `keep_chains = FALSE`;
- `keep_chains = TRUE` for exact backend;
- finite `bm`, `dm`, variance traces, and mixture weights;
- `dm_min <= dm <= dm_max`;
- `bm_min <= bm <= bm_max`;
- single-chain SD fields are zero;
- backend consistency checker passes.

BayesR-specific checks:

- `pi` and `pim` are length `K` per trait and sum approximately one;
- `ncomp` equals column sums of `comp_prob` within tolerance;
- `dm` equals row sum of non-null component probabilities;
- final component indices are in `0..K-1`.

Compatibility tests:

- `extract_stblr_finemap_loci(..., credible_sets = FALSE)` propagates
  `dm_sd` to `pip_sd`, `dm_min` to `pip_min`, `dm_max` to `pip_max`, and
  `bm_sd`/`bm_min`/`bm_max`;
- unsupported `updateLDswap = TRUE` errors clearly;
- scheduled mode errors or routes elsewhere until scheduled CSR BayesR exists.

Avoid expensive MCMC. Use small synthetic CSR fixtures and short chains.

## Suggested Implementation Prompt

Implement a new exact summary-statistics CSR BayesR backend only.

- Add `src/st_cpg_omp_csr_bayesr.cpp`.
- Use `src/st_cpg_omp_csr.cpp` as the implementation template.
- Adapt BayesR categorical component sampling and `sampleB` scaling from
  SBayesRC CSR code.
- Use `st_chain_utils.h` and `st_csr_common.h`.
- Support `nchains`, `keep_chains`, and `chain_seeds` from the start.
- Return standard `bm`/`dm` and standard chain summaries.
- Return BayesR `comp_prob`, `dm_component_mean`, and `ncomp`.
- Reject LD-swap.
- Add an internal `.stblr_csr_bayesr_experimental()` helper and a CSR BayesR
  formatter.
- Add formatter-level and tiny CSR tests.
- Do not implement scheduled CSR BayesR or LD-swap in this step.
