# Exact CSR BayesR Backend Design

## Executive Summary

Plain summary-statistics CSR BayesR is supported through an exact,
non-scheduled backend:

```text
src/st_cpg_omp_csr_bayesr.cpp
```

exporting:

```cpp
stblr_cpg_omp_csr_bayesr(...)
```

The implementation uses `src/st_cpg_omp_csr.cpp` as the primary template for
exact CSR likelihood updates, residual handling, chains, seed handling, and
standard output aggregation. Existing SBayesRC-style CSR files were used only
as references for categorical mixture math and component probability
accumulation. The BED BayesR backend defines the BayesR output contract:
standard `dm = P(component > 0)`, marker-by-component `comp_prob`, and
`dm_component_mean` for posterior mean component index. Full BED BayesR
slot layout, CPO/log-CPO diagnostics detail, and the original harmonization
design steps are archived in
`docs/dev/archive/stblr_bayesr_backend_design.md`; the current feature/slot
summary for `bed_bayesr` lives in
`docs/dev/stblr_backend_computation_inventory.md`.

The supported public R interfaces are `stblr_csr_bayesr()` and
`stblr_csr(method = "bayesr")`. Exact CSR BayesR now supports optional
active/null LD-swap/MH moves that relocate the full
`(component, b)` state from an active marker to a null LD friend. Scheduled CSR
BayesR remains future work.

The BayesR LD-swap/MH design is recorded separately in
`docs/dev/stblr_bayesr_ldswap_design.md`. The implemented first scope relocates
the full BayesR `(component, b)` state from an active marker to a null LD
friend, with the BayesC likelihood/proposal-ratio machinery reused only after
accounting for BayesR component priors. Active/active swaps and
annotation-specific swap prior terms are not implemented. Fixed global
`selection_s` supplies marker-specific variance terms for the active/null move.

## Implementation Status

Implemented in this pass:

- `src/st_cpg_omp_csr_bayesr.cpp`
- Rcpp export `stblr_cpg_omp_csr_bayesr(...)`
- internal formatter `.format_stblr_csr_bayesr_fit()`
- public helper `stblr_csr_bayesr()`
- high-level `stblr_csr(method = "bayesr")` dispatch
- compatibility alias `.stblr_csr_bayesr_experimental()`
- formatter and native smoke tests in `tests/testthat/test-bayesr-csr-backend.R`

The backend supports exact CSR updates, `nchains`, `chain_seeds`,
`keep_chains = TRUE`, standard `bm`/`dm` summaries, component probabilities,
posterior mean component index, variance traces, mixture-weight summaries, and
strict `updateE = TRUE` diagnostics. It also supports optional active/null
LD-swap/MH diagnostics. It does not implement scheduled CSR BayesR or
active/active BayesR LD-swap.

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
- `check_stblr_consistency()` for standard `dm`/`bm` and chain summary
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

Implemented file:

```text
src/st_cpg_omp_csr_bayesr.cpp
```

Implemented export:

```cpp
// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_bayesr(...)
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
  int updateE_start = 0,
  int updateE_every = 1,
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1
)
```

`updateLDswap` enables active/null full-state BayesR LD-swap. The tuning
arguments follow exact CSR BayesC conventions.

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
- `updateE_start >= 0`;
- `updateE_every >= 1`;
- `ld_swap_prob` and `ld_swap_r2` are in `[0, 1]`;
- `ld_swap_max_friends` is positive;
- `ld_swap_moves` is non-negative.

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

The implemented backend returns a named `Rcpp::List` instead of extending the
BayesC positional slot protocol. This keeps the BayesR path separate from
BayesC formatting and avoids collisions with LD-swap diagnostics.

Standard marker fields are marker-by-trait matrices:

- `bm`
- `dm`, standard `P(component > 0)`
- `bm_sd`, `bm_min`, `bm_max`
- `dm_sd`, `dm_min`, `dm_max`

BayesR-specific marker fields are:

- `component`, final zero-based component index, marker by trait
- `comp_prob`, list by trait; each element is marker by component
- `dm_component_mean`, posterior mean component index, marker by trait

The formatted component fields can be summarized without printing marker-level
matrices:

```r
summarise_components(fit)
```

This helper is an R-side diagnostic summary only. It does not implement
scheduled CSR BayesR, LD-swap, sampler changes, or native backend changes.

Variance and mixture summaries are:

- `vbs`, `vgs`, `ves`, `vle`, `vld`, trace by trait
- `covb`, `covg`, `cove`, `vb`, `vg`, `ve`
- `pi`, final mixture weights, trait by component
- `pim`, posterior mean mixture weights, trait by component
- `ncomp`, posterior mean component counts, trait by component
- `mixture_var`

Additional state fields are retained for diagnostics and formatter
compatibility:

- `wy`
- `r`
- `b`
- `marker_index`

`keep_chains = TRUE` is rejected in this first implementation. Compact
chain-level details can be added later as named fields such as `dm_chains`,
`bm_chains`, and `component_mean_chains`.

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

## LD-Swap Scope

Exact CSR BayesR supports active/null LD-swap. The move relocates the full
BayesR state `(component, b)` from an active marker to a null LD friend, updates
`r = X'y - X'Xb`, and accepts/rejects with the same likelihood/proposal-ratio
structure used by exact CSR BayesC. Under the plain CSR BayesR model without
fixed `selection_s`, global `pi` and global `mixture_var` make
component/effect prior terms cancel for this full-state relocation. With fixed
`selection_s`, the marker-specific prior variance terms enter the MH ratio.

Not implemented:

- active/active BayesR swaps;
- scheduled CSR BayesR swaps;
- annotation-specific swap prior terms.

## Residual Variance Update and Prior Scaling

This section records the inspection triggered by invalid residual variance
updates in the exact CSR BayesR backend when `updateE = TRUE`.

SBayesRC prior setup:

- The public SBayesRC wrapper is `stblr_csr_sbayesrc_generic()` in
  `R/stblr-csr-sbayesrc.R`.
- Its default mixture variance multipliers are `gamma = c(0, 0.01, 0.1, 1)`.
  These are scale factors, not absolute effect variances.
- Native SBayesRC kernels use active prior variance
  `vb_t * gamma[k]`, with `gamma[0] = 0` forcing the null component effect to
  exactly zero.
- The default active probability is `pi_init = 0.001`, distributed evenly over
  active components by `make_sbayesrc_alpha_init()` unless
  `active_comp_weights` is supplied.
- Marker probabilities in the overlapping-annotation SBayesRC path are
  produced by probit stick-breaking: `Phi(A_i alpha_j)`. With an intercept
  column, the intercept is initialized to encode the sparse baseline.
- `B` and `ssb_prior` are initialized as
  `(vy * h2) / (m * pi_vb_init)` and
  `((nub - 2) / nub) * (vy * h2) / (m * pi_prior_mean)`.
- `E` and `sse_prior` are initialized as `vy * (1 - h2)` and
  `((nue - 2) / nue) * vy * (1 - h2)`.
- Defaults are `nub = 4`, `nue = 4`, `h2 = 0.5`, `updateAlpha = TRUE`,
  `updateB = TRUE`, `updateE = TRUE`, and `adjE = 0.9`.
- `comp_init` defaults to all null components. `r_init` defaults to a rebuild
  from `b_init`, unless the caller explicitly supplies `use_r_init = TRUE`.
- `rebuild_r_before_updateE` exists but defaults to `FALSE`.

SBayesRC native update path:

- `src/st_sbayesrc_omp_csr.cpp`,
  `src/st_sbayesrc_omp_csr_annot.cpp`, and
  `src/st_sbayresrc_omp_csr.cpp` all call the shared
  `sampleE_ST_csr()` from `src/st_csr_common.h` when `updateE = TRUE`.
- The arguments are the current marker count, `nue`, mutable `ve_t`, full
  `b_t`, `wy_t`, current `r_t`, diagonal `sse_prior`, `yy_t`, `n[t]`, and the
  chain RNG.
- `sampleE_ST_csr()` uses the same identity as the BayesR path:
  `SSE = y'y - b'X'y - b'r`, with `r = X'y - X'Xb`.
- SBayesRC can optionally rebuild `r` immediately before `sampleE_ST_csr()`,
  but it does not do so by default.
- The same active effect posterior form is used in SBayesRC and plain CSR
  BayesR: component likelihood terms include
  `0.5 * log(vei / (vei + ww_i * vb * scale_k))` and the corresponding
  quadratic score term. The conditional effect draw uses
  `lhs = ww_i + vei / (vb * scale_k)`.
- `sampleB_*` updates use `sum b_i^2 / scale_k` over active components and do
  not divide by LD score or marker-specific annotation variance.
- Annotation effects in SBayesRC change mixture probabilities, not the
  component variance formula.

Comparison with plain CSR BayesR:

| Quantity | SBayesRC CSR | Plain CSR BayesR |
| --- | --- | --- |
| Mixture grid | `c(0, 0.01, 0.1, 1)` | `c(0, 0.01, 0.1, 1)` |
| Scale convention | `vb * gamma[k]` | `vb * mixture_var[k]` |
| Initial active probability | `0.001` | now `0.001`; previously `0.05` |
| Initial component state | all null unless supplied | all null unless supplied |
| Mixture probability model | marker-specific probit stick-breaking or class-specific rows | global trait/chain `pi` |
| Pi prior/update | annotation alpha update or class Dirichlet variants | global `Dirichlet(alpha + counts)` |
| Default alpha prior | encoded by sparse alpha init; class variants are caller-supplied | now centered on `pi` with strength `5e5`; previously flat `rep(1, K)` |
| `B`/`ssb_prior` scale | divided by `m * pi_*` | divided by `m * sum(pi[-1])` |
| `E`/`sse_prior` scale | `vy * (1 - h2)` convention | same convention |
| update order | markers, optional alpha/Pi, B, E, diagnostics | markers, B, E, Pi, diagnostics |
| `updateE` call | shared `sampleE_ST_csr()` | shared `sampleE_ST_csr()` with extra diagnostic precheck |
| `adjE` | `vei = ve + adjE * vg` after diagnostics | same |
| Null component | effect exactly zero | effect exactly zero |

The failing diagnostic had `nonzero_components=2894` out of `m=5000` at
iteration 3. That is not expected under the SBayesRC-style sparse prior:
`pi_init = 0.001` implies roughly 5 active markers initially for `m=5000`.
The earlier plain BayesR defaults implied 250 active markers
before any likelihood contribution, and the flat Dirichlet prior did little to
pull the global mixture weights back toward sparsity after dense early
assignments. A dense state can overfit the summary statistics enough that
`y'y - b'X'y - b'r` becomes negative, making the inverse-chi-square residual
scale invalid even when the SSE identity itself is correct.

Most likely causes evaluated:

- A. Mixture variance scale bug: not supported by code inspection. Both paths
  use `vb * scale_k`, and both default grids are scale factors.
- B. Prior probability bug: supported. The earlier wrapper used a much
  denser default `pi` and a flat Dirichlet prior.
- C. Effect variance bug: not supported by code inspection. `sampleB` and
  marker updates use the same `vb * scale_k` convention.
- D. Component likelihood bug: not supported by the inspected formulas; the
  determinant and quadratic terms are present and match SBayesRC.
- E. Initialization bug: partly mitigated by all-null component initialization,
  but the previous dense `pi` made early activation too permissive.
- F. updateE incompatibility: still possible on real data. Even with sparse
  priors, `updateE = TRUE` should fail strictly with diagnostics if a dataset
  enters an invalid residual-scale state.

Minimal fix applied:

- `stblr_csr_bayesr()` defaults to total active
  `pi = 0.001` instead of `0.05`.
- When `alpha` is omitted, it is now centered on the normalized `pi` with total
  strength `5e5`, matching the sparse CSR prior convention instead of using a
  flat `rep(1, K)`.

`updateE = TRUE` is enabled in `stblr_csr_bayesr()`. By default it
updates residual variance from zero-based iteration `0` and every iteration
after that (`updateE_start = NULL` maps to `0`, `updateE_every = 1`). These are
supported controls; callers can delay or thin residual-variance
updates with `updateE_start` and `updateE_every` when investigating a dataset.

The native backend rebuilds `r = X'y - X'Xb` immediately before each residual
variance update, verifies that null-component markers have exactly zero effect,
computes SSE diagnostics, and then calls the shared `sampleE_ST_csr()` only if
the residual scale is finite and positive. It does not clamp or silently skip
invalid SSE states. Successful fits return a lightweight
`updateE_diagnostics` matrix with per trait-chain update counts, minimum SSE,
minimum residual scale, maximum active-component count, maximum absolute
effect, and maximum absolute fitted quadratic term.

## R Wrapper Design

Use the public exact CSR BayesR helper directly:

```r
stblr_csr_bayesr(...)
```

The same exact backend is also available through the high-level CSR wrapper:

```r
stblr_csr(..., method = "bayesr")
```

`stblr_csr()` keeps BayesC as the default (`method = "bayesc"`). Under
`method = "bayesr"`, it performs only high-level argument dispatch and lets
`stblr_csr_bayesr()` own BayesR defaults and backend validation. BayesC-specific
prior arguments such as `pi_init`, `pi_prior_a`, and `pi_prior_b` are rejected
when explicitly supplied for BayesR; callers should use BayesR `pi` and
`alpha` instead. `scheduled = TRUE` errors clearly because scheduled CSR
BayesR is not implemented.

The internal helper should:

- accept `stats`, `ld_prefix`, priors, `mixture_var`, `pi`, `alpha`,
  `nchains`, `keep_chains`, `chain_seeds`, and MCMC controls;
- create `b_init`, `comp_init`, prior matrices, and marker/trait names using
  the same conventions as `stblr_csr()`;
- accept BayesR active/null LD-swap arguments and pass them to the CSR backend;
- support `keep_chains = TRUE` summaries and chain-level LD-swap diagnostics;
- call `stblr_cpg_omp_csr_bayesr()`;
- call a BayesR CSR formatter;
- set metadata:

```r
fit$input$model = "bayesr"
fit$input$backend = "csr_bayesr"
fit$input$scheduled = FALSE
fit$input$nchains = nchains
fit$input$keep_chains = FALSE
fit$input$mixture_var = mixture_var
fit$input$alpha = alpha
```

## Formatter Design

Add a CSR-capable BayesR formatter rather than forcing the BED formatter to
handle every layout immediately:

```r
.format_stblr_csr_bayesr_fit(...)
```

Implemented behavior:

1. Read named fields from the `Rcpp::List` returned by
   `stblr_cpg_omp_csr_bayesr()`.
2. Format standard marker fields and chain summaries as marker-by-trait
   matrices with marker and trait dimnames.
3. Format `comp_prob` as a named list by trait, each `m x K` with columns
   `component_0`, `component_1`, ...
4. Set `fit$dm <- 1 - comp_prob[[trait]][, "component_0"]`.
5. Expose `dm_component_mean`, `pi`, `pim`, `ncomp`, and `mixture_var`.
6. Validate `comp_prob` values are finite, within `[0, 1]`, and row sums are
   approximately one.

Longer term, `.format_stblr_bayesr_fit()` and `.format_stblr_csr_bayesr_fit()`
can share a small internal helper for component-probability formatting. Keep
BayesC `.format_stblr_fit()` stable.

Formatted CSR BayesR fits should pass:

```r
check_stblr_consistency(fit, require_chain_summaries = TRUE)
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
- `keep_chains = TRUE` rejects until compact chain output is implemented;
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
- `updateLDswap = TRUE` returns trait-level and, when `keep_chains = TRUE`,
  chain-level LD-swap diagnostics;
- scheduled mode errors or routes elsewhere until scheduled CSR BayesR exists.

Avoid expensive MCMC. Use small synthetic CSR fixtures and short chains.

## Suggested Next Implementation Prompt

Next work should stay outside the supported exact CSR BayesR surface unless it
is explicitly scoped.

- Implement scheduled CSR BayesR as a separate backend.
- Extend BayesR LD-swap only in separately scoped tasks, such as active/active
  swaps or marker-specific prior terms.
- Extend public docs or examples if real-data smoke tests identify missing
  usage guidance.
- Keep scheduled CSR BayesR and LD-swap extensions, such as active/active swaps
  and marker-specific prior terms, as separate future tasks.
