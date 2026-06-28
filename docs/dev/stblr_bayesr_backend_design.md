# ST-BLR BayesR Backend Design

## Executive Summary

BayesR support exists in two different forms:

- `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp` is a BED
  scheduled, multi-chain BayesR backend.
- `src/st_sbayesrc_omp_csr.cpp`, `src/st_sbayesrc_omp_csr_annot.cpp`, and
  `src/st_sbayresrc_omp_csr.cpp` are CSR SBayesRC-style samplers with
  BayesR-like mixture components, but they are annotation-oriented and do not
  follow the recently harmonized ST-BLR BayesC chain architecture.

The first harmonization step has been implemented for the existing BED BayesR
chain backend: standard `dm` now means `P(component > 0)`, the old posterior
mean component index is preserved separately, and standard `bm`/`dm` chain
summaries are returned. The backend remains lower-level and experimental
because no public R wrapper currently calls it.

Plain exact CSR BayesR now has a first internal experimental backend using
`src/st_cpg_omp_csr.cpp` as the architectural template and the SBayesRC CSR code
only as a source of mixture-update math. Scheduled CSR BayesR should come after
the exact CSR BayesR path is validated. LD-swap remains postponed for BayesR
until the Metropolis-Hastings state exchange is explicitly designed for both
effect and component assignment.

The detailed exact CSR BayesR design and implementation status are in
`docs/dev/stblr_csr_bayesr_design.md`. That path keeps CSR BayesR separate from
SBayesRC-style annotation models, uses an internal experimental R helper before
a public API, and preserves the standard chain-output convention used by the
harmonized BayesC backends.

## Current BayesR Backend Inventory

### C++ Files

Explicit BayesR/BayesR-like native files:

- `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp`
  - Packed BED, marker-scheduled, multi-chain BayesR backend.
  - Uses a K-component mixture prior with `c[0] = 0`, non-null `c[k] > 0`,
    component state `d[j]`, Dirichlet-updated `pi`, and per-component marker
    posterior probabilities.
- `src/st_sbayesrc_omp_csr.cpp`
  - CSR SBayesRC-style overlapping-annotation sampler with fixed gamma mixture
    components and learned annotation coefficients.
- `src/st_sbayesrc_omp_csr_annot.cpp`
  - CSR SBayesRC-like annotation-class comparison sampler.
- `src/st_sbayresrc_omp_csr.cpp`
  - CSR SBayesRC-style annotation-class sampler. The filename appears to have a
    transposed spelling, but it is exported as `stblr_cpg_omp_csr_sbayesrc_annot1`.

There is now a plain experimental `src/st_cpg_omp_csr_bayesr.cpp`, exported as
`stblr_cpg_omp_csr_bayesr()`. There is no scheduled
`src/st_cpg_omp_csr_scheduled_bayesr.cpp` today.

### R Wrappers

Rcpp-generated wrappers exist in `R/RcppExports.R` for:

- `stblr_cpg_omp_bed_marker_scheduled_chains_bayesr()`
- `stblr_cpg_omp_csr_sbayesrc()`
- `stblr_cpg_omp_csr_sbayesrc_annot()`
- `stblr_cpg_omp_csr_sbayesrc_annot1()`

Public user-facing R support exists for the SBayesRC generic CSR path:

- `stblr_csr_sbayesrc_generic()` in `R/stblr-csr-sbayesrc.R`
- SBayesRC helper functions in `R/sbayesrc-helpers.R`

No current user-facing R function calls
`stblr_cpg_omp_bed_marker_scheduled_chains_bayesr()`. The BED BayesR symbol is
exported through Rcpp but appears orphaned from public R workflows.

### Tests

Existing tests cover:

- SBayesRC helper functions and `format_sbayesrc_csr_fit()`
- BayesC BED scheduled chains
- BayesC CSR chain summaries and LD-swap
- backend consistency checks
- fine-mapping extraction with optional chain summaries

No test currently calls
`stblr_cpg_omp_bed_marker_scheduled_chains_bayesr()` directly, and no public BED
BayesR wrapper test exists.

### Compile Status

This inspection did not run a native compile because this is a design-only pass
and native compilation can update generated artifacts. The repository contains
`src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.o`, and the function is
present in both `R/RcppExports.R` and `src/RcppExports.cpp`, which suggests it
has compiled before. Current compile status should still be re-verified during
the implementation pass.

### Status Classification

`stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp` is best treated as
experimental, exported lower-level code that is currently unused or orphaned
from public R entry points. It now follows the standard marker-level `bm`/`dm`
and chain-summary convention, but it still does not have a public wrapper. The
SBayesRC CSR path is active public code, but it is not the same as a plain
BayesR ST-BLR backend.

## Relationship Between BayesC and BayesR Backends

BayesC uses one null/non-null indicator and one active variance component:

- `dm`: posterior inclusion probability.
- `bm`: posterior mean marker effect.
- `pi`: active-marker probability.

BayesR generalizes this to a zero component plus multiple non-zero normal
components:

- component `0`: null, marker effect fixed to zero.
- components `1..K-1`: non-null effect distributions with component-specific
  variance multipliers.
- `dm` should mean `P(component > 0)`, not posterior mean component index.
- `bm` remains the posterior mean marker effect.
- marker component posterior probabilities should be exposed separately as
  `comp_prob`, `pip_k`, or another explicitly component-specific object.

The existing SBayesRC formatter already uses `dm = P(component > 0)` and
`comp_prob` for marker-by-component posterior probabilities. Future plain
BayesR paths should follow that convention.

## Existing BED BayesR Scheduled Chains Status

The BED BayesR backend is exported as:

```cpp
stblr_cpg_omp_bed_marker_scheduled_chains_bayesr(...)
```

Its arguments mirror the BayesC BED scheduled chains backend, with BayesR
mixture additions:

- `pi`: initial K-component mixture probabilities.
- `c`: K-component variance multipliers, with `c[0] = 0`.
- `alpha`: Dirichlet prior shape vector for `pi`.

It supports:

- `nchains`
- `ncores`
- scheduled marker updates
- compact BED block reading
- `return_wy`
- `return_r`
- CPO/log-CPO diagnostics
- chain aggregation for core posterior summaries
- BayesR-specific `final_pi`, `mean_pi`, and marker component probabilities

It does not currently support:

- `keep_chains`
- `chain_seeds`

The job mapping is chain-major:

```cpp
job = chain * nt + trait
```

The seed formula is:

```cpp
seed + 1000003 * (trait + 1) + 9176 * (chain + 1)
```

This matches the default multi-chain seed formula used by the BayesC BED chain
backend, but it does not use `src/st_chain_utils.h`.

The result struct is `ChainResultBayesR`, analogous to
`ChainResultSTScheduled`. Per-chain `bm` is available at aggregation time as
`r.bm`. Per-chain component probabilities are available as `r.pip_k`. Standard
`r.dm` is now accumulated directly as `component > 0`; the old posterior mean
component index is preserved in `r.component_mean`.

Current return layout is 30 slots:

- slots `0:21` follow the BayesC BED scheduled chain base layout.
- slot `22` is BayesR-specific `pip_k`, flattened as `K * m` per trait.
- slots `23:28` are standard chain summaries:
  `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`, `dm_max`.
- slot `29` is BayesR-specific `component_mean`, the old posterior mean
  component index.

In R's one-based indexing, BayesR `pip_k` is `fit[[23]]`, chain summaries are
`fit[[24:29]]`, and posterior mean component index is `fit[[30]]`.

## Current R-side Status

The BED BayesR path remains experimental and lower-level. No public wrapper
routes user fits to `stblr_cpg_omp_bed_marker_scheduled_chains_bayesr()` yet.
The stable R-side path for tests and future wrappers is internal:

- `.format_stblr_bayesr_fit()` formats the raw Rcpp return.
- `.stblr_bed_marker_bayesr_experimental()` is an internal helper that calls
  the Rcpp backend, formats the result, and records minimal BayesR metadata.
- No new public API has been added for BED BayesR.

The formatted BayesR fit convention is:

- `dm`: standard non-null PIP, `P(component > 0)`.
- `bm`: posterior mean marker effect.
- `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, `bm_max`: standard
  chain-summary matrices with the same dimensions and dimnames as `dm`/`bm`.
- `comp_prob`: named list by trait; each element is a marker-by-component
  matrix. `component_0` is the null component, so `dm = 1 - component_0`.
- `dm_component_mean`: posterior mean component index, preserving the old
  BayesR-specific quantity separately from standard `dm`.

## Required Output Convention

For all BayesR ST-BLR fits, the standard fields should be:

- `fit$bm`: marker posterior mean effect.
- `fit$dm`: marker posterior probability of any non-zero component.
- `fit$bm_sd`, `fit$bm_min`, `fit$bm_max`: across-chain summaries of per-chain
  `bm`.
- `fit$dm_sd`, `fit$dm_min`, `fit$dm_max`: across-chain summaries of per-chain
  `dm`.

BayesR-specific fields should be separate from these standard fields:

- `fit$comp_prob`: marker-by-component posterior probabilities, preferably a
  list of matrices indexed by trait.
- `fit$pi` and `fit$pim`: final and posterior mean mixture weights with one
  column per mixture component.
- `fit$ncomp`: posterior mean number of markers per component, when tracked.
- `fit$component`: final component assignment, zero-based.
- `fit$gamma` or `fit$mixture_var`: component variance multipliers.

For the existing BED BayesR backend, rename or format slot `22` as
`comp_prob_raw` and expose a formatted `comp_prob` object rather than treating
the slot as an iteration trace named `pis`.

## BED BayesR Harmonization Design

Minimal alignment has been added without changing sampler math:

1. During aggregation, create `bm_sd_mat`, `bm_min_mat`, `bm_max_mat`,
   `dm_sd_mat`, `dm_min_mat`, and `dm_max_mat`, exactly as in
   `src/st_cpg_omp_individual_scheduled_chains.cpp`.
2. Use `r.bm` for effect summaries and per-chain `r.dm` for standard non-null
   PIP summaries.
3. Compute means first, then sample SD with denominator `nchains - 1` when
   `nchains > 1`.
4. For `nchains == 1`, SD remains zero and min/max equal the single-chain
   values.
5. Append slots zero-based `23:28` as:
   - `23`: `bm_sd`
   - `24`: `bm_min`
   - `25`: `bm_max`
   - `26`: `dm_sd`
   - `27`: `dm_min`
   - `28`: `dm_max`

`.format_stblr_fit()` already detects slots one-based `24:29` and names them
as the six standard chain summaries. The internal
`.format_stblr_bayesr_fit()` wrapper additionally formats BayesR component
probabilities from slot `23`, enforces standard `dm = 1 - P(component 0)`, and
exposes slot `30` as `dm_component_mean`.

`check_stblr_backend_consistency()` should work unchanged for BayesR fits if:

- `fit$dm` is marker-by-trait and in `[0, 1]`.
- `fit$bm` has the same dimensions and dimnames.
- `fit$input$nchains` is present.
- the six chain summary matrices are present for `nchains > 1`.

## Existing CSR BayesR Status

There is now a first plain exact summary-stat CSR BayesR backend:

- `src/st_cpg_omp_csr_bayesr.cpp`
- `stblr_cpg_omp_csr_bayesr(...)`
- internal formatter `.format_stblr_csr_bayesr_fit()`
- internal helper `.stblr_csr_bayesr_experimental()`

This path is experimental and internal. It supports exact CSR updates,
`nchains`, `chain_seeds`, standard chain summaries, `comp_prob`, and
`dm_component_mean`. It rejects `keep_chains = TRUE` and LD-swap. Scheduled CSR
BayesR remains future work.

The residual variance update investigation is documented in
`docs/dev/stblr_csr_bayesr_design.md` under "Residual Variance Update and Prior
Scaling". The experimental CSR BayesR wrapper now uses the same sparse active
probability convention as the SBayesRC/BayesC CSR helpers, but `updateE = TRUE`
remains disabled until real-data native diagnostics validate the residual
scale.

Other existing CSR mixture code is SBayesRC-style:

- `stblr_cpg_omp_csr_sbayesrc()` models marker-specific mixture probabilities
  from overlapping annotations through probit stick-breaking.
- `stblr_cpg_omp_csr_sbayesrc_annot()` and
  `stblr_cpg_omp_csr_sbayesrc_annot1()` use annotation-class-specific mixture
  probabilities.

Reusable ideas exist in these files:

- categorical component update by log posterior probability
- component variance multiplier handling
- `dm = component > 0`
- `comp_prob` accumulation
- Dirichlet updates for mixture weights in class-specific variants

Reusable infrastructure exists in:

- `src/st_csr_common.h` for CSR likelihood and variance helpers.
- `src/st_chain_utils.h` for chain seed, task, and thread helpers.

`src/cpg_samplers.h` only exposes BayesC-style sampler declarations and is not
a reusable BayesR kernel source today.

## Exact CSR BayesR Design

Implemented file name:

```text
src/st_cpg_omp_csr_bayesr.cpp
```

Implemented C++ export:

```cpp
stblr_cpg_omp_csr_bayesr(...)
```

Internal R-facing helper:

```r
.stblr_csr_bayesr_experimental(...)
```

The implementation starts from `src/st_cpg_omp_csr.cpp` for architecture, not
from the SBayesRC files. It uses the SBayesRC files for the BayesR component
update and output conventions.

Arguments to reuse from BayesC CSR:

- `wy`, `ww`, `yy`
- `b_init`
- `r_init`, `use_r_init`, `rebuild_r_before_updateE`
- `ld_prefix`
- `B`, `E`, `ssb_prior`, `sse_prior`
- `nub`, `nue`, `updateB`, `updateE`, `adjE`
- `n`, `nit`, `nburn`, `nthin`
- `ncores`, `seed`, `nchains`, `keep_chains`, `chain_seeds`

BayesR-specific arguments:

- `comp_init`: initial zero-based component states.
- `use_comp_init`
- `mixture_var` or `gamma`: component variance multipliers, with first value
  zero and all others positive.
- `pi`: initial K-component mixture weights.
- `alpha`: Dirichlet prior shape vector for global mixture weights.
- `updatePi`: whether to update global mixture weights.

Initial LD-swap support:

- Do not support LD-swap in the first exact CSR BayesR implementation.
- Reject `updateLDswap = TRUE` if the argument is present.

Standard outputs:

- base `bm`, `dm`, `wy`, `r`, `b`, `component`, marker index
- `vbs`, `vgs`, `ves`, `vle`, `vld`
- final and posterior mean mixture weights
- `bm_sd`, `bm_min`, `bm_max`
- `dm_sd`, `dm_min`, `dm_max`

BayesR-specific outputs:

- `comp_prob`: marker-by-component posterior probability per trait.
- `dm_component_mean`: posterior mean component index per marker and trait.
- `ncomp`: posterior mean component counts per trait.
- optional `component_trace` or `pi_trace` only if the storage cost and format
  are explicitly justified.

Chain aggregation:

- Use `src/st_chain_utils.h` with trait-major task mapping:
  `trait = task / nchains`, `chain = task % nchains`.
- Preserve single-chain seed behavior:
  `seed + 1000003 * (trait + 1)`.
- Use `chain_seeds[chain] + 1000003 * (trait + 1)` when explicit chain seeds
  are supplied.
- For `nchains > 1` without explicit seeds, use
  `seed + 1000003 * (trait + 1) + 9176 * (chain + 1)`.
- Average standard traces across chains as the BayesC CSR backend does.
- Aggregate `comp_prob` across chains by arithmetic mean.
- Aggregate `pi`/`pim` across chains by arithmetic mean.

Extractor compatibility:

- `extract_stblr_finemap_loci()` needs only `fit$dm`, `fit$bm`, marker names,
  and optional chain summary matrices. It should work unchanged if the BayesR
  fit follows the standard field convention.

## Scheduled CSR BayesR Design

Recommended file name:

```text
src/st_cpg_omp_csr_scheduled_bayesr.cpp
```

Recommended export:

```cpp
stblr_cpg_omp_csr_scheduled_bayesr(...)
```

This should be implemented after exact CSR BayesR. The scheduled backend should
share small BayesR-specific helpers with exact CSR BayesR where practical, but
should initially remain a separate backend to avoid broad refactors.

Required initial support:

- `nchains` from the start.
- standard chain summary outputs from the start.
- same BayesR mixture arguments and output convention as exact CSR BayesR.
- `keep_chains = FALSE` only, matching scheduled CSR BayesC for now.
- no LD-swap initially.

Scheduling controls should be interpreted against total non-null probability:

- Candidate activation should use `dm = P(component > 0)`.
- A marker is a likely null candidate when total non-null probability is below
  the scheduled threshold.
- Max non-null component probability can be tracked diagnostically, but should
  not replace total non-null probability as the default scheduler criterion.
- Expected effect contribution is more model-dependent and should be postponed.

## LD-Swap Considerations

For BayesC, an LD-swap move exchanges inclusion/effect state between correlated
markers. For BayesR, the state includes both effect and mixture component.

The key design question is whether a swap exchanges:

- only non-null status and effect,
- effect plus component assignment,
- or a proposal that changes component assignment with a separate acceptance
  ratio term.

The correct Metropolis-Hastings ratio depends on the component prior
probabilities, component variance multipliers, and any marker-specific mixture
probabilities. Because this is not a trivial extension of BayesC LD-swap,
BayesR LD-swap should be postponed until exact non-swap CSR BayesR is validated.

## Recommended Implementation Order

1. Existing BED BayesR chain output has been harmonized with standard `bm`/`dm`
   chain summaries and a small internal BayesR formatter for component
   probabilities.
2. Focused tests now cover formatted BED BayesR lower-level output, component
   probabilities, backend consistency, and extractor compatibility.
3. A first internal exact CSR BayesR backend now exists without LD-swap, using
   regular CSR BayesC as the architecture template and SBayesRC CSR code for
   mixture math.
4. Validate exact CSR BayesR with tiny native CSR fixture tests and explicit
   `chain_seeds`.
5. Decide whether compact `keep_chains` output and a public
   `stblr_csr_bayesr()` wrapper are needed.
6. Implement scheduled CSR BayesR after exact CSR BayesR is stable.
7. Design BayesR LD-swap later as a separate statistical-method change.

This corresponds to option E: create small BayesR-specific helpers and keep
backends separate initially. A broad BayesC/BayesR sampler refactor should wait
until the BayesR contracts and tests are stable.

## Testing Plan

Existing BED BayesR chains:

- Direct lower-level smoke test with `nchains = 1`.
- Direct lower-level smoke test with `nchains = 2`.
- Verify standard fields after formatting:
  `dm`, `bm`, `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, `bm_max`.
- Check dimensions, dimnames, finite values, and `min <= mean <= max`.
- Run `check_stblr_backend_consistency()` with `require_chain_summaries = TRUE`
  for multi-chain fits.

BayesR-specific outputs:

- `pi`/`pim` rows or columns sum to one, depending on chosen formatter shape.
- `comp_prob[[trait]]` has dimension `m x K`.
- Component probabilities are in `[0, 1]` and row sums are near one.
- `dm` equals the row sum of non-null component probabilities.
- `ncomp` equals column sums of `comp_prob`, within tolerance.

Future exact CSR BayesR:

- Small synthetic CSR fixture.
- Single-chain output.
- Multi-chain output.
- Chain summaries.
- `chain_seeds` determinism.
- `keep_chains` if supported from the start.
- extractor compatibility through `extract_stblr_finemap_loci()`.

Future scheduled CSR BayesR:

- scheduled multi-chain output.
- standard chain summaries.
- scheduled controls and candidate thresholds.
- explicit rejection of LD-swap and `keep_chains = TRUE` if unsupported.

## Suggested Next Implementation Prompt

The first BayesR harmonization step has been implemented:

- `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp` appends standard
  chain summary slots `23:28` using per-chain `r.bm` and standard non-null
  `r.dm`.
- R formatting exposes BayesR slot `22` as component probabilities and slot
  `29` as `dm_component_mean`.
- Focused tests cover formatted fields, component probabilities,
  `check_stblr_backend_consistency()`, and extractor compatibility.

Suggested next implementation prompt:

- Add a small experimental R wrapper for the BED BayesR backend if there is a
  clear user workflow for direct BED BayesR fits.
- Otherwise, implement exact CSR BayesR without LD-swap as a separate task.
