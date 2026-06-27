# ST-BLR CSR Multi-Chain Design

## Executive Summary

`stblr_csr()` now supports multiple independent MCMC chains per trait in the
regular summary-statistic CSR backend when `scheduled = FALSE` and in the
scheduled summary-statistic CSR backend when `scheduled = TRUE`. The scheduled
CSR backend remains LD-swap-free.

The repository also has useful multi-chain machinery in the packed-BED
scheduled backend. That implementation runs `nchains * nt` chain-trait jobs,
derives per-chain seeds, aggregates posterior summaries across chains, and
returns the existing 23-slot fit layout extended with the same SD/min/max
marker-summary fields as the regular CSR backend.

The recommended first implementation has been completed for the regular CSR
backend. The next step is not to copy that logic into every backend, but to
introduce small shared chain helpers and align chain-capable outputs in focused
follow-up work.

## Implementation Status

Implemented for the regular summary-statistic CSR backend:

- `stblr_csr(..., scheduled = FALSE)` accepts `nchains`, `keep_chains`, and
  vector `chain_seeds`.
- regular CSR C++ execution uses independent trait-by-chain OpenMP tasks.
- `nchains = 1` preserves the previous per-trait seed rule.
- `nchains > 1` aggregates `bm`, `dm`, final summaries, averaged traces, and
  LD-swap diagnostics across chains.
- multi-chain fits expose `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, and
  `bm_max`.
- `keep_chains = TRUE` returns compact per-chain `dm`, `bm`, and LD-swap
  diagnostics.

Implemented for the scheduled summary-statistic CSR backend:

- `stblr_csr(..., scheduled = TRUE)` accepts `nchains` and vector
  `chain_seeds`.
- scheduled CSR C++ execution uses independent trait-by-chain OpenMP tasks.
- `nchains = 1` preserves the previous per-trait seed rule.
- `nchains > 1` aggregates `bm`, `dm`, final summaries, and averaged traces
  across chains.
- scheduled multi-chain fits expose `dm_sd`, `dm_min`, `dm_max`, `bm_sd`,
  `bm_min`, and `bm_max`.
- `keep_chains = TRUE` remains unsupported for scheduled CSR.

Not implemented in this pass:

- LD-swap in the scheduled CSR backend;
- per-chain full trace retention or convergence diagnostics.

## Current Regular CSR Sampler

The exported C++ function is `stblr_cpg_omp_csr()` in
`src/st_cpg_omp_csr.cpp`. It is exported to R as `stblr_cpg_omp_csr()` in
`R/RcppExports.R` and registered as `_sblr_stblr_cpg_omp_csr` in
`src/RcppExports.cpp`.

The user-facing R function is `stblr_csr()` in `R/sparse_ld_bed_helper.R`.
When `scheduled = FALSE`, it builds common arguments and calls
`stblr_cpg_omp_csr()`.

The regular CSR sampler now:

- reads and builds one shared flat CSR LD object from `ld_prefix`;
- optionally builds LD-swap friend lists once;
- allocates one row per trait-chain task for `bm`, `dm`, final state, traces,
  and LD-swap diagnostics;
- runs `#pragma omp parallel for` over `task = 0 .. nt * nchains - 1`;
- maps `trait = task / nchains` and `chain = task % nchains`;
- sets `nthreads` from `ncores`, bounded by the number of tasks;
- preserves the old single-chain seed rule
  `seed + 1000003 * (t + 1)`;
- uses `seed + 1000003 * (t + 1) + 9176 * (chain + 1)` for default
  multi-chain seeds;
- uses `chain_seeds[chain] + 1000003 * (t + 1)` when chain seeds are supplied;
- aggregates posterior summaries across independent chains.

The main returned object is a `std::vector<std::vector<std::vector<double>>>`.
The R formatter `.format_stblr_fit()` names and reshapes it. Important base
slots:

- slot 1 / index 0: `bm`, posterior mean effects;
- slot 2 / index 1: `dm`, posterior inclusion probabilities;
- slots 8 to 10 / indices 7 to 9: `vbs`, `vgs`, `ves` traces;
- slots 21 to 22 / indices 20 to 21: `vle`, `vld` traces;
- slot 23 / index 22: LD-swap diagnostics for the regular CSR backend.

When `nchains > 1`, slots 24 to 29 / indices 23 to 28 contain `bm_sd`,
`bm_min`, `bm_max`, `dm_sd`, `dm_min`, and `dm_max`. When `keep_chains = TRUE`,
slots 30 to 32 / indices 29 to 31 contain compact chain-major `dm`, `bm`, and
LD-swap diagnostics.

`dm`, `bm`, `vbs`, `ves`, `pis`, `vgs`, `vld`, and `vle` are accumulated from
per-task local outputs after the OpenMP loop. `ld_swap` diagnostics are created
by summing per-task attempted and accepted counts within each trait, then
placed in result slot 23 with a four-value marker that `.format_stblr_fit()`
uses to recognize LD-swap output.

## Current Scheduled CSR Sampler

The exported scheduled summary-statistic CSR function is
`stblr_cpg_omp_csr_scheduled()` in `src/st_cpg_omp_csr_scheduled.cpp`. It is
exported to R as `stblr_cpg_omp_csr_scheduled()` in `R/RcppExports.R` and
registered as `_sblr_stblr_cpg_omp_csr_scheduled` in `src/RcppExports.cpp`.

`stblr_csr()` routes to it when `scheduled = TRUE`.

The scheduled CSR sampler:

- uses the same summary-statistic CSR inputs as the regular backend;
- adds marker scheduling controls such as `full_sweep_every`,
  `null_skip_base`, `candidate_threshold`, and LD-neighbor wakeup controls;
- parallelizes over `nt * nchains` trait-chain tasks;
- maps `trait = task / nchains` and `chain = task % nchains`;
- sets `nthreads = max(1, min(ncores, nt * nchains))`;
- preserves the old single-chain seed rule
  `seed + 1000003 * (t + 1)`;
- uses `seed + 1000003 * (t + 1) + 9176 * (chain + 1)` for default
  multi-chain seeds;
- uses `chain_seeds[chain] + 1000003 * (t + 1)` when chain seeds are supplied;
- returns chain SD/min/max summaries for `dm` and `bm`;
- returns the same base fit layout, with slot 23 used for the pi trace rather
  than LD-swap diagnostics.

The scheduled CSR sampler does not currently support LD-swap. The R wrapper
blocks this before calling C++:

```r
if (isTRUE(scheduled) && isTRUE(updateLDswap)) {
  stop("updateLDswap is currently implemented only for scheduled = FALSE.")
}
```

There is no `updateLDswap` argument in `stblr_cpg_omp_csr_scheduled()`, and the
scheduled C++ file does not contain LD-swap friend building or swap proposal
logic.

## Existing Chain-Related Code

The relevant existing multi-chain implementation is
`stblr_cpg_omp_bed_marker_scheduled_chains()` in
`src/st_cpg_omp_individual_scheduled_chains.cpp`, called by
`stblr_bed_marker()` when the packed-BED chains backend is selected.

That implementation already has useful machinery:

- an `nchains` argument;
- a `ChainResultSTScheduled` per chain-trait job;
- `njobs = nchains * nt`;
- job mapping `chain = job / nt`, `t = job % nt`;
- OpenMP parallelism over jobs with `nthreads = min(ncores, njobs)`;
- per-chain seeds:
  `seed + 1000003 * (t + 1) + 9176 * (chain + 1)`;
- posterior summaries averaged across chains before returning;
- averaged `bm`, `dm`, final `b/d`, variance traces, pi traces, and final
  variance/pi summaries;
- basic runtime diagnostics such as mean/max seconds and mean sample count.

What it does not currently provide:

- retained compact per-chain summaries in the R fit;
- trace concatenation across chains;
- convergence diagnostics such as R-hat or ESS;
- LD-swap diagnostics, because this BED scheduled chains backend is not the
  CSR LD-swap backend.

There are repeated commented historical copies in
`src/st_cpg_omp_individual_scheduled_chains.cpp`; only the active exported
function should be used as a template.

## Posterior-Only Fine Mapping Compatibility

`extract_stblr_finemap_loci()` in `R/extract-stblr-finemap-loci.R` already
looks for optional fields:

- `fit$dm_sd`, `fit$dm_min`, `fit$dm_max`;
- `fit$bm_sd`, `fit$bm_min`, `fit$bm_max`.

If present, it maps them into `pip_sd`, `pip_min`, `pip_max`, `bm_sd`,
`bm_min`, and `bm_max` in the extracted locus tables. If absent, those columns
are returned as `NA`.

No inspected CSR C++ backend currently emits these fields. The multi-chain CSR
implementation should therefore return:

- `fit$dm`: mean PIP across chains;
- `fit$bm`: mean posterior effect across chains;
- `fit$dm_sd`, `fit$dm_min`, `fit$dm_max`: chain summary of PIP;
- `fit$bm_sd`, `fit$bm_min`, `fit$bm_max`: chain summary of posterior effect;
- `fit$chains`: optional compact per-chain summaries when `keep_chains = TRUE`;
- `fit$ld_swap`: aggregated LD-swap diagnostics;
- `fit$ld_swap_chains`: optional per-chain LD-swap diagnostics when retained;
- `fit$input$nchains`.

This can be done in the R formatter after C++ returns enough chain summaries,
or directly through extra C++ return slots. A named `Rcpp::List` would be
cleaner, but it would be a larger compatibility change than extending the
existing slot protocol.

## Recommended API

Add these arguments to `stblr_csr()`:

```r
nchains = 1L
keep_chains = FALSE
chain_seeds = NULL
```

Behavior:

- `nchains = 1L` preserves the current output structure and random-number
  behavior as closely as possible.
- `nchains > 1L` runs independent chains for each trait in C++.
- `chain_seeds = NULL` derives deterministic chain seeds from `seed`.
- supplied `chain_seeds` should be an integer vector of length `nchains`, or
  a matrix/data frame with one seed per trait-chain pair if that flexibility
  is needed later.
- `keep_chains = FALSE` returns aggregated summaries only.
- `keep_chains = TRUE` returns compact per-chain `dm`, `bm`, final pi/variance
  summaries, diagnostics, and LD-swap diagnostics if applicable.

Parallelization should use a single OpenMP loop over total tasks:

```text
total tasks = number of traits * nchains
task = one independent genome-wide MCMC chain for one trait
```

Set `nthreads = max(1, min(ncores, nt * nchains))`. Avoid nested OpenMP.

For default seeds, preserve current single-chain seeds exactly:

```text
nchains = 1: seed + 1000003 * (t + 1)
nchains > 1: seed + 1000003 * (t + 1) + 9176 * (chain + 1)
```

This keeps existing reproducibility for `nchains = 1` while matching the BED
chains convention for additional chains.

## Implementation Options

### Option A: Use Scheduled Backend as Canonical Multi-Chain Backend

Pros:

- Scheduling code is already organized around a more selective marker update
  loop.
- There is a similar BED scheduled chains implementation to copy from.

Cons:

- The summary-stat scheduled CSR backend does not currently support `nchains`.
- It does not support LD-swap.
- `scheduled = FALSE, nchains > 1, updateLDswap = TRUE` is the immediate user
  target, and this option does not address it without porting LD-swap too.

### Option B: Refactor Single-Chain Core and Share Multi-Chain Wrappers

Pros:

- Best long-term structure.
- Avoids duplicating sampler state setup, posterior accumulation, seeding,
  diagnostics, and chain aggregation between regular and scheduled CSR.
- Makes it easier to add chain SD/min/max consistently.

Cons:

- Higher risk because it touches both C++ CSR samplers.
- Needs careful preservation of current `nchains = 1` random-number behavior.
- More work before the first useful multi-chain LD-swap implementation lands.

### Option C: Add Multi-Chain Task Loop to Regular CSR First

Pros:

- Directly supports the priority call:
  `scheduled = FALSE, nchains > 1, updateLDswap = TRUE`.
- Reuses existing regular CSR LD-swap code without porting it into scheduled
  mode first.
- Lower initial blast radius than refactoring both backends at once.
- Can still borrow the BED chains job loop, seeding, and aggregation pattern.

Cons:

- Some temporary duplication versus scheduled CSR.
- Requires extending the current 23-slot result format or adding a parallel
  chain-summary return convention.

Recommended choice after the completed CSR update: do not continue copying
chain logic into each backend. Use regular CSR as the canonical exact
summary-stat backend, use the BED scheduled chains backend as the canonical
individual-level chain backend, and introduce small shared seed/task/aggregation
helpers for scheduled CSR chain execution.

## LD-Swap Interaction

LD-swap is implemented in `src/st_cpg_omp_csr.cpp` only. The regular backend
contains:

- `build_ld_swap_friends_st_csr()`;
- `collect_ld_swap_candidates()`;
- `attempt_ld_swap_st_csr()`;
- per-trait counters for attempted and accepted swaps;
- result slot 23 LD-swap diagnostics;
- R formatting into `fit$ld_swap`.

The scheduled CSR backend lacks all of that and is blocked at the R level when
`scheduled = TRUE` and `updateLDswap = TRUE`.

For `scheduled = FALSE, nchains > 1, updateLDswap = TRUE`, each chain should
run LD-swap independently and maintain chain-level attempted/accepted counts.
Aggregate diagnostics should sum attempted and accepted counts across chains
within trait, then compute `accepted / attempted`. If `keep_chains = TRUE`,
retain one diagnostics row per trait-chain in `fit$ld_swap_chains`.

For `scheduled = TRUE, nchains > 1, updateLDswap = TRUE`, keep the existing
unsupported error initially. Supporting this later requires porting or sharing
the LD-swap helpers and deciding how swap proposals interact with scheduled
active/candidate marker lists after a swap.

## Output Structure

For compatibility, the first 23 return slots remain backend-specific base
slots. Chain-capable CSR fits add:

```text
fit$bm        matrix m x nt, mean across chains
fit$dm        matrix m x nt, mean across chains
fit$bm_sd     matrix m x nt, SD across chain-level bm
fit$bm_min    matrix m x nt, min across chain-level bm
fit$bm_max    matrix m x nt, max across chain-level bm
fit$dm_sd     matrix m x nt, SD across chain-level dm
fit$dm_min    matrix m x nt, min across chain-level dm
fit$dm_max    matrix m x nt, max across chain-level dm
fit$vbs/vgs/ves/vle/vld/pis averaged traces or documented mean traces
fit$ld_swap   per-trait aggregated diagnostics for regular CSR only
fit$chains    optional compact per-chain summaries for regular CSR only
fit$input$nchains
```

Do not concatenate full traces by default. Averaged traces preserve the
existing fit shape and keep memory bounded. If full per-chain traces become
necessary, gate them behind a separate option such as `keep_chain_traces`.

## Cleanup Needs

- The CSR regular and scheduled samplers have large duplicated single-chain
  setup and result-building logic.
- The 23-slot return format is overloaded: slot 23 means LD-swap diagnostics
  in the regular CSR backend and pi trace in the scheduled CSR backend.
- `.format_stblr_fit()` infers slot 23 meaning from shape, which works but will
  become fragile as more chain summaries are added.
- The BED chains backend proves that chain-task execution works, but its code
  is tied to packed BED data and does not directly solve CSR LD-swap.

## Testing Plan

Implemented tests cover:

- `nchains = 1` gives the same structure as current `stblr_csr()` output.
- `nchains = 1` preserves current seed behavior for a fixed small fixture.
- `nchains = 2` returns `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`,
  `bm_max`.
- `fit$dm` equals the mean of per-chain `dm` when regular CSR
  `keep_chains = TRUE`.
- `fit$bm` equals the mean of per-chain `bm` when regular CSR
  `keep_chains = TRUE`.
- `fit$input$nchains` is set.
- chain seeds are reproducible and distinct across trait-chain jobs.
- `extract_stblr_finemap_loci()` fills `pip_sd`, `pip_min`, and `pip_max`
  from a multi-chain fit.
- LD-swap diagnostics aggregate attempted and accepted counts across chains.
- `nchains > 1` with `updateLDswap = TRUE` works in `scheduled = FALSE`.
- current single-chain LD-swap tests continue to pass.
- `scheduled = TRUE, updateLDswap = TRUE` remains clearly blocked until that
  backend supports LD-swap.
- `scheduled = TRUE, nchains > 1` is implemented and tested.
- `scheduled = TRUE, keep_chains = TRUE` is clearly rejected.

## Suggested Next Codex Implementation Task

Design compact per-chain output for `stblr_csr(scheduled = TRUE)`.

Scope:

- decide whether scheduled compact output should include only `dm`/`bm` or
  also final pi/variance summaries;
- keep full traces out of the compact output unless explicitly requested;
- preserve the scheduled CSR pi trace in base slot 23;
- leave `scheduled = TRUE, updateLDswap = TRUE` unsupported.
