# ST-BLR CSR Multi-Chain Design

## Executive Summary

`stblr_csr()` currently runs one independent MCMC chain per trait in the
summary-statistic CSR backends. `ncores` controls OpenMP threads for the trait
loop, so with three traits and `ncores = 3` the observed output
`trait 0 used thread 0`, `trait 1 used thread 1`, and `trait 2 used thread 2`
is expected. It is not running multiple chains per trait.

The repository already has useful multi-chain machinery, but it is in the
packed-BED scheduled backend, not in the summary-statistic CSR backend used by
`stblr_csr()`. The BED implementation runs `nchains * nt` chain-trait jobs,
derives per-chain seeds, aggregates posterior summaries across chains, and
returns the existing 23-slot fit layout. The summary-statistic scheduled CSR
implementation does not currently expose `nchains` and still parallelizes only
over traits.

Recommended first implementation: add a multi-chain task loop to the regular
CSR backend, because that is the backend that already supports LD-swap and is
the path needed for `scheduled = FALSE, nchains > 1, updateLDswap = TRUE`.
While doing that, extract a single-chain CSR worker and shared aggregation
helpers so the scheduled CSR backend can reuse the same chain/task structure
later.

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

Not implemented in this pass:

- scheduled CSR multi-chain execution;
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

The regular CSR sampler:

- reads and builds one shared flat CSR LD object from `ld_prefix`;
- optionally builds LD-swap friend lists once;
- allocates one row per trait for `bm`, `dm`, final state, and traces;
- runs `#pragma omp parallel for` over `t = 0 .. nt - 1`;
- sets `nthreads = max(1, min(ncores, nt))`;
- seeds trait `t` with `seed + 1000003 * (t + 1)`;
- has no `chain`, `chains`, or `nchains` argument;
- has no current notion of independent chains per trait.

The main returned object is a 23-slot `std::vector<std::vector<std::vector<double>>>`.
The R formatter `.format_stblr_fit()` names and reshapes it. Important slots:

- slot 1 / index 0: `bm`, posterior mean effects;
- slot 2 / index 1: `dm`, posterior inclusion probabilities;
- slots 8 to 10 / indices 7 to 9: `vbs`, `vgs`, `ves` traces;
- slots 21 to 22 / indices 20 to 21: `vle`, `vld` traces;
- slot 23 / index 22: LD-swap diagnostics for the regular CSR backend.

`dm`, `bm`, `vbs`, `ves`, `pis`, `vgs`, `vld`, and `vle` are allocated in
`src/st_cpg_omp_csr.cpp` before the trait-parallel loop and filled from the
per-trait local accumulators inside that loop. `ld_swap` diagnostics are
created from `ld_swap_attempted_vec` and `ld_swap_accepted_vec`, then placed
in result slot 23 with a four-value marker that `.format_stblr_fit()` uses to
recognize LD-swap output.

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
- still parallelizes only over traits;
- sets `nthreads = max(1, min(ncores, nt))`;
- seeds trait `t` with `seed + 1000003 * (t + 1)`;
- does not expose `nchains`;
- does not return per-chain summaries;
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

- chain-level standard deviation, min, or max for `dm` or `bm`;
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
- Scheduled CSR remains single-chain until a follow-up.
- Requires extending the current 23-slot result format or adding a parallel
  chain-summary return convention.

Recommended choice: implement Option C first, but do it by extracting a small
regular CSR single-chain worker and a chain aggregation helper. That creates a
natural path to Option B without making the first implementation too broad.

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

For compatibility, `nchains = 1L` should keep the existing formatted fields.
For `nchains > 1L`, the recommended formatted object is:

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
fit$ld_swap   per-trait aggregated diagnostics
fit$chains    optional compact per-chain summaries
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

Add tests in a later implementation task:

- `nchains = 1` gives the same structure as current `stblr_csr()` output.
- `nchains = 1` preserves current seed behavior for a fixed small fixture.
- `nchains = 2` returns `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`,
  `bm_max`.
- `fit$dm` equals the mean of per-chain `dm` when `keep_chains = TRUE`.
- `fit$bm` equals the mean of per-chain `bm` when `keep_chains = TRUE`.
- `fit$input$nchains` is set.
- chain seeds are reproducible and distinct across trait-chain jobs.
- `extract_stblr_finemap_loci()` fills `pip_sd`, `pip_min`, and `pip_max`
  from a multi-chain fit.
- LD-swap diagnostics aggregate attempted and accepted counts across chains.
- `nchains > 1` with `updateLDswap = TRUE` works in `scheduled = FALSE`.
- current single-chain LD-swap tests continue to pass.
- `scheduled = TRUE, updateLDswap = TRUE` remains clearly blocked until that
  backend supports LD-swap.
- `scheduled = TRUE, nchains > 1` is either implemented and tested or clearly
  rejected with an informative error.

## Suggested Next Codex Implementation Task

Implement `nchains`, `keep_chains`, and `chain_seeds` for
`stblr_csr(scheduled = FALSE)` only.

Scope:

- update `stblr_csr()` argument validation and `fit$input`;
- add an exported regular CSR multi-chain C++ entry point or extend
  `stblr_cpg_omp_csr()` carefully;
- factor the regular CSR per-trait body into a single-chain worker returning
  `bm`, `dm`, final state, traces, variance/pi summaries, and LD-swap counts;
- run OpenMP over trait-chain jobs;
- aggregate means, SD/min/max for `bm` and `dm`, and LD-swap diagnostics;
- add optional compact `fit$chains` formatting;
- leave `scheduled = TRUE, updateLDswap = TRUE` unsupported;
- avoid changing `nchains = 1` behavior except for adding `fit$input$nchains`.
