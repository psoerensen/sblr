# ST-BLR Backend Harmonization Design

## Executive Summary

The regular summary-statistic CSR backend in `src/st_cpg_omp_csr.cpp` is now
the canonical exact ST-BLR summary-stat backend. It supports independent
trait-by-chain OpenMP jobs, deterministic chain seeds, optional per-chain
compact summaries, chain SD/min/max for `dm` and `bm`, and LD-swap diagnostics
aggregated across chains.

The related scheduled CSR and individual/BED backends now expose the same
marker-level chain stability summaries as regular CSR. Scheduled CSR remains
LD-swap-free, and the individual/BED chain support still lives in a separate
scheduled chains backend.

BayesR backend harmonization is tracked separately in
`docs/dev/stblr_bayesr_backend_design.md`. The existing BED BayesR scheduled
chains backend now uses standard `dm = P(component > 0)` and exposes standard
`bm`/`dm` chain summaries. It remains a lower-level experimental backend with
no public wrapper. Plain exact and scheduled CSR BayesR backends remain future
work.

Implementation status: scheduled CSR and individual/BED scheduled chains now
extend their existing 23-slot base results with CSR-compatible chain-stability
summaries in slots 24 to 29: `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`,
and `dm_max`. The R formatter exposes these as marker-by-trait matrices with
the same dimnames as `bm` and `dm`.

Recommended implementation scope for this pass is design only. The backends
are too divergent to safely add shared helpers or chain summaries everywhere
without first deciding the public return convention and reducing the large
amount of repeated historical code in the BED chain file.

## Current Backend Inventory

### Regular CSR Summary-Statistic Backend

File: `src/st_cpg_omp_csr.cpp`

Exported C++ function: `stblr_cpg_omp_csr()`

R export: `stblr_cpg_omp_csr()` in `R/RcppExports.R`

User-facing wrapper: `stblr_csr(..., scheduled = FALSE)` in
`R/sparse_ld_bed_helper.R`

This is the reference exact summary-statistic backend.

Chain arguments:

- `nchains`
- `keep_chains`
- `chain_seeds`

LD-swap arguments:

- `updateLDswap`
- `ld_swap_prob`
- `ld_swap_r2`
- `ld_swap_max_friends`
- `ld_swap_moves`

Seed handling:

- If `chain_seeds` is supplied, task seed is
  `chain_seeds[chain] + 1000003 * (trait + 1)`.
- If `nchains == 1` and no chain seeds are supplied, task seed remains
  `seed + 1000003 * (trait + 1)` to preserve the old single-chain behavior.
- If `nchains > 1` and no chain seeds are supplied, task seed is
  `seed + 1000003 * (trait + 1) + 9176 * (chain + 1)`.

OpenMP task structure:

- `ntasks = nt * nchains`.
- The active mapping is `trait = task / nchains`, `chain = task % nchains`.
- OpenMP runs `#pragma omp parallel for num_threads(nthreads) schedule(static)`
  over `task = 0 .. ntasks - 1`.
- `nthreads` is bounded by requested cores and task count.

Output layout:

- The base 23 slots preserve the existing formatter contract:
  - `result[0]`: `bm`
  - `result[1]`: `dm`
  - `result[2]`: `wy`
  - `result[3]`: `r`
  - `result[4]`: final `b`
  - `result[5]`: final `d`
  - `result[6]`: marker index
  - `result[7:9]`: `vbs`, `vgs`, `ves`
  - `result[10:15]`: covariance/final variance matrices
  - `result[16]`: final pi
  - `result[17]`: posterior mean pi
  - `result[18:19]`: diagnostics placeholders
  - `result[20:21]`: `vle`, `vld`
  - `result[22]`: LD-swap diagnostics with marker column
- If `nchains > 1`, slots 23 to 28 are added:
  - `bm_sd`, `bm_min`, `bm_max`
  - `dm_sd`, `dm_min`, `dm_max`
- If `keep_chains = TRUE`, slots 29 to 31 are added:
  - chain-major `dm`
  - chain-major `bm`
  - chain-major LD-swap diagnostics

Where chain summaries are produced:

- Per-task posterior summaries are stored in `bm_task` and `dm_task`.
- Per-trait means are accumulated across chains and divided by `nchains`.
- Min and max are computed during the chain aggregation loop.
- SD is computed after means are known using the sample SD denominator
  `nchains - 1`.
- `.format_stblr_fit()` maps the extra C++ slots to `fit$dm_sd`,
  `fit$dm_min`, `fit$dm_max`, `fit$bm_sd`, `fit$bm_min`, and `fit$bm_max`.

How `keep_chains` works:

- C++ returns compact chain-major vectors for `dm`, `bm`, and LD-swap
  diagnostics.
- `.format_stblr_fit()` reshapes those vectors into `fit$chains[[trait]]`, a
  list of `chain1`, `chain2`, etc.
- Each compact chain entry has marker-level `dm` and `bm`, and an `ld_swap`
  vector when diagnostics are present.
- `fit$ld_swap_chains[[trait]]` is also produced as a data frame.

LD-swap diagnostics:

- LD friends are built once from the shared CSR LD.
- Each chain tracks attempted and accepted swaps independently.
- Aggregated trait diagnostics sum attempted and accepted across chains.
- Acceptance rate is computed as `accepted / attempted`, or `0` when no
  attempts were made.

### Scheduled CSR Summary-Statistic Backend

File: `src/st_cpg_omp_csr_scheduled.cpp`

Exported C++ function: `stblr_cpg_omp_csr_scheduled()`

R export: `stblr_cpg_omp_csr_scheduled()` in `R/RcppExports.R`

User-facing wrapper: `stblr_csr(..., scheduled = TRUE)` in
`R/sparse_ld_bed_helper.R`

This is an approximate acceleration backend. It uses the same summary-statistic
CSR inputs but updates markers according to scheduling rules:

- full sweeps controlled by `full_sweep_every`;
- null marker skipping controlled by `null_skip_base` and `null_skip_max`;
- candidate marker refresh controlled by `candidate_threshold` and
  `candidate_lifetime`;
- optional LD-neighbor wakeup controlled by `wakeup_ld_neighbors`,
  `wakeup_diff_threshold`, and `wakeup_max_neighbors`.

Current chain and LD-swap status:

- It supports `nchains` in C++ using trait-major task mapping.
- The R wrapper permits `scheduled = TRUE, nchains > 1`.
- `keep_chains = TRUE` remains unsupported for scheduled CSR.
- It does not support LD-swap.
- The R wrapper rejects `scheduled = TRUE, updateLDswap = TRUE`.

Parallelization:

- Parallelizes over `nt * nchains` trait-chain tasks.
- For `nchains == 1` with no chain seeds, seeds trait `t` with
  `seed + 1000003 * (t + 1)`, preserving previous single-chain behavior.
- For default multi-chain runs, seeds task `(trait, chain)` with
  `seed + 1000003 * (trait + 1) + 9176 * (chain + 1)`.
- Uses `#pragma omp parallel for` over `task = 0 .. nt * nchains - 1`.

Output layout:

- Returns the same base fit shape used by the formatter, with 23 slots.
- Slot 22 is a pi trace for the scheduled backend, not LD-swap diagnostics.
- Slots 23 to 28 contain `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`,
  and `dm_max`.
- `.format_stblr_fit()` distinguishes LD-swap diagnostics from a pi trace by
  slot shape.

Duplicated code:

- BayesC marker update code is near-duplicated from regular CSR, with a
  scheduled update result containing `p1` and `diff`.
- Variance updates, pi updates, VLE/VLD computation, input validation, result
  construction, and per-trait sampler setup are largely duplicated.
- The shared CSR reader and residual rebuild already live in
  `src/st_csr_common.h`.

What could be shared:

- Result aggregation for chain means, SD/min/max, and LD-swap diagnostics.
- Common scalar update helpers where scheduled and exact code have identical
  math.
- A single-chain worker result type for CSR, if introduced carefully.

Recommendation:

Scheduled CSR now uses shared seed/task helpers for chain execution. Keep it
LD-swap-free until swap behavior with scheduled active/candidate marker lists
is designed.

### Individual/BED Sparse Backend

File: `src/st_cpg_omp_individual.cpp`

Main relevant exported function: `stblr_cpg_omp_bed_marker_sparse()`

This is an individual-level PLINK/BED backend. The file also contains older
block-oriented individual/BED implementations and commented historical code.

Current behavior:

- Supports multiple traits through columns of `y`.
- Parallelizes over traits.
- Seeds trait `t` with `seed + 1000003 * (trait + 1)`.
- Does not expose `nchains`.
- Does not return chain SD/min/max fields.
- Uses a 23-slot result for the sparse BED marker backend, including a pi trace
  in slot 22.
- Older block functions in the same file use 20-slot result layouts.

Recommendation:

Do not add `nchains` here now. The scheduled chains backend is already the
more mature individual/BED chain path. Extending this older sparse backend
would create another independent copy of chain aggregation logic.

### Individual/BED Scheduled Backend

File: `src/st_cpg_omp_individual_scheduled.cpp`

Exported C++ function: `stblr_cpg_omp_bed_marker_scheduled()`

This is the single-chain scheduled PLINK/BED backend.

Current behavior:

- Reads and packs BED data.
- Builds marker maps.
- Parallelizes over traits.
- Uses scheduled null-marker update logic similar to scheduled CSR.
- Seeds trait `t` with `seed + 1000003 * (trait + 1)`.
- Does not support chains.
- Returns a 23-slot layout with `vle`, `vld`, and a pi trace.
- Has extensive duplicated code with the scheduled chains backend.

Recommendation:

Treat this as the single-chain scheduled BED implementation. If the scheduled
chains backend is stable for single-chain operation too, consider routing
future user-facing scheduled BED calls through the chains backend with
`nchains = 1`, then gradually retire this as an internal/reference
implementation.

### Individual/BED Scheduled Chains Backend

File: `src/st_cpg_omp_individual_scheduled_chains.cpp`

Exported C++ function: `stblr_cpg_omp_bed_marker_scheduled_chains()`

This is the mature chain-capable PLINK/BED backend.

Argument names:

- `nchains`
- `ncores`
- `seed`
- plus scheduled controls, BED row/filtering controls, and output controls
  such as `return_wy`, `return_r`, `read_block_size`, and `progress_every`.

Seed formula:

- `seed + 1000003 * (trait + 1) + 9176 * (chain + 1)`.

OpenMP task structure:

- `njobs = nchains * nt`.
- Active mapping is `chain = job / nt`, `trait = job % nt`.
- OpenMP runs over jobs and stores each result in `job_results[job]`.
- This mapping differs from regular CSR, which uses `trait = task / nchains`
  and `chain = task % nchains`.

Aggregation:

- Per-job `ChainResultSTScheduled` objects contain `bm`, `dm`, final states,
  traces, final variances, final pi, mean pi, timing, and sample counts.
- Aggregation averages `bm`, `dm`, final state, traces, final variance/pi
  summaries, CPO diagnostics, and runtime diagnostics across chains.
- It computes CSR-compatible chain-stability summaries for marker-level `dm`
  and `bm`: sample SD for `nchains > 1`, zero SD for `nchains == 1`, and
  min/max across chain-level summaries.
- It does not return compact per-chain marker summaries.
- It does not have LD-swap diagnostics.

Differences from regular CSR chains:

- BED chain results are aggregated back into the same 23-slot base shape as
  single-chain scheduled BED, then extended with optional chain-summary slots
  24 to 29 following the regular CSR convention.
- Job mapping is chain-major rather than trait-major.
- It has a single-chain worker result struct, which is a useful architecture
  pattern, but its return convention is now behind the CSR chain convention.
- The file contains many commented historical copies; only the active top
  implementation should be used as a template.

Recommendation:

Align this backend to the regular CSR output convention in a focused follow-up:
add SD/min/max slots for `dm` and `bm`, and update the R formatter to surface
them for BED chain fits. Do not add compact `keep_chains` support until there
is a clear user-facing option for BED chains.

## Simplification Options

### Option A: Keep all five backends separate and add `nchains` to each

Pros:

- Direct and local changes.
- Each backend can preserve its current result layout.

Cons:

- Creates four or more copies of seed, task, aggregation, and diagnostics
  logic.
- Increases risk of divergent random-number and output semantics.
- Does not address the already large duplication between scheduled BED files.

Recommendation: reject.

### Option B: Use regular CSR and BED scheduled chains as canonical backends

Pros:

- Matches the current maturity of the codebase.
- Keeps exact summary-stat work centered in `st_cpg_omp_csr.cpp`.
- Uses the existing BED chain backend rather than extending older paths.
- Minimizes user-facing backend proliferation.

Cons:

- Some older backends remain present as internal/reference code.
- Requires wrapper decisions for when to route BED scheduled single-chain calls
  through the chains backend.

Recommendation: use as the main simplification direction.

### Option C: Refactor common chain machinery into shared helpers

Pros:

- Best long-term way to avoid duplicated chain aggregation.
- Gives one seed rule, one task mapping convention, one SD/min/max
  implementation, and one diagnostics aggregation convention.
- Makes scheduled CSR chain support easier later.

Cons:

- Requires touching active C++ files.
- Needs tests that prove single-chain reproducibility is preserved.
- Generic aggregation helpers may need both Armadillo matrix helpers and
  backend-specific result adapters.

Recommendation: do this incrementally after the design is accepted. Start with
header-only seed/task helpers and chain summary aggregation helpers that do not
change sampler math.

### Option D: Only update user-facing wrappers and leave older C++ internal

Pros:

- Lowest risk.
- Can make unsupported combinations fail clearly.
- Lets `extract_stblr_finemap_loci()` stay backend-agnostic by using optional
  fields when present.

Cons:

- Does not reduce native-code duplication.
- Cannot give scheduled CSR or BED chains the new SD/min/max fields by itself.

Recommendation: use together with Option B for short-term stability.

## Recommended Simplification

Use Option B plus a staged version of Option C.

1. Keep `st_cpg_omp_csr.cpp` as the canonical exact summary-stat backend.
2. Keep `st_cpg_omp_csr_scheduled.cpp` as approximate acceleration with
   multi-chain task execution and no LD-swap.
3. Keep LD-swap only in regular CSR for now.
4. Prefer `st_cpg_omp_individual_scheduled_chains.cpp` for BED chain support.
5. Avoid extending both `st_cpg_omp_individual.cpp` and
   `st_cpg_omp_individual_scheduled.cpp` with independent `nchains` logic.
6. Standardize chain-capable outputs around:
   - `dm`
   - `bm`
   - `dm_sd`, `dm_min`, `dm_max`
   - `bm_sd`, `bm_min`, `bm_max`
   - `input$nchains`
   - `input$keep_chains` where the backend has a keep option
7. Keep `extract_stblr_finemap_loci()` backend-agnostic by relying only on
   those optional fields.

## Shared Helper Candidates

A small header such as `src/st_chain_utils.h` is preferable to expanding
`src/st_csr_common.h`, because the helpers are not CSR-specific.

Recommended helpers:

```cpp
inline unsigned int stblr_chain_seed(int seed, int trait, int chain) {
 return static_cast<unsigned int>(
  seed + 1000003 * (trait + 1) + 9176 * (chain + 1)
 );
}

inline unsigned int stblr_seed_with_chain_base(int chain_seed, int trait) {
 return static_cast<unsigned int>(chain_seed + 1000003 * (trait + 1));
}
```

Task mapping should be standardized before use. The regular CSR mapping is
trait-major:

```cpp
inline int stblr_task_trait(int task, int nchains) {
 return task / nchains;
}

inline int stblr_task_chain(int task, int nchains) {
 return task % nchains;
}
```

If the BED chains backend is aligned later, use this mapping there too. If not,
document that BED remains chain-major internally.

Aggregation helpers:

- Mean across chain rows for an Armadillo matrix.
- Sample SD across chain rows, returning zero when `nchains == 1`.
- Per-marker min and max across chain rows.
- Sum diagnostics for attempted and accepted counts, then compute acceptance
  rate.

Keep the helpers small and boring. Do not introduce a generalized sampler
framework until two active backends are actually using the helper code.

## Implementation Scope Decision

Implemented scope:

- Regular CSR implements chains, optional compact chain output, and LD-swap.
- BED scheduled chains implements CSR-compatible chain summary fields.
- Scheduled CSR implements chains and CSR-compatible chain summary fields,
  while still rejecting LD-swap and compact per-chain output.

## Recommended Implementation Order

1. Add `src/st_chain_utils.h` with seed and task mapping helpers.
2. Use those helpers in regular CSR only, preserving behavior exactly.
3. Add focused C++ unit-style or R integration tests for seed reproducibility
   and chain aggregation.
4. Add chain SD/min/max aggregation to
   `st_cpg_omp_individual_scheduled_chains.cpp`.
5. Update `.format_stblr_fit()` so BED chain fits expose the same optional
   fields as CSR chain fits.
6. Consider routing scheduled BED single-chain calls through the chains backend
   with `nchains = 1` if tests show identical or acceptable behavior.
7. Scheduled CSR `nchains` uses the same trait-major task mapping and shared
   seed helpers as regular CSR.
8. Keep scheduled CSR LD-swap unsupported until swap behavior with scheduled
   active/candidate marker lists is designed.

## Testing Plan

CSR exact backend:

- `nchains = 1` keeps the existing output shape and seed behavior.
- `nchains = 3` returns `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`,
  `bm_max`.
- `keep_chains = TRUE` returns compact chain summaries whose row means equal
  `fit$dm` and `fit$bm`.
- LD-swap works with `nchains = 3`.
- LD-swap attempted and accepted counts aggregate across chains.
- `extract_stblr_finemap_loci()` maps optional fields into `pip_sd`,
  `pip_min`, `pip_max`, `bm_sd`, `bm_min`, and `bm_max`.

CSR scheduled backend:

- `scheduled = TRUE, nchains = 1` preserves the old single-chain seed rule.
- `scheduled = TRUE, nchains = 2` returns the same chain-summary fields as
  regular CSR.
- `scheduled = TRUE, updateLDswap = TRUE` still rejects clearly.
- `scheduled = TRUE, keep_chains = TRUE` rejects clearly until compact
  scheduled chain output is implemented.

Individual/BED chain backend:

- Chain-capable backend returns consistent `dm` and `bm` dimensions.
- After alignment, it returns `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`,
  and `bm_max`.
- Existing row-subset and scheduled BED tests continue to pass.
- `nchains = 1` through the chain backend remains valid.

API consistency:

- Invalid `nchains` errors.
- Invalid `keep_chains` errors where supported.
- Invalid chain seed inputs error where supported.
- Seed reproducibility is tested for `nchains = 1` and `nchains > 1`.

## Suggested Next Implementation Prompt

Implement Scope 2 for ST-BLR chain utilities.

Create `src/st_chain_utils.h` with deterministic seed helpers and a single
trait-major task mapping convention. Use it in `src/st_cpg_omp_csr.cpp` only,
without changing any function signatures or output layout. Add or update tests
that confirm `nchains = 1` output is unchanged for a fixed tiny CSR fixture and
that `nchains > 1` still returns the existing chain summary fields. Do not
touch scheduled CSR or BED backends in that pass.
