# ST-BLR CSR Fine-Mapping Implementation Overview

This note summarizes the current ST-BLR / CSR fine-mapping stack and separates
implemented C++ functionality from R-level prototype orchestration. It is based
on inspection of the current source tree, especially:

- `R/sparse_ld_bed_helper.R`
- `R/finemap-stblr-csr.R`
- `R/credible_sets.R`
- `src/st_cpg_omp_csr.cpp`
- `src/st_csr_common.h`
- `src/st_cpg_omp_csr_scheduled.cpp`
- `tests/testthat/test-finemap-stblr-csr.R`
- `tests/testthat/test-csr-ld-swap.R`
- `tests/testthat/test-credible-sets.R`

No source-code changes are implied here. The main conclusion is that global
CSR ST-BLR fitting and LD-swap proposals are implemented in C++. The preferred
fine-mapping workflow should treat the genome-wide posterior as the fitting
result and extract locus-level summaries from it. Local refitting remains
useful as reference or sensitivity code, but it should not be the default
workflow.

## Current Architecture

The implemented flow is:

```text
global sparse LD prefix + summary stats
    |
    v
stblr_csr()
    |
    +--> stblr_cpg_omp_csr()                  C++ base CSR sampler
    +--> stblr_cpg_omp_csr_scheduled()        C++ scheduled sampler
    |
    v
global fit with bm, dm, parameter traces, optional ld_swap diagnostics
    |
    v
make_credible_sets()
    |
    v
cs_global$locus_sets
    |
    v
extract_stblr_finemap_loci()
    |
    +--> R post-processing of fit$dm and fit$bm
    +--> optional regional credible sets from sparse LD
```

The local-refit reference path is:

```text
cs_global$locus_sets
    |
    v
finemap_stblr_csr()
    |
    +--> R loop over traits, loci, runs
    +--> R residualized local stats
    +--> R temporary CSR subset files
    +--> .stblr_run_local_csr()
            |
            v
        stblr_cpg_omp_csr()                  C++ sampler reused locally
```

The local fine-mapping model is not a distinct native backend. It reuses the
global C++ sampler on one local trait/locus/run at a time with fixed `ve`,
`vb`, and `pi`.

## Implemented in C++

### Base CSR ST-BLR sampler

`stblr_cpg_omp_csr()` in `src/st_cpg_omp_csr.cpp` is the central native
implementation for summary-statistics single-trait ST-BLR over one or more
traits. It:

- validates marker and trait dimensions;
- reads a disk-backed sparse LD prefix;
- converts stored LD correlations to the `X_i'X_j` scale using `stats$ww`;
- builds a symmetric in-memory CSR neighbor structure;
- reconstructs and incrementally maintains residual scores;
- runs BayesC-style marker updates;
- optionally updates `B`, `E`, and `pi`;
- runs traits in parallel with OpenMP where available; and
- returns posterior means, PIPs, traces, variance summaries, diagnostics, and
  LD-swap diagnostics.

The shared CSR reader and residual helpers are in `src/st_csr_common.h`.
The key contract is that the on-disk CSR stores sparse signed correlations,
while the sampler constructs the cross-product scale internally from `ww`.

### Scheduled CSR sampler

`stblr_cpg_omp_csr_scheduled()` in `src/st_cpg_omp_csr_scheduled.cpp` is a
separate C++ sampler path for scheduled sparse updates. The R wrapper rejects
`scheduled = TRUE` together with `updateLDswap = TRUE`, so LD-swap is currently
available only on the base non-scheduled CSR sampler.

### LD-swap / Strategy 3 proposals

LD-swap support is implemented in the base CSR sampler. The native code:

- builds high-LD friend lists from the already loaded CSR structure;
- chooses included markers with excluded LD friends as swap candidates;
- proposes moving an effect from an included marker to an excluded LD friend;
- computes the residual SSE before and after the move;
- applies an MH correction with the forward and reverse proposal probabilities;
- restores the old state if rejected; and
- records attempted, accepted, and acceptance-rate diagnostics per trait.

The R wrappers validate `updateLDswap`, `ld_swap_prob`, `ld_swap_r2`,
`ld_swap_max_friends`, and `ld_swap_moves`, pass them through to C++, and expose
the result as `fit$ld_swap`.

## Recommended Fine-Mapping Workflow

The recommended workflow is:

```r
fitMH <- stblr_csr(
  stats = stats,
  Glist = Glist,
  pi_init = 0.001,
  pi_prior_a = 1,
  pi_prior_b = 1,
  h2 = 0.3,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  ncores = 3,
  seed = 10,
  scheduled = FALSE,
  updateLDswap = TRUE,
  ld_swap_prob = 0.10,
  ld_swap_r2 = 0.05,
  ld_swap_moves = 5
)

cs_global <- make_credible_sets(
  fit = fitMH,
  Glist = Glist,
  trait = "D1",
  coverage = 0.95,
  min_r2 = 0.5,
  pip_cutoff = 0.001,
  locus_pip_cutoff = 0.01,
  max_locus_distance = 1e6
)

fm <- extract_stblr_finemap_loci(
  fit = fitMH,
  Glist = Glist,
  locus_sets = cs_global$locus_sets,
  trait = "D1",
  credible_sets = TRUE,
  cs_mode = "multi"
)
```

`extract_stblr_finemap_loci()` does not run MCMC. It uses genome-wide
posterior PIPs/effects from `fit$dm` and `fit$bm`, adds marker positions from
`Glist`, and optionally constructs regional credible sets from sparse LD.

## Implemented in R

### Global wrapper and formatting

`stblr_csr()` in `R/sparse_ld_bed_helper.R` is an R wrapper around the native
CSR samplers. It constructs priors, initializes sampler inputs, chooses between
base and scheduled native backends, formats the returned list into named fit
components, and records the input settings.

The wrapper is user-facing and exported. Its computational core is native code.

### Credible-set utilities

Credible-set construction is R-level post-processing in `R/credible_sets.R`:

- `make_credible_sets_from_ld()` builds single-signal LD/PIP credible sets from
  a dense regional LD matrix.
- `make_multisignal_credible_sets_from_ld()` builds approximate multi-signal
  credible sets from marginal PIPs and dense regional LD.
- `make_credible_sets()` extracts PIPs from a fitted ST-BLR object,
  defines or accepts loci, obtains regional LD either from dense input or by
  densifying sparse CSR slices, and returns `summary`, `sets`, `loci`, and
  `locus_sets`.

`locus_sets` is explicitly useful as input to `extract_stblr_finemap_loci()`
or, when local refitting is needed, `finemap_stblr_csr()`. These
credible sets are based on marginal PIPs and LD post-processing; the
multi-signal function is approximate and not a SuSiE-style per-effect
credible-set construction.

### Posterior locus extraction

`extract_stblr_finemap_loci()` in `R/extract-stblr-finemap-loci.R` is the
lightweight post-processing path. It:

- accepts predefined `locus_sets`, including `cs_global$locus_sets`;
- extracts marker PIPs from `fit$dm`;
- extracts marker effects from `fit$bm`;
- carries optional `fit$dm_sd`, `fit$dm_min`, `fit$dm_max`, `fit$bm_sd`,
  `fit$bm_min`, and `fit$bm_max` when present;
- adds chromosome and position information from `Glist`;
- computes locus lead markers, total PIP, and secondary PIP;
- optionally reads sparse LD and constructs single- or multi-signal credible
  sets; and
- returns a `"stblr_finemap"` object without a `runs` element.

This is now the main fine-mapping summary workflow because it reuses the
efficient genome-wide CSR fit and avoids repeated local sampler calls.

### Local refit fine mapping

`finemap_stblr_csr()` in `R/finemap-stblr-csr.R` is currently a working
R-level prototype. It:

- resolves marker and trait names;
- cleans supplied marker sets and aligns them to global marker order;
- resolves global `ve`, `vb`, and `pi` from posterior traces after burn-in,
  with optional overrides;
- reads the full CSR prefix into R once with `sparseLD_read_CSR()`;
- loops over requested traits, loci, and independent runs;
- builds local summary-statistics objects;
- optionally residualizes local `wy` against global posterior mean effects
  outside the locus;
- writes a temporary local CSR prefix per locus/trait;
- calls `.stblr_run_local_csr()` once per run;
- aggregates marker PIPs and posterior mean effects across runs;
- optionally constructs local credible sets from aggregated local PIPs; and
- returns a `"stblr_finemap"` object with summaries, marker tables, credible
  sets, loci, run summaries, and resolved parameters.

`.stblr_run_local_csr()` is a thin local runner. It calls the same native
`stblr_cpg_omp_csr()` backend with:

- one trait;
- local marker count;
- fixed `B = vb`, `E = ve`, and `pi = c(1 - pi, pi)`;
- `updateB = FALSE`;
- `updateE = FALSE`;
- `updatePi = FALSE`;
- `ncores = 1`; and
- optional LD-swap pass-through.

This is functional and should remain available as a local-refit/reference
method. It is not the preferred default because R repeatedly writes temporary
CSR subsets and crosses the R/C++ boundary for every run.

## Tests Covering the Current Behavior

Current tests exercise:

- credible-set construction, marker-name alignment, incomplete-set behavior,
  LD-neighborhood removal, and approximate multi-signal behavior;
- fine-mapping set cleaning, global parameter extraction after burn-in,
  aggregation of PIPs/effects across runs, and use of aggregated local PIPs for
  credible sets;
- posterior-only fine-mapping extraction without sparse LD;
- posterior-only fine-mapping extraction with a tiny sparse LD fixture;
- the guarantee that credible-set thresholds do not remove markers before local
  model fitting;
- CSR LD-swap argument validation;
- zero LD-swap diagnostics by default;
- LD-swap execution on a tiny CSR example; and
- LD-swap pass-through through `.stblr_run_local_csr()`.

The tests are useful regression coverage for the R prototype and the LD-swap
interface. They do not yet test a native local fine-mapping backend because no
such backend exists.

## Prototype Boundaries

The main prototype boundary is local refitting. If local refits remain
important after posterior-only extraction is validated, the following
operations are still R-level orchestration and are the best candidates for
native migration:

1. Looping over loci, traits, and runs.
2. Extracting local CSR subregions.
3. Writing temporary CSR subset files.
4. Re-reading those temporary subsets in the C++ sampler.
5. Constructing residualized local `wy`.
6. Aggregating run-level PIPs and effects.
7. Returning compact local run diagnostics.

The most expensive and least scalable part is the temporary CSR path:

```text
full CSR in R
    -> R subset extraction
    -> temporary CSR files
    -> C++ read_and_build_st_ld_csr()
    -> one local run
```

This should become an in-memory local subregion path.

Credible-set construction can stay in R initially. It is post-processing over
local dense LD and aggregated PIPs, and it is less tightly coupled to the MCMC
state than the local sampler itself. A later native implementation may be
useful for very large batches, but it is not the first bottleneck.

## Optional Native Local-Refit Backend

A possible later native target is a C++ entry point that runs local refits over
many loci and runs internally. This is no longer required for the main
fine-mapping workflow, but it could still be useful for sensitivity analyses.
Conceptually:

```text
stblr_cpg_omp_csr_finemap(
    global wy, ww, yy,
    global posterior mean effects,
    full CSR prefix or full in-memory CSR,
    list of locus marker indices,
    per-trait ve/vb/pi,
    runs/seeds/nit/nburn,
    use_residual,
    updateLDswap controls
)
```

The first version should preserve the current statistical behavior:

- one local model per trait/locus/run;
- fixed `ve`, `vb`, and `pi`;
- same BayesC marker update as `stblr_cpg_omp_csr()`;
- same residualization formula as `.stblr_residualize_wy_csr()`;
- same LD-swap proposal logic when enabled;
- same local PIP/effect summaries as `.stblr_aggregate_finemap_runs()`; and
- deterministic seed mapping compatible with the current R wrapper, or clearly
  documented if changed.

The main internal C++ building blocks should be:

- a reusable CSR loader for the full LD prefix;
- a function that builds an in-memory local `STLDCSR` view or copy from global
  marker indices;
- a local stats builder using global `wy`, `ww`, `yy`, optional global `bm`,
  and the local CSR;
- a single local-run sampler that reuses the existing marker update, residual,
  variance-summary, and LD-swap helpers; and
- a batch loop that accumulates sums and sums of squares for PIP/effect means
  and standard deviations without storing every full run unless requested.

The R wrapper can then become mostly validation, input alignment, argument
normalization, and output formatting.

## Suggested Local-Refit Migration Order

1. Factor reusable native helpers from `stblr_cpg_omp_csr()` without changing
   sampler behavior:
   - marker update;
   - residual rebuild;
   - LD-swap friend construction and proposal;
   - posterior accumulation.
2. Add an in-memory local CSR extraction helper from a loaded `STLDCSR` and a
   vector of global marker indices.
3. Add a native single-local-run function and compare it against
   `.stblr_run_local_csr()` on tiny deterministic examples.
4. Add a native batch function over runs for one locus/trait.
5. Extend the batch function over loci and traits.
6. Keep R credible-set post-processing unchanged until native local summaries
   are validated.
7. Only then consider native credible-set construction if profiling shows it is
   material.

## Open Design Points

- Whether the native backend should accept a full `ld_prefix` and read CSR
  itself, or accept an R-loaded CSR object. Reading once in C++ is likely
  cleaner for scale.
- Whether local CSR subregions should be copied into compact local indices or
  represented as views into global CSR. A compact copy is simpler and matches
  the current local sampler's indexing.
- Whether local runs should be parallelized across traits, loci, runs, or a
  combination. The current global sampler parallelizes across traits; local
  fine mapping may benefit more from locus/run parallelism.
- Whether to preserve exact current seed sequences. This matters for regression
  tests but may constrain parallel scheduling.
- How much run-level detail to return. The current R output stores compact run
  summaries; storing all local fits would be too large for genome-scale use.

## Practical Status

Implemented and reasonably stable:

- disk-backed sparse LD CSR infrastructure;
- global `stblr_csr()` wrapper;
- base C++ CSR sampler `stblr_cpg_omp_csr()`;
- scheduled C++ CSR sampler without LD-swap;
- optional LD-swap in the base CSR sampler;
- LD-swap diagnostics in formatted fits;
- R credible-set construction;
- R posterior fine-mapping extraction from `fit$dm` and `fit$bm`;
- R local fine mapping using repeated local C++ sampler calls as a reference
  method.

Still R-level prototype:

- local refit orchestration;
- local residualized stats construction;
- local CSR subsetting and temporary file writing;
- looping over loci/runs;
- aggregation of local runs;
- local credible-set post-processing.

Best candidates to move to C++ next:

- in-memory local CSR subsetting;
- local residualization;
- batched local sampler calls over loci/runs;
- run aggregation;
- LD-swap-enabled local fine mapping without temporary CSR files.
