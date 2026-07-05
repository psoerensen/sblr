# sblr large-scale validation checkpoint: limitations and future implementation notes

## Current status

The `sblr` package has now been migrated to the named `stblr_raw` backend-output schema (validated against the `stblr_raw_schema` contract). All active C++ BLR backends now return named raw objects, and the R layer converts these into a stable formatted fit object through the canonical `.as_stblr_fit()` formatter.

The old positional backend-output path has been removed or made unsupported. The current intended architecture is:

```text
C++ BLR backend
    ↓
stblr_raw object
    ↓
validated against stblr_raw_schema
    ↓
.as_stblr_fit() / canonical formatter
    ↓
stable formatted fit object
    ↓
summaries, plots, architecture summaries, credible sets, consistency checks
```

The full testthat suite passes:

```text
2626 tests, 0 failures, 0 errors
```

The test suite verifies field presence, formatter consistency, schema routing, LD-swap diagnostics, annotation/group/prior/component fields, single-chain and multi-chain marker summaries, and consistency checks on small fixtures.

No sampler math, RNG behaviour, OpenMP scheduling, BED decoding, LD-alignment logic, public arguments, defaults, or final formatted fit field names were intentionally changed during the migration.

## Stable formatted fit fields

The formatted fit object is intended to preserve stable public fields, including:

```text
bm, dm, wy, r, b, d,
vbs, vgs, ves, vle, vld,
pis, pi, pim,
covb, covg, cove,
input, chains,
ld_swap, ld_swap_chains,
selection_s*,
comp_prob,
dm_component_mean,
annotation/group/prior/component fields where relevant
```

Fields such as `ld_swap` and `ld_swap_chains` should be present in formatted fits. If diagnostics are unavailable or not applicable, they should be present with `NULL` values rather than omitted.

Marker-level chain-summary fields are also intended to be stable:

```text
bm_sd, bm_min, bm_max,
dm_sd, dm_min, dm_max
```

For single-chain fits, these summaries are represented as degenerate summaries:

```text
*_sd  = 0
*_min = corresponding mean
*_max = corresponding mean
```

This is a formatter convention only. It does not change sampler output.

## Important limitations before large-scale use

### 1. Tests are still mostly tiny-fixture tests

The full test suite passes, but the tests mostly exercise small synthetic or tiny package fixtures. They confirm that the package is structurally consistent, but they do not prove that all large-scale behaviours are optimal.

Large-scale validation is still needed for:

```text
memory use,
runtime,
OpenMP scaling,
long-chain stability,
large sparse-LD handling,
large BED input handling,
multi-trait behaviour,
annotation-heavy models,
BayesR/SBayesRC behaviour,
LD-swap diagnostics at scale,
posterior summary stability.
```

### 2. The schema is stable but still young

`stblr_raw` (validated against `stblr_raw_schema`) is now the canonical internal raw schema, but it is still a new package-level contract. Real data may reveal fields that need clearer naming, stronger validation, or better documentation.

Large-scale runs should therefore store:

```text
fit$schema / raw schema version,
model type,
backend name,
seed,
niter,
nburn,
nchains,
ncores,
number of markers,
number of traits,
number of nonzero LD entries,
annotation/group/prior metadata where relevant.
```

### 3. BED BayesR `pis` trace may still be a special case

One known gap from the migration work is that BED BayesR did not appear to emit a true per-iteration `pis` trace from the C++ sampler. The formatter should not invent this trace.

Therefore, if BED BayesR is used, confirm whether:

```text
fit$pis exists,
fit$pi exists,
fit$pim exists,
and whether the absence of fit$pis is expected for that backend.
```

If per-iteration `pis` is needed for BED BayesR, this should be implemented in the native sampler as a real tracked quantity.

### 4. Single-chain chain summaries are degenerate

For `nchains == 1`, fields such as `bm_sd`, `dm_sd`, `bm_min`, `bm_max`, `dm_min`, and `dm_max` are formatter-stabilised degenerate summaries.

This is useful for structural consistency, but it should not be interpreted as evidence of between-chain stability. For convergence assessment, use multiple chains.

### 5. Stack-imbalance warning was not reproduced

A previous R stack-imbalance warning was observed interactively after a learned annotation CSR run, but stress testing did not reproduce it. The current backend appears to avoid R API calls inside OpenMP worker regions.

Still, when running large jobs, it is worth monitoring for:

```text
stack imbalance warnings,
segmentation faults,
random R session instability,
warnings appearing after native-code calls,
differences between fresh-session and long-session behaviour.
```

If such warnings reappear, reproduce them in a fresh R session with a minimal script and record `sessionInfo()`.

### 6. No long MCMC validation yet

The current tests verify structure and small examples. They do not validate long-chain posterior behaviour on real data.

Large-scale runs should include sanity checks for:

```text
posterior means,
posterior variances,
inclusion probabilities,
variance components,
heritability estimates,
LD/LE decomposition,
annotation effects,
component probabilities,
chain summaries,
trace behaviour,
effective sample behaviour if available.
```

### 7. Memory use may become the main practical bottleneck

Large-scale BLR fits can become memory-heavy, especially when storing:

```text
posterior marker traces,
per-iteration b/d samples,
multi-trait arrays,
multi-chain outputs,
annotation/component traces,
LD-swap diagnostics,
large sparse-LD objects.
```

For real runs, it may be useful to implement or strengthen options for:

```text
posterior means only,
selective trace storage,
marker subset trace storage,
thinning,
component-only storage,
streaming summaries,
on-disk output,
checkpointing.
```

### 8. Formatter stability does not guarantee scientific calibration

The migration fixed the software contract. It does not itself validate model calibration.

Large-scale model validation should compare:

```text
known benchmark results,
expected SNP heritability ranges,
posterior inclusion probability distributions,
prediction accuracy,
LD/LE decomposition behaviour,
annotation enrichment behaviour,
BayesC versus BayesR behaviour,
summary-statistic versus BED-based behaviour where comparable.
```

## Recommended checks for first large-scale runs

For the first real-data run, start with a moderate chromosome or marker subset before running the full genome.

Recommended checks after each fit:

```r
check_stblr_consistency(fit)
names(fit)
str(fit, max.level = 1)
summary(colMeans(fit$vgs))
summary(colMeans(fit$ves))
summary(colMeans(fit$pis), if available)
```

Also check dimensions of important fields:

```r
dim(fit$bm)
dim(fit$dm)
dim(fit$vgs)
dim(fit$ves)
dim(fit$pis)
dim(fit$chains)
```

For multi-chain runs, check:

```r
dim(fit$bm_sd)
dim(fit$dm_sd)
summary(fit$bm_sd)
summary(fit$dm_sd)
```

For annotation models, check:

```r
names(fit$annotation)
names(fit$component)
dim(fit$comp_prob)
dim(fit$dm_component_mean)
```

For LD-swap diagnostics, check:

```r
fit$ld_swap
fit$ld_swap_chains
```

## Suggested future implementations

### 1. Formal schema validator

`.validate_stblr_raw()` is now a thin alias over `.is_stblr_raw()`. A deeper internal validator for raw `stblr_raw` objects (checked against `stblr_raw_schema`) before formatting could still be added. This should check:

```text
top-level schema blocks,
required field names,
field dimensions,
trait/marker consistency,
chain consistency,
annotation/group/prior/component consistency.
```

This would catch backend-output problems before they propagate into formatted fits.

### 2. Formal formatted-fit contract validator

Strengthen `check_stblr_consistency()` into the main contract validator for final fit objects.

It should distinguish:

```text
required fields,
optional fields,
present-but-NULL fields,
backend-specific fields,
single-chain degenerate summaries,
multi-chain genuine summaries.
```

### 3. Developer schema documentation

Update or create a developer note documenting:

```text
stblr_raw / stblr_raw_schema raw schema,
formatted fit contract,
backend-specific fields,
single-chain conventions,
LD-swap diagnostic conventions,
annotation/group/prior/component conventions.
```

This should live in `docs/dev/`.

### 4. Large-scale regression fixtures

Add one or more medium-sized regression fixtures that are larger than toy tests but still feasible in CI or manually.

For example:

```text
one small chromosome region,
one CSR fit,
one BED fit,
one annotation fit,
one BayesR or SBayesRC fit,
one multi-chain fit.
```

These should test structure and rough numerical sanity, not exact long-chain equality.

### 5. Backend inventory test

Maintain a backend inventory test that confirms every active backend:

```text
returns a valid stblr_raw object,
formats successfully,
passes check_stblr_consistency(),
contains expected public fields,
does not fall back to legacy positional output.
```

### 6. BED BayesR `pis` tracking

If needed, implement true per-iteration `pis` tracking for BED BayesR in C++.

This should not be patched in R. It should be emitted by the sampler if the trace is scientifically meaningful.

### 7. Storage-control options

Add or improve user-facing or internal options controlling how much posterior output is stored.

Potential options:

```text
store_b
store_d
store_pis
store_chains
store_marker_trace
store_component_trace
thin
save_every
posterior_means_only
```

This will become important for large-scale genome-wide runs.

### 8. Checkpointing and resumability

For long MCMC jobs, future support for checkpointing would be useful.

Potential checkpoint contents:

```text
current iteration,
current marker effects,
current inclusion indicators,
variance components,
pi values,
RNG state if reproducibility is required,
chain diagnostics,
partial posterior sums.
```

This is a larger design task and should not be mixed with formatter cleanup.

### 9. Better runtime and memory diagnostics

Large-scale runs should report or store:

```text
runtime per trait,
runtime per chain,
number of markers,
number of LD nonzeros,
memory estimates,
OpenMP thread count,
backend name,
scheduler type,
acceptance rates where relevant.
```

### 10. More explicit model metadata

The formatted fit should ideally contain enough metadata to identify:

```text
model family,
backend,
prior type,
BayesC/BayesR/SBayesRC,
CSR/BED,
annotation/group/prior/component mode,
scheduled/non-scheduled,
chains,
niter,
nburn,
ncores,
seed if available.
```

This will make downstream summaries and debugging much easier.

## Suggested immediate next step

Before full-scale production, run a moderate validation dataset and save:

```text
the fit object,
check_stblr_consistency(fit),
sessionInfo(),
git commit hash,
model call,
runtime,
memory use,
key posterior summaries,
warnings.
```

Only after this passes should the package be run on full genome-scale data.

