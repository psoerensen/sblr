# Public individual-level multivariate packed-BED contract

## 1. Purpose

`mtblr_bed()` fits public joint individual-level multivariate BayesC models
from one BED-backed genotype list.

## 2. Public function

```r
mtblr_bed(
  y, Glist, covar = NULL, chr = NULL, cls = NULL, rows = NULL,
  scale = TRUE, center = TRUE,
  residual_covariance = c("full", "diagonal"), method = "bayesC",
  trait_metadata = NULL, sets = NULL, block_size = 1000,
  beta = NULL, b = NULL, state = NULL, h2 = 0.5, pi = 0.001,
  models = NULL, pimodels = NULL, vg = NULL, vb = NULL, ve = NULL,
  ssb_prior = NULL, sse_prior = NULL, updateB = TRUE, updateE = TRUE,
  updatePi = TRUE, nub = 4, nue = 4, nit = 1000, nburn = 500,
  nthin = 1, seed = 1, memory_warning_gb = 8, verbose = FALSE
)
```

## 3. Numerical owner

Phase 17O `mtblr_bed_internal()` is the sole numerical route. The public
adapter prepares, calls it once, enriches its named raw result, and formats it.

## 4. Difference from stblr_bed

`mtblr_bed()` fits one joint MT likelihood with joint inclusion patterns and
effect/residual covariance. `stblr_bed()` fits existing separate scalar
trait-specific likelihoods.

## 5. Glist

Exactly one BED-backed `Glist` is required, with BED files, sample count,
individual and marker identities, selected columns, and allele frequencies.

## 6. Phenotype forms

`y` may be a numeric vector, numeric matrix, or all-numeric data frame. Vector
names and matrix/data-frame row names identify individuals. Missing trait names
become `T1`, `T2`, and so on; final trait names must be unique and nonempty.

## 7. Individual alignment

The established `stblr_bed()` preparation helper owns ID matching, explicit
rows, unmatched-ID warnings, duplicate rejection, and order. One selected row
set and its phenotype order are shared by every trait.

## 8. Complete-data policy

After alignment, `y` must have more than one row, at least one trait, finite
values, and positive finite variance in every column. Missingness, trait masks,
imputation, and automatic complete-case dropping are unsupported.

## 9. Phenotype centering

With `center=TRUE`, aligned column means are subtracted in R. With
`center=FALSE`, each aligned column must already meet
`abs(mean) <= 1e-10 * max(1, RMS)`. Means, variances, tolerances, and status are
recorded.

## 10. Phenotype scaling

Phenotypes are not scaled; their post-centering units are retained.

## 11. Covariates

Covariates are not fitted or projected. `covar` must be `NULL`; users may pass
phenotypes already adjusted for desired covariates.

## 12. Marker selection

The existing BED helper resolves `chr`, `cls`, default `rsidsLD`, physical BED
columns, names, and order. The adapter never changes physical marker order.

## 13. Frequencies

Selected `Glist$af` values are required, finite, and strictly inside `(0,1)`.
They are not computed, replaced, or clipped by the adapter.

## 14. Genotype scale

Only `scale=TRUE` is supported: standardized genotypes with missing genotype
codes mean-imputed after centering (therefore decoded as zero).

## 15. Models

Joint binary inclusion patterns and probabilities use `.mtblr_models()`:
unique patterns, a required null pattern, stable names, normalized
probabilities, full/restrictive defaults, and the existing pattern limit.

## 16. Sets

Explicit one-based marker sets use `.mtblr_sets()` and must be a disjoint
complete partition. Defaults come from the established BED helper's block/file
assignment without changing marker order.

## 17. Initialization

`beta`, `b`, and `state` accept marker-by-trait matrices or trait lists. Missing
`b` defaults to zero, missing `state` to `b != 0`, and missing `beta` to `b`.
State rows must be supplied model patterns and `b = state * beta`; inactive
latent values may remain nonzero.

## 18. Priors

Centered phenotype variances and `h2` define diagonal starting defaults.
`p_active` is one minus null-pattern probability. Defaults are
`vg=diag(vy*h2)`, `ve=diag(vy*(1-h2))`, `vb=vg/(m*p_active)`,
`ssb_prior=((nub-2)/nub)*vg/(m*p_active)`, and
`sse_prior=((nue-2)/nue)*ve`. Shared covariance validation requires finite SPD
matrices and degrees of freedom greater than `max(2, nt-1)`.

## 19. Full residual covariance

`residual_covariance="full"` is the default canonical same-individual
likelihood. Full SPD `ve` and `sse_prior` are supported.

## 20. Diagonal covariance

`"diagonal"` is an explicit reduction/testing model and requires exactly
diagonal `ve` and `sse_prior`; it is not an overlap approximation.

## 21. Execution

Execution is serial, one chain, BayesC method 4, with one Phase 17O native call
and one fit-local RNG owned by that route.

## 22. Raw output

The backend is `mt_bed_bayesc`, data level is `individual`, and output remains
named `mtblr_raw` schema version 1. Ordinary `wy` and `r` mean `X'Y` and final
`X'R`.

## 23. Formatted fit

The shared `.as_mtblr_fit()` creates `class = c("mtblr_fit", "list")`.
Phase 17P adds `bed_diagnostics`, `phenotype_preprocessing`, and
`memory_estimate`.

## 24. Diagnostics

MT BED diagnostics identify the Phase 17O owner/view, scale, covariance mode,
dimensions, workspace, Cholesky safeguards, E-update counts, and unsupported
sample residual, genetic-value, CPO, and LE/LD outputs.

## 25. Alignment metadata

Alignment records the individual policy, input/selected counts, row status and
rows, sample order, unmatched IDs, duplicate policy, complete phenotype,
centering, marker selection/order, by-construction orientation, standardized
scale, and covariate policy.

## 26. Memory estimate

The adapter reports component formulas and total GiB for packed genotypes,
phenotype/residual/effects/state/workspaces/maps/order/covariance/models,
marker summaries, and traces. It is an analytical working-memory estimate,
not measured RSS or measured peak RSS. Crossing `memory_warning_gb` warns but
does not stop execution.

## 27. Unsupported capabilities

Covariate fitting, missing phenotypes, CPO, LE/LD decomposition, individual
predictions, sample residual/genetic-value output, OpenMP, multichain,
BayesR/BayesRC, annotations, scheduling, and legacy eigen migration are absent.

## 28. Evolution

Adding these capabilities requires separate contracts for preprocessing and
evidence, observation masks, prediction definitions, chain ownership/seeds,
aggregation and diagnostics, retained-chain memory, and parallel dispatch.

## Phase 17Q multichain boundary

The public interface remains serial and unchanged. Phase 17Q specifies future
chain-only tasks, deterministic seeds, pooled summaries, primary-chain final
state, optional compact chains, and analytical worker-aware memory. These are
contracts for Phase 17R/17S, not current capabilities.

Phase 17R activates those contracts internally only. Public `mtblr_bed()` still
has no chain controls and continues to call the exact serial Phase 17O route.

## Phase 17S public multichain activation

The supported signature is:

```r
mtblr_bed(y, Glist, covar = NULL, chr = NULL, cls = NULL, rows = NULL,
  scale = TRUE, center = TRUE,
  residual_covariance = c("full", "diagonal"), method = "bayesC",
  trait_metadata = NULL, sets = NULL, block_size = 1000,
  beta = NULL, b = NULL, state = NULL, h2 = 0.5, pi = 0.001,
  models = NULL, pimodels = NULL, vg = NULL, vb = NULL, ve = NULL,
  ssb_prior = NULL, sse_prior = NULL,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE, nub = 4, nue = 4,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL, keep_chains = FALSE,
  memory_warning_gb = 8, verbose = FALSE)
```

One task is one complete joint-MT model. Defaults request one chain and one
worker, use native default seeds, and omit compact records. `seed` remains the
base seed; explicit signed 32-bit `chain_seeds` retain supplied order and may be
negative. Native default resolution is modulo-2^32 `seed + 9176*c`.

OpenMP is static and chain-level only; workers are capped by `nchains`. Without
OpenMP, a multi-core request warns once and runs serially. The package never
changes global BLAS settings, and one-thread BLAS is normally recommended for
concurrent package workers.

Marker, covariance, and probability posterior means pool retained samples.
Traces are iterationwise chain means. Final effects, binary state, residual
scores, covariance matrices, and probabilities come from primary chain 1;
binary states are not averaged. Stability fields are sample SD/min/max across
per-chain posterior means. They are not posterior standard deviations, R-hat,
ESS, or MCSE, and multiple chains do not establish convergence.

Compact retained records include marker means/final state, traces, covariance
summaries/finals, probabilities, seed, and diagnostics. They omit marker
residuals, `X'Y`, order, models, sets, phenotype and sample residual matrices,
genetic values, packed bytes, maps, and BED metadata.

The analytical memory estimate separates shared immutable bytes, private bytes
per worker, results per chain, pooled output, and optional compact output per
chain. It is an analytical upper-bound estimate, not measured RSS or measured
peak RSS.

## Phase 17T convergence boundary

Formal convergence diagnostics are not yet implemented. Future Tier 1 uses
original per-chain post-burn B/G/E diagonal traces independently of
`keep_chains`; pooled traces and marker stability are not R-hat or ESS inputs.
## Phase 17U status

An internal Tier 1 convergence engine exists, but public activation is deferred.
`mtblr_bed()` has no convergence formals, does not call the trace route, emits
no convergence warning, and returns no convergence fields in Phase 17U.

## Phase 17V convergence controls

The public signature adds `convergence=c("auto","none","core")` and
`convergence_control=NULL` after `keep_chains`. Automatic multichain fits
compute B/G/E-diagonal diagnostics; automatic single-chain fits report quiet
unavailability; `none` disables capture. Optional post-burn traces are separate
from compact chains. Memory and at most one advisory warning include requested
diagnostic scope. Passing diagnostics does not prove convergence.

## Phase 17W future extended diagnostics

Current modes remain `auto`, `none`, and `core`. A future explicit `extended`
mode will use nested failure-closed covariance, probability, and selected-marker
controls; `auto` remains Tier 1 only. Phase 17W changes no public surface.
## Phase 18 common fit boundary

MT BED retains its validated numerical and convergence routes, while its
formatted result now follows the family-neutral metadata ownership and naming
contract. Covariance mean/final states, probability mean/final/trace fields,
and chain-mean stability summaries use explicit names. The shared convergence
engine is a nonnumerical generalization of the Phase 17U implementation.
