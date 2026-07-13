# Unified BLR Framework: Phase 1 Report

## 1. Executive summary

Phase 1 introduced internal resolved R specifications, binding-neutral typed
C++ specifications and result vocabulary, an explicit CSR ownership contract,
a thin Rcpp specification-validation boundary, deterministic regression
fixtures, and reproducible runtime and memory baselines.

The existing production samplers and public execution routes were not changed.
In particular, unscheduled CSR BayesC remains the production implementation;
the new boundary validates and round-trips specifications but never invokes a
sampler.

## 2. Repository baseline

- Branch: `master`.
- Starting commit: `7bdea9fdf30a6faae6810456a7dff8808859dc2f`
  (`7bdea9f Plan for version 2`).
- Initial working tree: one pre-existing untracked file,
  `docs/dev/blr_framework_phase1_report.md`, containing the superseded compiler
  blocker report. No tracked changes were present.
- Initial `git diff --check`: passed.
- Recent commits: `7bdea9f`, `afb2c2e`, `1b3e8d6`, `adac7c6`, and `56be6dc`.

The environment uses R 4.4.1 on Windows 11 with Rtools44. The package uses the
GNU C++17 standard selected by this R toolchain. Native compilation succeeded.
The initial compiled load took approximately 300 seconds, which was a
successful long-running compilation rather than a timeout or compiler failure.
The full baseline suite passed:

```text
FAIL 0
WARN 0
SKIP 0
PASS 3181
```

## 3. Files inspected

The following design and repository files were inspected before editing:

- `AGENTS.md`;
- `docs/dev/blr_framework_implementation_plan.md`;
- `docs/dev/blr_model_capability_matrix.md`;
- `docs/dev/blr_reduction_test_matrix.md`;
- `docs/dev/blr_framework_phase0_audit.md`;
- `docs/dev/stblr_raw_schema.md`;
- `docs/dev/stblr_backend_naming.md`;
- `docs/dev/stblr_backend_computation_inventory.md`;
- `docs/dev/archive/stblr_csr_multichain_design.md`;
- `DESCRIPTION`, `src/Makevars`, and `src/Makevars.win`;
- existing CSR interface, raw-schema, backend-consistency, multichain,
  `selection_s`, and LD-swap tests;
- the public CSR wrapper and relevant generated Rcpp registration files.

Protected sampler files were inspected only as needed to identify current
contracts and fixtures. They were not edited.

## 4. Files changed

- `R/blr-model-spec.R`: adds internal R constructors, validation, and internal
  wrappers for the C++ validation bridge.
- `src/blr_spec.h`: adds binding-neutral typed specifications and validation.
- `src/blr_result.h`: adds the typed result vocabulary and dimension
  validation.
- `src/blr_csr_contract.h`: adds the shared, read-only CSR resource contract.
- `src/blr_phase1_rcpp.cpp`: adds the thin Rcpp conversion and validation
  boundary; it contains no sampler call.
- `R/RcppExports.R` and `src/RcppExports.cpp`: generated registrations for two
  internal Rcpp functions. `NAMESPACE` was not changed.
- `tests/testthat/test-blr-framework-phase1.R`: adds specification, source,
  dimension, ownership, and production-regression tests.
- `tools/benchmarks/blr_phase1_csr_bayesc.R`: adds the reproducible tiny and
  moderate baseline driver.
- `docs/dev/blr_framework_phase1_report.md`: replaces the pre-existing,
  incorrect untracked blocker report with this validated Phase 1 record.

## 5. R resolved specification

The internal `blr_resolved_spec` is built from explicit `data`, `model`,
`mcmc`, `output`, and `execution` sections plus schema name and version. The
first supported combination is CSR data, independent traits, scalar BayesC,
binary state, global binary probability, unit scale, scalar-independent trait
and residual covariance, and the unscheduled `stblr_cpg_omp_csr` backend.

Validation checks dimensions, canonical marker and trait identifiers, sample
sizes, MCMC controls, explicit chain-seed length, output flags, schema values,
operator/backend tags, and the CSR resource contract. Errors identify the
invalid field. Unsupported models, operators, scheduled BayesC, and
multivariate covariance policies are rejected.

The specification stores metadata and identifiers, not LD values or large
statistics arrays. Marker and trait order are preserved. It is internal,
unexported, is not attached to fit objects, and is deliberately not connected
to the current public route.

## 6. Typed C++ specification

`sblr::core` defines enum classes for data representation, design, scaling,
kernel, model family, state space, probability and scale policies, covariance
policies, operator, and backend. `DataSpec`, `ModelSpec`, `McmcControl`,
`OutputSpec`, `ExecutionSpec`, and `ResolvedSpec` encode the supported contract.

Validation uses standard C++ exceptions. The core headers use only standard
library types and contain no R, Rcpp, Armadillo, Python, or binding-specific
types or calls. They compile independently with the repository's GNU C++17
configuration and do not require newer language features.

## 7. Typed result vocabulary

The result vocabulary comprises `MarkerSummary`, `StateSummary`,
`TraceSummary`, `VarianceSummary`, `DiagnosticsSummary`, `ChainSummary`,
`OptionalResultDimensions`, and `BlrResult`. It is a contract only; no existing
sampler returns these structures in Phase 1.

Validated canonical dimensions are:

- marker effects and marker PIP: markers by traits;
- traces: retained samples by parameter dimension;
- trait covariance: traits by traits;
- optional component probability: markers by components by traits;
- optional pattern probability: markers by patterns.

Optional future result families are represented only by dimension metadata,
avoiding speculative large containers. The result types are independent of R
list layout.

## 8. CSR ownership contract

`CsrResourceSpec` records the resource identifier, marker count, shared
read-only status, absence of per-chain data payload, and the requirement that
storage outlive all chain executions. Validation rejects an empty resource,
zero marker count, mutable/non-shared storage, per-chain data ownership, or an
insufficient lifetime.

The source contract states that chain state owns only mutable effects,
residuals, parameters, accumulators, RNG state, and workspace. Marker order and
scaling remain data-contract properties. The current CSR reader and operator
remain in their existing implementation; no second reader was added and no
storage was moved.

## 9. Rcpp conversion boundary

The internal bridge converts an R resolved specification to the typed C++
structures, validates it, and returns a normalized R representation. The round
trip preserves enum meanings, marker and trait counts and order, all MCMC
values, optional explicit chain seeds, output flags, CSR ownership fields, and
operator/backend tags.

A separate internal bridge exercises typed result-dimension validation. Both
functions are registered only through internal generated R wrappers. Neither
function invokes, dispatches, or references a sampler.

## 10. Reference fixtures

The regression tests use the existing tiny CSR construction path with three
markers and one trait. They exercise fixed seeds, explicit chain seeds,
one-chain and multichain execution, one and two cores where supported, default
and fixed-disabled `selection_s`, and disabled LD-swap.

The protected identities are exact repeated-call equality, one-core/two-core
equality, explicit-chain-seed behavior, marker and trait order, stable raw and
formatted dimensions, absence of Phase 1 fields from public fits, and an
unchanged `stblr_raw_v1` schema. All reference checks passed through the
existing public and backend paths.

## 11. Performance and memory baseline

Commands:

```text
Rscript tools/benchmarks/blr_phase1_csr_bayesc.R --fixture=tiny --nrep=3 --chains=all --cores=all --output=all --peak=TRUE
Rscript tools/benchmarks/blr_phase1_csr_bayesc.R --fixture=moderate --nrep=3 --chains=all --cores=all --output=all --peak=TRUE
```

The tiny fixture has 3 markers, 1 trait, 8 sampling iterations, 2 burn-in
iterations, and 8 retained samples. Timings ranged from 0.00 to 0.59 seconds
and were frequently below clock resolution. Sampled whole-process peak RSS
ranged from approximately 142.4 MB to 143.9 MB. This fixture is useful for
correctness, not performance inference.

The moderate fixture has 2,000 markers, 2 traits, 80 sampling iterations, 20
burn-in iterations, and 80 retained samples. Across one or two chains, one or
two cores, and minimal or ordinary output, individual timings ranged from 0.20
to 2.89 seconds. Configuration means ranged from 0.21 to 1.46 seconds. The
first configuration showed substantial warm-up/order variability
(0.65--2.89 seconds), so these short runs do not support claims about relative
speed. Sampled whole-process peak RSS ranged from approximately 144.7 MB to
149.5 MB.

Peak memory was measured in a child R process using installed optional
`processx` and `ps` packages, polling RSS every 10 milliseconds. It includes R,
package loading, fixture construction, and retained output; it is not a precise
sampler-only instantaneous peak. No required benchmark dependency was added.
The environment was R 4.4.1 on Windows 11, with 12 logical processors reported,
Rcpp 1.1.1, RcppArmadillo 15.2.3-1, pkgload 1.4.1, and testthat 3.3.2.

## 12. Test results

Baseline validation:

- `Rcpp::compileAttributes(".")`: passed with no generated changes.
- `pkgload::load_all(".", compile = TRUE)`: passed after approximately 300
  seconds; compilation emitted existing native compiler warnings but no
  compiler error.
- focused CSR/schema/consistency/selection/LD-swap tests: passed.
- full baseline suite: 3,181 passed, 0 failed, 0 warned, 0 skipped (45.0
  seconds).

Post-change validation:

- standalone GNU C++17 syntax checks for all three binding-neutral headers:
  passed;
- final `Rcpp::compileAttributes(".")`: passed;
- final `pkgload::load_all(".", compile = TRUE)`: passed after 225.4 seconds,
  a successful long-running native build;
- new Phase 1 tests: 102 passed, 0 failed, 0 warned, 0 skipped (2.6 seconds);
- combined focused Phase 1, CSR BayesC, multichain, raw-schema,
  backend-consistency, `selection_s`, LD-swap, and public-interface tests:
  1,222 passed, 0 failed, 0 warned, 0 skipped (17.4 seconds);
- full post-change suite: 3,283 passed, 0 failed, 0 warned, 0 skipped (40.2
  seconds).

R reported that the installed testthat package was built under R 4.4.3. This
external package-version notice was outside testthat's suite counts; the suite
itself reported zero warnings.

## 13. Public behavior statement

No production sampler source changed. No public execution route, function
argument, default, result field, `stblr_raw_v1` field, marker traversal order,
burn-in or thinning rule, seed formula, RNG call order, OpenMP scheduling,
chain aggregation, residual or variance update, LD-swap behavior, or
`selection_s` behavior changed. `R/sparse_ld_bed_helper.R` and all protected
sampler files remain untouched.

## 14. Deviations and blockers

There are no remaining Phase 1 blockers.

During development, the initial moderate benchmark invocation exited nonzero
with `stblr_cpg_omp_csr: inconsistent trait dimensions.` The full command and
diagnostic were captured. The cause was confined to the new benchmark fixture:
it supplied a per-trait sample-size vector where the current wrapper expects a
shared scalar. Correcting that fixture made the same benchmark matrix pass;
production code was not changed.

An early regression assertion expected a changed seed to alter marker means,
but the original weak tiny scores produced an all-zero state. The fixture was
strengthened, after which deterministic repeated-seed and explicit-seed checks
passed. This was a test-fixture correction, not a sampler change.

The short benchmark timings and sampled whole-process RSS are baselines with
the limitations stated above, not precision performance comparisons.

## 15. Recommended Phase 2 boundary

Migrate the existing unscheduled CSR BayesC kernel into the new typed
boundaries while preserving its mathematics, RNG ordering, speed, and memory
use.

## 16. Readiness marker

PHASE 1 COMPLETE — CONTRACTS AND REGRESSION BASELINES VALIDATED
