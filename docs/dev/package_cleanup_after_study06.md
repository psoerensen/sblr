# Package cleanup after Study 06

## Scope and provenance

This cleanup establishes the auditable `sblr` 0.2.0 baseline after the
continuous-alpha sampler-development endpoint. It does not change a likelihood,
prior, posterior target, sampler transition, residual policy, public default,
or formatted result contract.

| Item | Recorded value |
|---|---|
| Starting source HEAD | `3c7b97f3de76a6e19a7e82bf73b2b2b10bf83d34` |
| Branch | `master` |
| Source package version | `0.2.0` |
| Installed package before cleanup | `0.2.0` at `C:/Users/au223366/AppData/Local/R/win-library/4.4/sblr` |
| R and toolchain | R 4.4.1 UCRT; GCC/G++ 13.2.0; OpenMP enabled |
| BLAS/LAPACK | R 4.4.1 Windows built-in BLAS/LAPACK |
| Read-only `sblrbench` HEAD | `21ef774ec084d11932e132186ccc085054b9aa56` |
| Scientific endpoint | PMA-R3; ordinary production sampler retained |

Both repositories were clean at the precondition gate. `sblrbench` was read
only throughout this task.

## File classification

The audit used the following durable categories. The table groups files with a
common role; detailed Study 06 evidence remains in the individual documents.

| Category | Principal files | Cleanup disposition |
|---|---|---|
| Production | `R/stblr-csr-sbayesrc.R`, `R/stblr-bed-bayesrc-internal.R`, `R/stblr-block-eigen.R`, annotation and block native owners under `src/` | Preserved; no scientific code or defaults changed |
| Internal production helper | `R/sbayesrc-helpers.R` and canonical raw-to-fit/schema helpers | Preserved; not conflated with reference samplers |
| Supported localization feature | public LD-swap controls and their native/tests | Preserved as established functionality |
| Development reference | `R/bayesrc-coordinated-transition-reference.R`, `R/sbayesrc-block-px-reference.R`, `R/sbayesrc-block-particle-reference.R`, `R/sbayesrc-particle-marginal-reference.R` | Retained as non-exported, non-production mathematical references with explicit file headers |
| Mathematical reference/oracle tests | coordinated, PX, particle-Gibbs, particle-marginal, coupling-ratio and pairwise tests under `tests/testthat/` | Retained to protect exact posterior identities; names and helpers identify reference status |
| Production contract tests | BayesRC annotation, block residual, trace, schema, formatter, LD-swap and ordinary RNG tests | Preserved |
| Study 06-specific tooling | `research/sbayesrc/tools/` | Retained as historical/development tooling, excluded from source builds, and never required by installed-package tests |
| Development documentation | Study 06 decision records and method derivations under `docs/dev/` | Preserved and linked from the endpoint navigation page |
| Obsolete/dead | None demonstrated | Nothing deleted |

## Installed-package test dependency

`tests/testthat/test-coupling-tempering-offline-ratios.R` sourced
`research/sbayesrc/tools/study06_partial_exchange_feasibility.R`. The research script is a large,
historical Study 06 audit and is intentionally not installed with the package,
so installed-package testing failed before any assertions ran.

The test now uses a small independent oracle in
`tests/testthat/helper-coupling-tempering-offline-ratios.R`. The helper contains
only the generic probability and exchange-ratio mathematics needed by the
test. It deliberately does not call production transition code, preserving the
scientific value of the independent check. The historical tool remains intact,
and no Study 06 logic was promoted into the public package API.

## Public API and ordinary RNG audit

`NAMESPACE`, roxygen exports, the public `stblr_*` entry points, and block/BED
wrappers expose no PX, particle, particle-marginal, coordinated, tempering,
partial-exchange, or pair-Gibbs sampler option. The reference functions are
dot-prefixed and unexported. Internal diagnostic switches remain disabled and
are not advertised as supported interfaces. Established LD-swap options are
unchanged.

The cleanup changes test support, comments, NEWS, and developer documentation
only. Production R/native code and function signatures are unchanged; no new
branch or RNG draw is introduced. Ordinary default draws therefore retain the
existing RNG contract.

## NEWS and build-content audit

`NEWS.md` now separates user-visible 0.2.0 changes from one concise
development/reference entry. Reference-method validation is not presented as a
supported sampler feature.

Repository-local `tools/` is explicitly excluded by `.Rbuildignore`: those
scripts remain versioned scientific provenance but are neither installed code
nor source-package test dependencies. The source-build audit must confirm that
the test helper and its test are shipped, `tools/` and ignored local Study 06
evidence are absent, and no compiled or check artifacts remain in the working
tree.

## Validation record

The final validation results are recorded here after running the cleanup gate:

| Check | Result |
|---|---|
| Focused offline-ratio test | Passed in source tree and from the installed package (13 expectations) |
| Reference and ordinary-RNG focused tests | Passed; explicit disabled tempering, block 1/1 scheduling, and disabled PX preserve ordinary output |
| `devtools::test()` | Passed with 0 failures, one opt-in skip, and the known covariance warning |
| `R CMD build .` and content inspection | Passed; 255 entries, helper/test present, `tools/`, local results, and compiled artifacts absent |
| `R CMD check --no-manual --as-cran` | Passed: 0 errors, 0 warnings, 5 notes; installed tests reported 4,792 passes, 3 expected skips, and 1 known warning |
| Installed version and production smokes | Installed `0.2.0`; installed BED BayesR, annotation backends, and block-eigen test smokes passed |
| `git diff --check` | Passed at the final scope review |

The five check notes were: new-submission metadata, 6.0 MB installed size,
unavailable clock verification, missing pandoc plus the existing top-level
`CLAUDE.md`, and the existing non-portable `-march=native` compilation flag.
The installed-package test skips were the CRAN-gated particle reference, the
source-tree-only architecture assertion, and the opt-in fresh-process
reproducibility check. The single warning is the known MT covariance diagnostic
warning for a deliberately non-positive fixture.

## Decision

**CLEAN-R1 — clean baseline established.** The package tests and installed
source check pass, the historical `tools/` dependency is removed, production
and reference roles are explicit, and the ordinary scientific implementation
is unchanged. This cleanup commit should become the baseline for separate
annotation-selection development after review. No Study 06 scientific fit was
run during validation.
