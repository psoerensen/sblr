# Phase 18 — unified STBLR/MTBLR alignment and cleanup

## 1. Executive summary

Phase 18 aligns the public STBLR and MTBLR boundaries around canonical family,
model, and operator identifiers; common chain/convergence controls; explicit
fit names; one scalar convergence engine; and one metadata/memory vocabulary.
The final correction establishes six compositional scientific models and one
common five-trace trait variance decomposition without duplicating S kernels or
changing the convergence mathematics.

## 2. Repository baseline

The clean baseline was `230640cb330fe8ef89fd1ac20c6a5241658a3e07`
on `master`, using R 4.4.1 UCRT, GCC/G++/Fortran 13.2, Make 4.4.1,
OpenMP with 12 reported processors, and the reference R BLAS/LAPACK. The
baseline source suite reported 5,602 expectations: 5,600 passed and two
established opt-in cases skipped. Package checking had zero errors and zero
warnings; the three pre-existing notes concerned the Phase 17C source path,
installed size, and legacy scalar `std::cout` output.

## 3. Architecture before Phase 18

The package mixed method casing, overlapping wrappers, family-specific output
aliases, different convergence availability, and an absent public scalar
block-eigen entry. The detailed inventory is in `blr_unified_architecture.md`.

## 4. Final model/operator matrix

Canonical public routes are STBLR BayesC/SBayesC/BayesR/SBayesR CSR and block
eigen, SBayesRC annotation CSR/block eigen, BayesC/BayesR/BayesRC BED, the
specialized annotation CSR policies, and MTBLR BayesC CSR/block eigen/BED.
Unsupported cells fail before numerical execution.

## 5. Naming conventions

Machine identifiers are `stblr`/`mtblr`,
`bayesc`/`sbayesc`/`bayesr`/`sbayesr`/`bayesrc`/`sbayesrc`, and
`csr`/`block_eigen`/`packed_bed`/`dense_reference`. S models reuse existing
kernels with explicit `maf_s`; probability policies are `global`,
`fixed_marker`, `group`, `learned_logistic`, and
`annotation_probit_stick`.

## 6. Public API before and after

The canonical exports are `stblr_csr`, `stblr_block_eigen`, `stblr_bed`,
`stblr_csr_annot`, `mtblr_csr`, `mtblr_block_eigen`, and `mtblr_bed`.
Redundant model wrappers, the dense `sblr` dispatcher, the experimental BED
marker entry, and the legacy convergence screen are no longer exported.

## 7. Public method values

Only exact lowercase method values are accepted. No case-normalization aliases
remain at canonical public boundaries.

## 8. Common argument semantics

Canonical routes use `nit`, `nburn`, `nthin`, `seed`, `nchains`, `ncores`,
`chain_seeds`, `keep_chains`, `convergence`, `convergence_control`,
`memory_warning_gb`, and `verbose` with shared meanings and validation.

## 9. Thread terminology

`ncores` requests concurrent logical MCMC tasks; `nthreads` is reserved for
operator preparation, decoding, or other non-chain work.

## 10. Seed and logical-task policy

An ST task is trait × chain; an MT task is a complete joint chain. Seeds are
resolved before dispatch and never depend on worker assignment.

## 11. Compact-chain structure

`fit$chains` is a deterministic flat list. Records identify family, model,
operator, trait (NA for MT), chain, seed, retained count, and supported compact
posterior state. `keep_chains=FALSE` returns NULL.

## 12. Common fit layout

All canonical fits expose family/model/operator plus input, data, diagnostics,
convergence, convergence traces, chains, and memory estimate. Scientific
marker and trace fields remain directly accessible when defined.

## 13. Probability naming cleanup

`pi`, `pim`, and `pis` became `pi_final`, `pi_mean`, and `pi_trace`.
`comp_prob` became `component_probabilities`; old aliases are absent.

## 14. Covariance naming cleanup

MT posterior means use `cov_b_mean`, `cov_g_mean`, and `cov_e_mean`; primary
final states use `cov_b_final`, `cov_g_final`, and `cov_e_final`.

## 15. Chain-summary naming cleanup

Marker stability fields are `bm_chain_mean_sd|min|max` and
`dm_chain_mean_sd|min|max`, explicitly distinguishing chain-mean variation
from posterior uncertainty.

## 16. Formatter consolidation

`.as_stblr_fit()` and `.as_mtblr_fit()` remain the family formatters.
`.blr_finalize_fit()` owns common metadata, class, and public-name cleanup.

## 17. Shared convergence engine

The dependency-free `.blr_convergence_*` engine is shared. The Phase 17U rank,
folded R-hat, ESS, posterior-SD, and mean-MCSE mathematics were not changed.

## 18. Trace-bundle contract

Adapters produce `blr_convergence_trace_bundle` version 1 arrays ordered by
iteration, logical chain, and scalar descriptor with family/model/operator and
identity metadata.

## 19. STBLR convergence integration

Canonical ST routes diagnose the available unpooled post-burn `vbs`, `vgs`,
`ves`, `vle`, and `vld` traces independently of posterior-summary thinning and
compact-chain retention.

The same names now define MT core convergence. CSR, block eigen, and packed BED
capture chain-private `vle` from current effective effects and the represented
operator diagonal, then set `vld = vgs - vle` at the same completed-iteration
checkpoint as `vbs`/`vgs`/`ves`. Public traces are iteration × trait and
convergence arrays are iteration × chain × quantity. Full MT covariance means
and final states remain separate `cov_*_mean`/`cov_*_final` matrices.

## 20. MT CSR multichain/convergence alignment

`mtblr_csr_chains_raw_internal()` reads and constructs each owned CSR operator
once, creates immutable views shared by all logical chains, and dispatches one
complete joint model per chain with static OpenMP scheduling. Mutable sampler
state and RNG are chain-private. Results return in chain order; the R adapter
pools retained draws, retains chain 1 final state, and computes convergence
from the unpooled native chain traces without rerunning a chain.

## 21. MT block-eigen multichain/convergence alignment

`mtblr_block_eigen_chains_raw_internal()` parses each descriptor, reads its BED
provenance, reconstructs each owned operator, and transforms scores once before
the same deterministic logical-chain dispatcher runs. It uses the same seed,
pooling, primary-state, convergence, output, and memory contracts as MT CSR.

## 22. ST block-eigen public route

`stblr_block_eigen()` exposes only validated BayesC, BayesR, and SBayesRC
routes with one-based blocks, provenance, filter metadata, and common controls.

## 23. Convergence warning behavior

Warnings are one main-thread aggregated advisory per fit at most. Auto
single-chain unavailability is quiet; details remain in `fit$convergence`.

## 24. Memory accounting

Analytical estimates separate shared operator data, private worker state,
per-chain results, convergence capture/workspace, retained chains/traces, and
formatted output. They are not measured RSS or peak RSS.

## 25. CSR/block-eigen reductions

Permanent owners execute BayesC, BayesR, and SBayesRC fits against unfiltered
CSR and block-eigen representations. BayesC agrees within `1e-7`, BayesR
within `1e-4` for scalar fields and `1e-8` for component probabilities, and
SBayesRC agrees within `1e-12`. The BayesC convergence summaries agree within
`1e-12`.

## 26. BED/summary reductions

The executable standardized-data null-state BayesC reduction agrees within
`1e-12` for marker means, inclusion probabilities, variance traces, and final
probability. Joint full residual covariance is not identified by marginal
summary inputs and remains explicitly non-comparable.

## 27. One-trait ST/MT reductions

The executable matched null-state one-trait BayesC reduction agrees within
`1e-12` for `bm`, `dm`, and `vgs`, and exactly for `pi_final`. ST and MT
`vbs`/`ves` intentionally differ because the scalar and joint covariance prior
parameterizations are not identical; the owner asserts that distinction.

## 28. Public API protection

Unsupported method/operator combinations fail explicitly and experimental
low-level routes are not namespace exports.

## 29. Reproducibility

Logical seeds are invariant to worker count. Deterministic task order and
chain-private RNG preserve repeated and serial/OpenMP equality where native
parallel dispatch is available.

## 30. Test ownership

Permanent owners are `test-blr-unified-public-contract.R`,
`test-blr-unified-convergence.R`, `test-blr-operator-reductions.R`, and
`test-blr-unified-reproducibility.R`; historical scientific owners remain.

## 31. Architecture audit

`blr_phase18_unified_architecture_audit.R` checks exports, capability and
naming contracts, formatter/convergence ownership, one-time MT operator
preparation, deterministic native dispatch and seeds, scheduled trace capture,
retention/thinning independence, executable reductions, memory ownership, and
raw-schema stability. All guards pass.

## 32. Mutation sensitivity

`blr_phase18_mutation_sensitivity.R` guards 60 failure modes, including alias
reintroduction, seed/trace mistakes, duplicate mathematics, per-chain operator
reconstruction, scheduled trace loss, memory terminology, missing executable
reductions, and scope leakage. All 60 guards pass.

## 33. Workflow and documentation migration

Maintained workflows use lowercase methods, canonical routes, explicit output
names, and the common fit layout. `07_unified_operators_workflow.R` provides a
single operator/family map.

## 34. Focused tests

The initial incomplete validation exposed the four recovery groups: serial
R-level MT chain orchestration, absent scheduled-CSR task traces,
classification-only reductions, and historical assertions tied to removed
aliases and fields. Recovery used focused owners after each production change,
then migrated the affected historical files individually. The unified public,
convergence, reproducibility, and executable reduction owners now pass, as do
the affected BayesR, annotation, raw-schema, LD-swap, BED, CSR, and Phase 17
scientific owners.

The historical migration touched 62 tracked test/fixture/helper files. The
failures were classified as removed exports, renamed arguments and output
fields, canonical nesting, flat compact-chain ordering, additive convergence
objects or warnings, and lowercase model/operator identifiers. No failure was
resolved by restoring a removed alias. Focused recovery owners and the final
fast tier completed without failures or warnings.

## 35. Fast tier

The updated fast filter includes public contracts, shared convergence, small
integrations, executable reductions, naming guards, and package checking. It
completed with 2,834 passed expectations, zero failures, zero warnings, and
one established opt-in peak-RSS skip.

## 36. Full source suite

The first post-implementation run was intentionally retained in this report:
it completed in 159.6 seconds and exposed more than 400 historical/API
migration failures plus 30 advisory warnings. Subsequent recovery runs exposed
the remaining Phase 9A full-object comparisons, annotation-chain wrapper and
nesting assumptions, 56 LD-swap/diagnostic API assertions, and six finemap
fixtures using the old chain-summary names. These failures were classified and
migrated rather than hidden or skipped. The definitive final source run completed
with 5,695 passed expectations, zero failures, zero warnings,
zero errors, and two established opt-in skips (`SBLR_RUN_EXTENDED_REPRODUCIBILITY`
and `SBLR_RUN_PEAK_RSS`).

## 37. Installed tests

Portable Phase 18 owners are included in built-package testing; source-only
hash and architecture assertions use narrow source-checkout guards. The final
built-tarball check ran the installed `testthat.R` suite successfully.

## 38. Package check

The initial built-tarball check completed with **1 ERROR and 3 NOTEs** because
the installed suite still contained 12 historical API failures and two test
warnings (`4,595` passes and `110` source-only or opt-in skips). After migration,
two sandboxed helper-driven checks encountered Windows process-pipe or temporary
installation errors; equivalent commands outside that restricted execution
context completed normally. The definitive built-tarball `R CMD check` then
completed in 701.1 seconds with **0 ERRORs, 0 WARNINGs, and 3 classified NOTEs**:
the established long fixture path, installed package size (5.4 MiB, 4.3 MiB
under `libs`), and pre-existing native `std::cout` symbols. Examples, installed
tests, documentation, loading, unloading, and all other check stages passed.

## 39. Diff and artifact hygiene

The protected Phase 17B fixture and tracked reference tarball have the exact
starting-commit blob hashes. Compiler/DLL, diagnostic check-directory, and
source-tarball artifacts produced during validation were removed with targets
resolved inside the repository. Final `git diff --check`, whitespace/EOL, and
artifact results are recorded in the completion handoff.

## 40. Deviations and blockers

The four recovery blocker groups are resolved. Fixed/learned/group annotation
raw-field inventory tests explicitly disable convergence because they own raw
payload structure rather than trace capture; the dedicated convergence owner
continues to exercise genuine task-private traces. The intermediate sandbox-only
process-pipe and temporary-installation failures were execution-environment
limitations; the definitive unrestricted built-tarball check passed. No
implementation blocker remains.

Before Phase 19 was committed, the temporary Phase 18 interpretation of the
`s` prefix as a MAF-S scale marker was corrected. The final convention uses
the prefix exclusively for summary-statistics data; `selection_s` independently
controls MAF scaling. Existing numerical kernels were reused unchanged, and
new fits record semantic version 2, prior kernel, data level, and effect-scale
policy separately.

## 41. Recommended next phase

> implement MTBLR BayesR as one coherent pattern-by-mixture model across the aligned CSR, block-eigen, and packed-BED operators, reusing the Phase 18 chain, convergence, output, naming, memory, and warning infrastructure.

## 42. Readiness marker

PHASE 18 COMPLETE — STBLR AND MTBLR ARCHITECTURES ALIGNED
