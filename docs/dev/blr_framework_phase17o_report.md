# Unified BLR Framework Phase 17O report

## 1. Executive summary

The internal serial one-chain individual-level multivariate BayesC packed-BED
core is active. It is reachable only through `mtblr_bed_internal()`; no public
adapter or export was added.

## 2. Repository baseline

Work began on clean `master` at
`524d7ca65082e85892bbf680be6d430c86c8788d`. The baseline used R 4.4.1 UCRT,
GCC/G++/Fortran 13.2.0, LAPACK 3.12.0, and the R default BLAS. Hosted CI status
was not visible from the local checkout. Baseline source tests reported 4669
passes, zero failures and two established opt-in skips. The built package had
zero errors, zero warnings, and the two established notes.

## 3. Phase 17N verification

| Phase 17N decision | Phase 17O implementation | Validation owner |
|---|---|---|
| `PackedBedMatrix` | one adapter-local move-only owner | Phase 17O architecture test/audit |
| `BedPackedGenotypeView` | one immutable borrowed view | Phase 17O architecture test/audit |
| standardized genotypes | one `MtBedMarkerMap` per marker | BED execution and identity tests |
| complete centered `Y` | central native validation | rejection tests |
| sample residual | `arma::mat n x nt` | dense-X identities |
| full `E` canonical | binding-neutral inverse-Wishart update | full-E tests |
| diagonal reduction | traitwise inverse-chi-square update | dense summary reductions |
| common MT result path | `MtDefaultCoreResult` and existing adapters | raw-schema tests |

## 4. Packed owner and view

The binding constructs exactly one
`PackedBedMatrix owner = read_bedfiles_to_packed_matrix(...)`, then one
`BedPackedGenotypeView<PackedBedMatrix>` borrowing `owner.data`. The owner is
not moved after view construction and outlives maps, preparation, MCMC,
finalization, and raw conversion.

## 5. Marker maps

For frequency `p`, BED codes map to
`(2-2p)/sqrt(2p(1-p)), 0, (1-2p)/sqrt(2p(1-p)), (0-2p)/sqrt(2p(1-p))`.
`xx` is the increasing-sample sum of squared decoded values. Frequencies and
positive finite `xx` are validated.

## 6. Data and state types

`MtBedDataView` borrows the packed view, maps, phenotype, marker `wy`, and
marker order. `MtBedInitialState` owns latent effects, effective effects,
binary states, `B`, `E`, and `pi`. `MtBedCoreDiagnostics` owns only numerical
counts. The posterior result is `MtDefaultCoreResult`.

## 7. Central validation

The core validates packed dimensions, maps, finite centered phenotypes,
marker order, unique binary models with a null pattern, normalized
probabilities, disjoint complete zero-based sets, SPD covariance/prior
matrices, degrees of freedom, diagonal-mode restrictions, consistent initial
latent/effective/state values, method, residual policy, and MCMC controls.
The binding separately validates paths, `cls`, one-based rows, frequencies,
phenotype dimensions, and unique nonempty trait names before reading.

## 8. Preparation order

The fixed order is BED read, view construction, marker maps, marker `X'Y`,
stable marker order, central validation, sample-residual rebuild, and finally
fit-local RNG construction. BED files are closed by the reader before MCMC.

## 9. Sample residual

The core owns column-major `arma::mat R`. It starts from `Y` and subtracts
decoded `x_j b_j'` in input marker order. Marker updates apply
`R -= x_j (b_new-b_old)'`.

## 10. Marker kernel

For `s=x_j'R+w b_j`, `Omega=E^-1`, `P=B^-1`, and model mask `D_k`, the
production kernel uses
`C_k=P+w D_k Omega D_k`, `rhs_k=D_k Omega s`, and
`log(pi_k)-log|C_k|/2+rhs_k'C_k^-1 rhs_k/2`. All models retain the full
latent dimension; zero probabilities have negative-infinite weight.

## 11. Marker sampling

Each visit decodes once, computes scores in trait order, evaluates the shared
kernel, draws one uniform in model order, draws `nt` normals in trait order,
solves `mean+L^-T z`, masks effective effects, updates `R`, and stores latent,
effective, and binary state. Null patterns still draw all latent coordinates.

## 12. Full E

Full mode forms the symmetrized `R'R`, adds `nue*sse_prior`, and draws
`IW(nue+n, S_post)` through the single binding-neutral `std::mt19937`
implementation in `blr_mt_covariance_rng.h`.

## 13. Diagonal E

Diagonal mode draws traits in increasing order using
`(sum(R_t^2)+nue*sse_prior_tt)/ChiSq(n+nue)` and zeros off-diagonals. This is
the deterministic corrected-summary reduction mode, not the canonical
same-individual likelihood.

## 14. B update reuse

For every set, the core calls existing `sampleBset()` then
`sampleB_latent()`, inverts the resulting `B`, and updates that set. After all
sets it calls existing global `sampleB()`. Disabled updates preserve supplied
`B`.

## 15. Pi and retention

`cmodel` starts at one per model, marker selections increment it, and existing
`samplePi()` runs after marker updates. Marker means are post-burn thinned;
covariance and probability accumulators preserve the corrected MT ownership
rules; traces have length `nit+nburn`.

## 16. Genetic covariance

Every iteration computes `U=Y-R` and the symmetrized `G=U'U/n`. No additional
centering or `n-1` denominator is used.

## 17. Final marker outputs

Preparation retains `wy=X'Y`. After the final state the reusable workspace
reconstructs `r=X'R_final` in marker, sample, and trait order. The sample-space
residual and genetic values are not returned.

## 18. Core result and finalization

The new core fills `MtDefaultCoreResult`, then calls exactly
`finalize_mt_default_result()`, `make_mt_default_legacy_result()`, and the
shared `mtblr_legacy_to_raw()`.

## 19. Native signature

The implemented signature is the Phase 17N signature documented in
`blr_mt_bed_internal_contract.md`, ending in `residual_covariance`, `nit`,
`nburn`, `nthin`, `seed`, and `method=4`.

## 20. Raw schema

The route returns `mtblr_raw` version 1 with backend `mt_bed_bayesc`, data
level `individual`, standard marker/trace/variance/model namespaces, and
`diagnostics$mt_bed`. Diagnostics name the owner/view/scale/workspace,
dimensions, covariance mode, jitter, E-update counts, and unsupported output
classes.

## 21. Deterministic marker oracle

Two- and three-trait full-`E`, diagonal-`E`, null, single-trait, shared,
correlated covariance, unequal-probability, and zero-probability cases agree
with the independent Phase 17N R oracle to `1e-12`.

## 22. Diagonal summary reductions

Fixed controls, B-only, E-only, pi-only, all updates, multiple sets, three
traits, and one trait reduce to corrected dense MT. Continuous fields meet
`1e-10`; states and marker order are exact.

## 23. Full-E scientific identities

Two-trait, three-trait, multiple-set, zero/nonzero initialization, fixed, and
all-update configurations satisfy `R=Y-Xb`, marker `r=X'R`, `G=U'U/n`, binary
states, normalized probabilities, and SPD final `B`/`E`. Full mode permits
off-diagonal residual covariance.

## 24. Initialization

Zero, active nonzero, and inactive latent/nonzero configurations are accepted
when consistent. Nonbinary states, absent patterns, inactive effective
effects, active latent/effective disagreement, dimension errors, and
nonfinite values are rejected. Residuals are always rebuilt internally.

## 25. Reproducibility

Repeated calls with one seed are exact, intervening summary fits do not alter
the result, and a different seed changes stochastic output. The core uses no
R RNG, random device, worker seed, OpenMP, or global RNG.

## 26. Ownership and memory

The audit reports one owner, one view, one RNG, zero per-chain packed bytes,
zero MCMC-time BED reads/file handles, and one fit-lifetime decoded workspace.
Analytical bytes are `m*stride` packed, `40m` maps, `8nnt` each for `Y` and
`R`, `8n` workspace, `8mnt` each latent/effective effects, and `4mnt` states.

## 27. Runtime

The Phase 17O benchmark records small/moderate diagonal/full total call
signals, object sizes, and analytical storage. Binding-internal preparation,
MCMC, and final-reconstruction times are honestly marked not separately
observable; tiny timings are not performance claims.

## 28. Summary-MT protection

Dense Phase 17C, CSR Phase 17I/17J, and block-eigen Phase 17L/17M permanent
tests pass unchanged. Summary likelihood representation seams were not edited.

## 29. Scalar BED protection

Canonical packed-BED BayesC, BayesR, and BayesRC code was not edited and its
permanent owners pass unchanged.

## 30. Research-route protection

The new headers and binding do not call `mtblr_eigen()`.

## 31. Public-interface protection

There is no `mtblr_bed()` function or NAMESPACE export. No existing public
signature or default changed.

## 32. Generated wrappers

`Rcpp::compileAttributes()` adds exactly the internal
`mtblr_bed_internal()` and `mtblr_bed_marker_contract_internal()`
registrations in `R/RcppExports.R` and `src/RcppExports.cpp`.

## 33. Installed-check behavior

Numerical execution, validation, marker-oracle, full-E identity, reduction,
raw-schema, and reproducibility tests are portable. Only source-text
architecture checks skip narrowly without repository source.

## 34. Mutation sensitivity

The Phase 17O mutation audit detects all 31 required ownership, decoding,
conditional, update-order, output, reuse, scope, and protection mutations.

## 35. Tests and CI

The fast workflow filter includes `17o`. Focused and exact fast-tier tests
passed. The final source suite reported 4805 expectations: 4803 passed, zero
failed, zero warnings, zero errors, and two established opt-in skips. The
built-package test process completed successfully.

## 36. Package check

The built-tarball package check completed with zero errors and zero warnings.
It reported three established notes in two legacy categories: the long Phase
17C fixture path, plus installed-size and existing `std::cout` diagnostics.
No Phase 17O note was introduced.

## 37. Diff hygiene

The final audit checks generated wrappers, NAMESPACE, existing fixtures,
line endings, compiled artifacts, tarballs, check directories, and temporary
files.

## 38. Deviations and blockers

No statistical or architectural blocker remains. Timing components inside the
single native binding are not separately instrumented; the benchmark labels
those component columns unavailable instead of estimating them.

## 39. Recommended next phase

> expose a public `mtblr_bed()` adapter that reuses the validated Phase 17O internal route and existing Glist/BED alignment, while defining phenotype centering evidence, covariate policy, initialization defaults, full-versus-diagonal residual-covariance controls, trait/marker metadata, diagnostics, memory warnings, and public examples without changing the Phase 17O numerical implementation.

## 40. Readiness marker

PHASE 17O COMPLETE — INTERNAL INDIVIDUAL-LEVEL MT PACKED-BED CORE ACTIVE


R CMD check completed with 0 errors, 0 warnings, and 3 notes.

The long Phase 17C fixture-path note and compiled std::cout note are unchanged.
The installed-package-size note is new in Phase 17O and is classified as a
non-functional packaging-size note. It does not indicate a numerical,
interface, schema, test, or portability failure.
