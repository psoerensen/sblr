# BLR Phase 2 provider/operator checkpoint

## Status and scope

**Status:** `READY FOR INDEPENDENT PHASE 2 VERIFICATION`

Phase 2 implements shared data ownership, marker alignment, operator-resource,
and likelihood-provider adapters. It does not change a sampler, posterior,
native signature, RNG stream, seed, schedule, retention rule, or output
meaning. The current MT covariance hybrid remains legacy and is not routed to
`blr_raw` version 2 as sampled $V_b$.

## Executable source trace

| Representation | Existing owner/view and operation path | Phase 2 decision |
|---|---|---|
| Complete or sparse CSR | `SparseLdCsrStorage` and `SparseLdCsrView` in `src/blr_sparse_ld_csr.h`; `CsrOperator` in `src/st_ld_operator.h`; construction in `src/st_cpg_omp_csr.cpp` and `src/st_cpg_omp_csr_bayesr.cpp` | Reuse unchanged. The view is immutable and the operator already supplies diagonal, residual-score update, rebuild, quadratic, and residual operations. |
| Packed BED | Move-only `PackedBedMatrix` in `src/packed_bed.h`; templated `BedPackedGenotypeView` in `src/blr_bed_family_types.h`; scheduled ST contexts construct borrowed views after one owner is prepared | Reuse unchanged. BED stays sample-space and is not materialized as a dense cross-product. |
| Block eigen | `BlockEigenStorage` and `BlockEigenView` in `src/blr_block_eigen.h`; `BlockEigenDispatchOperator` in `src/st_ld_operator.h`; preparation in `src/st_block_eigen_execution.h` | Reuse unchanged. Dense reconstructed and retained-coordinate residual policies remain distinct. |
| Dense cross-product | No separate maintained public ST hot loop; dense matrices are reference/reconstruction objects in existing reduction tests | Implement an R qualification adapter only. Do not introduce a duplicate native algebra stack. |

`PackedBedMatrix` deletes its copy constructor. `SparseLdCsrView`,
`BlockEigenView`, and `BedPackedGenotypeView` borrow immutable storage.
Residual scores, marker effects, fitted values, RNGs, workspaces, accumulators,
and diagnostics remain execution-task state and are not fields of an operator
resource.

## Implemented R contract

`R/blr-provider-operators.R` supplies small internal constructors and validators
for:

- one ordered global marker universe with aligned alleles, coding, and stable
  one-based indices;
- immutable operator-resource references and serializable legacy view
  descriptors;
- likelihood providers with ordered nonempty trait sets, resource IDs,
  unique resource-local-to-global maps, sufficient statistics, sample sizes,
  likelihood/error regimes, scales, populations, and provenance;
- provider collections that own resources once and validate every provider
  reference and alignment;
- explicit provider presence matrices, where an absent marker has no
  likelihood contribution;
- dense, complete/thresholded CSR, full-rank block-eigen, and retained-rank
  block-eigen reference actions used only for deterministic qualification.

`validate_blr_resolved_spec()` now invokes the same marker/resource/provider
validators. Maintained ST wrappers still forward identical native arguments in
identical order. Their legacy native operator construction is unchanged; the
resolved specification records the shared Phase 2 contract around it.

Several providers may reference one resource. The collection owns the one
storage reference, while providers contain only `operator_resource_id`. For
independent providers estimating one effect, incompatible declared effect
scales fail before numerical evaluation. Allele, coding, order, score-shape,
sample-size, and map mismatches also fail before evaluation.

## Representation contracts

Complete aligned dense and complete CSR resources represent the same declared
cross-product. Thresholded CSR is marked as an approximation.

For block $b$,

$$
C_b \approx Q_b\Lambda_bQ_b^\top.
$$

A full-rank block-eigen resource is exact for its declared block-diagonal
operator. It equals an original dense cross-product only if omitted cross-block
entries are zero or the likelihood explicitly excludes them. A retained-rank
resource equals its reconstructed retained operator, not the original dense
operator.

The packed-BED resource describes selected marker/sample ownership, coding,
centering, standardization, and provenance. Existing native decoding remains
authoritative and unchanged.

## Qualification boundary

| Capability | Status after Phase 2 |
|---|---|
| Shared R marker/resource/provider construction and validation | Implemented and active in maintained Phase 1 ST resolved specifications |
| Native CSR, BED, and block-eigen immutable views | Implemented before Phase 2, source-traced, adopted as the common boundary, and unchanged |
| Dense/CSR/block-eigen R operator reductions | Qualification-only deterministic reference |
| Multiple independent providers and heterogeneous marker subsets | Structurally represented and independently evaluable; not yet accumulated in one production posterior |
| Common-sample multi-trait BED provider | Structurally represented as one coupled provider; not connected to an MT sampler |
| Correct Cheng MT transition or sampled full $V_e$ | Deferred |
| Unified scheduler, seed contract, and retention | Deferred to Phase 3 |

## Permanent evidence

`tests/testthat/test-blr-provider-operators.R` verifies:

- identity, subset, permutation, missing-marker, duplicate-map, range, allele,
  coding, and score-order contracts;
- one shared immutable storage reference used by multiple providers;
- dense versus complete CSR equality;
- full-rank block eigen versus its declared block-diagonal operator;
- retained-rank self-consistency and deliberate difference from the original
  dense operator;
- block-order, provider-order, and provider-local marker-order invariance;
- compatible versus incompatible shared-effect scales;
- one non-factorized two-trait common-sample BED provider;
- preservation of legacy seed- and retention-contract version 0 in maintained
  resolved specifications.

Existing operator, BED decoding, block-eigen, public-interface, and trajectory
tests remain the executable oracle for native behavior. No Phase 2 test changes
a frozen scientific expectation.

## Qualification results

Checkpoint measurements on 2026-08-14 were:

- 41 Phase 2 provider/operator expectations passed;
- the frozen Phase 0 contract fixture passed 60 expectations, and the Phase 1
  contract suite passed unchanged;
- the independent research provider/operator reference passed, including
  heterogeneous-provider and retained-rank checks;
- native compilation, package loading, source-package building, generated-
  interface comparison, and private `Rcpp::compileAttributes()` comparison
  passed; generated wrappers remained unchanged;
- the generated-interface audit passed with 74 wrappers, 74 registrations, and
  31 exports, and the architecture audit passed 19 of 19 checks;
- the complete ordinary source suite had no Phase 2 failure. Nine frozen MD5
  expectations in `test-bayesc-engine-extraction.R` and
  `test-logvar-block-eigen.R` differed under the current compiler build and
  reproduced identically in a clean-HEAD copy using the same DLL. They were not
  updated;
- the documentation audit retained its known unrelated historical-link false
  positive in `sbayesrc_s_em_phase5c.md`.

These are checkpoint measurements, not permanent test-count invariants.

## Deferred decisions

Phase 3 owns unified logical-task scheduling and seed derivation. Later summary
likelihood phases own accumulation of heterogeneous providers and overlap-aware
error models. The first corrected joint-MT vertical slice owns the Cheng state,
$V_b$, and staged $V_e$ policies. Phase 2 deliberately does not make any of
those decisions executable.
