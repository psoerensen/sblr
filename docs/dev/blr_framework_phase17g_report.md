# Unified BLR Framework Phase 17G report

## 1. Executive summary

The corrected typed default `mtblr()` route is now the sole supported public
multivariate sampler. Three redundant or defective CPG implementations and the
unused `mtblr_hybrid()` prototype were removed. `mtblr_eigen()` and
`mtblr_cpg_omp_csr()` remain only as explicitly unsupported, native/internal
research evidence. No CSR or block-eigen sampler was implemented.

## 2. Repository baseline

- Branch: `master`.
- Starting and Phase 17F2 commit: `bb418f0a929a39e07c3dd2fa2876a54d03b69a76`.
- Initial status: clean; `git diff --check` passed.
- Hosted CI: no combined hosted result was visible locally, so none is claimed.
- Toolchain: R 4.4.1 UCRT; GCC 13.2.0; Rtools44.
- Phase 17F2 baseline: 52 files, 363 blocks, 4,286 passes, zero failures and
  warnings, two opt-in skips, and 31.86 seconds of recorded test time. The
  pre-edit validation rerun also passed; its wall time was 70.5 seconds.

## 3. Previous multivariate inventory

The previous surface comprised `mtblr()`, `mtblr_cpg()`, `mtblr_cpg_arma()`,
`mtblr_cpg_omp()`, `mtblr_eigen()`, `mtblr_hybrid()`, and
`mtblr_cpg_omp_csr()`.

## 4. Dependency audit

`mt_cpg.cpp`, `mt_cpg_arma.cpp`, and `mt_cpg_omp.cpp` each owned their exported
sampler and local marker-update helpers. Their only R callers were the retired
`sblr()` dispatch branches; their other references were generated wrappers,
historical documents, tests, and exploratory tooling.

`cpg_samplers.cpp/.h` formerly held seven helpers. Full-repository call search
showed that only `samplePi_cpg()` has a retained caller:
`mtblr_cpg_omp_csr()`. The other six helpers were used only by the deleted CPG
routes and were removed. Default `mtblr()` uses its own helpers in `mtblr.cpp`;
the retained eigen implementation shares those default-file helpers. Scalar
and packed-BED implementations do not call the removed helpers.

## 5. Public routing before and after

Before:

```text
sblr(algorithm)
  -> mtblr | mtblr_cpg | mtblr_cpg_arma | mtblr_cpg_omp | mtblr_eigen
```

After:

```text
sblr(algorithm = "default")
  -> mtblr
  -> typed core
  -> typed finalizer
  -> legacy positional adapter
  -> unchanged R formatting
```

Non-default values fail before summary preparation; there is no fallback or
silent remapping.

## 6. Retired CPG routes

- `mtblr_cpg()`: removed as a redundant historical progression with defective
  update/retention semantics and no permanent references.
- `mtblr_cpg_arma()`: removed as a redundant optimization branch with the same
  unsupported legacy contract.
- `mtblr_cpg_omp()`: removed because its statistical contract was redundant and
  worker identity entered seed construction.

The three production files and their generated wrappers/registrations are gone.
Git history and historical reports retain the evidence.

## 7. Hybrid retirement

`mtblr_hybrid()` had no R caller, fixture owner, workflow, or unique supported
capability. Its complete function region and generated wrapper were removed.
The adjacent corrected default and retained eigen functions were not rewritten.

## 8. Retained internal research routes

| Route | Reachability | Status | Restriction |
|---|---|---|---|
| `mtblr_eigen()` | generated native/internal wrapper only | unsupported research | temporary evidence until a shared scalar/MT block-eigen representation exists |
| `mtblr_cpg_omp_csr()` | generated native/internal wrapper only | unsupported research | local one-LD `LDCSR`; not canonical and not a public-interface basis |

Neither route is selected by `sblr()` or exported through `NAMESPACE`.

## 9. Rcpp export cleanup

`Rcpp::compileAttributes()` removed wrappers and registrations for
`mtblr_cpg`, `mtblr_cpg_arma`, `mtblr_cpg_omp`, and `mtblr_hybrid`. It retained
exactly the default, eigen-research, and CSR-research MT wrappers. MT generated
wrappers fell from seven to three. `NAMESPACE` did not change.

## 10. Public interface

The `algorithm` argument remains with default `"default"` for signature
stability. Any other value raises: “Only algorithm = \"default\" is supported;
experimental multivariate backends were retired.” The method arguments,
hidden R-derived seed, native default signature, 20 positions, names,
dimensions, orientations, classes, and conditional correlations are unchanged.

## 11. Corrected scientific contract

Phase 17C remains authoritative: update controls, `it >= nburn`, relative
thinning, accumulator-specific denominators, normalized updated `pim`, fixed
disabled state, and one fit-local RNG all pass unchanged. Phase 17E typed-core
and Phase 17F typed-finalization tests pass without duplicated references.

## 12. Shared CSR implications

The retained local `LDCSR` is explicitly experimental and noncanonical. Future
MT CSR must reuse the exact canonical scalar CSR storage, ownership, marker
ordering, allele orientation, and validation vocabulary.

## 13. Trait-specific LD requirement

The future bundle requires one LD operator per trait/study. Values and sparsity
patterns may differ. Row offsets and indices may be shared only when truly
identical; independent structures must not be forced into a union pattern.

## 14. Shared block-eigen implications

Future MT block-eigen execution must reuse the scalar representation and allow
per-trait/study blocks, eigenvectors, eigenvalues, ranks, tolerances, and
reference metadata. The retained eigen route was not migrated or modified.

## 15. Tests

- Focused Phase 17A/17C/17E/17F route, reference, core, and finalizer tests:
  passed after wrapper regeneration.
- Permanent fast filter: passed in 56.18 seconds; one expected peak-RSS skip.
- Complete ordinary suite: 52 files, 364 blocks, 4,294 passing expectations,
  zero failures, zero warnings, and two opt-in skips.
- Package compilation/load: passed.
- `devtools::check()` could not start `R CMD build` because Windows denied
  `processx` write-pipe creation (`system error 5`). A direct `R CMD check`
  reached DESCRIPTION validation and stopped because this research snapshot
  lacks explicit `Author` and `Maintainer` fields. No compile, test, or Phase
  17G source diagnostic failed.

## 16. Reference preservation

Phase 17C remains the sole current owner. All three corrected native/raw and all
three formatted fixtures passed with exact structure and numerical tolerance
`1e-12`. Phase 17B and Phase 17C RDS files have no Git diff and were not
regenerated.

## 17. Quantitative simplification

| Measure | Before | After | Change |
|---|---:|---:|---:|
| compiled `.cpp` sources | 26 | 23 | -3 |
| active MT native entries | 7 | 3 | -4 |
| public MT algorithm branches | 5 | 1 | -4 |
| generated MT wrappers | 7 | 3 | -4 |
| deleted parallel CPG source lines | 2,123 | 0 | -2,123 |
| `mtblr_hybrid()` source lines | 423 | 0 | -423 |
| ordinary test files | 52 | 52 | 0 |
| ordinary test blocks | 363 | 364 | +1 |
| ordinary expectations | 4,286 | 4,294 | +8 (0.19%) |

The total diff removes more than 3,300 lines after helper, wrapper, example,
and hybrid cleanup. Current test-recorded time was 113.61 seconds (114.78
seconds wall time); this is not directly comparable with the Phase 17F2
31.86-second recorded baseline because the current Windows debug build and host
load were materially slower. There is no numerical performance claim.

## 18. Production protection

`blr_mt_default_core_impl.h` and `blr_mt_default_finalize_impl.h` are
byte-unchanged. Corrected behavior is protected primarily by deterministic
references. Scalar CSR, packed-BED, and block-eigen numerical sources are
unchanged. Active-source wrapper and research-file hashes were removed from old
phase tests in accordance with the Phase 17F2 hash policy; canonical numerical
references remain active.

## 19. Deviations and blockers

No implementation blocker remains. Package checking is limited by the Windows
`processx` access-denied condition and, when invoked directly, the repository's
pre-existing missing DESCRIPTION `Author`/`Maintainer` fields. Compilation,
package loading, focused tests, fast tests, and the complete suite all
succeeded. Hosted CI status was not available.

## 20. Recommended next phase

Audit and formalize the exact canonical scalar CSR storage, ownership,
marker-order, allele-orientation, and validation contracts for reuse by a
trait-specific multivariate LD-operator bundle, without yet implementing the
multivariate CSR sampler.

## 21. Readiness marker

PHASE 17G COMPLETE — PUBLIC MULTIVARIATE ROUTES CONSOLIDATED
