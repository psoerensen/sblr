# Unified BLR Framework Phase 17H report

## 1. Executive summary

One binding-neutral sparse-LD CSR storage and borrowed-view contract is active.
All scalar CSR backends share the owner definition, and ordinary canonical CSR
BayesC actively uses the borrowed view without changing statistical work.

## 2. Repository baseline

Branch `master` started clean at `481dfb65dab9a5e03a028b88c7959ea7101f44b7`
(Phase 17G). Hosted CI status was not visible locally. R 4.4.1/Rtools44 compiled
and loaded the package. Baseline ordinary tests passed: 52 files, 364 blocks,
4,294 expectations, zero failures/warnings, two opt-in skips; test-recorded
runtime was approximately 31.9 seconds.

## 3. Previous CSR architecture

`STLDCSR` owned buffers, `CsrOperator` copied storage and diagonal,
`CsrBayesCDataView` repeated neutral pointer fields, and `CsrResourceSpec`
described one scalar read-only resource shared across chains.

## 4. Complete backend inventory

| Backend | Owner/view after 17H | Status |
|---|---|---|
| ordinary CSR BayesC | shared owner + active shared view | canonical |
| scheduled CSR BayesC | shared owner alias | canonical, view migration later |
| CSR BayesR | shared owner + `CsrOperator` | canonical |
| CSR SBayesRC | shared owner + `CsrOperator` | canonical |
| fixed-prior CSR BayesC | shared owner alias | canonical |
| group CSR BayesC | shared owner alias | canonical |
| learned-annotation CSR BayesC | shared owner alias | canonical |
| block-eigen scalar | separate eigen representation | unchanged |
| dense corrected MT | no CSR | authoritative public MT |
| `mtblr_cpg_omp_csr` | local `LDCSR` | unsupported research |
| `mtblr_eigen` | legacy eigen | unsupported research |

## 5. Disk contract

Files, types, zero-based indexing, metadata, raw-correlation values, and marker
row order are specified in `blr_sparse_ld_csr_contract.md`.

## 6. Owning storage

`sblr::core::SparseLdCsrStorage` owns `uint64` row pointers, `int` columns,
`float` pre-scaled off-diagonal values, and builder diagnostics. `STLDCSR` is
the sole compatibility alias, not another definition.

## 7. Borrowed view

`SparseLdCsrView` borrows immutable row pointers, columns, values, and the
separate double-precision diagonal for one marker domain. It has no R, trait,
sample, chain, RNG, or result fields.

## 8. Scalar BayesC integration

Disk builder -> `SparseLdCsrStorage` -> existing `CsrOperator` owner ->
`CsrOperator::view()` -> `CsrBayesCDataView::ld` -> centralized validation ->
unchanged ordinary BayesC update loop -> canonical raw converter.

## 9. Other scalar backends

All use the shared owner through one alias. Their model-specific numerical
sources were not migrated to the borrowed view; this avoids broad churn.

## 10. Validation

Mandatory view checks cover pointers, dimensions, monotonic offsets, terminal
nnz, columns, off-diagonal policy, finiteness, and positive diagonal. Symmetry,
duplicate, and sort checks remain builder/deep-validation concerns.

## 11. Ownership and lifetime

The adapter owner outlives all chain tasks; the view is read-only; chains own
only effects, residuals, state, RNG, traces, and workspaces. Moved-from storage
is never referenced.

## 12. Copy audit

Phase 17H adds zero O(nnz) copies and zero per-chain copies. The existing
`CsrOperator(const STLDCSR&, const arma::rowvec&)` still makes one O(nnz) and
one diagonal compatibility copy. Optimizing it is deferred.

## 13. Marker order

CSR disk row order and R summary order are conventionally aligned. Explicit
ID matching exists for selection-S MAF alignment; synthetic native IDs and
ordinary paths leave gaps requiring Phase 17I validation.

## 14. Allele orientation and scaling

Genotypes are standardized `(dosage-2p)/sqrt(2p(1-p))`; `wy=x'y`, `ww=x'x`,
and internal `xij=x_i'x_j`. BED-derived orientation is by construction.
External effect/other alleles, strand, panel, ancestry, and study metadata are
not comprehensively represented or validated.

## 15. Trait-specific LD compatibility

The same view supports fully shared operators, shared row/index patterns with
trait-specific values/diagonals, and fully independent operators. No MT bundle
or sampler was introduced.

## 16. Resource metadata

Current `CsrResourceSpec` remains correct for one scalar resource shared across
chains. Phase 17I needs a small external per-trait resource list with explicit
sharing mode; `ResolvedSpec` was not generalized prematurely.

## 17. Experimental MT CSR comparison

The local research type uses the same widths and pre-scaled float convention,
but has a separate builder, weaker metadata/validation, and one-shared-LD
assumption. It remains noncanonical, unsupported, and not publicly routed.

## 18. Public interfaces

No R argument, default, native signature, raw schema, formatted schema,
wrapper, or namespace behavior changed.

## 19. Numerical protection

Phase 3 ordinary CSR BayesC passed its permanent exact references after the
qualification-only integration. Phase 6 BayesR, Phase 8 SBayesRC, Phases
9C/9E/9G, Phase 10D scheduled BayesC, Phase 17C corrected dense MT, Phase 17E
typed core, and Phase 17F finalization owners all passed. No fixture was
regenerated.

## 20. Performance and memory

The audit script reports exact component byte formulas. Views are fixed-size;
per-chain CSR bytes are zero; no MCMC-time I/O occurs. Calculated storage is
not peak RSS.

## 21. Tests and CI

Phase 17H adds one focused architectural file (52 expectations) to the fast
tier. Focused Phase 17H, all permanent numerical owners, and fast CI passed.
The complete ordinary suite passed (53 files, 369 blocks, 4,346 expectations,
zero failures/warnings, two opt-in skips; 40.1 seconds wall time in the summary
run). Extended process isolation could not start because `processx` was denied
permission to create a Windows pipe.
The attempted test-local dynamic validator probe encountered the known Windows
MSYS signal-pipe limitation; permanent tests instead cover every mandatory
rejection branch structurally, while package compilation verifies the code.

## 22. Deviations and blockers

Owning storage retains the established `ptr`/`idx`/`xij` field names to avoid
broad numerical churn. The existing adapter copy remains. Direct test-time
dynamic compilation is unavailable in this Windows process environment; no
test-only Rcpp export or wrapper churn was introduced. `devtools::check()`
likewise stopped before package build when `processx` could not create the
`Rcmd.exe` write pipe (system error 5); compilation/loading and ordinary tests
passed independently.

## 23. Recommended next phase

Implement an internal trait-specific MT BayesC CSR core using the corrected
typed dense MT statistical contract and one canonical `SparseLdCsrView` per
trait/study, with the dense route retained as the numerical oracle and without
yet changing the public `sblr()` interface.

## 24. Readiness marker

PHASE 17H COMPLETE — SHARED SCALAR/MT CSR CONTRACT ESTABLISHED
