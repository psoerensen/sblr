# Unified BLR Framework Phase 17A Report

## 1. Executive summary

Every remaining active non-packed-BED backend was inventoried, classified, and
prioritized without production migration. The seven scalar CSR families retain
canonical architecture. The three block-eigen routes remain internal research
implementations. The public multivariate `sblr()` family remains supported but
uses a legacy architecture; its default `mtblr()` route is the highest-value
migration candidate. The OpenMP CPG variant has a worker-sensitive RNG risk and
must be corrected before migration. No numerical source, route, export, wrapper,
schema, or `NAMESPACE` entry changed.

## 2. Repository baseline

| Item | Result |
|---|---|
| Branch | `master` |
| Starting commit | `8454320` (`Clarify experimental packed-BED BayesC routes`) |
| Phase 16A commit | `8454320` |
| Initial status | clean |
| Compiler/toolchain | Rtools 4.4, GCC/MinGW-w64 with OpenMP |
| Compile/load | passed |
| Baseline full suite | 5,651 passed, 0 failed, 0 warnings, 8 skipped |

## 3. Repository-wide backend inventory

| Family/route | Representation | Binding | Status | Tests/references |
|---|---|---|---|---|
| ordinary CSR BayesC | sparse CSR LD | public | canonical | strong permanent references |
| scheduled CSR BayesC | sparse CSR LD | public | canonical | strong permanent references |
| CSR BayesR | sparse CSR LD | public | canonical | strong permanent references |
| CSR SBayesRC | sparse CSR LD plus annotations | public | canonical | strong permanent references |
| fixed-prior CSR BayesC | sparse CSR LD plus fixed priors | public | canonical | strong permanent references |
| group CSR BayesC | sparse CSR LD plus groups | public | canonical | strong permanent references |
| learned-annotation CSR BayesC | sparse CSR LD plus annotations | public | canonical | strong permanent references |
| block-eigen BayesC | CSR plus fit-local BED-derived dense blocks | internal | experimental | smoke test only |
| block-eigen BayesR | CSR plus fit-local BED-derived dense blocks | internal | experimental | smoke test only |
| block-eigen SBayesRC | CSR plus fit-local BED-derived dense blocks | internal | experimental | smoke test only |
| `mtblr()` | dense individual data, multivariate | public through `sblr()` | supported legacy | smoke test only |
| `mtblr_cpg()` | dense individual data, multivariate | public through `sblr()` | experimental legacy | smoke test only |
| `mtblr_cpg_arma()` | dense individual data, multivariate | public through `sblr()` | experimental legacy | smoke test only |
| `mtblr_cpg_omp()` | dense individual data, multivariate | public through `sblr()` | experimental with RNG risk | smoke test only |
| `mtblr_eigen()` | eigen representation, multivariate | public through `sblr()` | experimental legacy | smoke test only |
| `mtblr_hybrid()` | dense hybrid multivariate | native-only | unused candidate | no deterministic reference |
| `mtblr_cpg_omp_csr()` | CSR individual-level multivariate | native-only research | experimental | no deterministic reference |

Headers and implementation headers serving those routes are active production
or shared helpers; preparation utilities such as BED checks, sparse-LD readers,
and sampler distributions are not independent backend routes.

## 4. R route inventory

| R route | Reachability | Native destination | Formatting |
|---|---|---|---|
| `stblr_csr()` | supported public | ordinary, scheduled, or BayesR CSR entry | `stblr_raw_v1` formatter |
| `stblr_csr_annot()` | supported public | CSR SBayesRC entry | `stblr_raw_v1` formatter |
| `stblr_csr_prior_annot()` | supported public | fixed-prior CSR entry | `stblr_raw_v1` formatter |
| `stblr_csr_group_annot()` | supported public | group CSR entry | `stblr_raw_v1` formatter |
| `stblr_csr_learn_annot()` | supported public | learned-annotation CSR entry | `stblr_raw_v1` formatter |
| `.stblr_csr_*_block_eigen()` | internal research | three block-eigen entries | standard scalar formatter path |
| `sblr(algorithm=...)` | supported public legacy | five multivariate entries | positional result named/formatted in R |

No separate active non-BED dense scalar wrapper was found. Packed-BED routes
remain protected controls and are outside this audit.

## 5. Native export inventory

All 17 inventoried native symbols have generated R wrappers. The seven scalar
CSR exports have active public callers. The three block-eigen exports have
internal callers. Five multivariate exports are selected by public `sblr()`.
`mtblr_hybrid` and `mtblr_cpg_omp_csr` have no active package R caller and are
therefore native-only; the latter retains research-script evidence.

## 6. Usage evidence

Scalar CSR routes have public R routing, extensive test coverage, developer
documentation, fixtures, and canonical benchmarks. Block-eigen routes appear in
internal R helpers and `test-stblr-block-eigen.R`, but not public exports or
normal workflows. The multivariate family is exposed by `sblr()` and appears in
tests and examples; coverage is smoke-oriented. Draft/research material refers
to `mtblr_cpg_omp_csr()`. No current package workflow calls `mtblr_hybrid()`.
Historical phase reports were not counted as current usage.

## 7. Capability verification

The capability metadata was confirmed for the canonical CSR families and
corrected at the status level for block-eigen and legacy multivariate routes.
CSR routes are scalar and differ by mixture, scheduler, prior, group, and
annotation contracts. Block-eigen is a representation/operator specialization,
not a distinct public model family. Multivariate routes model correlated traits
and therefore remain scientifically distinct from scalar CSR and packed-BED.

## 8. Architecture maturity

| Route group | Context/core boundary | Aggregate/converter | Ownership | Classification |
|---|---|---|---|---|
| seven scalar CSR routes | typed, binding-neutral | singular typed paths | borrowed immutable inputs | canonical architecture |
| three block-eigen routes | reuse canonical CSR numerical cores through operator | scalar CSR presentation | fit-local blocks | audited but noncanonical |
| public/default multivariate | monolithic binding/numerics | positional R formatting | incompletely explicit | legacy architecture |
| multivariate variants | monolithic progression variants | positional or native-only | incompletely explicit | legacy architecture |

## 9. RNG ownership and reproducibility

| Backend | Engine/identity | Risk | Evidence |
|---|---|---|---|
| canonical scalar CSR | logical-chain engines and seed mapping | logical-chain safe | strong permanent references and fresh-process tests |
| block-eigen | inherits canonical CSR core | logical-chain safe within internal route | smoke tests; frozen references absent |
| `mtblr()`, `mtblr_cpg()`, `mtblr_cpg_arma()`, `mtblr_eigen()` | one native engine from an R-generated seed | single-call/single-engine; explicit replay unavailable publicly | smoke test only |
| `mtblr_cpg_omp()` | iteration plus `omp_get_thread_num()` seed contribution | worker-sensitive risk | source-level reproduction of worker-dependent formula |
| `mtblr_hybrid()`, `mtblr_cpg_omp_csr()` | native single-engine paths | unknown across processes/fits | no deterministic reference |

No production defect was changed. The `mtblr_cpg_omp()` issue can change draws
with worker assignment and is a correctness prerequisite for any migration.

## 10. Threading and dispatch

Canonical CSR dispatch has deterministic logical tasks, explicit result order,
and tested worker independence. Block-eigen uses the same scalar execution
contracts after fit-local operator construction. Multivariate OpenMP variants
parallelize lower-level work rather than typed trait-chain tasks; notably CPG
OpenMP binds stochastic identity to the worker. Exception and aggregation
boundaries are not uniformly typed in the legacy multivariate implementations.

## 11. Data ownership and memory

| Family | Ownership | Scaling assessment |
|---|---|---|
| scalar CSR | borrowed immutable sparse inputs, chain-local state | explicit and safe |
| block-eigen | fit-local dense eigen blocks plus CSR/BED preparation | potentially copy-heavy; valuable only when block structure pays off |
| dense multivariate | dense genotype/phenotype and covariance workspaces | potential scaling risk; ownership not expressed as typed views |
| CSR multivariate | sparse representation with legacy workspaces | potentially valuable, but ownership and copy frequency require measurement |

Completed-fit and peak-RSS evidence exists for canonical workstreams, not for
the active multivariate candidates. Peak RSS must be established before their
migration.

## 12. I/O behavior

Scalar CSR data are prepared before MCMC and numerical cores do not perform
hidden file reads. Block-eigen helpers read/prepare BED information and build
dense blocks before numerical execution; the internal numerical core does not
reread BED per MCMC step. Multivariate dense routes consume in-memory objects.
No evidence of MCMC-time file reopening was found, but page-cache and temporary
dense-block costs are not benchmarked for block-eigen.

## 13. Schema maturity

| Family | Schema status |
|---|---|
| canonical scalar CSR | `stblr_raw_v1` compliant with stable formatted fits |
| internal block-eigen | uses scalar raw/format conventions but lacks a public support contract |
| public multivariate `sblr()` | legacy but stable positional native output named and reshaped in R |
| native-only multivariate | backend-specific/undocumented |

The multivariate route does not yet follow the intended named raw-backend
boundary. No field, order, class, dimname, or actual-`NULL` behavior changed.

## 14. Deterministic-reference coverage

| Backend group | Raw/formatted fixtures | Fresh/core/worker coverage | Strength |
|---|---|---|---|
| seven canonical CSR families | permanent phase-specific fixture families | same/fresh/core-order as applicable | strong permanent references |
| block-eigen | none frozen | deterministic tiny smoke only | smoke test only |
| five public multivariate routes | none frozen | ordinary smoke coverage | smoke test only |
| two native-only multivariate routes | none | none | no deterministic reference |

Reference capture, replayable seed policy, recursive first-difference reporting,
and raw/formatted schema decisions are minimum prerequisites for multivariate
migration or retirement decisions.

## 15. Scientific uniqueness and duplication

Scalar CSR variants are scientifically distinct through mixture, scheduling,
fixed-prior, group, or learned-annotation behavior. Block-eigen is
representation-specific and performance-specialized; it may reproduce scalar
models but can retain computational research value. Multivariate modeling is
scientifically unique. The multiple dense multivariate sources appear to be
progression variants rather than separately supported scientific contracts.
`mtblr_hybrid()` is an apparently superseded duplicate; the CSR multivariate
variant may have representation-specific value.

## 16. Block-eigen family

The three internal routes accept prepared scalar inputs plus block/eigen
information, construct fit-local dense blocks, and invoke canonical CSR BayesC,
BayesR, or SBayesRC numerical policies. They are internal research routes with
validation and smoke tests, no public namespace exposure, no permanent fixture,
and no benchmark/memory baseline. They remain **retain experimental (P2)**.
Before migration, establish frozen references, route contracts, dense-block
ownership/RSS, and a workload showing computational value over plain CSR.

## 17. Multivariate family

`mtblr()` is the authoritative default selected by public `sblr()` and is the
P1 migration candidate. `mtblr_cpg()`, `mtblr_cpg_arma()`, and `mtblr_eigen()`
are retained experimental P3 progression/specialization variants.
`mtblr_cpg_omp()` is retained but classified P0: correct its worker-sensitive
RNG before migration. `mtblr_hybrid()` is a P4 retire/remove candidate because
it is native-only with no caller or unique tested contract.
`mtblr_cpg_omp_csr()` remains native-only research P2 because sparse
multivariate execution may be valuable, but it needs references and ownership
evidence. The public wrapper's hidden R-generated seed prevents user-controlled
exact replay and should be part of the next boundary audit, not changed here.

## 18. Other scalar/non-packed-BED routes

No additional active scalar model route was discovered. BED readers, sparse-LD
construction, genotype cross-products, distributions, and chain utilities are
shared preparation or numerical helpers. Commented/historical fragments are
not active routes and remain historical inactive.

## 19. CI coverage

Canonical CSR families compile in normal package CI and have targeted framework,
schema, and reference tests. Block-eigen has tiny package tests but no frozen
reference or benchmark job. Public multivariate routes compile and receive smoke
coverage, but lack deterministic fixtures, fresh-process matrices, and memory
baselines. Native-only multivariate routes lack CI execution evidence. Phase 17A
adds only fast source/status audit coverage.

## 20. Complete disposition matrix

| Backend/route | Reachability | Value/maturity | RNG/reference | Disposition | Priority |
|---|---|---|---|---|---|
| ordinary/scheduled CSR BayesC | supported public | canonical | safe/strong | canonical — no action | P2 |
| CSR BayesR/SBayesRC | supported public | canonical | safe/strong | canonical — no action | P2 |
| fixed/group/learned-annotation CSR BayesC | supported public | canonical distinct policies | safe/strong | canonical — no action | P2 |
| three block-eigen routes | internal research | performance-specialized | safe/smoke | retain experimental | P2 |
| `mtblr()` | supported public legacy | scientifically unique, authoritative | weak references | retain supported — migrate | P1 |
| `mtblr_cpg()`, `mtblr_cpg_arma()`, `mtblr_eigen()` | supported public legacy algorithms | progression variants | smoke only | retain experimental | P3 |
| `mtblr_cpg_omp()` | supported public legacy algorithm | threaded progression variant | worker-sensitive | correct defect before migration | P0 |
| `mtblr_hybrid()` | native-only | no current use/contract | no references | retire/remove candidate | P4 |
| `mtblr_cpg_omp_csr()` | native-only research | representation-specific value | no references | retain experimental | P2 |

## 21. Migration-readiness matrix

| Candidate | Route/value | Fixtures/RNG | Seams/benchmark | Readiness |
|---|---|---|---|---|
| public default `mtblr()` | clear, actively public, unique multivariate value | reference work and replay policy required | monolithic seams; no baseline | reference work required first |
| `mtblr_cpg_omp()` | public variant | worker-sensitive RNG defect | legacy seams | correctness correction required first |
| block-eigen family | internal, performance-specialized | fixture work required | operator seam exists; RSS absent | reference work required first |
| `mtblr_cpg_omp_csr()` | internal research, sparse value | no references | unclear ownership/aggregation | needs evidence |
| `mtblr_hybrid()` | native-only, no demonstrated value | none | no support contract | not worth migration |

## 22. Correctness and maintenance risks

| Source | Evidence | Consequence | Recommended action |
|---|---|---|---|
| `src/mt_cpg_omp.cpp` | seed includes `omp_get_thread_num()` | RNG/worker assignment can alter posterior draws | bounded RNG correction before any migration |
| `R/interface_mtblr.R` | random hidden seed and positional result naming | exact user replay and raw-schema validation are absent | establish references and route/schema contract |
| block-eigen helpers | dense fit-local blocks without RSS baseline | possible memory regression/page-cache sensitivity | measure ownership, peak RSS, and representative workload |
| native-only MT exports | wrappers without package callers | maintenance and stale-code risk | prove research use or retire in a later bounded phase |

## 23. Protected canonical backends

Canonical packed-BED and scalar CSR numerical sources were hash-protected and
remained byte-identical. Permanent packed-BED references were 9/9 raw and 9/9
formatted exact. Permanent scalar CSR reference families were 30/30 raw and
30/30 formatted exact: ordinary BayesC 6/6, scheduled BayesC 3/3, BayesR 6/6,
SBayesRC 6/6, and fixed-prior/group/learned-annotation BayesC 3/3 each.
Block-eigen and multivariate sources also remained byte-identical during this
audit; no expected fixture was regenerated.

## 24. Public API and schema

Public R functions, arguments, defaults, route selection, native signatures,
generated wrappers, `NAMESPACE`, `stblr_raw_v1`, formatted fields, field order,
dimnames, classes, and actual R `NULL` behavior are unchanged.

## 25. Tests

Phase 17A adds route-inventory, reachability, capability/status, maturity, RNG,
reference-strength, block-eigen, multivariate, public-boundary, and protected
source tests. Phase 17A focused validation passed 47 expectations. The enabled
canonical Phase 11B/13A/14A/14E fresh-process selection passed 211 expectations
with no failure, warning, or skip. The full suite passed 5,698 expectations with
zero failures, zero test warnings, and eight intentional opt-in skips. Native
compile/load and both inventory scripts passed. R emitted only the non-failing
notice that `testthat` was built under R 4.4.3 while validation used R 4.4.1.

## 26. Recommended next phase

Establish deterministic raw/formatted references and a migration-boundary audit
for the authoritative public `sblr(algorithm = "default")` / `mtblr()` backend
while classifying duplicate multivariate implementations for retention or
retirement.

## 27. Readiness marker

PHASE 17A COMPLETE — REMAINING BACKEND ROUTES AUDITED AND PRIORITIZED
