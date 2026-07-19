# Unified BLR Framework Phase 17B Report

## 1. Executive summary

The authoritative public multivariate route, `sblr(algorithm = "default")` to
native `mtblr()`, is fully audited and frozen without production migration.
Three frozen native positional and three formatted references now cover
updated covariance/probability controls, nominally fixed controls, correlated
traits, explicit marker sets, and explicit state patterns. The route is
single-chain, single-engine, in-memory, and independent of OpenMP worker identity.

The audit found a public correctness defect: `updateB = FALSE does not keep B
fixed`, because two marker-covariance updates execute unconditionally inside the
set loop. It also confirmed legacy posterior denominators and strict burn-in
comparisons that require an explicit compatibility decision. Separately, the
native sparse execution branch of public experimental `mtblr_cpg_omp()` remains
worker-sensitive. These findings are frozen as evidence, not corrected here.

## 2. Repository baseline

| Item | Result |
|---|---|
| Branch | `master` |
| Starting/Phase 17A commit | `dc429e1` (`Audit remaining backend routes and priorities`) |
| Initial status | clean |
| Toolchain | R 4.4.1, Rtools44 GCC/MinGW-w64, C++17/OpenMP |
| Baseline compile/load | passed |
| Baseline full suite | 5,698 passed, 0 failed, 0 warnings, 8 skipped |

## 3. Public route and call graph

```text
exported sblr()
  -> validate/prepare yy, Xy, XX, n, models, priors, covariances and sets
  -> algorithm == "default"
  -> generate seed with sample.int(.Machine$integer.max, 1)
  -> .Call("_sblr_mtblr", ...)
  -> generated native wrapper mtblr()
  -> src/mtblr.cpp:4050 authoritative execution
  -> one Gibbs loop and inline summaries
  -> 20-position native std::vector result
  -> R/interface_mtblr.R assigns names/orientations
  -> 18 BayesC fields plus conditional rb/rg/re correlations
  -> ordinary named R list
```

There is exactly one implementation behind `algorithm = "default",` with no
dimension-, environment-, or error-based fallback. Alternative implementations
are selected only by explicit public algorithm values. `sblr()` is exported;
the native generated wrapper is internal. There is no separate formatter or
`stblr_raw_v1` validation on this legacy path.

## 4. Multivariate implementation inventory

| Implementation | Source | Reachability | Representation/threading | Role |
|---|---|---|---|---|
| `mtblr()` | `src/mtblr.cpp` | authoritative public default | dense/list summary statistics, serial | authoritative public |
| `mtblr_cpg()` | `src/mt_cpg.cpp` | experimental public | dense CPG, serial | algorithmic progression variant |
| `mtblr_cpg_arma()` | `src/mt_cpg_arma.cpp` | experimental public | Armadillo dense, serial | representation/progression variant |
| `mtblr_cpg_omp()` | `src/mt_cpg_omp.cpp` | experimental public | Armadillo plus OpenMP | threaded variant; correction required |
| `mtblr_eigen()` | `src/mtblr.cpp` | experimental public | eigen input convention, serial | representation variant |
| `mtblr_cpg_omp_csr()` | `src/mt_cpg_omp_csr.cpp` | native-only research | file-backed CSR preparation | internal representation research |
| `mtblr_hybrid()` | `src/mtblr.cpp` | native-only | dense hybrid/prototype | retirement candidate |

No additional active multivariate numerical route was discovered.

## 5. Statistical model

The public route accepts trait-specific summary statistics: `wy[t][i]` is the
marker/trait cross-product, `ww[t][i]` the diagonal marker cross-product, and
`XXvalues[t][i][j]` the dense contribution used to update summary residuals.
There are `nt` traits and `m` markers. Effective effects `b[t][i]`, latent
effects `beta[t][i]`, binary inclusion states `d[t][i]`, and residual summaries
`r[t][i]` are trait-major.

Each marker selects one joint trait pattern from `models[k][t]`. The default
patterns are the `2^nt` binary combinations in R `expand.grid` order. Conditional
on a pattern, a latent `nt`-vector is drawn from a multivariate normal with
precision based on marker covariance precision `Bi`, residual precision `Ei`,
and trait-specific `ww`. Effective effects are the latent draw masked by the
selected pattern. A marker has one joint state and pattern-specific trait
indicators.

`B` is the latent marker-effect covariance. The active route updates it using
both set-specific scaled-inverse-chi-square/correlation reconstruction and an
inverse-Wishart latent-effect update, then optionally performs another global
scaled-inverse-chi-square/correlation update. `E` is residual covariance, but
the called overload updates only its diagonal; initialized off-diagonal entries
remain unchanged. `G` is derived genetic covariance. `pi` contains global joint
pattern probabilities.

Trait-specific sample sizes and summary inputs are accepted. Cross-trait sample
overlap is not explicitly represented. Missing values have no explicit policy.
Effects default to zero, heritability to 0.5, and default state mass to `1-pi`
on the all-zero pattern with remaining mass split equally.

## 6. Trait, marker, state, and covariance ordering

| Dimension | Input/internal/output contract |
|---|---|
| Traits | list order of `Xy`/`XX`; zero-based native `t`; marker-matrix columns, trace columns, and covariance rows/columns preserve that order |
| Markers | input vector/matrix order; zero-based native `i`; public rows restore input names |
| Update order | descending `max_t (wy/ww)^2`, restricted by zero-based sets; `o` stores the zero-based permutation |
| States | default two-trait order `00`, `10`, `01`, `11`; native indices zero-based; `d` is binary |
| Covariances | trait-by-trait Armadillo column-major internally; R reconstruction retains trait order |

`pi` and `pim` names are underscore-joined state patterns. No triangle-packed
covariance representation or assignment base conversion is used.

## 7. Priors and update order

Defaults are `vy = diag(yy/(n-1))`, `vg = diag(vy)*h2`,
`ve = diag(vy)*(1-h2)`, and `vb = diag(vy)*h2/(m*pi)`. Marker and residual
scales are `((nub-2)/nub)*diag(vg)/(m*pi)` and
`((nue-2)/nue)*diag(ve)`; both degrees of freedom default to four.

One iteration executes:

1. reset model counts to one;
2. per set, call `sampleBset()` and `sampleB_latent()`;
3. invert `B`, traverse all scored markers, and call
   `sampleBetaCPG_Mt_latent()` for members of the set;
4. accumulate retained marker summaries under `it > nburn && it % nthin == 0`;
5. optionally sample/accumulate `pi`;
6. optionally sample/accumulate global `B`;
7. compute/accumulate `G`;
8. optionally sample the diagonal of `E`, invert it, and accumulate it;
9. increment the retained count.

Traces have `nit + nburn` rows. Covariance and probability accumulators divide
by `nit`, although strict `it > nburn` yields `nit-1` unthinned contributions;
marker summaries divide by `nsamples`. These legacy semantics are frozen.

The control is defective: `updateB = FALSE does not keep B fixed`, because
step 2 is unconditional. The fixed-control fixture shows final `B` changes while
`E` and `pi` remain fixed.

## 8. RNG ownership

| Engine/distribution | Owner/lifetime | Classification |
|---|---|---|
| `std::mt19937 gen(seed)` | `mtblr()` stack, one fit | single-chain safe |
| uniform state draw | marker helper call | fit-local via borrowed engine |
| normal latent/Wishart draws | helper call | fit-local |
| gamma probability draws | probability update | fit-local |
| chi-square covariance draws | covariance helper call | fit-local |

R draws the native seed from its global RNG. Replay requires resetting R's seed;
there is no public explicit seed argument. Default execution has no worker
contribution, `static` RNG, `thread_local` RNG, or cached distribution state.

## 9. `mtblr_cpg_omp()` worker-sensitive RNG risk

The sparse native branch was run with identical input/seed in isolated one- and
two-thread processes. Its first difference was:

```text
path: root$1$1
index: 1
one thread: 0.0614042
two threads: 0.05991059
```

The source uses `seed + 100000 * it + omp_get_thread_num()`. Static scheduling
therefore changes marker-to-engine assignment. The public dense wrapper's full
indices currently trigger serial fallback, but the exported native symbol can
reach the risky sparse branch. Disposition: **correct RNG before retention** or
retire/disable it; it is not a migration oracle.

## 10. Threading and dispatch

`mtblr()`, `mtblr_cpg()`, `mtblr_cpg_arma()`, and `mtblr_eigen()` are serial.
The default has no chains, task map, worker ID, reduction, or core control;
changing `OMP_NUM_THREADS` leaves frozen default outputs exact.
`mtblr_cpg_omp()` parallelizes markers within graph-color sets using static
scheduling, worker-local engines, and a critical count reduction.
`mtblr_cpg_omp_csr()` is separately threaded/representation-specific. Failures
cross the Rcpp boundary; there is no typed failure result.

## 11. Data representation

Default R conversion turns each dense `XX` matrix into data-frame/list columns
and creates full zero-based indices for every marker. Nested vectors store
summary statistics and marker state; Armadillo handles covariances. The Armadillo
variant moves more state to dense matrices. Eigen changes `ww` and scales
residual quantities. CSR consumes a file prefix. Hybrid is a dense prototype.

## 12. Ownership and memory

| Object | Dimensions | Owner/copy risk |
|---|---|---|
| R `XX` and native `XXvalues` | `nt*m*m` | immutable R input plus full native conversion; dominant O(nt m²) duplication |
| `XXindices` | `m*m` | native copied full indices |
| `wy`, `ww`, `b`, `beta`, `r`, `d` | `nt*m` each | fit-owned nested vectors, several mutable copies |
| `B`, `E`, `G`, inverses | `nt*nt` | fit-local Armadillo matrices |
| posterior marker arrays | `nt*m` | fit-local accumulators |
| traces | `nt*(nit+nburn)` | fit-local accumulators |
| helper workspaces | trait/model sized | repeated marker-update allocation |
| positional result | 20 slots | duplicates summaries during conversion |

There are no chain or thread copies on default execution. Dense LD conversion
and helper allocation dominate scaling risk.

## 13. I/O behavior

Default execution is entirely in-memory, opens no file, and performs no
MCMC-time I/O. R-to-C++ conversion duplicates dense data. The native-only CSR
variant prepares files per fit and needs a separate audit. Page-cache behavior
does not apply to the default route.

## 14. Legacy raw/formatted schema

| Position | Field | Meaning/shape |
|---:|---|---|
| 1--7 | `bm`, `dm`, `wy`, `r`, `b`, `d`, `o` | marker by trait matrices |
| 8--10 | `vbs`, `vgs`, `ves` | iteration by trait traces |
| 11--13 | `covb`, `covg`, `cove` | posterior covariance summaries |
| 14--16 | `vb`, `vg`, `ve` | final covariance matrices |
| 17--18 | `pi`, `pim` | final and legacy mean pattern probabilities |
| 19--20 | `pitrait`, `pimarker` | BayesR-only positions, discarded for BayesC |

BayesC formatting retains 1--18, names/orients them, and conditionally appends
`rb`, `rg`, `re` only for positive covariance diagonals. Absent correlation
fields are omitted, not actual `NULL`. There is no class, timing, failure object,
or schema metadata. This is a legacy positional schema, not `stblr_raw_v1`.

## 15. Frozen references

Three configurations in `blr_phase17b_mt_default` cover: all updates; correlated
initial covariance with nominally fixed controls; and explicit marker sets plus
four explicit joint patterns. Each has two traits, four named markers, explicit
metadata and seed. **Raw references: 3/3 within `1e-12`. Formatted references:
3/3 within `1e-12`.** Types, dimensions, names, state/order fields, and schema
are exact. Armadillo covariance/solve results vary by a few ulps across repeated
executions, so bitwise comparison is not justified. Generation is manual
maintenance tooling only.

## 16. Reproducibility

| Comparison | Result |
|---|---|
| A; A | numerically equivalent within `1e-12` with exact structure |
| A; B; A | numerically equivalent within `1e-12` with exact structure |
| fresh versus reused process | numerically equivalent within `1e-12` |
| `OMP_NUM_THREADS=1;2` default | within `1e-12`; public core control unsupported |
| one/multiple chains | unsupported; default is single-chain |
| intervening multivariate fit | no effect with R seed reset |
| intervening canonical fits | canonical controls remain exact; resetting R seed isolates default fit |

The lack of a public seed argument leaves replay dependent on R's global RNG.

## 17. Scientific identities

Frozen outputs are finite, marker-by-trait dimensions/names are stable, final
states are binary, and `dm` is bounded. `B`, `G`, `E` and summaries are symmetric
and positive semidefinite within tolerance. Covariance dimnames retain trait
order and trace length is `nit+nburn`. The fixed fixture proves `E` and `pi`
stay fixed while `B` changes, preserving defect evidence.

## 18. Alternative implementation comparisons

| Variant | Classification |
|---|---|
| `mtblr_cpg()` | algorithmic progression; related model, different update/RNG path |
| `mtblr_cpg_arma()` | Armadillo progression/representation variant |
| `mtblr_cpg_omp()` | threaded variant with unsafe sparse worker RNG |
| `mtblr_eigen()` | eigen representation with distinct scaling |
| `mtblr_cpg_omp_csr()` | native-only sparse representation research |
| `mtblr_hybrid()` | historical prototype without caller/reference |

No alternative is an exact numerical duplicate suitable as an oracle.

## 19. Variant dispositions

| Implementation | Final disposition |
|---|---|
| `mtblr()` | retain as authoritative public frozen legacy |
| `mtblr_cpg()` | retain experimental |
| `mtblr_cpg_arma()` | retain experimental |
| `mtblr_cpg_omp()` | correct RNG before retention |
| `mtblr_eigen()` | retain experimental |
| `mtblr_cpg_omp_csr()` | retain internal representation research |
| `mtblr_hybrid()` | retire/remove candidate |

## 20. Typed audit contracts

`src/blr_mt_default_audit_types.h` defines audit-only `MtDataSpec`,
`MtModelSpec`, `MtCovariancePriorSpec`, `MtExecutionAuditSpec`,
`MtOwnershipAuditSpec`, `MtChainResultVocabulary`, and
`MtExecutionResultVocabulary`. They use standard C++ only, invoke no sampler,
consume no RNG, contain no Rcpp/SEXP or file handle, and are not production
dependencies.

## 21. Future extraction seam

The migration unit should be complete single-fit numerical execution, not a
chain/thread abstraction. In `src/mtblr.cpp`: native entry/decoding starts at
4050; initialization ends at 4176; engine and MCMC begin at 4177--4181; marker
updates are 4185--4236; probability/covariance updates are 4251--4295; execution
ends at 4300; summary/result construction is 4302--4392; return is 4393. Binding
naming/orientation begins at `R/interface_mtblr.R:252`. Extraction must preserve
set/marker order, strict burn-in, helper order, engine, draws, and positional
result. It must wait for the `updateB` correction decision.

## 22. Future aggregation and converter boundary

Marker/covariance/probability accumulators currently live inside the sampler;
division and 20-slot copying follow the loop. R interleaves orientation, naming,
truncation, and correlation calculation. A future typed aggregate should own all
numeric summaries/finals/traces and retained counts. One numerical aggregator
should finalize denominators. One named converter should own names, orientation,
legacy/public field order, correlations, classes, and schema metadata. No such
boundary is activated here.

## 23. Performance, memory, and I/O baseline

| Workload | Five times (s) | Mean/median | Completed-fit RSS |
|---|---|---|---|
| 2 traits, 4 markers, updated | 0.09, 0.02, 0.00, 0.02, 0.00 | 0.026/0.02 | 132.45 MiB |
| 2 traits, 4 markers, nominal fixed | 0.00, 0.01, 0.02, 0.01, 0.00 | 0.008/0.01 | 132.46 MiB |
| 3 traits, 80 markers, updated | 0.03, 0.03, 0.04, 0.03, 0.01 | 0.028/0.03 | 121.52 MiB |

Tiny and moderate `XX` objects are approximately 0.0025 and 0.1776 MiB.
Completed-fit RSS is not peak RSS; peak RSS was not sampled. Dense matrices may
dominate memory. Tiny timings are regression signals, and no cross-implementation
speed ranking is made. Default execution has no I/O/page-cache dependency.

## 24. Public API and schema

Arguments, defaults, algorithm values, routes, signatures, positional output,
public names/order/orientation, conditional correlation fields, and dimnames are
unchanged. The legacy schema is documented and frozen, not redesigned.

## 25. Protected backends

All multivariate variants, block-eigen files, canonical scalar CSR and packed-BED
sources, generated wrappers, and `NAMESPACE` remain byte-identical to `dc429e1`.
Canonical reference families remain protected controls.

## 26. Tests

Phase 17B covers route singularity, 3/3 raw and 3/3 formatted references,
same/fresh-process behavior, thread-environment behavior, scientific identities,
the update-control defect, variant dispositions, P0 risk evidence, audit types,
and protected hashes. Ordinary focused validation passed 126 expectations with
one intentional fresh-process skip; the enabled fresh-process run passed
127/127 with no skip. The full suite passed 5,824 expectations with zero
failures, zero test warnings, and nine intentional opt-in skips. Native
compile/load, the expected-risk audit, and the benchmark passed. R emitted only
the non-failing notice that `testthat` was built under R 4.4.3.

## 27. Risks and blockers

| Source | Consequence | Recommended action |
|---|---|---|
| unconditional `sampleBset/sampleB_latent` | `updateB=FALSE` changes `B` | bounded public correction before migration |
| strict burn-in plus `/nit` summaries | systematic legacy scaling ambiguity | decide/correct with explicit reference transition |
| hidden R seed | replay depends on global RNG | explicit resolved-seed contract in future boundary |
| CPG OpenMP worker seed | thread-dependent posterior | correct before retention or retire |
| dense `XX` conversion | O(nt m²) duplication | peak-RSS/ownership protection |
| positional conditional schema | absent rather than `NULL` fields | deliberate future schema/converter contract |

## 28. Recommended next phase

Correct the authoritative public `mtblr()` marker-covariance update-control
defect so `updateB = FALSE` is honored, decide and test the strict-burn-in
posterior-denominator policy, and establish post-correction deterministic
references before any mechanical extraction.

## 29. Readiness marker

PHASE 17B COMPLETE — PUBLIC MULTIVARIATE CONTRACT AND MIGRATION AUDIT READY

## Post-review hardening

The Phase 17B route assertion now uses the Phase 12A zero-safe source-count
helper and directly distinguishes zero, one, and duplicate matches. The fast
GitHub Actions gate includes ordinary Phase 17B tests, while extended validation
sets `SBLR_RUN_PHASE17B_FRESH="true"` for the existing full suite. Hosted CI
execution is not claimed by this local review.

Current capability metadata now classifies default `mtblr()` as the
authoritative supported public legacy implementation, explicitly noncanonical,
with its `updateB=FALSE` defect and denominator-policy decision still blocking
extraction. Alternative implementations retain their separate dispositions.

Fixture metadata now states `structure_exact_numeric_tolerance`, numerical
tolerance `1e-12`, and `structure_exact=TRUE`. Names, types, dimensions,
dimnames, classes, order, and positional schema remain exact; only numerical
comparisons use the narrow tolerance. The three existing RDS files received
metadata-only rewrites, with their `raw` and `fit` objects proven identical
before and after rewriting.

A direct permanent runtime test freezes the current legacy probability summary:
strict zero-based `it > nburn` accumulation contributes `nit - 1` probability
vectors and final `pim` divides by `nit`, so `sum(pim)=(nit-1)/nit`. Source-level
assertions retain both statements. This is evidence for Phase 17C, not the
desired future policy. No numerical production behavior, route, API, wrapper,
or schema changed, and Phase 17C has not started.

Local post-review validation passed 147 ordinary Phase 17B expectations with
one intended fresh-process skip, 148/148 with fresh-process validation enabled,
48 Phase 12A workflow-structure expectations with its existing peak-RSS skip,
and the full suite at 5,845 passes, zero failures, zero test warnings, and nine
intentional opt-in skips. Native compile/load and `git diff --check` passed.
These are local results; hosted GitHub Actions execution is not claimed.
