# Unified BLR Framework Phase 16A Report

## 1. Executive summary

Every noncanonical packed-BED BayesC route now has an explicit disposition. The
scheduled single-chain route is retained unchanged as an explicitly experimental
deterministic/reference implementation. The sparse route is retained unchanged
as an explicitly experimental performance-policy implementation because its
`null_update_prob` scheduler is genuinely distinct. Neither is a fallback from
canonical `stblr_bed(method = "bayesc")`; no experimental numerical migration was
started. The commented duplicate implementations in the single-chain source are
classified as historical audit artifacts rather than active routes.

## 2. Repository baseline

Branch `master` began clean at commit `de0f5a7` (`Consolidate packed-BED common
infrastructure`), which is also the Phase 15B commit. `git diff --check` passed.
R 4.4.1 and Windows Rtools44/GCC were used. `compileAttributes()`, a fresh native
build, and the pre-edit full suite succeeded; the Phase 15B checkpoint contained
5,621 passing expectations and eight opt-in skips. Canonical references began at
3/3 raw and 3/3 formatted for each of BayesC, BayesR, and BayesRC.

## 3. Complete route inventory

| Route | R wrapper | Native export/source | Chains/scheduler | Aggregation/converter | Status |
|---|---|---|---|---|---|
| scheduled multichain | `stblr_bed(method="bayesc")`; lower-level multichain branch | `stblr_cpg_omp_bed_marker_scheduled_chains()` / `st_cpg_omp_individual_scheduled_chains.cpp` | logical trait-chain tasks; adaptive scheduled | typed native aggregate and named converter | canonical public |
| scheduled single-chain | `stblr_bed_marker(backend="auto"|"scheduled", nchains=1)` for one BED file | `stblr_cpg_omp_bed_marker_scheduled()` / `st_cpg_omp_individual_scheduled.cpp` | trait-parallel historical adaptive scheduler | route-local summaries/conversion | experimental public lower-level |
| sparse marker scheduler | `stblr_bed_marker(backend="sparse")` | `stblr_cpg_omp_bed_marker_sparse()` / `st_cpg_omp_individual.cpp` | trait-parallel stochastic null updates/candidates | route-local summaries/conversion | experimental public lower-level |
| commented scheduled duplicates | none | commented blocks in `st_cpg_omp_individual_scheduled.cpp` | none | none | historical inactive audit artifact |

There is one active implementation for each named export and no old/new selector
inside any route.

## 4. Reachability

The canonical call graph is `stblr_bed()` -> BayesC preparation ->
`stblr_cpg_omp_bed_marker_scheduled_chains()` -> typed context/core/results ->
aggregation -> converter. It never calls either experimental export.

The lower-level exported `stblr_bed_marker()` selects the historical single-chain
export only for a one-file, one-chain `auto`/`scheduled` request; multichromosome
or multichain requests select the canonical multichain export. Its explicit
`backend="sparse"` branch selects the sparse export and rejects more than one BED
file. Both native exports also remain internally callable through generated Rcpp
wrappers. No other package dispatcher reaches them.

## 5. Capability comparison

| Feature | Canonical multichain | Experimental single-chain | Experimental sparse |
|---|---|---|---|
| chains | one or more | exactly one | exactly one |
| BED files | one or more | one in current dispatcher | one |
| scheduler | deterministic null skip/candidate policy | historical deterministic null skip/candidate policy | stochastic `null_update_prob` plus candidate policy |
| RNG identity | shared logical task/seed contract | historical trait seed; chain-owned engine | trait-owned engine |
| typed core/result/aggregate | yes | no | no |
| chain summaries | canonical degenerate/multichain | no canonical chain vocabulary | no canonical chain vocabulary |
| public schema | `stblr_raw_v1` canonical | formatted through raw compatibility path | formatted through raw compatibility path |
| unique support value | supported public implementation | deterministic historical oracle | research scheduler/performance policy |

The single-chain route is a functional subset, but its preserved RNG mapping makes
it useful as a historical oracle. The sparse route does not define a different
BayesC posterior; it is execution-policy distinct and may exhibit different
finite-run behavior and costs.

## 6. Deterministic comparisons

The permanent Phase 11B single-chain fixture remains exact. A matched canonical
one-chain call is intentionally not bit-identical: the established first sampled
difference is marker posterior effect `raw$marker$bm`, caused by historical seed
mapping rather than a model difference. The sparse route is same-process exact
under repeated matched calls. It is not required to equal the canonical route,
because stochastic null visitation changes RNG consumption and finite-run state.

## 7. Usage evidence

`stblr_bed_marker()` is exported, documented, tested for row handling, and used in
the research examples and workflow notes. Phase 11B fixtures directly protect the
single-chain route. Sparse calls occur in experimental/example workflows but not
in the canonical `stblr_bed()` dispatcher or ordinary canonical fixtures. The two
native symbols occur in generated wrappers and historical architecture reports.

## 8. Final disposition matrix

| Route | Before | Unique capability | Reachability/usage | Coverage | Final disposition | Action |
|---|---|---|---|---|---|---|
| scheduled single-chain | noncanonical | historical seed/execution oracle | exported lower-level; fixtures/workflows | exact Phase 11B reference | retain as explicitly experimental | clarify support and add permanent Phase 16A checks |
| sparse marker scheduler | noncanonical | stochastic null-update scheduler | exported lower-level; research examples | repeated-call deterministic check | retain as explicitly experimental | clarify support and add permanent Phase 16A checks |
| commented historical duplicates | inactive | none at runtime | unreachable | source inventory | historical audit artifact | retain classified; no compiled symbol |

## 9. Removed routes

None. Removal was not justified: it would discard either a permanent deterministic
oracle or a distinct research scheduling policy and would unnecessarily disrupt
an intentionally exported lower-level research interface. Therefore no C++ file,
native symbol, wrapper, dispatch branch, fixture, or benchmark was removed.

## 10. Retained experimental routes

Both active routes are unsupported for general use and are not canonical. The
single-chain route exists for deterministic regression/research comparisons. The
sparse route exists only to study sparse visitation policy. They retain their
existing inputs, outputs, RNG ownership, genotype decoding, and limitations; no
canonical schema or reproducibility guarantee is inferred beyond their permanent
tests.

## 11. Deferred migrations

No migration is recommended. If research establishes that the sparse scheduler
provides material value, a separate deterministic-reference and typed-boundary
audit would be required. The historical single-chain route should remain an
oracle, not be migrated into a second canonical implementation.

## 12. Public API

`stblr_bed(method=c("bayesc","bayesr","bayesrc"))` remains the supported public
surface. `stblr_bed_marker()` remains an exported lower-level research interface;
its `auto`/`scheduled` one-file/one-chain path and `sparse` option are explicitly
experimental. Arguments, defaults, method values, and return schema are unchanged.

## 13. Native exports

No export was removed. `Rcpp::compileAttributes()` produced no intentional wrapper
change. Both experimental exports and the canonical multichain export remain
synchronized in `R/RcppExports.R` and `src/RcppExports.cpp`; `NAMESPACE` is
unchanged.

## 14. Canonical BayesC protection

The canonical per-chain core, task/seed helpers, aggregation, converter, and
adapter numerical behavior are byte-identical to the starting commit. Canonical
BayesC retains 3/3 raw and 3/3 formatted exact references plus same-process,
fresh-process, core-order, chain-count, and worker-independent reproducibility.

## 15. BayesR/BayesRC protection

BayesR and BayesRC production numerical sources remain byte-identical. Each
retains 3/3 raw and 3/3 formatted exact references. The fixed-alpha
BayesRC-to-fixed-pi-BayesR reduction remains exact.

## 16. Documentation and workflow cleanup

The implementation plan, capability matrix, reduction matrix, computation
inventory, backend naming note, and `stblr_bed_marker()` parameter documentation
now distinguish supported canonical routing from experimental lower-level routing.
Existing workflows remain intentionally experimental and no canonical workflow
was redirected.

## 17. Performance, memory, and I/O

`tools/benchmarks/blr_phase16a_experimental_bayesc_disposition.R` records five
tiny repetitions per retained route, completed-fit RSS, BED size, and disposition
rationale. These are capability/disposition signals, not canonical speed rankings.
Completed-fit RSS is not peak RSS and page-cache effects apply. Both routes decode
fit-local packed storage before sampling; no new decoding, copying, or MCMC-time
I/O was introduced. On the 6-sample/2-marker, 7-byte BED fixture, all five timed
repetitions rounded to 0.00 seconds for each route; completed-fit RSS was 32.6 MiB
for scheduled single-chain and 32.7 MiB for sparse. Peak RSS was not sampled.

## 18. Tests

Phase 16A adds route inventory, reachability, distinction, RNG/ownership,
single-chain exact-reference, sparse repeated-call, public-support, and protected
source tests. Final validation recorded:

- Phase 16A focused: 56 expectations across eight tests, all passed
- canonical references: 9/9 raw and 9/9 formatted exact
- enabled canonical fresh-process selection: 157 expectations at the canonical
  checkpoint; the Phase 11B, 13A, 14A, and 14E opt-in processes were rerun and
  passed in Phase 16A
- full suite: 5,651 passes, zero failures, zero warnings, eight opt-in skips

## 19. Risks and limitations

The sparse route lacks a frozen raw/formatted fixture and remains intentionally
noncanonical; its protection is deterministic repetition plus static policy
checks. The single-chain file contains large commented historical blocks, retained
as classified audit history to avoid unrelated source churn. Experimental routes
do not receive the typed ownership/aggregation guarantees of the canonical route.
`git diff --check` passes. The tracked content and ignore-space-at-EOL statistics
are identical; Git-for-Windows reports only its normal future LF-to-CRLF checkout
notices, with no line-ending-only churn.

## 20. Recommended next phase

Audit the remaining active non-packed-BED backend families and prioritize only
those with real public or research use, beginning with a repository-wide route
and usage inventory rather than immediate migration.

## 21. Readiness marker

PHASE 16A COMPLETE — EXPERIMENTAL PACKED-BED BAYESC ROUTES DISPOSED
