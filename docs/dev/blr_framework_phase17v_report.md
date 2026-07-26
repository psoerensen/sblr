# Phase 17V implementation report

## 1. Executive summary

Public Tier 1 MT BED convergence diagnostics are active through `mtblr_bed()`.

## 2. Repository baseline

Work began from clean `master` at `e31799abfdfe08327e4e47a4d28bbfa5810037d4`.
The host uses R 4.4.1 UCRT, GCC/GFortran 13.2, OpenMP, and R BLAS/LAPACK.
Baseline source and installed tests passed; check had zero errors/warnings and
the three established notes.

## 3. Phase 17U verification

The trace bundle, extractor, diagnostic mathematics, native route, wrappers,
and registration are unchanged.

## 4. Public signature

`convergence=c("auto","none","core")` and `convergence_control=NULL` follow
`keep_chains`; every prior formal retains its order and default.

## 5. Convergence-mode validation

`match.arg()` accepts exactly `auto`, `none`, and `core`.

## 6. Control validation

The unique named controls are `warn`, three positive finite thresholds, and
`keep_traces`; defaults are TRUE, 1.01, 100, 0.05, and FALSE.

## 7. Route selection

Trace capture is required when traces are retained or diagnostics are requested
with at least two chains. Otherwise the ordinary chains route is selected.

## 8. Exactly-once execution

One function is selected before one `do.call()`; routes are never combined and
chains are never rerun.

## 9. Disabled diagnostics

`none` uses the ordinary route and returns a validated zero-row
`not_requested` object.

## 10. Automatic single-chain behavior

`auto` with one chain uses the ordinary route and returns quiet Tier 1
`unavailable_single_chain` rows.

## 11. Automatic multichain behavior

`auto` with at least two chains captures post-burn Tier 1 traces and computes
the Phase 17U summaries.

## 12. Explicit core behavior

`core` computes when possible; single-chain requests retain validated
unavailability and warn once when enabled. Trace retention still selects the
trace route.

## 13. Raw integration

The additive object is stored at `raw$diagnostics$convergence`; raw version 1
is retained and raw objects without this field remain valid.

## 14. Formatter integration

The sole `.as_mtblr_fit()` conditionally carries the convergence object.

## 15. Public fit output

Every public MT BED fit has `fit$convergence` and a stable present-but-NULL
`fit$convergence_traces` unless traces were requested.

## 16. Convergence warnings

One deterministic main-thread advisory reports status, counts, extrema,
constant-chain mismatches, chain advisory, and `fit$convergence` location.

## 17. Warning suppression

`warn=FALSE`, `none`, and automatic single-chain unavailability emit no
convergence warning; candidate text remains stored when applicable.

## 18. Threshold customization

R-hat, per-chain ESS, and relative-MCSE thresholds are validated and passed to
the unchanged engine.

## 19. Memory integration

Trace capture, one-quantity workspace, summary output, and retained traces are
added to ordinary requested and execution analytical totals.

## 20. Memory-warning text

The pre-execution warning identifies mode, capture, retention, diagnostic GiB,
and reiterates that estimates are not measured RSS or peak RSS.

## 21. Diagnostic scope

Scope is Tier 1 post-burn per-chain B/G/E diagonals only, never pooled traces.

## 22. R-hat semantics

Rank-normalized split and folded R-hat are reported with their maximum.

## 23. ESS semantics

Bulk, q05, q95, tail-minimum, and raw-scale mean ESS follow Phase 17U.

## 24. MCSE semantics

Only posterior-mean MCSE and its posterior-SD ratio are supported.

## 25. Fixed quantities

Disabled B/E updates produce `not_updated`; zero placeholders are not diagnosed.

## 26. Chain and draw sufficiency

One chain is unavailable; R-hat needs four post-burn draws and ESS/MCSE six;
four/five draws can be partial. Fewer than four chains is advisory only.

## 27. Optional traces

Retained arrays are iteration by chain by `3*nt`, exclude burn-in, and receive
no additional thinning. Tier 2 and marker traces remain unsupported.

## 28. keep_chains independence

Compact-chain retention is statistically and structurally independent of
convergence summaries and traces.

## 29. Default numerical reduction

Default single-chain ordinary numerical and sampler fields equal Phase 17S;
only convergence and convergence-memory metadata are additive.

## 30. Public/internal equality

Public computed summaries and traces reduce exactly to the explicit Phase 17U
native-plus-R sequence after timing normalization.

## 31. Warning tests

Deterministic threshold, partial, unavailable, suppression, none, and quiet-auto
cases enforce one aggregated warning at most.

## 32. Reproducibility

Repeated, serial/OpenMP, explicit/default seed, trace/chain-retention, fresh,
and intervening-fit owners remain exact after timing normalization.

## 33. Public API protection

CSR, block-eigen, scalar BED, and `sblr()` formals are unchanged; only
`mtblr_bed()` gains the two requested controls.

## 34. Native protection

All Phase 17U native files retain their starting hashes.

## 35. Wrapper and registration protection

`Rcpp::compileAttributes()` produces no wrapper or registration change.

## 36. Mutation sensitivity

All 50 activation, routing, warning, memory, output, and protection guards pass.

## 37. Benchmark

The 180-case public regression grid completed in 18 seconds. It covered full
and diagonal covariance, `none`/`auto`/`core`, one/two/four chains and requested
cores, both compact-chain settings, and valid trace-retention settings. It
recorded route, dispatch/total time, ordinary/diagnostic memory, object sizes,
and warning status without a speedup claim.

## 38. Existing-route protection

Permanent Phase 15/17, packed-BED, summary-MT, schema, BED, and backend owners
remain active.

## 39. Installed-check behavior

Mode/control, routing, diagnostics, warnings, memory, reductions, and formatting
tests are portable; only source/hash assertions skip without a checkout.

## 40. Generated documentation

Roxygen updates only `man/mtblr_bed.Rd`; NAMESPACE is unchanged.

## 41. Tests and CI

Phase 17V contributes 72 focused expectations. The exact fast filter passed
2,665 expectations (2,664 successes, zero failures/warnings/errors, and the
established peak-RSS opt-in skip). The final full source suite passed 5,523
expectations: 5,521 successes, zero failures, zero warnings, zero errors, and
the two established opt-in skips. Phase 17V is included in the exact fast
workflow filter.

## 42. Package check

Built-tarball installation and installed tests passed. `R CMD check` completed
with zero errors, zero warnings, and the three established notes: the long
Phase 17C fixture path, installed size (5.2 MB; `libs` 4.2 MB), and legacy
scalar-backend `std::cout` references. No new note was introduced.

## 43. Diff hygiene

Final checks compare ordinary and ignore-EOL statistics and protect fixtures,
native sources, wrappers, registration, NAMESPACE, DESCRIPTION, and artifacts.

## 44. Deviations and blockers

No contract deviation or implementation blocker remains.

## 45. Recommended next phase

> audit and formalize extended MT BED convergence diagnostics for lower-triangular B/G/E covariance entries, null versus active model-probability mass, optional pattern-specific probabilities, and explicitly selected marker effects and inclusion states, including trace ownership, memory controls, output scope, warning policy, and computational feasibility, without yet implementing Tier 2 or Tier 3 traces.

## 46. Readiness marker

PHASE 17V COMPLETE — PUBLIC MT BED TIER 1 CONVERGENCE DIAGNOSTICS ACTIVE
