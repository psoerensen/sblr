# Unified BLR Framework Phase 13B Report

## 1. Executive summary

The deterministic packed-BED BayesR per-chain numerical execution was
mechanically extracted into one guarded implementation header without changing
Phase 13A behavior.

## 2. Repository baseline

Branch `master`; starting/Phase 13A commit `582ba76` (`Audit packed-BED BayesR
migration boundary`); initial tree and `git diff --check` clean. R 4.4.1,
Rtools44, GNU C++17/OpenMP. Phase 13A baseline: 4,998 full-suite passes, no
failures/warnings, four opt-in skips; enabled Phase 13A fresh matrix 61/61.

## 3. Production call graph

`stblr_bed(method="bayesr")` -> `.fit_stblr_bed_bayesr()` -> native export ->
BED validation/fit-local decoding -> OpenMP trait-chain dispatch -> extracted
`run_one_bayesr_chain()` -> inline aggregation -> inline raw-schema construction
-> canonical raw-to-fit formatter.

## 4. Extraction target

`run_one_bayesr_chain()` was the sole active per-chain implementation. It owns
chain RNG/sampler/scheduler state and returns `ChainResultBayesR`; dispatch,
genotype ownership, aggregation and conversion were already separable.

## 5. Extraction seam

The moved block began at `struct ChainResultBayesR` (Phase 13A source line 46)
and ended immediately before `// [[Rcpp::export]]` (through line 629, including
the stable separating blank line). The adapter now includes
`blr_bed_bayesr_core_impl.h` at that exact seam.

## 6. Files changed

The BayesR source replaces the contiguous block with one include. The new core
header owns the moved block plus guard/comment. Phase 13A tests were adjusted
only for the intentional location/hash change. Phase 13B adds structural/
reference tests and this report; plan/matrix mark migration in progress.

## 7. Lines mechanically moved

```text
MECHANICAL_LINES=584
HEADER_BODY_LINES=584
IDENTICAL=TRUE
```

Comparison against commit `582ba76` normalized only header guard/comment.
There are no changed numerical lines.

## 8. Extracted execution content

The header contains `ChainResultBayesR`, marker component sampling, local
Dirichlet update, adaptive skip helper, chain seed/RNG construction, effects,
residuals, assignments, pi/variance state, scheduler buckets/lists, the sole
MCMC loop, full/adaptive traversal, packed-genotype access, posterior/component
accumulation, CPO, retained counts, timing/failure and finalization.

## 9. Component semantics preservation

Component 0 remains the scale-zero null; components 1..K-1 remain zero-mean
normals with variance `vb*c[k]`. Input/prior/output order, inverse-CDF order,
conditional normal draw, `dm=1-P(null)`, final/mean pi and component state are
unchanged.

## 10. Scheduler preservation

Initial jitter distribution and draw sites, null due scheduling, active/
candidate initialization, reservations, adaptive skip formula, candidate expiry,
50-iteration compaction, iteration-zero/periodic full sweeps and exact
active-candidate-due order are byte-identical. Skipped-marker RNG/genotype
behavior is unchanged.

## 11. RNG ownership preservation

One chain-local mt19937, uniform, normal and uniform-integer jitter distribution
remain constructed after the resolved trait-chain seed. Variable-shape Gamma
draws remain local. No static/thread-local, worker-owned or fit-persistent state
was introduced.

## 12. Genotype and I/O preservation

BED validation, blocked file decoding, marker maps and fit-owned packed storage
remain adapter-side. The header borrows prepared storage; it has no file open,
fseek/fread, path ownership, additional decode or per-chain genotype copy.

## 13. Existing aggregation and conversion

OpenMP dispatch, task-result storage/failures, cross-chain means/SD/min/max,
component/trace/final/CPO/timing aggregation and all R raw construction remain
inline and mechanically unchanged in the adapter.

## 14. Exact references

Phase 13A raw references: 3/3 exact. Formatted references: 3/3 exact. Phase 11A
BayesR audit reference remains unchanged. No fixture was regenerated.

## 15. Reproducibility

Repeated A, A-B-A, normalized 1-2-2-1, fresh/reused process, intervening BayesC,
intervening BayesRC and different-chain-count sequences remain exact.

## 16. Component identities

Final/mean pi and per-marker component probabilities sum to one; dm remains one
minus null probability; states remain in range; scale identity and fixed-pi
normalization behavior remain unchanged.

## 17. Reductions and nonreductions

The documented BayesR/BayesC nonreduction and full-sweep different-skip-base
nonreduction remain unchanged, including scheduler-initialization jitter as the
latter's cause.

## 18. Protected backends

BayesRC, canonical/experimental BED BayesC, canonical CSR BayesR/other CSR,
block-eigen, multivariate, wrappers, NAMESPACE, APIs and schemas are unchanged.

## 19. Tests

Phase 13A and new Phase 13B structural/reference tests, enabled fresh-process
matrix, BayesR focused protections and full-suite results passed. Phase 13A/13B
focused: 100 passes and one opt-in skip; enabled fresh-process: 104/104. Full
suite: 5,041 passes, zero failures/warnings and four opt-in skips. The Phase 13A
fresh path was run separately and passed.

## 20. Deviations and blockers

None. Conditional existing progress Rcout remains inside the lexical helper by
design; binding neutrality is Phase 13C work and moving it now would violate
mechanical identity.

## 21. Recommended Phase 13C

> replace the lexically dependent packed-BED BayesR chain include with an explicit typed per-chain execution context, callable numerical core, and typed per-chain result while retaining task dispatch, aggregation, and R conversion until all Phase 13A references pass again.

## 22. Readiness marker

PHASE 13B COMPLETE — PACKED-BED BAYESR PER-CHAIN EXECUTION BLOCK EXTRACTED
