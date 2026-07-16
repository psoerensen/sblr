# Unified BLR Framework Phase 9F3 Report

## 1. Executive summary

Learned-annotation CSR BayesC migration is complete with one typed numerical core, one active marker loop, one named binding-layer result converter, and one wrapper-level multichain aggregation path. Statistical behavior and public objects are unchanged.

## 2. Repository baseline

- Branch and starting commit: `master`, `ab5c78d` (`Activate typed learned-annotation BayesC execution boundary`).
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1 UCRT and Rtools44 GCC 13.2.0.
- Baseline: 4,294 passes, 0 failures, warnings, or skips.
- Phase 9A comparable annotation benchmark (2,000 markers, 100 retained): 1 chain/1 core mean 0.483 s, 2 chains/1 core 0.850 s, retained 2 chains/2 cores 0.887 s; completed-fit RSS 128.6--131.8 MiB. This is completed-fit process RSS, not peak memory.

## 3. Phase 9F2 structure inventory

| Item | Classification | Resolution |
|---|---|---|
| typed context, active policy, callable core, typed result | retain permanently | unchanged |
| `cpg_annot_raw_v1()` inline binding converter | centralize/rename now | renamed to the stable typed-result converter |
| marker/trace/trait/annotation conversion helpers | retain permanently | sole converter uses them; they contain binding-only shape conversion |
| chain payload conversion helper | retain permanently | preserves optional chain schema and order |
| raw vector aliases used by aggregation | defer with justification | wrapper aggregation remains deliberately native and schema-preserving |
| duplicate converter fragments | remove now | no second top-level converter existed; both return branches now call the same converter |
| Phase 9F1/9F2 inline-converter assertions | rename for stable use | now protect the named converter |
| obsolete learned-source MD5 expectations | remove now | already removed in Phase 9F2; structural and frozen-reference tests are permanent |
| context/result fields | retain permanently | all are consumed by validation, execution, conversion, aggregation, or diagnostics |
| `std::cout` diagnostics | retain permanently | binding-neutral and outside worker marker loops |

## 4. Files changed

- `src/st_cpg_omp_csr_annot.cpp`: renamed and typed the sole converter; adapted both post-aggregation return branches.
- Phase 9F1/9F2 tests: replaced the obsolete inline-converter name assertion.
- `tests/testthat/test-blr-framework-phase9f3.R`: permanent architecture, exact-reference, reproducibility, and protection checks.
- `tools/benchmarks/blr_phase9f3_csr_learned_annotation_bayesc.R`: post-migration timing and completed-fit RSS baseline.
- Framework plan and capability matrix: mark the backend migrated and ready for canonicalization.
- This report records migration closure.

## 5. Final execution path

The public learned-annotation R route performs existing validation/alignment, prepares annotation policy and CSR state, constructs `CsrLearnedAnnotationBayesCExecutionContext`, calls `run_csr_learned_annotation_bayesc()`, aggregates chains once at wrapper level, passes a typed result to the named converter, validates `stblr_raw_v1`, and uses the unchanged canonical raw-to-fit formatter.

## 6. Centralized result converter

`stblr_csr_learned_annotation_bayesc_result_to_raw(const CsrLearnedAnnotationBayesCExecutionResult&, ...)` is the sole top-level binding converter. It preserves top-level ordering, R types, dimensions, classes, real `NULL`, marker/trait/annotation/coefficient/chain order, schema metadata, marker summaries, variance/global outputs, annotation summaries, diagnostics, LD-swap, selection, and optional chain payloads. It performs no sampling, proposals, links, clipping, or aggregation.

## 7. Wrapper-level multichain aggregation

Each core invocation produces one native raw payload. The wrapper aggregates marker and coefficient summaries, effective values, LD-swap counts, and optional chain payloads exactly once in input chain order. Conversion happens only after aggregation. `keep_chains` continues to control slots 30--34. Moving this boundary would risk schema, coefficient diagnostic, and chain-order changes, so it remains outside the numerical core.

## 8. Native adapter boundary

The adapter is limited to R decoding, alignment, native preparation, typed-context construction, core invocation, existing wrapper aggregation, one converter call, and exception translation. It contains no MCMC/marker loop, RNG draw, proposal/acceptance calculation, link transformation, clipping, global update, residual/variance update, or posterior formula.

## 9. Numerical core

The guarded implementation header remains unchanged in Phase 9F3. It contains the single callable numerical implementation and the sole active marker loop and remains free of Rcpp and Python binding types.

## 10. Learned-annotation policy preservation

Column-major annotation layout/order, no implicit native intercept, coefficient dimensions/order/initialization, centered-logistic inclusion probabilities, exponential multipliers, priors, random-walk proposals, update frequency, bounds and clipping order, proposal/acceptance counters, global `pi`/`B` interactions, and unsupported cases are unchanged.

## 11. Logging decision

The Phase 9F2 `std::cout` diagnostics are retained. They occur before or after numerical worker regions or in post-region trait reporting, never in the marker hot loop; they consume no RNG, change no control flow or returned state, and preserve diagnostic content without R API calls from workers.

## 12. Exact frozen references

All permanent Phase 9A fixtures pass exactly: learned-annotation raw 3/3 and formatted 3/3.

## 13. Reproducibility

Exact checks cover repeated calls, one/two cores, reversed `1,2,2,1` ordering, intervening group fits, fixture chain seeds and retained/dropped chains, fixed/probability-only/multiplier-only/both learning modes, explicit update frequency and bounds, and disabled LD-swap. Only documented execution metadata is normalized.

## 14. Public API and schema

Arguments, defaults, native exported signatures, routing, `NAMESPACE`, generated wrappers, `stblr_raw_v1`, formatted fit, dimensions, names/classes, and actual `NULL` fields are unchanged.

## 15. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior BayesC, group BayesC, block-eigen, BED, scheduled, and multivariate source protections remain exact against `ab5c78d`; their focused and full-suite tests remain active.

## 16. Performance and memory

Commands: `Rscript tools/benchmarks/blr_phase9a_annotation_backends.R` and `Rscript tools/benchmarks/blr_phase9f3_csr_learned_annotation_bayesc.R`. The post-migration 2,000-marker both-learning rows were 0.43--0.49 s (mean 0.460) for 1 chain/1 core, 0.89--0.94 s (mean 0.922) for 2 chains/1 core, and 0.87--0.99 s (mean 0.936) for retained 2 chains/2 cores. Probability-only mean was 0.406 s, multiplier-only 0.400 s, and fixed coefficients 0.286 s. Completed-fit RSS was 128.1--137.5 MiB. Five timed repetitions followed one warm-up per primary configuration. Compared with Phase 9A, comparable means differ by roughly -5%, +8%, and +6%, within the observed Windows timing variability; RSS is comparable. No unexplained material regression is present. RSS was sampled after completed fits and is not interval peak memory.

## 17. Test results

- Baseline full suite: 4,294 passes.
- Phase 9F3: 33 passes after correcting one over-broad source-scan assertion; no product defect was found.
- Phase 9F1--9F3 focused suite: 105 passes, 0 failures, warnings, or skips.
- Final fresh-build/full suite: 4,327 passes, 0 failures, warnings, or skips.

## 18. Deviations and blockers

The first sandboxed focused rebuild could not discover Rtools and returned no compiler diagnostic. A fresh process with explicit Rtools PATH compiled successfully. One initial architecture assertion scanned native helper definitions rather than the exported adapter body; its scope was corrected. The benchmark uses completed-fit RSS because interval process sampling is unavailable. No blocker remains.

## 19. Recommended next phase

> canonicalize and stabilize learned-annotation CSR BayesC, remove remaining migration-only wording or aliases, retain permanent exact fixtures, establish the post-migration baseline as canonical, and then begin the scheduled-CSR contract and RNG audit.

## 20. Readiness marker

PHASE 9F3 COMPLETE — LEARNED-ANNOTATION BAYESC MIGRATED WITH BEHAVIOR PRESERVED
