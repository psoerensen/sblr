# BLR permanent test-contract ownership

This document assigns one primary executable owner to every supported BLR contract. Supporting tests may exercise a route for a different purpose, but must not repeat its complete reference family or process matrix.

## Test tiers

- **Tier 1 — fast:** canonical reference owners, focused units, public/schema contracts, scientific identities, and essential architecture boundaries. No child processes, thread-process matrices, benchmarks, or peak-memory work.
- **Tier 2 — complete ordinary:** Tier 1 plus inexpensive supported and experimental integration tests. It has no enabled fresh-process or peak-memory work.
- **Tier 3 — extended/manual:** `test-blr-extended-reproducibility.R`, peak RSS, experimental RNG-risk tooling, and manually invoked benchmarks.

## Numerical references

| Contract | Permanent owner | Fixture/support | Tier | Status | Superseded duplicates removed |
|---|---|---|---|---|---|
| ordinary CSR BayesC | `test-blr-framework-phase3.R` | Phase 2 CSR fixtures | 1 | current | Phases 1–2 |
| ordinary CSR BayesR | `test-blr-framework-phase6.R` | Phase 5A fixtures | 1 | current | Phases 5A–5B |
| ordinary CSR SBayesRC | `test-blr-framework-phase8.R` | Phase 7A fixtures | 1 | current | Phases 7A–7B3 |
| fixed-prior CSR BayesC | `test-blr-framework-phase9c.R` | Phase 9 annotation fixtures | 1 | current | Phases 9B1–9B3 |
| group CSR BayesC | `test-blr-framework-phase9e.R` | Phase 9 annotation fixtures | 1 | current | Phases 9D1–9D3 |
| learned-annotation CSR BayesC | `test-blr-framework-phase9g.R` | Phase 9 annotation fixtures | 1 | current | Phases 9F1–9F3 |
| scheduled CSR BayesC | `test-blr-framework-phase10d.R` | Phase 10B fixtures | 1 | current | Phases 10A–10C3 |
| scheduled packed-BED BayesC | `test-blr-framework-phase11d.R` | Phase 11B fixtures | 1 | current | Phases 11A–11C3 |
| packed-BED BayesR | `test-blr-framework-phase13e.R` | Phase 13A fixtures | 1 | current | Phases 13A–13D |
| packed-BED BayesRC | `test-blr-framework-phase14e.R` | Phase 14A fixtures | 1 | current | Phases 14A–14D |
| corrected dense MT BayesC | `test-blr-framework-phase17c.R` | Phase 17C corrected fixtures | 1 | current | Phases 17D–17F repetitions |

## Public and scientific contracts

| Contract | Permanent owner | Tier | Notes/action |
|---|---|---|---|
| R arguments/defaults and routing | CSR/BED interface tests; Phase 17C for the sole MT default route | 1 | retired MT algorithms fail early; behavioral public calls preferred |
| native signatures | interface and generated-wrapper tests | 1 | wrapper hash only where deliberate immutability applies |
| raw schemas and actual NULL semantics | `test-stblr-raw-schema.R` | 1 | sole raw-schema owner |
| formatted schemas, names, dimensions | backend interface tests | 1 | model-specific format owners |
| marker/trait ordering | model reference owner | 1 | included in exact-structure reference comparison |
| fixed update controls | model reference owner; Phase 17C for MT | 1 | behavioral |
| retained-count policy | Phase 4 scalar unit; Phase 17C MT | 1 | behavioral |
| probability normalization | model reference owner | 1 | behavioral |
| covariance symmetry/PSD and orientation | Phase 17C for MT; model-specific owners elsewhere | 1 | behavioral |
| binary states and bounded inclusion means | model reference owner | 1 | behavioral |
| residual consistency | backend consistency tests | 1 | behavioral |

## Reproducibility and architecture

| Contract | Permanent owner | Tier | Notes/action |
|---|---|---|---|
| same-process repeated calls | each canonical model reference owner | 1 | one focused comparison per model |
| intervening fits | model-specific owner only where RNG policy differs | 2 | broad repeated matrices removed |
| fresh process and thread environment | `test-blr-extended-reproducibility.R` | 3 | representative scalar, scheduled, BED, and MT routes |
| logical-chain and seed mapping | Phase 4 and Phase 15B | 1 | direct unit contracts, not process matrices |
| one canonical core/converter/aggregation path | final canonical architecture owner per family | 1 | narrow source assertions |
| binding-neutral MT types/core | Phase 17E | 1 | forbids Rcpp/SEXP and adapter Gibbs/RNG work |
| one MT numerical finalizer | Phase 17F | 1 | owns division/adaptation boundary |
| shared CSR naming requirements | Phase 17F | 1 | documentation contract only |
| shared CSR storage/view and validation | `test-blr-framework-phase17h.R` | 1 | architecture only; Phase 3 owns ordinary BayesC numerics |
| internal trait-specific MT CSR reduction | `test-blr-framework-phase17i.R` | 1 | shared/trait-specific/independent operators; Phase 17C remains dense oracle owner |

## Historical evidence and hashes

| Contract | Permanent owner | Fixture | Tier | Policy |
|---|---|---|---|---|
| defective MT fixed-B output | Phase 17B | historical config 2 RDS | 1 | read frozen object; never run current sampler |
| legacy MT denominator | Phase 17B | historical config 1 RDS | 1 | read frozen object; never run current sampler |
| Phase 17B fixture integrity | Phase 17B | three historical RDS files | 1 | permanent MD5/SHA integrity |
| corrected MT fixtures | Phase 17C | three corrected RDS files | 1 | exact structure, numeric tolerance 1e-12 |
| historical route dispositions | Phase 16A/17A documentation checks | reports/inventory | 2 | no repeated numerical execution |

Frozen binary fixtures and generated wrappers may retain hashes. Actively maintained numerical source files are protected by deterministic references, behavioral identities, and narrow architectural assertions instead of repeated whole-file hashes. Function-region hashes are temporary only while an untouched legacy function shares a translation unit with an active implementation; they should disappear when that legacy route is separately disposed.
| Public MT CSR normalization, biological alignment, named raw schema, and fit formatting | `test-blr-framework-phase17j.R` | ordinary fast/full | current |
| Source/installed test-context and path contract | `test-blr-framework-phase17j2.R` | ordinary fast/full | source assertions narrow-skip only |
| Canonical block-filtered storage/view and filter mathematics | `test-blr-framework-phase17k.R` | generated tiny BED | ordinary and installed; source structure blocks skip narrowly |
# Phase 17L ownership

`test-blr-framework-phase17l.R` owns MT block-eigen contracts, internal execution reductions, descriptor validation, transformed-summary identities, shared-core structure, and public/research exclusion. Phase 17K continues to own filtering mathematics and scalar block-eigen behavior.
