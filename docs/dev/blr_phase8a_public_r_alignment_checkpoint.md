# Phase 8A public R alignment checkpoint

## Status and boundary

Phase 8A simplifies the public R framework without changing any native
sampler, posterior target, seed, schedule, retained index, or convergence
capture. Raw-v2 remains authoritative. Formatted `stblr_fit` and `mtblr_fit`
objects retain the Phase 0--7 fields; the new access functions provide one
canonical route to those fields and to retained raw-v2 arrays.

The cleanup removes the unreachable pre-Phase-5 MT covariance-hybrid R stack,
the standalone stale capability matrix, and obsolete tests that existed only
for the removed hybrid. Current capability remains owned by executable public
dispatch, resolved specifications, raw schemas, and focused route tests.

## Before and after architecture

Before Phase 8A, the supported MT wrappers coexisted with approximately 1,500
lines of unreachable BED/CSR covariance-hybrid preparation, separate BayesR
and BayesRC normalization helpers, and a legacy raw-v1 multi-chain formatter.
Scientific quantities were also read independently by credible-set,
fine-mapping, component-summary, and workflow helpers.

After Phase 8A the maintained R path is:

```text
thin public ST/MT adapter
  -> shared exact-name and execution-control normalization
  -> model-specific resolved specification and provider preparation
  -> unchanged native sampler
  -> validated blr_raw v2
  -> established shared formatter
  -> extract_posterior() / summarise_posterior() / extract_diagnostics()
```

Provider-specific adapters and model-specific prior resolution remain separate
where their scientific meanings differ.

## Canonical public inputs

| Scientific concept | Canonical argument | ST | MT | Meaning | Validation owner |
| --- | --- | --- | --- | --- | --- |
| model | `method` | yes | yes | Declared supported model for the operator | public adapter and resolved model |
| post-burn iterations | `nit` | yes | yes | Sampling iterations after burn-in under the maintained route contract | shared execution controls |
| burn-in | `nburn` | yes | yes | Completed burn-in iterations | shared execution controls |
| thinning | `nthin` | yes | yes | Retention interval | shared execution controls |
| chains | `nchains` | yes | yes | Logical chains | shared execution controls |
| seed | `seed` | yes | yes | Master seed under seed-contract v1 | shared execution controls |
| per-chain seed bases | `chain_seeds` | yes | yes | Optional chain seed bases | shared execution controls |
| workers | `ncores` | yes | yes | Requested chain workers | shared execution controls |
| compact chain records | `keep_chains` | yes | yes | Retain formatted chain-final records | shared execution controls |
| convergence mode | `convergence` | yes | yes | Declared convergence capture policy | public execution boundary |
| convergence controls | `convergence_control` | yes | yes | Extended ST controls; `NULL` on current MT routes | public execution boundary |
| trace retention | `keep_traces` | model dependent | yes | Retain contracted convergence traces | execution contract |
| memory ceiling | `memory_limit_bytes` | route dependent | yes | Incremental fit-allocation ceiling | R preflight and native guards |
| marker IDs | `global_marker_ids` or provider-specific BED metadata | route dependent | yes | Stable global marker order | provider/resolved-spec boundary |
| traits | `trait_ids` | implicit or route dependent | yes | Stable trait order | provider/resolved-spec boundary |
| marker scale | `marker_multipliers` | model dependent | BayesR | Positive marker-specific variance multiplier | model prior resolver |
| positive scales | `component_scales` | BayesR | BayesR | Ordered positive mixture scales | model prior resolver |
| activity prior | `activity_pattern_dirichlet_prior` | no | yes | Dirichlet prior over MT activity patterns | MT resolved specification |
| marker covariance prior | `marker_covariance_prior_degrees_of_freedom`, `marker_covariance_prior_scale` | no | yes | Inverse-Wishart prior for full $V_b$ | MT resolved specification |

ST inclusion probabilities, MT activity-pattern probabilities, ST scalar
variances, MT covariance matrices, and summary-provider residual scales remain
deliberately distinct.

## Canonical outputs

| Scientific quantity | Raw-v2 field | Formatted field | ST shape | MT shape | Canonical access |
| --- | --- | --- | --- | --- | --- |
| realised-effect mean | `posterior$realised_effect_mean` | `bm` | marker x trait | marker x trait | `extract_posterior(fit, "realised_effects")` |
| trait PIP | `posterior$pips` | `dm` | marker x trait | marker x trait | `extract_posterior(fit, "pips")` |
| retained realised effects | `draws$realised_effects` | retained in raw-v2 | draw x chain x marker x trait | same | state `"retained"` |
| final realised effects | `final$realised_effects` | `b` convenience view | chain x marker x trait in raw | same | state `"final"` |
| component probabilities | `posterior$component_probabilities` | `component_probabilities` | trait-indexed marker x component | marker x component for MT pattern-scale | `"component_probabilities"` |
| activity-pattern probabilities | `posterior$activity_pattern_probabilities` | same | unavailable | marker x pattern | `"activity_pattern_probabilities"` |
| all-active probability | `posterior$pleiotropic_probabilities` | same | unavailable | marker | `"pleiotropic_probabilities"` |
| effect covariance | `posterior$marker_covariance_mean` and `draws$marker_covariance` | `cov_b_mean` and retained raw | scalar/size-one trait convention where applicable | trait x trait | `"effect_covariance"` |
| genetic covariance | authoritative stored posterior-derived value where defined | `cov_g_mean` | `NULL` for models that do not define or retain it | trait x trait, including 1 x 1 | `"genetic_covariance"` (posterior only) |
| residual covariance | `posterior$residual_covariance_mean` and `draws$residual_covariance` | `cov_e_mean` and retained raw | scalar/size-one trait convention where applicable | trait x trait or `NULL` | `"residual_covariance"` |
| predictions | `draws$predictions` where available | `predictions` | observation x trait | observation x trait or `NULL` | `"predictions"` |
| execution diagnostics | `diagnostics` | `diagnostics` | named namespace | named namespace | `extract_diagnostics()` |

Extraction preserves every declared dimension, including one trait and one
chain. Scientifically unavailable values return `NULL`; the fixed raw-v2
schema continues to require its established present-but-`NULL` fields.

## Access responsibilities

- The formatter converts validated raw-v2 to the existing stable fit fields.
  It does not define a second posterior.
- `extract_posterior()` retrieves a stored posterior, retained, or final
  quantity. It supports explicit axis selection and never summarizes or drops
  dimensions.
- `summarise_posterior(..., quantity = ...)` summarizes the retained draw and
  chain axes from raw-v2. It does not remove burn-in again. The established
  global ST trace-summary mode remains available when `quantity` is omitted.
- `extract_diagnostics()` owns diagnostic retrieval and returns explicit
  `sampler`, `execution`, `convergence`, `providers`, and `provenance`
  namespaces. It assembles existing raw-v2 fields—including task seeds and
  IDs, retained/convergence indices, worker/scheduler state, provider/resource
  maps, and provenance—without recomputation or posterior-array copies.

Credible-set, fine-mapping, component-summary, and maintained workflow helpers
now obtain effects and PIPs through the canonical extractor.

## Validation ownership

Public wrappers own exact spelling and basic boundary types. Resolved specs own
scientific and cross-field consistency. Providers own resources, maps, coding,
and ordering. Native code retains immediate allocation and numerical guards.
Raw-v2 validation owns returned scientific axes, identities, and unavailable
fields. Phase 8A removes only redundant or unreachable R validation; it does
not remove an independent native safety guard.

## Removed code and metrics

The following were removed rather than wrapped:

- the legacy MT BED covariance-hybrid wrapper and formatter;
- the legacy MT CSR normalization and raw-v1 formatter;
- legacy MT BayesR/BayesRC specification helpers;
- the legacy MT summary multi-chain combiner;
- the stale standalone capability matrix;
- obsolete tests whose only subject was that removed stack;
- duplicate direct PIP/effect/diagnostic access helpers.

The stale executable-tool audit classified the two Study 06 cache drivers and
the particle-reference driver as obsolete research launchers; their scientific
record remains in the explicitly research-only Study 06 developer documents.
The old block-residual validation launcher was already non-executable at clean
HEAD and is superseded by the maintained block-residual qualification record
and tests. The old model benchmark targeted removed MT BayesRC/hybrid controls.
Those five launchers were deleted rather than given compatibility wrappers;
current architecture and generated-interface audits remain maintained.

For the reproducible scope of public MT adapters, unified/raw access,
summaries, credible sets, and fine-mapping, the cleanup changes 8,383 R source
lines in 13 files to 6,199 lines in 11 files, and 155 top-level functions to
113. Compatibility/alias line hits fall from 62 to 49. Directly relevant tests
fall from 1,136 lines in seven files to 770 lines in six files. Package exports
increase from 31 to 33 only because `extract_posterior()` and
`extract_diagnostics()` establish the canonical access routes.

The extractor/summarizer/diagnostic implementation count uses this explicit
rule: count top-level functions in the scoped R files whose primary
responsibility is scientific extraction, posterior summarization, or
diagnostic validation/access. It is 10 before and 11 after. The clean-checkpoint
functions are `.blr_validate_compact_transition_diagnostics()`,
`summarise_posterior()`, `summarise_components()`, `.stblr_extract_pip()`,
`.extract_sparseLD_region_dense()`, `.stblr_extract_trait_pip()`,
`.stblr_get_global_parameter()`, `.stblr_extract_stat_vector()`,
`.stblr_extract_stat_scalar()`, and `.stblr_extract_bm()`. The Phase 8A
functions are `.blr_validate_compact_transition_diagnostics()`,
`summarise_posterior()`, `summarise_components()`, `.stblr_extract_pip()`,
`.extract_sparseLD_region_dense()`, `.stblr_get_global_parameter()`,
`.stblr_extract_stat_vector()`, `.stblr_extract_stat_scalar()`,
`extract_posterior()`, `extract_diagnostics()`, and
`.blr_summarise_retained()`. The count increases because the two public access
functions make responsibilities explicit; duplicated caller implementations
and their LOC still decrease materially.

## Scientific-neutrality evidence

Clean-checkpoint outputs were frozen before editing for ST CSR BayesC, ST CSR
BayesR, fixed-marker annotation, MT BED BayesC, MT CSR BayesR, and a two-chain
retained-rank MT block-eigen fit. The comparison includes posterior fields,
retained draws, final states, covariance and probability quantities, task
seeds, retained indices, convergence, diagnostics, and canonical formatted
fields. All six frozen objects are bitwise identical after normalizing only the
disposable temporary BED-file path. Phase 8A changes no native source or
generated native interface.

## Deferred work and Phase 8B handoff

Phase 8A does not unify scientifically different prior namespaces, introduce a
generic provider class, or rewrite legacy ST trace storage. A separate Phase
8B in `sblrbench` should use `extract_posterior()`, quantity-mode
`summarise_posterior()`, and `extract_diagnostics()`; remove direct `$bm`,
`$dm`, and diagnostic traversal; consolidate duplicate study helpers; and move
heavy validation to the benchmark repository. No `sblrbench` file is changed
by this checkpoint.
