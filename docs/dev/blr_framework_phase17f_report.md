# BLR framework Phase 17F report

## 1. Executive summary

Phase 17F moved all posterior numerical finalization for the authoritative public default `mtblr()` route into one typed, binding-neutral function. `MtDefaultFinalResult` now owns finalized posterior means, traces, final state, and retained-count metadata. The public native adapter only constructs typed inputs, calls the core and finalizer, and copies values into the unchanged legacy 20-position result. Numerical execution, the public R interface, and all schemas are unchanged.

The phase also establishes shared ST/MT naming, canonical scalar-CSR reuse, per-trait LD operator, block-eigen reuse, marker alignment, and study/overlap requirements. It implements none of those future operators.

## 2. Repository baseline

- Branch: `master`.
- Starting and Phase 17E commit: `3a40a0b` (`Activate typed public multivariate core`).
- Initial status: clean.
- Visible hosted CI: no combined hosted status was available locally; no hosted-CI success is claimed.
- Toolchain: R 4.4.1, Rtools44 GNU C++17, Rcpp, RcppArmadillo, and OpenMP on Windows.
- Baseline compilation succeeded.
- Baseline full suite: 6,282 passed, 0 failures, 0 warnings, 12 skips.
- Phase 17E benchmark baseline: updated mean/median 0.032/0.010 s; fixed 0/0 s; explicit-set 0.008/0.010 s; moderate 0.064/0.040 s. Completed-fit RSS ranged from 121.4180 to 124.1328 MiB.

## 3. Public call graph

```text
sblr(algorithm = "default")
-> mtblr() R validation and dense-summary preparation
-> unchanged _sblr_mtblr native wrapper
-> thin public mtblr() adapter
-> run_mt_default_core(...)
-> MtDefaultCoreResult
-> finalize_mt_default_result(...)
-> MtDefaultFinalResult
-> unchanged 20-position positional mapping
-> unchanged R naming, orientation, and conditional correlations
```

## 4. Previous inline finalization

Phase 17E left posterior divisions in the public native adapter after the typed core returned. That boundary mixed numerical posterior means with positional allocation and copying. Phase 17F supersedes only that inline numerical-finalization portion; the core and positional compatibility boundary remain separate.

## 5. Finalization audit

The audited numerical work comprised six active division statements: marker-effect and inclusion accumulators by `marker_retained_count`; marker, genetic, and residual covariance accumulators by their respective counts; and pattern-probability accumulators by `pi_retained_count`. Four count-dependent branches preserved zero summaries when B, E, or pi updates were disabled. Legacy `pitrait` and `pimarker` accumulators are unsupported and remain zero.

Schema work remains in `mtblr.cpp`: `result.resize(20)`, per-position shaping, copying finalized values, copying borrowed `wy`, mapping unsupported positions 19–20, and returning the positional object. The line-level seam is the call to `finalize_mt_default_result(std::move(core_result))`; everything after its field aliases is compatibility allocation and copying, with no posterior arithmetic.

## 6. Typed finalized result

`MtDefaultFinalResult` owns:

- dimensions: `nt`, `m`, `nmodels`;
- retained counts: marker, B, G, E, and pi counts;
- finalized marker values: `bm`, `dm`, `r`, `b`, `d`, and `marker_order`;
- traces: `vbs`, `vgs`, and `ves`;
- posterior covariance means: `covb`, `covg`, and `cove`;
- final covariance matrices: `vb`, `vg`, and `ve`;
- probability results: `pi_final` and `pi_mean`;
- unsupported legacy values: `pitrait` and `pimarker`.

The descriptive probability names deliberately distinguish final probabilities from posterior means. The legacy adapter maps them to positions `pi` and `pim`. The type contains no R binding, schema, name, class, file, or operator field.

## 7. Finalization function

```cpp
inline MtDefaultFinalResult finalize_mt_default_result(
    MtDefaultCoreResult core_result);
```

It consumes the owning core result by value, moves reusable buffers, performs the existing divisions and disabled-update branches, and returns one binding-neutral finalized result. It contains no RNG, sampler, Rcpp construction, positional allocation, or I/O.

## 8. Arithmetic preservation

```text
PHASE17E_FINALIZATION_STATEMENTS=6
PHASE17F_FINALIZATION_STATEMENTS=6
ARITHMETIC_ORDER_EQUIVALENT=TRUE
DENOMINATOR_POLICY_EQUIVALENT=TRUE
DISABLED_UPDATE_POLICY_EQUIVALENT=TRUE
```

Loop nesting, trait/marker/state order, container types, and division order are unchanged. The only nonqualification changes are ownership moves and removal of meaningless division of permanently zero unsupported `pitrait`/`pimarker` values; their observable zeros are unchanged.

## 9. Responsibility boundaries

- `MtDefaultCoreResult`: unnormalized posterior accumulators, traces, final state, and contribution counts.
- `MtDefaultFinalResult`: finalized posterior means, traces, final state, and the same counts.
- Legacy adapter: positions, sizes, compatibility copying, borrowed `wy`, and unsupported legacy fields.

Neither typed result owns R-facing schema concerns.

## 10. Ownership and copies

The core result is transferred by value into the finalizer. Marker arrays, states, order, traces, accumulators, final matrices, probabilities, and unsupported fields are moved where safe. `pi_mean` requires one small O(K) finalized allocation. Moved-from core fields are not read afterward. The finalized result owns its data until the legacy adapter performs unavoidable compatibility copies into nested result positions. Borrowed `wy` is copied only into position 3. No new dense O(nt x m^2) copy or duplicate MCMC workspace was introduced.

## 11. Legacy positional schema

Positions remain exactly: 1 `bm`, 2 `dm`, 3 `wy`, 4 `r`, 5 `b`, 6 `d`, 7 `o`, 8 `vbs`, 9 `vgs`, 10 `ves`, 11 `covb`, 12 `covg`, 13 `cove`, 14 `vb`, 15 `vg`, 16 `ve`, 17 `pi`, 18 `pim`, 19 `pitrait`, and 20 `pimarker`. Types, dimensions, trait-major orientation, marker/state/iteration order, and unsupported zero semantics are unchanged.

## 12. Public R interface

`R/interface_mtblr.R`, generated wrappers, native signature, `NAMESPACE`, public arguments/defaults, hidden seed behavior, `.Call` target, public names, orientation, conditional correlation fields, and omission behavior are byte-identical.

## 13. Shared ST/MT naming plan

| Scientific concept | Current scalar name | Current MT legacy name | Phase 17F typed name | Planned canonical shared name | Compatibility decision |
|---|---|---|---|---|---|
| sample size | `n` | `n` | `data.n` | `sample_size`/`sample_sizes` | singular for ST, ordered per trait/study for MT |
| marker IDs | marker names | marker names | future metadata | `marker_ids` | shared identity convention |
| trait IDs | trait names | trait names | future metadata | `trait_ids` | MT ordered labels; ST one trait |
| summary cross-products | `wy` | `wy` | `data.wy` | `wy` | same scientific meaning |
| marker diagonal | `ww` | `ww` | `data.ww` | `ww` | same meaning, per trait for MT |
| sparse LD | `SparseLdCsrView` | experimental CSR variant | future operator bundle | `SparseLdCsrView` | reuse scalar representation exactly |
| dense LD | dense summary input | `XXvalues`/`XXindices` | borrowed data view | `ld_operator` vocabulary | preserve legacy names until schema phase |
| initial marker effects | initial `b` | `b` | `initial_state.b` | `initial_effects` | explicit initial-state meaning |
| posterior mean effects | `bm` | `bm` | `final_result.bm` | `effect_mean`/public `bm` | public compatibility retained |
| posterior inclusion | `dm` | `dm` | `final_result.dm` | `inclusion_mean`/public `dm` | shared meaning |
| final effects | `b` | `b` | `final_result.b` | `effect_final` | distinguish final from mean |
| final states | `d` | `d` | `final_result.d` | `state_final` | shared concept, model-specific state space |
| effect variance/covariance | `vb` | `B`/`vb` | `vb`, `covb` | effect variance/covariance | retain scalar-vs-matrix distinction |
| genetic variance/covariance | `vg` | `G`/`vg` | `vg`, `covg` | genetic variance/covariance | retain scalar-vs-matrix distinction |
| residual variance/covariance | `ve` | `E`/`ve` | `ve`, `cove` | residual variance/covariance | retain scalar-vs-matrix distinction |
| mixture/pattern probabilities | `pi`/`pis`/`pim` | `pi`/`pim` | `pi_final`/`pi_mean` | final/mean/trace names explicit | do not conflate traces and means |
| iterations | `nit` | `nit` | `execution.nit` | `iterations` | shared execution vocabulary |
| burn-in | `nburn` | `nburn` | `execution.nburn` | `burnin` | shared execution vocabulary |
| thinning | `nthin` | `nthin` | `execution.nthin` | `thinning` | shared execution vocabulary |
| input metadata | `input` | legacy R fields | future metadata | `input` | preserve public schema until deliberate migration |
| diagnostics | backend fields | legacy route-specific | none in numerical results | `diagnostics` | binding layer owns presentation |

Names are shared only where meanings agree. Scalar variance and multivariate covariance, and final/mean/trace probabilities, remain explicitly distinct.

## 14. Shared CSR requirement

Future canonical MT execution must reuse the existing canonical scalar `SparseLdCsrView`; no MT-specific CSR format is permitted. An MT operator bundle will hold one CSR view per trait/study. This is a design contract only; Phase 17F adds no CSR execution or representation.

## 15. Trait-specific LD support

With shared structure, immutable row offsets, column indices, and marker order may be shared while LD values, `ww`, `wy`, and `n` remain trait-specific. With independent structure, every trait may have its own offsets, indices, and values. Different ancestries, cohorts, populations, panels, and thresholds must not be forced into an artificial union sparsity pattern.

## 16. Shared block-eigen requirement

Future canonical MT block-eigen execution must reuse the scalar representation and conventions, with one operator per trait/study. Blocks, eigenvectors, eigenvalues, retained rank, tolerance, and reference metadata may differ. Decompositions may be shared only for genuinely identical LD matrices. No block-eigen code changed.

## 17. Marker and allele alignment

Marker row i in every trait-specific summary and LD operator must represent the same canonical marker ID and effect-allele orientation. The R validation/alignment layer owns ID matching, trait/study ordering, allele orientation, duplicates, missing-marker policy, and panel metadata. The first canonical MT CSR route should use explicit marker intersection unless a scientifically specified union-with-mask likelihood is separately established.

## 18. Study, ancestry, and overlap metadata

Future metadata must record trait, study, population, ancestry, LD reference, sample size, marker set, and sample-overlap policy. Independent, known-overlap, and unknown-overlap studies must be distinguished explicitly. Residual covariance `E` is not silently treated as a GWAS sample-overlap model. Phase 17F implements no overlap likelihood.

## 19. Phase 17B historical fixtures

All three historical RDS files remain byte-identical. SHA-256 values are:

- config 1: `82AF2F814C48BC6E5E4B8D7F748DC26BB7E2BC7F058D0D784D8F4508FDC874C9`
- config 2: `CF32B7C68BA943855BE131E0E0BB4C24AB76749DA6C8747A78870E601EED869C`
- config 3: `589E3A634BAB335C06E6F1D47063F28A0ACFE3B1C62F1EBC3A57A52DED2213D9`

## 20. Phase 17C corrected references

Raw references passed 3/3 and formatted references passed 3/3 with exact structure and numerical tolerance 1e-12. Corrected fixture SHA-256 values remain:

- config 1: `2AD9543841591A58BE91AD568E9E6A2E7CF424376EA1866D24ED791A6C64C806`
- config 2: `4F19291B02BE4AA82201E63BFC8701BCD3579B9E2AE8467BFD0938128AF91F3A`
- config 3: `8F7BD269392839D15F3CBE1A02DD2FFAE781C83E9EA70DA8ECC12C2B79FB7001`

## 21. Reproducibility

Same-process A;A and A;B;A comparisons, OMP thread environment 1 versus 2, intervening canonical CSR fit, and intervening packed-BED fit all passed with exact structure and tolerance 1e-12. The opt-in Phase 17F fresh-process comparison passed against corrected formatted configuration 1. Multiple chains remain unsupported.

## 22. Scientific identities

Finite outputs, marker/trait dimensions, names/order, binary states, bounded `dm`, covariance symmetry and positive semidefiniteness within tolerance, trace lengths/orientation, fixed B/E/pi controls, normalized updated `pim`, retained counts, zero disabled summaries, and covariance-diagonal/trace-mean identities all passed.

## 23. Alternative multivariate protection

`mt_cpg.cpp`, `mt_cpg_arma.cpp`, `mt_cpg_omp.cpp`, and `mt_cpg_omp_csr.cpp` are unchanged. Protected `mtblr_hybrid()` and `mtblr_eigen()` regions are unchanged. No alternative route calls the finalizer. The `mtblr_cpg_omp()` worker-sensitive RNG risk classification remains active.

## 24. Scalar, packed-BED, and block-eigen protection

Canonical scalar CSR numerical sources/references, canonical packed-BED sources/references and shared infrastructure, and block-eigen numerical sources remain unchanged. Generated wrappers and `NAMESPACE` remain unchanged.

## 25. Performance, memory, copies, and I/O

Phase 17F matched-workload measurements were:

| Workload | timings (s) | mean | median | completed-fit RSS MiB |
|---|---|---:|---:|---:|
| updated covariance/probability | 0.11, 0.01, 0.00, 0.02, 0.02 | 0.032 | 0.020 | 121.0430 |
| fixed B/E/pi | 0.00, 0.00, 0.00, 0.00, 0.01 | 0.002 | 0.000 | 121.1562 |
| explicit marker sets | 0.00, 0.00, 0.00, 0.00, 0.00 | 0.000 | 0.000 | 121.1562 |
| moderate dense | 0.03, 0.03, 0.06, 0.04, 0.05 | 0.042 | 0.040 | 121.1016 |

Compared with Phase 17E, these tiny-workload differences show no unexplained material regression. Completed-fit RSS is not peak RSS; sampled peak RSS was not available. Dense XX remains O(nt x m^2). The core-to-final move adds no such copy; legacy nested output copying remains the compatibility cost. There is no MCMC-time I/O, and tiny timings are regression signals rather than performance rankings.

## 26. CI coverage

The fast workflow filter includes ordinary Phase 17F tests while retaining prior coverage. Extended CI sets `SBLR_RUN_PHASE17F_FRESH=true`; the fresh test is opt-in locally. Benchmarks are not in the fast gate. Workflow files were inspected locally, but hosted CI success is not claimed.

## 27. Tests

- Phase 17B historical, Phase 17C corrected, Phase 17D extraction, and Phase 17E typed-core focused tests passed.
- Phase 17F ordinary: 180 passed, 0 failures, 0 warnings, 1 expected fresh-process skip.
- Phase 17F fresh-process enabled: 181 passed, 0 failures, 0 warnings, 0 skips.
- Corrected raw: 3/3; corrected formatted: 3/3.
- Alternative and canonical protected-backend checks passed.
- Full suite: 6,462 passed, 0 failures, 0 warnings, 13 opt-in skips.

## 28. Deviations and blockers

No blockers remain. The finalizer omits division of permanently unsupported zero `pitrait`/`pimarker` accumulators; their observable values and legacy positions are unchanged. Peak RSS was not measured, and hosted CI status was not visible. No public schema, chain aggregation, named converter, generic operator, CSR, eigen, or overlap implementation was introduced.

## 29. Recommended next phase

Audit and formalize reuse of the canonical scalar `SparseLdCsrView`, marker-alignment contract, and validation vocabulary for a trait-specific multivariate LD-operator bundle, supporting shared CSR structure with trait-specific values and fully independent CSR structures, without yet changing the public default dense route.

## 30. Readiness marker

PHASE 17F COMPLETE — TYPED MULTIVARIATE FINALIZATION AND SHARED NAMING CONTRACT ACTIVE
