# Unified BLR Framework Phase 17C Report

## 1. Executive summary

The authoritative supported public legacy `mtblr()` route now honors the
marker-covariance update control and uses a coherent retained-sample policy.
This correction changed neither its public route nor its legacy positional and
formatted schemas, and it did not begin numerical extraction or migration.

## 2. Repository baseline

Phase 17C began clean on `master` at `19a64e0` (`Harden Phase 17B review
coverage`), the committed Phase 17B review-fix baseline. R 4.4.1 and Rtools 4.4
with g++/OpenMP compiled successfully. Baseline focused Phase 17B validation
passed 148 expectations with one opt-in skip; the baseline full suite passed
5,845 expectations with zero failures, zero warnings, and nine skips.

## 3. Defects corrected

The public loop previously ran `sampleBset()` and `sampleB_latent()` even when
`updateB=FALSE`. It also excluded `it=nburn`, thinned on absolute iteration
indices, and divided several accumulators by `nit` or one unrelated marker
count. These were respectively an update-control defect, an off-by-one error,
an incorrectly anchored thinning policy, and denominator/count mismatches.

## 4. Files changed

- `src/mtblr.cpp`: only public `mtblr()` guards, retained conditions, counters,
  denominators, and a zero-retained guard.
- `tests/testthat/test-blr-framework-phase17b.R`: historical fixtures are now
  tested as immutable pre-correction evidence rather than current production.
- `tests/testthat/test-blr-framework-phase17c.R`: corrected references,
  identities, reproducibility, source contracts, and protected-source audits.
- `tests/testthat/fixtures/blr_phase17c_mt_default_corrected/`: three separate
  corrected fixture pairs and their helper.
- `tools/fixtures/generate_blr_phase17c_mt_default_corrected_fixtures.R`: manual
  corrected-only generator.
- `tools/benchmarks/blr_phase17c_mt_default_corrected.R`: corrected benchmark.
- both framework workflows: ordinary/fresh Phase 17C coverage.
- implementation plan, capability/reduction matrices, computation inventory,
  naming guide, and this report: current corrected status.
- Phase 17A and Phase 9D1 tests: obsolete whole-file `mtblr.cpp` hashes were
  removed; Phase 17C protects the untouched alternative function regions.

## 5. Corrected `updateB` contract

The set-local `sampleBset()` and latent `sampleB_latent()` calls share one
`if (updateB)` guard. The later global `sampleB()` remains under its established
`updateB && method==4` guard. `Bi` is recomputed from the current `B` after the
guard, so fixed runs use the inverse of the supplied matrix. Runtime tests prove
the final `B` equals the supplied matrix exactly when updates are disabled.

## 6. Corrected iteration contract

The loop remains `it=0,...,nburn+nit-1`. Burn-in is `it < nburn`; post-burn is
`it >= nburn`. Marker summaries retain when
`(it - nburn) % nthin == 0`, so the first post-burn iteration is eligible and
the marker count is `ceiling(nit/nthin)`.

## 7. Accumulator-specific counts

| Accumulator | Condition | Count | Denominator | Disabled behavior |
|---|---|---:|---|---|
| `bm`, `dm` | post-burn and relative thinning | `ceiling(nit/nthin)` | `marker_retained_count` | valid configurations require a positive count |
| `covb` | `updateB && it >= nburn` | `nit` or 0 | `covb_retained_count` | zero matrix |
| `covg` | `it >= nburn` | `nit` | `covg_retained_count` | not disabled |
| `cove` | `updateE && it >= nburn` | `nit` or 0 | `cove_retained_count` | zero matrix |
| `pim` | `updatePi && it >= nburn` | `nit` or 0 | `pi_retained_count` | zero vector |

Finalization divides only for positive counts. Unsupported native positions
18--19 retain their existing zero semantics and `nit` denominator.

## 8. RNG-order preservation

With `updateB=TRUE`, the set-local and global call order is unchanged and no
draw site moved. Phase 17B versus Phase 17C comparisons preserve `wy`, `r`,
`b`, `d`, `o`, all variance traces, final `B/G/E`, and final `pi` within
`1e-12`. Fixed-`B` trajectories intentionally change because obsolete draws
are skipped; no dummy draws are consumed.

## 9. Historical Phase 17B fixtures

The historical RDS files remain byte-identical:

| Fixture | SHA-256 |
|---|---|
| config-1 | `82AF2F814C48BC6E5E4B8D7F748DC26BB7E2BC7F058D0D784D8F4508FDC874C9` |
| config-2 | `CF32B7C68BA943855BE131E0E0BB4C24AB76749DA6C8747A78870E601EED869C` |
| config-3 | `589E3A634BAB335C06E6F1D47063F28A0ACFE3B1C62F1EBC3A57A52DED2213D9` |

They retain the defective fixed-`B` trajectory and `(nit-1)/nit` probability
normalization as permanent historical evidence.

## 10. Corrected Phase 17C fixtures

Three new raw/formatted pairs use exact structure and numerical tolerance
`1e-12`. Counts are respectively: config 1 `9/18/18/18/18`, config 2
`14/0/14/0/0`, and config 3 `10/20/20/20/20` for marker/covb/covg/cove/pi.
Their SHA-256 hashes are `2AD9543841591A58BE91AD568E9E6A2E7CF424376EA1866D24ED791A6C64C806`,
`4F19291B02BE4AA82201E63BFC8701BCD3579B9E2AE8467BFD0938128AF91F3A`, and
`8F7BD269392839D15F3CBE1A02DD2FFAE781C83E9EA70DA8ECC12C2B79FB7001`.

## 11. Legacy-versus-corrected differences

For updated configurations, differences are confined to retained summaries:
`bm`, `dm`, covariance posterior means, `pim`, and formatted correlations
derived from corrected covariance means. The fixed-control configuration also
changes downstream states/traces because `B` is genuinely fixed. Its final
`B`, `E`, and `pi` exactly equal supplied values; disabled posterior summaries
remain zero.

## 12. Public schema

The 20-position native result, public BayesC field names, dimensions,
orientations, conditional correlation fields, arguments, route, and formatter
are unchanged. No `stblr_raw_v1` conversion was introduced.

## 13. Reproducibility

`A;A`, `A;B;A`, OMP thread environment 1 versus 2, intervening canonical CSR,
intervening canonical packed-BED, and fresh versus reused process comparisons
all pass with exact structure and numerical tolerance `1e-12`. Multiple chains
remain unsupported. The hidden R-derived fit seed contract is unchanged.

## 14. Scientific identities

Tests confirm finite outputs, stable marker/trait order, binary states, bounded
`dm`, symmetric positive-semidefinite covariance matrices, correct trace
lengths, fixed `B/E/pi`, finite corrected means, retained-count identities, and
`sum(pim)=1` whenever `updatePi=TRUE`. Covariance diagonals independently match
post-burn trace means using the corresponding accumulator schedule.

## 15. Alternative multivariate protection

`mt_cpg.cpp`, `mt_cpg_arma.cpp`, `mt_cpg_omp.cpp`, and `mt_cpg_omp_csr.cpp`
remain at their Phase 17B hashes. Normalized function-region hashes prove
`mtblr_hybrid()` (`03ed109628b874223db109d2ec654827`) and `mtblr_eigen()`
(`4e5c38ede3345de10a684ab38470bf7b`) are unchanged. The OpenMP CPG worker-seed
P0 risk remains explicitly tested.

## 16. Scalar and packed-BED protection

Canonical CSR, packed-BED, and block-eigen protected source hashes and their
permanent reference tests pass. Generated wrappers and `NAMESPACE` retain their
baseline hashes.

## 17. Performance, memory, and I/O

Five repetitions produced: updated tiny `0.19,0.02,0.00,0.01,0.01` s (mean
0.046, median 0.01; completed-fit RSS 131.88 MiB); corrected fixed `0.01,0,0,0,0`
s (mean 0.002, median 0; 120.77 MiB); explicit sets `0.01,0.02,0.01,0.02,0.02`
s (mean 0.016, median 0.02; 120.77 MiB); and moderate 80-marker/3-trait dense
`0.06,0.03,0.03,0.06,0.04` s (mean 0.044, median 0.04; 121.27 MiB).
Phase 17B matched tiny means were 0.026 s updated and 0.008 s fixed. These tiny
timings are regression signals, not speed claims. Completed-fit RSS is not peak
RSS; peak RSS was not measured. Dense `XX` remains `O(nt * m^2)`. All inputs
are in memory and there is no MCMC-time disk I/O or page-cache dependency.

## 18. CI coverage

The fast filter includes ordinary Phase 17C tests. Extended CI sets
`SBLR_RUN_PHASE17C_FRESH: "true"`; benchmarks are not in the fast gate.

## 19. Tests

- Phase 17B historical evidence: 148 passed, one opt-in skip locally.
- Phase 17C ordinary: 184 passed, one opt-in fresh-process skip.
- Phase 17C fresh enabled: 185 passed, no failure.
- Corrected raw references: 3/3; corrected formatted references: 3/3.
- Update-control, retained-count, identities, alternatives, and canonical
  protections: passed.
- Full suite: 6,030 passed, zero failures, zero warnings, ten opt-in skips.

## 20. Deviations and blockers

No correctness blocker remains for the two scoped policies. Numerical tolerance
remains `1e-12` because Armadillo solves can differ by a few ulps. Completed-fit
RSS is not peak RSS. Experimental `mtblr_cpg_omp()` retains its separate P0 RNG
risk; it was neither changed nor corrected here.

## 21. Recommended next phase

Mechanically extract the corrected authoritative public `mtblr()` single-fit
numerical execution into one guarded implementation boundary while preserving
Phase 17C references, covariance/state update order, fit-local RNG ownership,
legacy positional output, and public formatting.

## 22. Readiness marker

PHASE 17C COMPLETE — PUBLIC MULTIVARIATE CORRECTNESS CONTRACT ESTABLISHED
