# Retained block-eigen performance validation

The retained operator changes per-marker work from a packed dense block sweep of
order `m_b` to a factor-column update of order `k_b`. Its immutable payload is
`sum_b 4 m_b k_b` bytes for float factors, plus transformed scores, diagonals,
mapping, and diagnostics. Each chain adds `sum_b 8 k_b` bytes for double reduced
residuals. These quantities are reported by the native builder.

## Earlier illustrative R-level benchmark

An earlier version of `tools/validation/benchmark_low_rank_operator.R` was run
on the validation host for block sizes 250, 500, 1000, and 2000. The deliberately
truncated case used `k/m = 0.25`. These timings are R-level operator
microbenchmarks and are not native MCMC throughput measurements. The illustrative
loop remains available only when
`SBLR_RUN_ILLUSTRATIVE_R_BENCHMARK=true` is set explicitly.

| m | k | sparse sweep (s) | packed/dense sweep (s) | full-rank factor sweep (s) | 25% factor sweep (s) | eigendecomposition (s) |
|---:|---:|---:|---:|---:|---:|---:|
| 250 | 63 | 0.0005 | <0.0005 | 0.0005 | <0.0005 | 0.03 |
| 500 | 125 | 0.0005 | 0.0015 | 0.0030 | 0.0005 | 0.04 |
| 1000 | 250 | 0.0030 | 0.0070 | 0.0075 | 0.0020 | 0.16 |
| 2000 | 500 | 0.0030 | 0.0405 | 0.0425 | 0.0125 | 0.95 |

The corresponding immutable payload estimates were:

| m | sparse bytes | packed float bytes | full-rank factor bytes | 25% factor bytes | packed chain bytes | 25% reduced chain bytes |
|---:|---:|---:|---:|---:|---:|---:|
| 250 | 7,992 | 125,500 | 254,000 | 65,504 | 2,000 | 504 |
| 500 | 15,992 | 501,000 | 1,008,000 | 255,000 | 4,000 | 1,000 |
| 1000 | 31,992 | 2,002,000 | 4,016,000 | 1,010,000 | 8,000 | 2,000 |
| 2000 | 63,992 | 8,004,000 | 16,032,000 | 4,020,000 | 16,000 | 4,000 |

## Native retained-factor hot-path benchmark

The updated benchmark calls development-only native entry points and excludes
factor construction, eigendecomposition, rank selection, transformed-score
construction, R object construction, and file I/O from measured regions. It
uses one 1,000-marker block, one requested OpenMP thread, a fixed seed, two
warm-ups, and seven measured repetitions. Times below are medians in seconds;
parentheses contain minima.

The validation host was Windows 10 x64 with 12 logical processors. The package
was built with GCC 13.2.0, OpenMP enabled, and effective optimized flags
`-O2 -Wall -mfpmath=sse -msse2 -mstackrealign` (the command also contained the
earlier `-O3 -march=native`, overridden by the later `-O2`).

| Operation | k=1000 | k=750 | k=500 | k=250 |
|:--|--:|--:|--:|--:|
| 1,000 `corrected_rhs` calls | 0.000167 (0.000166) | 0.000131 (0.000130) | 0.000089 (0.000089) | 0.000049 (0.000048) |
| 1,000 `apply_difference` calls | 0.000143 (0.000142) | 0.000122 (0.000122) | 0.000081 (0.000080) | 0.000047 (0.000047) |
| one `fitted_quadratic` | 0.00000150 (0.00000149) | 0.00000110 (0.00000109) | 0.00000074 (0.00000071) | 0.00000031 (0.00000031) |
| one direct `quadratic_form` reference | 0.002187 (0.002057) | 0.001639 (0.001541) | 0.001026 (0.001025) | 0.000516 (0.000506) |
| one complete residual rebuild | 0.001259 (0.001231) | 0.000930 (0.000922) | 0.000617 (0.000614) | 0.000307 (0.000307) |
| one BayesC marker sweep | 0.000460 (0.000456) | 0.000462 (0.000420) | 0.000355 (0.000354) | 0.000294 (0.000293) |
| one BayesR marker sweep | 0.000838 (0.000796) | 0.000838 (0.000766) | 0.000735 (0.000705) | 0.000652 (0.000647) |

The complete BayesC Gibbs measurement enables marker updates, the BayesC
mixture-probability update, effect-variance update, residual-variance update,
and genetic-variance calculation. The BayesR measurement enables marker and
component updates, component-probability update, effect-variance update,
residual-variance update, and genetic-variance calculation. Its direct
reference repeats the pre-correction arithmetic (direct effect-based
quadratics and residual reconstruction); it is not a separately built
pre-patch binary.

| Complete iteration | k=1000 | k=750 | k=500 | k=250 |
|:--|--:|--:|--:|--:|
| BayesC optimized (s) | 0.000748 (0.000461) | 0.000427 (0.000418) | 0.000358 (0.000357) | 0.000297 (0.000295) |
| BayesC direct reference (s) | 0.001711 (0.001671) | 0.001352 (0.001331) | 0.000964 (0.000960) | 0.000599 (0.000597) |
| BayesC relative median speedup | 2.29x | 3.17x | 2.70x | 2.02x |
| BayesR optimized (s) | 0.000798 (0.000789) | 0.000767 (0.000765) | 0.000714 (0.000710) | 0.000793 (0.000756) |
| BayesR direct reference (s) | 0.006218 (0.005729) | 0.004507 (0.004428) | 0.003141 (0.003129) | 0.001933 (0.001924) |
| BayesR relative median speedup | 7.79x | 5.88x | 4.40x | 2.44x |

The principal ordinary-iteration saving is state-aware evaluation of the
fitted quadratic as `sum((w - r)^2)`, which is order `sum_b k_b`, instead of a
fresh factor multiplication of order `sum_b m_b k_b`. Periodic rebuilding is
still required to measure and repair accumulated residual drift, but its
default frequency is once per 100 completed retained-low-rank iterations plus
one final check.

The evidence supports only the conditional claim: retained factors help when
effective rank is materially smaller than marker count. A full-rank float factor
uses about twice the immutable storage of a packed float dense block and had no
speed advantage here. CSR remains preferable for genuinely sparse blocks, and
packed dense storage can remain preferable for small, dense, nearly full-rank
blocks. This task does not add automatic hybrid selection.

## Baseline test isolation

At clean commit `96487b3194fc1f8c6789060da5f2e2a0eea89974`, the packed-BED
fixed-alpha BayesRC versus fixed-pi BayesR reduction test already has its three
current failures: `bm[1, 1]` is 0.1050150 versus 0.1440519, `dm[1, 1]` is 0.3
versus 0.4, and the first marker's component probabilities are
`c(0.7, 0.1, 0.1, 0.1)` versus `c(0.6, 0.2, 0.2, 0.0)`. These baseline
failures are unrelated to the retained-low-rank hot path; this optimization
does not change packed-BED behavior, BayesR prior calibration, or BayesRC.
