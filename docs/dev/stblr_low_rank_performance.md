# Retained block-eigen performance validation

The retained operator changes per-marker work from a packed dense block sweep of
order `m_b` to a factor-column update of order `k_b`. Its immutable payload is
`sum_b 4 m_b k_b` bytes for float factors, plus transformed scores, diagonals,
mapping, and diagnostics. Each chain adds `sum_b 8 k_b` bytes for double reduced
residuals. These quantities are reported by the native builder.

`tools/validation/benchmark_low_rank_operator.R` was run on the validation host
for block sizes 250, 500, 1000, and 2000. The deliberately truncated case used
`k/m = 0.25`; timings are R-level operator microbenchmarks and are not presented
as native MCMC throughput claims.

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

The evidence supports only the conditional claim: retained factors help when
effective rank is materially smaller than marker count. A full-rank float factor
uses about twice the immutable storage of a packed float dense block and had no
speed advantage here. CSR remains preferable for genuinely sparse blocks, and
packed dense storage can remain preferable for small, dense, nearly full-rank
blocks. This task does not add automatic hybrid selection.
