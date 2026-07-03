# C++ Code Review: sblr Package

## Context

This is a read-only inspection of the sblr R package's C++ code. No files have been modified.
The goal is to document suspected bugs, numerical risks, and R/C++ interface issues with suggested tests.

---

## 1. Critical Bugs

### 1.1 Double sampling of B within one MCMC iteration — `mt_cpg.cpp` and `mt_cpg_arma.cpp`

**Severity: Critical (silent wrong posterior).**

When `method == 4` and `updateB = TRUE`, the marker-variance matrix B is sampled **twice** per iteration:

- **First** inside the `sets` loop (to obtain Bi for marker updates):
  - `mt_cpg.cpp:313` — `sampleB_cpg(nt, m, nub, B, beta, ssb_prior, gen);`
  - `mt_cpg_arma.cpp:302` — `sampleB_cpg_arma(nt, m, nub, B, beta_mat, ssb_prior_mat, gen);`
- **Second** in the `updateB` block after the marker loop:
  - `mt_cpg.cpp:362` — `sampleB_cpg(nt, m, nub, B, beta, ssb_prior, gen);`
  - `mt_cpg_arma.cpp:355` — `sampleB_cpg_arma(nt, m, nub, B, beta_mat, ssb_prior_mat, gen);`

After the second sampling, `Bi` is **not recomputed**. The next iteration then uses a stale `Bi` that does not correspond to the B that was just stored. The developer left a comment in `mt_cpg_arma.cpp:298-301` flagging this:

```cpp
// NOTE:
// This preserves your current logic:
// B is sampled once per set.
// If this was not intended, move it outside this loop.
```

**Additional issue:** Inside the `sets` loop, `sampleB_cpg` is called once per set, each time using **all m markers** (not just those in the current set). This means B is re-sampled `sets.size()` times during marker updates with inconsistent partial effects, which is statistically wrong.

**Suggested tests:**
- Run with `sets = list(1:m)` (one set, all markers) and compare to a build with the set-loop B-sampling removed. Posteriors should match but won't if the second B sample is wrong.
- Run with `updateB = TRUE` and check whether `vbs` traces look overdispersed.

---

### 1.2 Off-by-one in burn-in boundary — `mt_cpg.cpp` and `mt_cpg_arma.cpp`

**Severity: High (systematically discards one post-burn-in sample).**

All posterior accumulation conditions use `it > nburn` (strictly greater):

| File | Lines | Condition |
|---|---|---|
| `mt_cpg.cpp` | 344, 356, 368, 380, 393 | `it > nburn` |
| `mt_cpg_arma.cpp` | 333, 347, 361, 375, 394, 399 | `it > nburn` |

The loop runs `it = 0 ... nit+nburn-1`. The first post-burn-in iteration is `it == nburn`, but with `it > nburn` it is **excluded**. This means the first usable sample is from `it == nburn + 1`, so at most `nit - 1` samples are actually collected when `nthin == 1`. If `nit == 1`, **zero** samples are collected and `nsamples == 0`.

The condition should be `it >= nburn`.

**Suggested test:**
- Run `mtblr_cpg(... nit=1, nburn=0, nthin=1, ...)`. With `it > nburn` the single draw at `it=0` is excluded (`0 > 0` is FALSE), `nsamples = 0`, and R sees NaN in `bm`. With `it >= nburn` it works correctly.

---

### 1.3 Division by zero when `nsamples == 0` — `mt_cpg.cpp:453`, `mt_cpg_arma.cpp:436`

**Severity: High (crash or NaN output).**

Both files divide accumulated statistics by `nsamples` at the output stage:

```cpp
result[0][t][i] = bm[t][i] / nsamples;   // posterior mean of b
result[1][t][i] = dm[t][i] / nsamples;   // posterior inclusion prob
```

Also `cvbm/nsamples`, `cvgm/nsamples`, `cvem/nsamples`, `pis[i]/nsamples` in the same block.

If `nsamples == 0` (which happens with the off-by-one above when `nit == 1, nburn == 0`), all these produce `NaN` or `Inf` without an error. There is no guard.

**Suggested test:**
- Call `mtblr_cpg(nit=1, nburn=0, nthin=1, ...)` and check the returned list for NaN.

---

## 2. Medium-Priority Bugs

### 2.1 Unchecked Cholesky for the selected-model draw

**Severity: Medium (silent garbage on near-singular matrices).**

In both `mt_cpg.cpp:152` and `mt_cpg_arma.cpp:110`, after model selection the Cholesky is called again to draw the posterior:

```cpp
arma::chol(L, C, "lower");   // return value discarded
```

The model-loop Cholesky was checked (`!arma::chol(L, C)` set `loglik[k] = -INFINITY`), but numerically the matrix could degenerate between the two calls (e.g. due to a slightly different code path). If this Cholesky fails, `L` retains its previous value and the draw is silently wrong.

**Suggested test:**
- Intercept cases where `ww[t][i]` is very small, forcing near-zero `C(t,t)`, and verify the sampler does not produce silent NaN.

---

### 2.2 Unguarded division by chi-squared draw in `sampleE_cpg` / `sampleE_cpg_arma`

**Severity: Medium (possible infinite variance sample).**

`sampleE_cpg` (`cpg_samplers.cpp:234-237`) and `sampleE_cpg_arma` (`cpg_samplers.cpp:65-68`):

```cpp
std::chi_squared_distribution<double> rchisq(n[t] + nue);
double chi2 = rchisq(gen);
E(t, t) = Se(t, t) / chi2;   // no guard
```

There is no `std::max(chi2, 1e-300)` guard, unlike `sampleE_ST_csr` (`st_csr_common.h:340`) which correctly writes:

```cpp
const double chi2 = std::max(rchisq(gen), 1e-300);
```

Additionally, neither function checks that `Se(t,t) > 0`. Accumulated floating-point error in the residuals can cause `sse` to go slightly negative, giving a negative scale. Division of a negative scale by a small positive chi2 produces a large negative variance, which then propagates to the next Cholesky (which would fail or give garbage).

**Suggested test:**
- Run many iterations with small `n` and small `nue`; check that `E(t,t)` is always finite and positive after sampling.

---

### 2.3 Potential division by zero in x2 marker-ranking computation

**Severity: Medium (NaN propagates into sort order).**

`mt_cpg.cpp:255` and `mt_cpg_arma.cpp:249`:

```cpp
x2[t][i] = (wy[t][i] / ww[t][i]) * (wy[t][i] / ww[t][i]);
```

If `ww[t][i] == 0` (a monomorphic marker or one with zero variance), this produces `Inf` or `NaN`. `std::sort` with `NaN` comparisons is undefined behavior.

**Suggested test:**
- Pass a dataset where one marker has `ww == 0` for at least one trait; check whether `sort` hangs or crashes.

---

## 3. Numerical and Memory-Safety Issues

### 3.1 `indptr` integer overflow in `csr_ld.cpp`

**File: `csr_ld.cpp:89`.**

```cpp
csr.indptr[i + 1] = static_cast<int>(nnz);
```

`nnz` is `size_t`. For a genome-wide LD matrix with millions of markers and large windows, `nnz` can easily exceed `INT_MAX` (~2.1 billion), causing silent wraparound. The old CSR struct uses `std::vector<int>` for `indptr`, which cannot represent values > 2^31. The newer `STLDCSR` in `st_csr_common.h` correctly uses `uint64_t` throughout.

**Suggested test:**
- Construct a case where total LD entries exceed 2^31 (or mock `nnz` with a large value) and verify `indptr` values are not negative.

---

### 3.2 Diagonal excluded when `ld_threshold >= 1.0` in `csr_ld.cpp`

**File: `csr_ld.cpp:71, 77`.**

```cpp
if (k == msize) ld = 1.0f;           // force diagonal to 1
if (std::abs(ld) > ld_threshold)     // strict greater-than
```

When `ld_threshold == 1.0` (or higher), `std::abs(1.0) > 1.0` is `false`, so the diagonal entry is **excluded** even though it was forced to 1. Any sampler relying on this CSR having a diagonal would compute wrong residuals.

**Suggested test:**
- Call `readLD_to_CSR_R(threshold = 1.0)` and verify the returned matrix has a diagonal entry for each marker.

---

### 3.3 Duplicate LD entries produce corrupted symmetric CSR

**File: `st_csr_common.h:148-243`.**

`read_and_build_st_ld_csr` builds a symmetric CSR by expanding a one-sided LD file. It counts degrees by iterating all non-diagonal pairs and adding 1 to both `degree[i]` and `degree[j]`. If the input file contains **both** `(i,j)` and `(j,i)`, each pair is counted and inserted **twice**, doubling the off-diagonal X'X values silently. No check verifies the input is strictly one-sided (upper or lower triangular).

**Suggested test:**
- Construct a binary LD file containing both `(i,j)` and `(j,i)` entries and pass it to a sampler; compare the resulting posteriors to one using a correct one-sided file.

---

### 3.4 Missing 64-byte alignment on MinGW/Rtools Windows

**File: `packed_bed.h:27-30`.**

```cpp
#elif defined(_WIN32)
    return std::malloc(nbytes_aligned);   // NOT 64-byte aligned
```

On Windows with non-MSVC compilers (Rtools uses GCC/MinGW), `_aligned_malloc` is unavailable and `std::malloc` is used as a fallback. `std::malloc` guarantees only fundamental alignment (typically 8 or 16 bytes), not the 64 bytes needed for AVX-512 SIMD. If the caller assumes 64-byte alignment for vectorized operations, this is undefined behavior. The `posix_memalign` path (Linux/macOS) and the `_aligned_malloc` path (MSVC) are both correct.

**Suggested test:**
- On Windows with Rtools, run `bed_pairwise_xtx_check` with a large BED file and compare results to the Linux build.

---

## 4. R/C++ Interface Issues

### 4.1 `cls` indices stored 1-based in `rcpp_to_int_list`; subtracted at use site in `packed_bed.h`

**Files: `packed_bed_matrix.cpp:66`, `packed_bed.h:171`.**

`rcpp_to_int_list` validates `x[k] > 0` and stores the raw 1-based R integer. The conversion to 0-based happens later at:

```cpp
const long long offset = 3LL + static_cast<long long>(cls - 1) * static_cast<long long>(nbytes_bed);
```

This is consistent and correct, but the two-step convention (store 1-based, convert at use) is fragile. If any new code path uses the stored index without subtracting 1, it will silently read the wrong marker column.

---

### 4.2 Genetic covariance normalization inconsistency in `computeG_cpg_arma`

**File: `cpg_samplers.cpp:92`.**

Diagonal entries use `n[t]` as divisor: `G(t,t) = b_t' XX b_t / n[t]`.
Off-diagonal entries use `sqrt(n[t1] * n[t2])`: `G(t1,t2) = 0.5*(b_t1' XX b_t2 + b_t2' XX b_t1) / sqrt(n[t1]*n[t2])`.

When `n[t1] != n[t2]`, the off-diagonal is in different units from the diagonal, so `G` is not a valid covariance matrix and cannot be used for genetic correlation computation (`G[t1,t2] / sqrt(G[t1,t1] * G[t2,t2])`). Same pattern in `computeG_cpg` (`cpg_samplers.cpp:275`).

**Suggested test:**
- Run with two traits of identical effects but different n; verify that the computed genetic correlation (`cvgm[1][2] / sqrt(cvgm[1][1] * cvgm[2][2])`) equals 1.0.

---

## 5. Summary Table

| # | File | Location | Severity | Category |
|---|---|---|---|---|
| 1.1 | `mt_cpg.cpp`, `mt_cpg_arma.cpp` | lines 313/362, 302/355 | Critical | Double B sampling, stale Bi |
| 1.2 | `mt_cpg.cpp`, `mt_cpg_arma.cpp` | multiple `it > nburn` | High | Off-by-one in burn-in |
| 1.3 | `mt_cpg.cpp`, `mt_cpg_arma.cpp` | result assembly | High | Division by zero (nsamples) |
| 2.1 | `mt_cpg.cpp`, `mt_cpg_arma.cpp` | lines 152, 110 | Medium | Unchecked Cholesky return |
| 2.2 | `cpg_samplers.cpp` | lines 67, 236 | Medium | Unguarded chi2 division |
| 2.3 | `mt_cpg.cpp`, `mt_cpg_arma.cpp` | lines 255, 249 | Medium | Division by zero in sort key |
| 3.1 | `csr_ld.cpp` | line 89 | Medium | int overflow in indptr |
| 3.2 | `csr_ld.cpp` | line 77 | Medium | Diagonal excluded at threshold=1 |
| 3.3 | `st_csr_common.h` | lines 148-243 | Medium | Duplicate LD entries corrupt CSR |
| 3.4 | `packed_bed.h` | lines 27-30 | Low | No 64-byte alignment on MinGW |
| 4.1 | `packed_bed_matrix.cpp`, `packed_bed.h` | lines 66, 171 | Low | Two-step 1→0 base convention |
| 4.2 | `cpg_samplers.cpp` | lines 92, 275 | Low | G normalization mismatch |

---

## 6. Verification Approach (no tests written)

These tests are described only. Do not implement them unless explicitly requested.

1. **Burn-in boundary**: `nit=1, nburn=0, nthin=1` — with the current code `nsamples` should be 0; result `bm[t][i]` will be NaN.
2. **Double B sampling**: compare posteriors from a single large set vs the same run split into two sets with `updateB=TRUE`. They should agree but won't due to bug 1.1.
3. **nsamples guard**: assert `nsamples > 0` before division.
4. **Chi2 floor**: run 100,000 E-sampling steps with `n=2, nue=1`; check no `E(t,t)` is infinite.
5. **x2 stability**: pass `ww[t][i]=0` for one marker; check no crash or NaN propagation.
6. **CSR indptr overflow**: mock `nnz = 2^31 + 1` and verify overflow behavior of `static_cast<int>`.
7. **LD threshold diagonal**: `readLD_to_CSR_R(threshold=1.0)` should include diagonal; currently does not.
8. **Genetic correlation**: two identical-effect traits with different n should give genetic correlation 1.0.
