# BLR Framework Phase 17D Report

## 1. Executive summary

The corrected authoritative public `mtblr()` single-fit numerical setup and
Gibbs execution were mechanically moved into one guarded lexical implementation
header. Phase 17C numerical behavior, RNG order, retained-count policy, legacy
20-position finalization, and public R formatting remain unchanged. No typed
production migration, aggregation extraction, or converter extraction began.

## 2. Repository baseline

Phase 17D began on `master` at clean commit `48d9498` (`Correct public
multivariate update and retention semantics`), which is also the Phase 17C
commit. No combined hosted Actions status was visible locally, so hosted CI is
not claimed. The Windows Rtools44 toolchain compiled and loaded the package.
The Phase 17C baseline was 6,030 passes, zero failures, zero warnings, and ten
opt-in skips. Its matched benchmark baseline is recorded in section 21.

## 3. Public call graph

`sblr(algorithm = "default")` calls `mtblr()`, which retains its R validation,
dense summary preparation, hidden R-derived seed, generated `.Call` wrapper,
and native signature. Native `mtblr()` includes the one lexical execution block,
then constructs the legacy 20-position result inline. `R/interface_mtblr.R`
continues to name, orient, and conditionally add correlations to the public fit.

## 4. Extraction target

The complete single-fit numerical block was selected because it is the smallest
existing boundary that preserves initialization, the one Gibbs loop, fit-local
RNG ownership, helper call order, mutable numerical state, and the corrected
retention contract without inventing production types.

## 5. Exact extraction seam

The moved block begins at `// Define local variables`, immediately after the
opening brace of public native `mtblr()`. It ends after the
`marker_retained_count <= 0.0` guard and immediately before `// Summarize
results`. That finalization comment, the 20-position allocation and filling,
and the return remain in `src/mtblr.cpp`; `mtblr_hybrid()` still follows it.

## 6. Lexical dependency inventory

Function inputs used lexically are `wy`, `ww`, `yy`, `b`, `XXvalues`,
`XXindices`, `sets`, `B`, `E`, `ssb_prior`, `sse_prior`, `models`, `pi`, `nub`,
`nue`, `updateB`, `updateE`, `updatePi`, `n`, `nit`, `nburn`, `nthin`, `seed`,
and `method`. Existing translation-unit helpers remain `sampleBset()`,
`sampleB_latent()`, `sampleBetaCPG_Mt_latent()`, `samplePi()`, `sampleB()`,
`computeG()`, and `sampleE()`.

Values produced for inline finalization include `nt`, `m`, `nmodels`, the five
retained counters, `bm`, `dm`, `wy`, `r`, `b`, `d`, `order`, `vbs`, `vgs`,
`ves`, `cvbm`, `cvgm`, `cvem`, `B`, `G`, `E`, `pi`, `pis`, `pistrait`, and
`pismarker`. The header relies on the translation unit's standard-library and
Armadillo includes and owns no binding metadata.

## 7. Files changed

- `src/mtblr.cpp`: replaced the moved numerical block with one lexical include.
- `src/blr_mt_default_core_impl.h`: added the guarded mechanical block.
- Phase 17B/17C/17D tests: redirected structural checks, preserved historical
  assertions, and added extraction, reference, reproducibility, and protection
  coverage.
- Fast and extended workflows: added ordinary and opt-in fresh Phase 17D tests.
- Phase 17D benchmark: reused the matched Phase 17C workloads and convention.
- Implementation plan, capability matrix, and this report: recorded the
  noncanonical lexical boundary and validation evidence.

No production R file, generated wrapper, `NAMESPACE`, fixture, or alternative
backend source changed.

## 8. Lines mechanically moved

```text
MECHANICAL_LINES=240
HEADER_BODY_LINES=240
IDENTICAL=TRUE
```

After excluding the include guard and implementation-detail comment, the header
body is text-identical to the Phase 17C block. Its normalized MD5 is
`7e8ea9e4812ce57a701416f8896a97cc`. The unchanged inline finalization region's
normalized MD5 is `01e41f91932d420df012cfc2e9b7b20d`.

## 9. Extracted numerical content

The header contains dimensions, five contribution counters, all fit-local
state and workspaces, residual and covariance initialization, marker scoring
and descending update order, one fit-local engine, the complete Gibbs loop,
set/state/effect/residual/covariance/probability updates, posterior
accumulation, relative thinning, counter increments, and the retained-marker
guard. It contains no result construction, R formatting, file I/O, export,
worker identity, Python, or alternative path.

## 10. Correctness-contract preservation

Set-local `sampleBset()` and `sampleB_latent()` remain jointly guarded by
`updateB`; the global `sampleB()` remains guarded by `updateB && method == 4`.
Post-burn-in remains `it >= nburn`, marker thinning remains
`(it - nburn) % nthin == 0`, and marker, B, G, E, and probability summaries
retain their accumulator-specific denominators. Disabled B/E/probability
summary conventions and the positive marker-retained guard are unchanged.

## 11. RNG preservation

Exactly one `std::mt19937 gen(seed)` remains in the public path at the same
location relative to initialization and helper calls. All draw sites, helper
order, marker and set traversal, and covariance/state sampling order are
unchanged. There is no public worker-ID contribution, static stochastic state,
or `thread_local` stochastic state. The unused historical `random_device`
declaration remains mechanically unchanged.

## 12. Marker/state/covariance preservation

Input trait order, output marker order, descending marker-update permutation,
explicit set traversal, joint trait-pattern state order, latent and effective
effect handling, residual changes, and B/G/E calculation and update order are
identical to Phase 17C. No vector expression, allocation, or formula was
rewritten.

## 13. Legacy finalization

The 20-position allocation, resizing, marker means, trace copies, covariance
means, final B/G/E, final and mean probabilities, legacy positions 19-20, and
return remain inline and unchanged in `src/mtblr.cpp`. No typed result or
aggregation object was introduced.

## 14. R interface and schema

`R/interface_mtblr.R`, generated R/C++ wrappers, `NAMESPACE`, native signature,
hidden seed behavior, algorithm route, positional schema, names, orientations,
conditional `rb`/`rg`/`re`, and omission semantics are byte-identical to the
Phase 17C baseline.

## 15. Phase 17B historical fixtures

The three RDS files remain byte-identical with SHA-256 values:

- config 1: `82AF2F814C48BC6E5E4B8D7F748DC26BB7E2BC7F058D0D784D8F4508FDC874C9`
- config 2: `CF32B7C68BA943855BE131E0E0BB4C24AB76749DA6C8747A78870E601EED869C`
- config 3: `589E3A634BAB335C06E6F1D47063F28A0ACFE3B1C62F1EBC3A57A52DED2213D9`

## 16. Phase 17C corrected references

All 3/3 corrected raw and 3/3 corrected formatted references pass with exact
structure and numerical tolerance `1e-12`. The committed fixture SHA-256 values
remain `2AD954...C806`, `4F1929...F3A`, and `8F7BD2...3001`; no fixture was
regenerated or modified.

## 17. Reproducibility

`A;A`, `A;B;A`, OMP thread environments 1 versus 2, intervening canonical CSR,
intervening canonical packed-BED, and fresh versus reused-process comparisons
remain exact in structure and equal within `1e-12`. The Phase 17D opt-in fresh
process comparison against corrected formatted configuration 1 passes.

## 18. Scientific identities

Finite-value, marker/trait dimension and ordering, binary state, bounded `dm`,
symmetric positive-semidefinite covariance, trace length and orientation,
fixed B/E/pi, retained-count, disabled-summary, and `sum(pim) = 1` identities
all pass.

## 19. Alternative multivariate implementations

`mt_cpg.cpp`, `mt_cpg_arma.cpp`, `mt_cpg_omp.cpp`, and `mt_cpg_omp_csr.cpp` are
unchanged. Normalized protected regions for `mtblr_hybrid()` and
`mtblr_eigen()` remain `03ed109628b874223db109d2ec654827` and
`4e5c38ede3345de10a684ab38470bf7b`. No alternative includes the new header;
the `mtblr_cpg_omp()` worker-sensitive RNG risk test remains active.

## 20. Scalar, packed-BED, and block-eigen protection

Canonical scalar CSR, canonical packed-BED, shared packed-BED family, and
block-eigen numerical source hashes and permanent references remain unchanged.
Generated wrappers and `NAMESPACE` are unchanged.

## 21. Performance, memory, and I/O

The Phase 17C matched baseline was: updated mean 0.048 s/median 0.01/RSS 120.58
MiB; fixed-B 0.010/0.01/120.63; explicit sets 0.012/0.01/120.63; moderate dense
0.052/0.05/121.52. A representative Phase 17D run produced updated
`0.47,0.01,0.01,0.03,0.03` s (mean 0.110, median 0.03, RSS 127.98 MiB), fixed-B
`0.03,0.03,0.02,0.03,0.02` (0.026, 0.03, 122.71), explicit sets
`0.04,0.03,0.03,0.03,0.03` (0.032, 0.03, 119.51), and moderate dense
`0.20,0.22,0.19,0.20,0.19` (0.200, 0.20, 120.72).

These short Windows debug-build timings are regression signals, not speed
claims; the identical extracted body provides no source-level change in work or
allocation, while system/toolchain load makes the observed timing comparison
noisy. Completed-fit RSS is not peak RSS; peak RSS was not sampled. Dense `XX`
remains `O(nt * m^2)`. All inputs remain in memory and there is no MCMC-time I/O.

## 22. CI coverage

The fast workflow filter includes ordinary Phase 17D tests. The extended
workflow retains all existing opt-in variables and adds
`SBLR_RUN_PHASE17D_FRESH: "true"`. Benchmarks are not part of the fast gate.
Hosted CI success is not claimed.

## 23. Tests

- Phase 17B historical tests: passed with the normal opt-in skip.
- Phase 17C corrected tests: passed with the normal opt-in skip.
- Phase 17D ordinary: 90 passed and one opt-in fresh-process skip.
- Phase 17D fresh-process enabled: 91 passed, zero failures.
- Corrected raw/formatted references: 3/3 and 3/3.
- Alternative and canonical protection tests: passed.
- Full suite: 6,120 passed, zero failures, zero warnings, and eleven opt-in skips.

## 24. Deviations and blockers

One full-suite invocation reached the 15-minute command harness timeout without
a test or compiler diagnostic; a clean rerun with a longer allowance produced
the result above. Tiny benchmark timings varied materially under the Windows
debug/load environment, but the executed 240-line body is mechanically
identical and completed-fit RSS and I/O behavior show no systematic regression.
No numerical, reference, schema, ownership, or extraction blocker remains.

## 25. Recommended next phase

Replace the lexically dependent corrected `mtblr()` execution include with
explicit typed multivariate data, model, covariance-prior, execution-context,
and result contracts plus one callable binding-neutral MT BayesC core, while
retaining legacy finalization and R formatting until all Phase 17C references
pass again.

## 26. Readiness marker

PHASE 17D COMPLETE — CORRECTED PUBLIC MULTIVARIATE EXECUTION BLOCK EXTRACTED
