# Internal individual-level multivariate packed-BED contract

## 1. Purpose

Phase 17N established this implementation contract. Phase 17O now implements
the internal, individual-level, serial multivariate BayesC route without
changing its physical genotype, likelihood, state, covariance, ownership,
memory, output, or reduction decisions.

## 2. Public status

Phase 17P exposes `mtblr_bed()` as the supported public adapter. It prepares
and validates public Glist/phenotype inputs and invokes the unchanged Phase 17O
`mtblr_bed_internal()` exactly once. The deterministic
`mtblr_bed_marker_contract_internal()` remains internal.

## 3. Existing scalar BED family

`stblr_bed(method="bayesc")`, `"bayesr"`, and `"bayesrc"` accept phenotype
matrices but execute one scalar likelihood for each trait and chain. They do not
fit a joint multivariate likelihood and their scalar marker conditionals are
not the basis of the future MT conditional.

| family | owner | reader | traversal | statistical status |
|---|---|---|---|---|
| scheduled BayesC | `FastPackedBedMatrix` | `read_bedfiles_to_fast_packed_matrix_blocked()` | packed lookup plus `BedPackedGenotypeView` | trait-specific |
| BayesR | `FastPackedBedMatrixBR` | `br_read_bed_blocked()` | shared BayesR packed utilities plus `BedPackedGenotypeView` | trait-specific |
| BayesRC | `FastPackedBedMatrixBR` | `br_read_bed_blocked()` | BayesR utilities plus `BedBayesRCPackedGenotypeView` | trait-specific |

The readers run before chain execution. The chain cores borrow immutable
storage; they do not open BED files during MCMC.

## 4. Existing summary MT family

The dense `sblr()` route, `mtblr_csr()`, and `mtblr_block_eigen()` share the
corrected Phase 17C/17I joint-pattern model, latent/effective effects, update
order, result types, finalization, `mtblr_raw` version 1, and fit formatter.
Phase 17O must reuse those statistical and result contracts. It must not use
legacy `mtblr_eigen()`.

## 5. Packed owners

| owner | definition/reader | allocation and alignment | stride | copy/move | consumers and lifetime |
|---|---|---|---|---|---|
| `PackedBedMatrix` | `packed_bed.h`; `read_bedfiles_to_packed_matrix()` | raw pointer, allocation rounded to 64 bytes; POSIX uses `posix_memalign`, MSVC uses `_aligned_malloc`, MinGW `_WIN32` currently uses `malloc` and therefore does not prove 64-byte address alignment | `round_up(ceil(n/4),64)` | copy deleted, move enabled | block-eigen and packed inspection; owner spans construction/use |
| `FastPackedBedMatrix` | scheduled BayesC binding; blocked reader | `std::vector<uint8_t>` zero-initialized | `ceil(n/4)` | vector value owner, normally moved/retained once | all BayesC traits/chains borrow one owner |
| `FastPackedBedMatrixBR` | `st_bed_bayesr_common.h`; `br_read_bed_blocked()` | `std::vector<uint8_t>` zero-initialized | `ceil(n/4)` | vector value owner, normally moved/retained once | BayesR and BayesRC traits/chains borrow one owner |

All are SNP-major. All preserve BED-file order, `cls` order within file, and
selected-row order. The blocked vector readers read markers in blocks but seek
to every requested marker. `PackedBedMatrix` reads one requested marker at a
time. File handles are closed before returning the owner.

The owners are logically equivalent for selected samples and markers, but they
are not representation-identical: stride, address alignment, padding beyond
`bytes_per_marker`, allocation, and exception machinery differ. Phase 17N does
not consolidate them.

Complete implementation inventory:

| owner type | reader/source | consumer | allocation/alignment | stride/bytes per marker | rows/markers | move/copy | borrowed view | lookup/standardization/missing | `x'x` / `x'e` | residual update/rebuild | threading and I/O | lifetime/tests |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `PackedBedMatrix` | `read_bedfiles_to_packed_matrix()` in `packed_bed.h` | block-eigen builder and packed inspection | raw pointer; 64-byte allocation intent with documented MinGW caveat | rounded-to-64 stride / `ceil(n/4)` logical bytes | supplied row and file/`cls` order | move-only | direct owner today; common view selected for Phase 17O | `get_bed_code`; caller maps; missing excluded from AF and mapped by caller | packed block construction and inspection helpers | block-eigen marker-space rebuild, no scalar sample residual | optional OpenMP AF/build work; files opened/read/closed before use | one construction owner; Phase 17K and packed-matrix tests |
| `FastPackedBedMatrix` | blocked reader in scheduled BayesC binding | BayesC trait/chain cores and aggregate | `vector<uint8_t>`, no address-alignment promise | compact `ceil(n/4)` / same | supplied row and file/`cls` order | vector value type but one fit owner, borrowed by chains | `BedPackedGenotypeView` | `fast_get_bed_code_from_row`; four-value map; missing zero or `2p` | `marker_xx_from_packed_scheduled_chains()` / `bed_marker_dot_residual_scheduled_chains()` | packed rank-one update; rebuild through `bed_xb_from_b_scheduled_chains()` | blocked pre-fit I/O; AF/maps may use OpenMP; chain threading outside view | owner spans all chains; Phase 11/15 and canonical BayesC references |
| `FastPackedBedMatrixBR` | `br_read_bed_blocked()` in `st_bed_bayesr_common.h` | BayesR and BayesRC trait/chain cores | `vector<uint8_t>`, no address-alignment promise | compact `ceil(n/4)` / same | supplied row and file/`cls` order | vector value type but one fit owner, borrowed by chains | common view for BayesR; narrower BayesRC view | `br_get_bed_code`; four-value map; missing zero or `2p` | `br_marker_xx()` / `br_dot_residual()` | `br_update_residual()`; rebuild through `br_xb()` | blocked pre-fit I/O; AF/maps may use OpenMP; chain threading outside view | owner spans all chains; Phase 13/14/15 and BayesR/BayesRC references |

## 6. Packed views

`BedPackedGenotypeView<PackedGenotype>` contains an immutable storage
reference, packed pointer/size, marker/sample counts, bytes per marker, and
stride. BayesC and BayesR use it. BayesRC uses the narrower
`BedBayesRCPackedGenotypeView`, which retains the owner reference and dimensions
but obtains row pointers through the owner.

Phase 17O reuses `BedPackedGenotypeView` unchanged. Its numerical code shall use
`packed_markers + marker*stride`; the `storage` reference is a lifetime anchor
and validation aid, not a required hot-loop dependency. The view is immutable,
non-owning, Rcpp-free, SEXP-free, file-handle-free, RNG-free, and
sampler-state-free.

## 7. Physical BED format

Only PLINK BED with header bytes `0x6c 0x1b 0x01` is accepted. The third byte
requires SNP-major mode. A marker occupies `ceil(n_bed/4)` physical bytes, four
samples per byte, two bits per sample, least-significant pair first:

| bits/code | logical dosage |
|---|---:|
| `00` / 0 | 2 |
| `01` / 1 | missing |
| `10` / 2 | 1 |
| `11` / 3 | 0 |

In the final partial byte only the first `n mod 4` pairs are samples. Source BED
padding bits are not required to be zero. Selected-row repacking starts from a
zero buffer, so its unused pairs are zero. `PackedBedMatrix` also zeroes stride
padding. Every active loop checks `sample<n` (or `jbase+lane<n`), so neither
partial-byte nor stride padding is treated as a sample.

## 8. Marker and row selection

BED files are traversed in the supplied order. Requested `cls` entries are
one-based physical marker columns and are emitted in their supplied order.
Selected rows are one-based at the R boundary, converted once to zero-based,
and emitted in supplied order; they are never sorted. With no row selection,
physical sample order is retained.

## 9. Genotype standardization

For supplied or computed frequency `p_j`, `scale=TRUE` maps observed dosage
`g_ij` to

```text
x_ij = (g_ij - 2 p_j) / sqrt(2 p_j (1-p_j)).
```

Phase 17O requires finite `p_j` strictly inside `(0,1)` and supports only this
standardized scale. It does not silently zero an invalid/monomorphic marker.

For the existing `scale=FALSE` scalar paths the map is dosage 2 to 2, dosage 1
to 1, dosage 0 to 0, and missing to `2p_j`.

Frequencies are the frequency of the dosage counted by BED code 0 (the
two-copy coded allele). Repository code does not establish that this is always
the minor, effect, reference, or alternate allele, so no stronger biological
name is assigned. Scalar BayesC/BayesR/BayesRC use supplied `Glist$af` after
marker selection when present; otherwise they calculate
`sum(nonmissing dosage)/(2*n_nonmissing)` from the selected packed rows.
`PackedBedMatrix` exposes the same calculation.

## 10. Missing genotype handling

With `scale=TRUE`, missing maps to zero: mean imputation after centering. With
`scale=FALSE`, missing maps to `2p_j`: mean dosage imputation. Frequency
calculation excludes missing genotypes. Phase 17O uses standardized
missing-to-zero behavior only.

## 11. Owner recommendation

Phase 17O selects option 1: use the existing move-only `PackedBedMatrix` and
`read_bedfiles_to_packed_matrix()` as the one fit-owned genotype owner, then
construct the existing `BedPackedGenotypeView`. This avoids touching the
protected scalar readers and avoids a fourth owner. The cost is non-blocked
per-marker reading and 64-byte-rounded strides; that is acceptable for the
first internal serial route and must be benchmarked before a public phase.

The MinGW address-alignment caveat above means Phase 17O must not issue aligned
SIMD loads merely from the stride contract. It may use unaligned loads or
verify address alignment separately.

## 12. Phenotype contract

Phase 17O accepts an owned `arma::mat Y` with `n>1`, `nt>0`, finite values,
unique nonempty column names at the R adapter, and exactly the shared selected
BED rows in exactly the same order for every trait. There is one genotype view.
Phenotypes are pre-adjusted for covariates and each column must be centered
before native execution.

Column scaling is not required and is not performed. This preserves phenotype
units in `B`, `E`, and genetic covariance while permitting users to scale
externally.

Current `.make_bed_marker_data()` converts vectors to one named trait, uses
matrix columns as traits, matches phenotype row names to unique `Glist$ids`,
warns and drops unmatched phenotype IDs, rejects duplicate phenotype IDs,
duplicate Glist IDs, duplicate selected rows, invalid rows, and empty matches,
and preserves phenotype/row order. Without names it requires `nrow(y)==Glist$n`.
It currently supplies missing trait names as `T1...`, but does not reject
duplicate trait names, nonfinite values, constant columns, uncentered values,
zero-column matrices, or all nonfinite columns explicitly; downstream failure
may be indirect. Those are not supported claims.

## 13. Covariates

`stblr_bed(covar=...)` stops and tells users to pass pre-adjusted phenotypes.
Phase 17O likewise has no covariate argument and requires pre-adjustment.

Before Phase 17P can expose covariates, it must specify design alignment,
rank/alias handling, intercept and centering policy, whether projection is
fixed or effects are sampled jointly, uncertainty propagation, missing
covariates, priors, outputs, and prediction semantics.

## 14. Missing phenotypes

Phase 17O supports a complete finite phenotype matrix only. It does not use
complete-case deletion internally, trait masks, pattern-specific samples, or
latent phenotype imputation.

Trait-specific missingness would change marker `x'x`, scores, effective sample
sizes, cross-trait likelihood terms, residual cross-products, `B` and `E`
updates, and predictive diagnostics. It therefore requires a later explicit
masked-likelihood contract rather than an adapter convenience.

## 15. Statistical model

For standardized `X` (`n x m`) and centered `Y` (`n x nt`):

```text
Y = X B_eff + R,
R_i ~ N_nt(0,E),
beta_j ~ N_nt(0,B),
D_j = diag(model-pattern_j),
b_j = D_j beta_j.
```

`beta_j` is the full latent effect and `b_j` is the effective marker effect.
The model prior is the existing unique binary pattern matrix with a null
pattern and probabilities `pi`. Sets are the existing disjoint complete
partition. Marker effects and `B` follow the corrected Phase 17C semantics.

## 16. Sample-space residual

The canonical state is column-major `arma::mat R(n,nt)`:

```text
R = Y - X B_eff.
```

It is always rebuilt from `Y`, `X`, and initial effective effects. No external
sample residual is accepted. For marker `j`:

```text
w_j = x_j'x_j
s_j = x_j'R + w_j b_j
R <- R - x_j (b_j,new - b_j,old)'.
```

`arma::mat` gives contiguous trait columns for scores, direct per-sample
cross-trait rows for full-`E` work, BLAS-compatible matrix products, and
unambiguous Rcpp conversion. A vector-of-vectors layout is rejected.

## 17. Full-E marker conditional

Let `Omega=E^{-1}`, `P=B^{-1}`, and for pattern `k` let `D_k` be its diagonal
mask. Removing the current marker contribution gives the likelihood kernel

```text
exp(beta_j' D_k Omega s_j
    - 0.5 w_j beta_j' D_k Omega D_k beta_j).
```

Combining it with the full latent prior gives

```text
C_k   = P + w_j D_k Omega D_k
rhs_k = D_k Omega s_j
V_k   = C_k^{-1}
m_k   = V_k rhs_k
```

and, up to the common model-independent term `0.5 log|P|`,

```text
log q_k = log(pi_k)
          - 0.5 log|C_k|
          + 0.5 rhs_k' C_k^{-1} rhs_k.
```

The omitted term is common because every pattern retains the same
`nt`-dimensional latent prior. For the null pattern, `D=0`, `C=P`, `rhs=0`;
after restoring the common term its log weight is exactly `log(pi_null)`.

After categorical selection:

```text
beta_j,new = m_k + L_k^{-T} z,  z ~ N(0,I),  L_k L_k' = C_k
b_j,new    = D_k beta_j,new
R          = R - x_j (b_j,new-b_j,old)'.
```

Inactive latent coordinates remain sampled. Because `B` may correlate latent
coordinates, their posterior can change when active coordinates are observed;
they never enter `b_j` or the residual directly. They do enter the established
latent `B` update. SPD failure after the established bounded Cholesky safeguard
is an error.

## 18. Diagonal-E reduction

For `E=diag(e_1,...,e_nt)`, `Omega` is diagonal and

```text
C_k = P + diag(w_j d_kt/e_t)
rhs_kt = d_kt s_jt/e_t.
```

This is the active corrected summary-MT marker conditional. Exact numerical
identity requires the same standardized `X`, `X'X`, `X'Y`, marker-space
residual state `X'R`, `B`, diagonal `E`, patterns, probabilities, sets, marker
order, Cholesky safeguards, and RNG state. The shared individual design also
requires one common `w_j`; summary views must expose that same diagonal for all
traits.

This is the permanent Phase 17O oracle. Diagonal `E` is a deterministic
reduction mode, not the canonical same-individual likelihood.

## 19. B update

Phase 17O preserves the corrected Phase 17C update control and order:

1. when `updateB`, each set invokes `sampleBset()` on active effective effects;
2. `sampleB_latent()` then draws
   `B ~ IW(nub+m, ssb_prior + sum_j beta_j beta_j')`, including all latent
   coordinates, including inactive ones;
3. marker updates for that set use the resulting `B^{-1}`;
4. after all marker updates, the established global `sampleB()` updates
   trait diagonals from active effective effects and reconstructs shrunk
   shared-marker correlations with its existing SPD floor.

When `updateB=FALSE`, none of these execute and supplied `B` remains exact.
This contract is intentionally the corrected shared MT behavior, not the
scalar BED variance update. Phase 17O must call shared helpers or a mechanically
identical extracted boundary; it must not reinterpret the update.

## 20. E update

Diagonal reduction mode uses, for each trait,

```text
SSE_t = sum_i R_it^2
E_tt = (SSE_t + nue*sse_prior_tt) / ChiSq(n+nue),
E_ts = 0 for t != s.
```

Full canonical mode uses

```text
S_E = 0.5*(R'R + (R'R)')
S_post = nue*sse_prior + S_E
E | R ~ IW(nue+n, S_post).
```

`sse_prior` must be symmetric positive definite, `nue>nt-1`, and
`S_post` must pass a symmetric-positive-definite Cholesky check. Symmetrization
addresses roundoff; invalid inputs are rejected. A narrowly bounded diagonal
jitter may be attempted only under the same documented deterministic policy as
the marker conditional and must be reported diagnostically.

Existing `rwishart()`/`rinvwishart()` use `std::mt19937`, Armadillo, SPD
inversion, and Bartlett draws and are mathematically suitable. They currently
live in `mtblr.cpp`; Phase 17O should place a binding-neutral, single
implementation in an internal header with reduction tests before use.

Full `E` is canonical for complete same-individual data. Diagonal `E` is the
test/reduction mode. A future public interface should default to full `E`.

## 21. Genetic covariance

Genetic values are

```text
U = X B_eff
G = U'U / n.
```

This matches scalar BED `sum(g_i^2)/n` on the diagonal and the same-design
summary identity `B_eff'X'X B_eff/n`. No `n-1` denominator or additional
centering is introduced; centered phenotypes and the standardized genotype
contract define the working scale. `G` is symmetrized after multiplication.

## 22. LE/LD components

Phase 17O initially retains only `vgs`, `ves`, and `vbs` plus their covariance
matrices. It does not report `vle`/`vld`. Traitwise scalar LE formulas do not
define cross-trait LE covariance, and inventing one would make `vld=G-LE`
ambiguous. A later phase may add a separately derived covariance decomposition.

## 23. Marker order and sets

Sets are an explicit disjoint complete partition. Every iteration performs a
full deterministic sweep; there is no adaptive null scheduling. Phase 17O uses
the corrected summary-MT marginal-score order computed once from `X'Y` and
`x'x`, with input marker index as the deterministic tie break. For each set it
scans that stable global order and updates members of the set. This preserves
the Phase 17C RNG/reduction contract.

## 24. Marker decoding strategy

Phase 17O decodes one marker once per visit into one fit-owned reusable
`arma::vec`/double workspace of length `n`. Scores for all traits, full-`E`
algebra, and residual updates reuse it. The workspace is allocated once in core
state, overwritten at each marker, never shared across fits, and maps missing
standardized genotypes to zero.

Direct packed traversal is competitive for `nt=1` but repeats lookup for each
trait and residual pass. The double workspace costs `8n` bytes and gives stable
operation order and BLAS-compatible precision for small/moderate `nt`.

## 25. RNG

Execution is serial and one-chain. One fit-local `std::mt19937` is constructed
from the explicit seed after BED reading, phenotype validation, marker-map
construction, and workspace allocation, immediately before core-state
initialization. The core owns the engine. No OpenMP, worker seed, R RNG, static
distribution state, or global cache is permitted.

A later multichain executor must give each logical chain a separate owner of
state, residual, workspace, and RNG while allowing immutable genotype sharing;
that executor is outside Phase 17O.

## 26. Ownership

```text
R/Rcpp adapter
  owns converted inputs, names, metadata
  owns one PackedBedMatrix
    borrowed by one BedPackedGenotypeView
  transfers/owns Y for the fit
core state
  owns R, decoded marker workspace, effective effects, latent effects,
  inclusion states, B, E, pi, marker maps/order, and one RNG
result
  owns retained summaries, traces, final states, and final covariances
binding result conversion
  owns raw names and metadata
```

There is one genotype owner, zero per-marker or per-chain genotype copies, zero
MCMC-time BED reads or file handles, zero mutable statics, and zero global
caches.

## 27. Memory

Let `q=ceil(n/4)`, `stride=64*ceil(q/64)`, `K` be pattern count, and `L` the
retained/trace length:

| object | analytical bytes |
|---|---:|
| packed owner | `m*stride` |
| phenotype `Y` | `8*n*nt` |
| sample residual `R` | `8*n*nt` |
| effective effects | `8*m*nt` |
| latent effects | `8*m*nt` |
| inclusion state | `4*m*nt` |
| decoded marker | `8*n` |
| four-code map plus `x'x` | `40*m` |
| marker order plus set label | `8*m` |
| six covariance work matrices | `48*nt^2` |
| per-model dense work upper bound | `8*K*(nt^2+2nt+2)` |
| minimum `B/E/G/pi` traces | `8*L*(3nt+K)` |

Transient BED reading additionally uses `ceil(n_bed/4)` for the unblocked owner
or block-sized buffers in the existing blocked readers. The analytical total
is not measured RSS or peak RSS. Completed-fit RSS includes R, allocator,
Armadillo, Rcpp, package, and output overhead; peak RSS also includes
construction and conversion temporaries.

## 28. Runtime complexity

With pattern active count `a_k`:

| operation | cost |
|---|---:|
| BED construction | `O(m*ceil(n_bed/4) + mn)` with selected-row repacking |
| frequency calculation | `O(mn)` |
| marker maps / `x'x` | `O(mn)` |
| initial `XB_eff` | `O(mn*nt)` |
| one marker decode | `O(n)` |
| one marker score | `O(n*nt)` |
| one residual update | `O(n*nt)` |
| model weights | dense contract `O(K*nt^3)`; active-subspace optimization could be `O(sum_k a_k^3)` only after an exact reduction |
| full sweep | `O(m*(n*nt + K*nt^3))` |
| corrected B updates | `O(m*nt^2 + nt^3)` per invoked global/latent update, plus set-local scans |
| diagonal E | `O(n*nt)` |
| full E | `O(n*nt^2 + nt^3)` |
| genetic values/covariance | `O(mn*nt + n*nt^2)` |
| final marker scores | `O(mn*nt)` |

No complexity claim implies a measured performance result.

## 29. Output schema

Phase 17O reuses `mtblr_raw` version 1. Common marker fields retain their
meaning, `meta$data_level="individual"`, and
`meta$backend="mt_bed_bayesc"`. No schema version is required because the
ordinary marker fields remain posterior/effective marker summaries and
marker-space sufficient statistics.

## 30. Marker-space outputs

Before MCMC, compute and retain

```text
wy_jt = x_j'y_t
```

in input marker order using the decoded double workspace and increasing sample
order. Do not maintain marker residuals during MCMC. After the final state,
recompute

```text
r_jt = x_j'R_t
```

once, in the same marker/sample/trait order, for raw output. Both are
`m x nt` doubles. This adds `8*m*nt` retained bytes for `wy` and transient/final
`r`, avoiding an `O(m*nt)` marker-residual state during sampling.

## 31. Sample-space outputs

The minimal internal result owns only what is needed for finalization and
marker-space reconstruction. Final sample residuals and genetic values may
remain core temporaries and are not added to `mtblr_raw` version 1. Posterior
mean individual predictions, fitted values, and predictive diagnostics are
deferred.

## 32. CPO

CPO is unsupported and absent in Phase 17O. Scalar traitwise log-CPO is not a
joint multivariate predictive density. Joint leave-one-individual-out CPO,
traitwise marginal CPO, and conditional-trait CPO require distinct definitions
and validation.

## 33. Phase 17O implementation

The active internal binding is:

```cpp
Rcpp::List mtblr_bed_internal(
  Rcpp::CharacterVector bed_files,
  int n_bed,
  Rcpp::List cls,
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::NumericVector af,
  Rcpp::NumericMatrix Y,
  std::vector<std::vector<double>> beta_init,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<int>> state_init,
  const std::vector<std::vector<int>>& sets,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<std::vector<int>> models,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  std::string residual_covariance,
  int nit,
  int nburn,
  int nthin,
  int seed,
  int method = 4
);
```

Only `method=4`, `residual_covariance` equal to `"full"` or reduction
`"diagonal"`, standardized genotypes, complete centered phenotypes, serial
one-chain execution, and explicit models/sets are accepted. The implementation
uses one `PackedBedMatrix`, one existing borrowed view, one decoded-marker
workspace, sample-space residuals, the verified full-`E` conditional, corrected
MT B/pi/retention semantics, `MtDefaultCoreResult`-compatible finalization, and
`mtblr_raw` version 1. Phase 17P adds only public R metadata and formatting.

## 34. Phase 17P public adapter and later requirements

Phase 17P defines Glist/phenotype ID alignment,
centering evidence and metadata, covariate policy, duplicate/constant/nonfinite
validation, missingness policy, allele/frequency/scale provenance, initialization
defaults, full-versus-diagonal residual covariance controls, diagnostics,
individual prediction outputs, memory warnings, and installed public examples.
Covariates, missingness, predictions, CPO, LE/LD, parallelism, and multichain
execution remain later work.

Phase 17N itself was audit-only and implemented no individual-level MT sampler.
Phase 17O is the first implementation and remains internal-only.

## Phase 17Q concurrency seam

The Phase 17O core is safe for independent concurrent calls after preparation
is separated from its Rcpp binding: its shared inputs are immutable and each
call owns state, residual, workspace, diagnostics, and RNG. Phase 17Q does not
perform that refactor or add OpenMP.

Phase 17R performs the narrow preparation refactor. The serial route and new
internal chains route share one stationary packed owner and fresh immutable
view. Each dispatched unchanged core owns its residual, workspace, RNG, and
state; aggregation occurs only after deterministic failure inspection.

## Phase 17S public dispatch note

The Phase 17O serial function remains the exact maintenance reference. Public
execution now uses the unchanged Phase 17R chains adapter for every chain count.

## Phase 17T convergence seam

The numerical core remains unchanged. Existing `vbs/vgs/ves` contain every
iteration; future diagnostics select only `nburn+1` through `nburn+nit`, and
`nthin` does not thin diagnostic traces.
## Phase 17U internal convergence seam

`mtblr_bed_convergence_trace_internal()` shares the Phase 17R preparation,
dispatch, aggregation, finalizer, legacy adapter, and raw converter. Before
typed chain results are discarded it extracts only post-burn B/G/E diagonal
traces into an Rcpp-free bundle. It neither changes nor reruns the Phase 17O
core, and the ordinary internal chains route remains numerically exact.

Phase 17V activates this unchanged route conditionally from public R. Native
types, extraction, execution, aggregation, wrappers, and registration remain
unchanged.
