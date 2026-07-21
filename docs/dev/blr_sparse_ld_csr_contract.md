# Canonical sparse-LD CSR contract

## 1. Purpose

`sblr::core::SparseLdCsrStorage` and `SparseLdCsrView` are the one numerical
CSR vocabulary for scalar ST-BLR and future trait/study-specific MT-BLR. They
describe storage only; sampler, chain, trait, result, and R metadata are out of
scope.

## 2. Disk format

For prefix `P`, the scalar builder reads `P.row_ptr.u64.bin` (`uint64`),
`P.col_idx.u32.0based.bin` (`uint32`, zero based), `P.values.f32.bin` (`float`
raw correlations), and `P.meta.txt`. Mandatory metadata keys are `n_variants`
and `nnz`. Disk row order is marker order. Canonical input is one triangular
off-diagonal representation; the builder does not infer or reorder markers.

## 3. Internal storage

The owning fields are `vector<uint64_t> ptr`, `vector<int> idx`, and
`vector<float> xij`, plus `input_nnz`, `symmetric_nnz`, `max_abs_r`, and
`max_abs_xij`. The legacy payload names are intentionally retained to avoid
broad numerical-source churn. There is no second canonical owner definition;
`STLDCSR` is an alias.

## 4. Borrowed view

`SparseLdCsrView` contains `marker_count`, immutable `row_ptr` plus size,
immutable `column_index`, immutable `offdiag_xij` plus count, and a pointer to
the immutable `arma::rowvec` diagonal. It owns nothing. Storage and diagonal
must outlive every view and all tasks using it.

## 5. Value meaning

Disk values are correlations `r_ij`. Internal values are
`xij = r_ij * sqrt(ww_i * ww_j)`, narrowed to `float`. `wy` and `ww` retain
their existing double precision.

## 6. Diagonal policy

Diagonal `ww = X'X` values are separate double-precision data. Disk diagonal
entries are ignored and internal CSR rows contain off-diagonal values only.

## 7. Symmetry

Each input off-diagonal is inserted in both directions. Upper- or
lower-triangular input works. Already symmetric or repeated input creates
duplicate internal entries and is therefore outside the canonical input
contract. Columns are not sorted; zero values are retained; empty rows are
accepted and can receive mirrored entries.

## 8. Marker order

CSR row `i`, summary marker `i`, effect position `i`, output row `i`, marker ID
`i`, and effect-allele orientation `i` must denote the same marker. Native
synthetic IDs do not prove this invariant; R preparation currently owns it.

## 9. Allele orientation

`make_summary_stats()` and BED-derived construction establish orientation by
construction. External summaries and CSR files require caller alignment.
`Glist$rsidsLD` to `Glist$rsids` is explicitly matched for `selection_s`, but
ordinary CSR fitting does not universally validate effect allele, other allele,
strand, or reference-panel metadata.

## 10. Validation

Builder-time checks cover metadata, file sizes, indices, finite correlations,
positive finite `ww`, and fill counts. Mandatory view-time checks cover a
positive marker count, pointer presence and size, zero/nondecreasing row
pointers, terminal `nnz`, pointer presence for nonempty payloads, in-range
off-diagonal columns, finite values, and a finite positive diagonal of length
`m`. Symmetry, duplicate detection, and sorting are optional deep checks, not
new runtime requirements, because the current builder does not guarantee the
latter two for arbitrary input.

## 11. Ownership and copying

The builder owns storage. `CsrOperator` currently makes one pre-existing
compatibility copy of storage and diagonal. Its view borrows that owner;
chain tasks share it read-only and make no CSR copy. No disk access occurs
during MCMC. Removing the adapter copy is a separate optimization.

## 12. Scalar use

All scalar CSR backends share the owning definition through `STLDCSR`.
Ordinary canonical CSR BayesC actively composes `SparseLdCsrView`. Scheduled
BayesC, BayesR, SBayesRC, fixed-prior, group, and learned-annotation routes
remain on narrow aliases/adapters pending model-specific mechanical migration.

## 13. Future MT use

A future bundle can hold one view per trait/study. Views may share every
pointer, share row/index pointers with trait-specific values and diagonals, or
use fully independent structures. The view itself encodes no sharing policy.

## 14. Study and ancestry metadata

External metadata must eventually include trait ID, study ID, population,
ancestry, LD reference, sample size, marker IDs, alleles, marker-set policy,
and overlap policy. These do not belong in the numerical view.

## 15. Sample overlap

Independent, known-overlap, and unknown-overlap studies require an explicit
future likelihood policy. Residual covariance alone is not declared to model
GWAS sample overlap.

## 16. Noncanonical representations

The retained internal `mt_cpg_omp_csr.cpp` `LDCSR` is unsupported research: it
assumes one shared LD object and is not a future canonical type. Builder-side
`SparseLDCSR` in `sparse_ld_bed_core.cpp` is a transient disk-construction
buffer, not an execution storage contract.

## 17. Evolution policy

Changes extend this shared contract or use explicit compatibility adapters;
they must not introduce parallel scalar and MT file/storage formats. Marker
alignment and trait-resource bundles are Phase 17I concerns.

## Phase 17I activation

`MtSparseLdBundleView` now composes one canonical view per trait/study for the
internal MT BayesC route. Fully shared pointers, shared patterns with distinct
values, and independent patterns are supported. Biological alignment and a
public route remain Phase 17J responsibilities.
Phase 17J public resolution distinguishes a fully shared pre-scaled operator,
one shared raw-correlation reference rebuilt with trait-specific diagonals, and
trait-specific reference resources. The binary CSR format is unchanged.
