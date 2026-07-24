# Internal trait-specific MT BayesC CSR contract

## 1. Purpose

Provide internal trait/study-specific sparse-LD execution for corrected MT
BayesC using the canonical scalar `SparseLdCsrView`.

## 2. Public status

`mtblr_csr_internal()` is maintenance-only. `sblr()` remains dense/default and
no namespace export or public CSR algorithm exists.

## 3. Statistical target

The target is exactly the corrected Phase 17C dense contract: joint trait
patterns, covariance policies, RNG order, update controls, and retained counts.

## 4. Shared statistical core

Dense and CSR wrappers instantiate `run_mt_bayesc_core_impl<DataView>()`.
There is one Gibbs loop and one active latent marker conditional.

## 5. Data contract

`MtCsrDataView` borrows `wy`, `yy`, `n`, and `MtSparseLdBundleView`, which holds
one canonical CSR view per trait/study.

## 6. Fully shared mode

One prefix builds one storage and diagonal owner. All trait diagonals must be
exactly identical; `nt` views borrow the same buffers.

## 7. Trait-specific mode

`nt` prefixes build `nt` owners. Values, diagonals, sparsity, sample sizes, and
references may differ, while marker count and positional domain remain common.

## 8. Marker-domain invariant

Trait CSR row, summary marker, joint effect position, output row, canonical ID,
and effect-allele orientation at position `i` must denote the same marker.

## 9. Scaling

`wy=x'y`, diagonal `ww=x'x`, and CSR `xij=r_ij*sqrt(ww_i*ww_j)` use the
standardized-genotype scale.

## 10. Ownership

The adapter fully constructs reserved storage and diagonal owners before views.
Owners remain stable through core execution; no chain or global cache exists.

## 11. Disk I/O

All canonical CSR reads, validation, scaling, and symmetrization finish before
the Gibbs core begins. There is no MCMC-time I/O.

## 12. Output

Both routes return `MtDefaultCoreResult`, use `finalize_mt_default_result()`,
and use one binding-neutral legacy 20-position adapter.

## 13. Dense reduction

The oracle places the diagonal first and float-rounded CSR entries in stored
row order. Shared, updated, initialized, set, three-trait, trait-value, and
independent-pattern reductions pass at `1e-12`.

## 14. Study and ancestry use

A trait dimension may denote phenotype, study, ancestry-specific phenotype, or
cohort. Each can use distinct `n`, `ww`, LD, sparsity, and reference panel.

## 15. Sample overlap

Residual covariance is not declared to model GWAS sample overlap. No new
cross-study sampling-error likelihood is introduced.

## 16. Memory and copying

Shared mode has one storage/diagonal owner and `nt` views; trait-specific mode
has `nt` owners/views. Per-chain CSR bytes are zero. Dense oracle storage is
test-only.

## 17. Public Phase 17J requirements

Phase 17J must design R/Glist ID and allele alignment, study/ancestry/reference
metadata, explicit marker intersection, naming, validation, API, and stable
named MT raw/fit conversion without changing this core.
Phase 17J exposes this validated numerical route through `mtblr_csr()`. Public
marker/allele/scale/resource validation remains R-owned; the internal core
continues to accept an already aligned positional domain.
# Parallel Phase 17L representation

The canonical CSR route remains unchanged. Phase 17L adds a sibling block-eigen data view that uses the same representation seam and shared statistical core, result, finalizer, and legacy adapter.

# Phase 17M raw conversion protection

CSR and public block-eigen raw boundaries now share one representation-neutral
legacy-to-raw converter. The CSR backend name, namespaces, dimensions,
diagnostics, alignment enrichment, and formatted output remain unchanged.
