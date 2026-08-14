# BLR architecture

The package has two statistical families and three operator families:

| Family | Logical task | Operators |
|---|---|---|
| `stblr` | trait × chain | `csr`, `block_eigen`, `packed_bed` |
| `mtblr` | one complete joint chain | `csr`, `block_eigen`, `packed_bed` |

The seven canonical fitting functions are `stblr_csr()`,
`stblr_csr_annot()`, `stblr_block_eigen()`, `stblr_bed()`, `mtblr_csr()`,
`mtblr_block_eigen()`, and `mtblr_bed()`.

Current annotation-informed effect-variance LV backends are available only for
ST CSR and retained block eigen through `annotation_model = "log_variance"`
with BayesC or BayesR. BED-LV and MT-LV are proposed work, not current routes.
Current and proposed annotation-provider ownership is audited in
[the architecture audit](annotation_prior_architecture_audit.md),
[capability matrix](annotation_prior_architecture_matrix.md), and
[implementation plan](annotation_prior_architecture_implementation_plan.md).

Prepared operator and annotation data are immutable and owned once per fit.
Sampler state, RNG, residuals, workspaces, accumulators, and diagnostic traces
are logical-task private. Seed resolution depends on logical task identity, not
worker assignment. OpenMP workers do not construct R or Rcpp objects.

Native kernels remain separate where sample-space and summary-statistics
likelihoods genuinely differ. The R layer owns public validation, canonical
metadata, formatting, convergence planning, warning aggregation, and memory
preflight.

## Scalar operator roles

The scalar SBayesR/SBayesRC routes expose distinct likelihood contracts rather
than interchangeable storage formats:

| Route | Contract |
|---|---|
| packed BED | Individual-level reference for the supplied genotypes and selected samples. |
| complete CSR | Summary-statistics reference when the complete aligned cross-product operator is computationally feasible. |
| thresholded or windowed CSR | Explicit approximation. Positive definiteness and internal residual identities do not establish fidelity to the source likelihood. |
| retained block eigen | Canonical scalable SBayesRC summary-statistics route, using a projected low-rank likelihood and one global projected residual variance. |

The retained route is closest to the original SBayesRC/GCTB eigen-LD
likelihood, but it is not identical to GCTB's block-specific residual-variance
implementation. Hard-sparse CSR is retained as a public route; its metadata
must identify the approximation and must not imply equivalence to complete
CSR. See `sbayesrc_reference_crosswalk.md` and
`study06_sbayesrc_stabilization.md` for the evidence and qualifications.

Raw backends return named `stblr_raw` or `mtblr_raw` schema-version-1 objects.
Validation precedes one canonical family formatter; positional output and
legacy schema fallback are not supported. Formatted fits use model-semantics
version 2. See `blr_output_schema.md` for field ownership and
`blr_convergence_contract.md` for observational trace capture.

## Approved Phase 0 target architecture

The preceding paragraphs describe current implementation. The approved target
is the shared architecture in
[the unified framework design](blr_unified_framework_design.md), especially
Sections 24--30. It is not yet implemented.

The target distinguishes analysis mode from execution policy:

- `single_trait`, `independent_traits`, and `joint_multitrait` identify
  posterior factorization and statistical coupling;
- `serial`/`parallel` plus `none`, `chains`, `traits`, or `trait_chains`
  identify execution only;
- traits inside a joint chain are never independent tasks.

ST and MT reuse one resolved-specification, global-marker, immutable operator-
resource, likelihood-provider, scheduler, RNG, retention, diagnostics,
provenance, raw-result, and formatted-fit infrastructure. Model-specific state,
probability, covariance, residual, and marker-transition policies remain
separate where posterior targets differ. A second MT BED/CSR/eigen/result stack
is rejected.

An operator resource owns reusable BED, dense, CSR, or block-eigen storage and
provenance. A likelihood provider references a resource and owns a nonempty
ordered trait set, sufficient statistics, local-to-global map, sample size,
scale, population, and residual/error regime. Several providers may reference
one resource. A common-sample full-$V_e$ MT likelihood uses a multi-trait
provider and is not factorized into singleton providers; independent summary
providers may have different $X_t^\top X_t$ and marker coverage.

A full-rank block-eigen representation is exact for its declared block-
diagonal operator. It equals an original dense cross-product only when omitted
cross-block entries are zero or explicitly excluded by the declared
likelihood. A retained-rank representation equals its reconstructed
approximation, not the original dense operator.

The current MT covariance hybrid is replacement-only evidence. It must not be
reused in the target statistical policy or reported as authoritative sampled
$V_b$. Phase 1 introduces contracts around current kernels; later phases
replace the joint MT transition at an explicit scientific checkpoint.
