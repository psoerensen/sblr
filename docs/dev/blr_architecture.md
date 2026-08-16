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

## Approved unified architecture and Phase 3 execution bridge

The preceding paragraphs describe the native implementation. The approved
target is the shared architecture in
[the unified framework design](blr_unified_framework_design.md), especially
Sections 24--30. Phase 1 implements the R-side resolved specification and
validated raw-v2 bridge for eligible maintained ST fits. Phase 2 implements
the common R global-marker map, immutable operator-resource, and likelihood-
provider adapters used by those resolved specifications. It also adopts the
existing immutable native views as the shared native operator boundary:
`SparseLdCsrView`, `BlockEigenView`, and `BedPackedGenotypeView`. Phase 3 adds
one R logical-task plan and one native execution contract around the qualified
existing ST schedulers. It does not replace representation-specific hot loops
or their arithmetic. Phase 4a adds one internal qualification-only joint-MT
policy: two common-sample traits, one packed-BED resource/provider, Cheng
MT-BayesC$\Pi$, fixed full $V_e$, and an authoritative sampled $V_b$. Phase 4b
extends that same route with an optional authoritative sampled full
inverse-Wishart $V_e$; fixed mode remains the Phase 4a zero-RNG policy. It
composes the shared Phase 1--3 infrastructure and existing immutable BED view
rather than creating a second MT stack. Public legacy MT routes remain
unchanged. See
[the Phase 4a checkpoint](blr_phase4a_cheng_mt_bayesc_checkpoint.md) and
[the Phase 4b checkpoint](blr_phase4b_sampled_residual_covariance_checkpoint.md).

Phase 1 public wrappers inspect original call spellings before R partial
matching can accept a scientific-control prefix. Resolved specifications then
validate global marker/allele order, registered likelihood and model policies,
provider/resource alignment, declared priors, and compute bounds before native
dispatch. Compatibility identifiers select migration/formatting behavior only;
they never weaken scientific validation.

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

The R collection owns each resource once; providers retain only its resource
ID and their own ordered map and sufficient statistics. Production legacy
wrappers currently resolve this metadata before dispatch, while their existing
native kernels keep constructing the same representation-specific immutable
owner/view in the same place and order. Dense/CSR/eigen vector evaluation in
the Phase 2 R adapter is a qualification oracle, not a replacement sampler
kernel. A common-sample provider may own several ordered traits and reference
one BED resource, but that representation is not connected to the legacy MT
sampler. See [the Phase 2 checkpoint](blr_phase2_provider_operator_checkpoint.md).

A full-rank block-eigen representation is exact for its declared block-
diagonal operator. It equals an original dense cross-product only when omitted
cross-block entries are zero or explicitly excluded by the declared
likelihood. A retained-rank representation equals its reconstructed
approximation, not the original dense operator.

The current MT covariance hybrid is replacement-only evidence. It is rejected
by the Phase 1 converter and must not be
reused in the target statistical policy or reported as authoritative sampled
$V_b$. Later phases replace the joint MT transition at an explicit scientific
checkpoint.

## Phase 3 execution boundary

Qualified ordinary CSR, packed-BED, block-eigen, fixed-marker annotation, and
learned-logistic annotation routes resolve canonical logical tasks, exact
uint32 task seeds, and retained post-burn indices before native dispatch.
Existing static OpenMP loops consume those immutable descriptors; task-local
RNGs and mutable states remain private. Sampler-local diagnostics distinguish
requested cores, configured workers, actual team size, and the worker ID of
each canonical task without consuming RNG.

Scheduled CSR, log-variance annotation, group annotation, BayesRC/SBayesRC,
and current MT samplers remain version 0. The joint-MT task topology can be
constructed and tested without calling the legacy MT covariance sampler.
Provider, block, set, and region traversal are not logical execution modes.
See [the Phase 3 checkpoint](blr_phase3_execution_checkpoint.md).
