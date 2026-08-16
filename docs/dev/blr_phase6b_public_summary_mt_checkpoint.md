# Phase 6B public independent-summary Cheng MT checkpoint

## Status and boundary

Phase 6B promotes the independently qualified Phase 6A implementation through
the public `mtblr_csr()` and `mtblr_block_eigen()` interfaces. Both functions
resolve public descriptors into the same Phase 2 provider/resource collection,
invoke the single Phase 6A native kernel, validate `blr_raw` version 2, and use
the shared MT formatter. No legacy covariance-hybrid adapter is involved.

The public boundary supports feasible general $T\geq2$, complete binary Cheng
activity patterns, sampled Dirichlet activity-pattern probabilities, and one
authoritative inverse-Wishart $V_b$. Each provider owns exactly one trait, a
score vector, a fixed positive residual scale, sample-size metadata, one local
marker map, and either a CSR or block-eigen operator. Multiple providers for
one trait are allowed only under an explicit independent-summary declaration.

Sample overlap, a sampled full $V_e$, sampled provider scales, MT-BayesR,
MT-BayesRC, annotations, regional covariance, templates, factors, and marker
scheduling remain deferred.

## Likelihood and transition

For provider $p$ assigned to trait $t(p)$, the public functions use

$$
\log L_p=-\frac{1}{2\phi_p}
\left(
\boldsymbol\alpha_{p,t(p)}^\top C_p\boldsymbol\alpha_{p,t(p)}
-2\boldsymbol\alpha_{p,t(p)}^\top\mathbf s_p
\right)+\text{constant}.
$$

Independent-provider terms add. Missing provider markers add no information.
The marker score and diagonal information are

$$
h_{jt}=\sum_{p:t(p)=t,\,j\in p}\frac{r_{p,j}^{(-j)}}{\phi_p},
\qquad
d_{jt}=\sum_{p:t(p)=t,\,j\in p}\frac{C_{p,jj}}{\phi_p}.
$$

Pattern draws, null collapse, conditional completion, the Dirichlet update,
and the single inverse-Wishart $V_b$ update are unchanged from Phase 6A.

## Public descriptors

Both interfaces accept a named `providers` list, a named operator-resource
descriptor list, global marker/allele metadata, ordered trait IDs, explicit
$V_b$ initialization and degrees-of-freedom/scale prior, and the common Phase 3
MCMC, seed, chain, convergence, trace, and memory controls. Public argument
names are exact; obsolete hybrid arguments are not translated.

CSR descriptors contain the maintained explicit CSR payload or a
`sparseLD_read_CSR()` payload plus a declared cross-product diagonal. CSR
approximation semantics remain explicit and production execution is sparse.

Block-eigen descriptors contain provider-local blocks and their retained
eigenpairs. Full-rank blocks are exact relative to the declared block-diagonal
operator. Retained-rank blocks target the retained reconstruction. Omitted
cross-block terms are never reconstructed and no global eigendecomposition is
formed.

## Output and provenance

Both public functions return the standard `mtblr_fit`. The validated raw-v2
object remains attached through the established output contract. Marker,
trait, pattern, chain, draw, provider, and covariance axes remain named and
non-dropped. PIPs, markerwise pattern probabilities, all-active probabilities,
$\boldsymbol\Pi$, authoritative $V_b$, provider metadata, approximation
metadata, seeds, retention/convergence indices, and compact pattern diagnostics
retain their Phase 6A meanings.

Predictions and full residual-covariance fields are explicitly `NULL`; no
summary quadratic is relabelled as SSE or true genomic variance.

## Legacy cleanup and validation

The standalone legacy MT CSR translation unit was removed. Legacy summary
wrappers and native registrations were removed; historical implementations
interleaved with maintained BED source remain compiled only as unregistered,
unreachable source to avoid risky unrelated BED surgery. Generated interfaces
contain only the corrected Phase 6A summary entry points.

Focused public tests establish CSR and block-eigen public/internal identity,
heterogeneous maps, retained-rank provenance, serial/parallel identity,
raw-v2/formatter consistency, overlap rejection, unsupported-method rejection,
and legacy-entry isolation. See
[`../../examples/workflows/mt_summary_cheng_workflow.R`](../../examples/workflows/mt_summary_cheng_workflow.R)
for a small public workflow.

## Deferred work

Overlap-aware likelihoods, mixed CSR/block-eigen public collections, public
provider classes, MT mixture models, annotations, regional covariance, and
marker scheduling require separate qualification and promotion.
