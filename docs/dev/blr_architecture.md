# BLR architecture

The package has two statistical families and three operator families:

| Family | Logical task | Operators |
|---|---|---|
| `stblr` | trait × chain | `csr`, `block_eigen`, `packed_bed` |
| `mtblr` | one complete joint chain | `csr`, `block_eigen`, `packed_bed` |

The seven canonical fitting functions are `stblr_csr()`,
`stblr_csr_annot()`, `stblr_block_eigen()`, `stblr_bed()`, `mtblr_csr()`,
`mtblr_block_eigen()`, and `mtblr_bed()`.

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
