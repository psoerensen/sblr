# sblr Note Strategy

## Purpose

This file defines a documentation structure for `sblr` that makes the main
conceptual path readable without deleting the derivations, assumptions,
implementation sketches, cautions, or research ideas in the existing Quarto
drafts.

The primary source material is:

- `docs/sblr_theory.qmd`
- `docs/annotation_informed_priors_stblr_mtblr.qmd`
- `examples/workflows/sparse_ld_bed_workflow.R`
- `examples/workflows/annotation_based_models.R`

The drafts should remain unchanged until their content has been migrated and
checked. When the new notes are written, copies of the drafts should be
retained under `docs/notes/drafts/` as an explicit archive. The archive is not
a substitute for mapping content into maintained notes. The workflow scripts
remain executable source documents and should be linked rather than copied
wholesale into the notes.

## Documentation Principles

1. **Preserve technical depth.** Do not shorten the documentation by deleting
   equations, derivations, assumptions, implementation details, or cautions.
   Move dense material into technical appendices and link to it from the main
   notes.
2. **Provide a readable main path.** Main notes should explain what problem is
   being solved, the model choices, the package interface, and how to interpret
   outputs. They should contain enough equations to define each model, but
   defer full derivations and implementation sketches to appendices.
3. **Separate model axes explicitly.** Every model discussion should identify:
   data representation, trait structure, prior structure, and whether
   annotation parameters are fixed, estimated in a separate empirical-Bayes
   stage, or learned jointly.
4. **Distinguish implemented interfaces from research proposals.** Current
   public wrappers should be labeled as implemented. Draft ideas such as
   annotation-aware MTBLR sharing-model priors, sparse annotation storage, and
   warm-started regional MTBLR should be retained but labeled as proposed or
   experimental where they are not public package interfaces.
5. **Keep package help and notes complementary.** Rd help should remain the
   argument-level reference. Notes should explain model meaning, assumptions,
   workflow choices, diagnostics, and relationships among wrappers.
6. **Keep examples compact and reproducible.** Main notes should show short
   workflow excerpts and output interpretation. Long-running scripts remain in
   `examples/workflows/`, with links and explicit runtime cautions.
7. **State inference limits.** Annotation enrichment, genetic sharing, and
   posterior association are not causal identification. Preserve the draft's
   causal-interpretation cautions and downstream-analysis guidance.
8. **Treat workflows as executable architecture documents.** The workflow
   scripts show how package layers compose, which ordering and scaling
   contracts matter, and how comparable outputs are interpreted. Notes should
   explain those choices and link to the scripts rather than treating them as
   incidental examples.

## Classification Framework

Use the following four-axis classification consistently in the notes.

| Axis | Categories | Documentation rule |
|---|---|---|
| Data | Individual-level BED; summary statistics plus dense LD; summary statistics plus disk-backed CSR LD | Introduce representations before models and state what sufficient statistics are required. |
| Traits | Single-trait; parallel single-trait; multi-trait with covariance/sharing | Do not describe parallel ST-BLR as a full joint multi-trait model. |
| Prior | Plain/global; fixed marker-specific annotation prior; group prior; learned annotation prior; SBayesRC-style mixture prior | Give each prior family its own definition and comparison table. |
| Estimation | Fixed/external prior; empirical-Bayes multi-stage prior; directly learned annotation effect; proposed hierarchical extension | Label whether uncertainty in annotation parameters is propagated. |

### Required Distinctions

- **Summary-statistics versus individual-level models:** individual-level BED
  workflows read genotype data directly; summary-statistics workflows operate
  on `wy`, `ww`, `yy`, LD/CSR structures, and possibly cross-trait `SSY` and
  overlap matrices.
- **Single-trait versus multi-trait:** single-trait models use scalar marker
  effects and variances; multi-trait models use effect vectors, covariance
  matrices, and potentially trait-sharing models.
- **Plain versus annotation-informed priors:** plain models use global prior
  parameters; annotation-informed models allow priors to vary by marker,
  group, or annotation profile.
- **Empirical-Bayes versus direct learned annotation models:** empirical-Bayes
  models estimate annotation quantities in a prior stage and treat them as
  fixed later; direct learned models update annotation effects during the
  genomic sampler.
- **Fixed, group, and SBayesRC-style priors:** fixed priors pass marker-specific
  inclusion/variance values; group priors assign markers to one group;
  SBayesRC-style priors use overlapping continuous annotations to determine
  mixture-component probabilities through probit stick-breaking.
- **Workflow scripts versus conceptual notes:** conceptual notes define stable
  models, contracts, and choices; workflow scripts demonstrate the exact
  sequencing of current wrappers. The notes should cite and explain scripts,
  but should not duplicate every runtime-heavy call or present scripts as
  formal tests.

## Proposed File Tree

```text
docs/
  notes/
    NOTE_STRATEGY.md
    index.qmd
    model_overview.qmd
    data_representations.qmd
    sparse_ld_and_bed_workflows.qmd
    single_trait_models.qmd
    multi_trait_models.qmd
    annotation_informed_models.qmd
    workflows_and_outputs.qmd
    technical_summary_statistics.qmd
    technical_multitrait_overlap.qmd
    technical_sparse_ld_csr.qmd
    technical_annotation_priors.qmd
    technical_sbayesrc_annotations.qmd
    drafts/
      sblr_theory.qmd
      annotation_informed_priors_stblr_mtblr.qmd
```

The two files under `drafts/` should be archival copies created only when the
new notes are written. The current source drafts under `docs/` must not be
deleted or overwritten during migration.

## Main Explanatory Notes

### `index.qmd`

**Purpose:** Provide the entry point, terminology, reading paths, and model
selection map.

**Audience:** New package users, collaborators, and readers deciding which
model or workflow applies.

**Detail level:** Short and navigational. No long derivations.

**Sections:**

1. What `sblr` is and its experimental status
2. Four documentation axes: data, traits, priors, estimation
3. Choose a workflow
4. Choose a model
5. Main notes versus technical appendices
6. Runtime, convergence, and reproducibility cautions
7. Public API and example-script links

**Concepts/functions/workflows:** `sblr()`, `stblr_bed_marker()`, `stblr_csr()`,
the four annotation wrappers, `examples/workflows/annotation_based_models.R`,
and the other clean workflows listed in `examples/README.md`.

**Draft reuse:** Motivation paragraphs from both drafts, rewritten as a concise
navigation overview. No technical content is removed; details are linked to
the relevant notes and appendices.

**Appendix links:** All five technical appendices.

### `model_overview.qmd`

**Purpose:** Define the common Bayesian linear regression model and show how
the package's model families differ.

**Audience:** Readers with basic Bayesian regression knowledge who need a
coherent package-level model map.

**Detail level:** Moderate. Include defining equations and comparison tables,
but not full marker-update derivations.

**Sections:**

1. Common model: \(y = Xb + e\)
2. Summary-statistic form: `wy`, `ww`, LD, residual score `r`
3. Sparse mixture priors and posterior inclusion probabilities
4. Single-trait, parallel ST-BLR, and joint MTBLR
5. Global versus marker/group/annotation-specific priors
6. Implemented models versus proposed extensions
7. Model comparison matrix

**Equations/concepts:** \(wy=X^\top y\), \(r=X^\top y-X^\top Xb\),
BayesC-style spike-and-slab, posterior means `bm`, PIPs `dm`, variance
components, effect covariance, trait-sharing models.

**Functions/workflows:** `sblr()`, `stblr_csr()`, `stblr_bed_marker()`,
`stblr_csr_prior_annot()`, `stblr_csr_learn_annot()`,
`stblr_csr_group_annot()`, `stblr_csr_sbayesrc_generic()`.

**Draft reuse:** Annotation draft sections "Motivation", "Summary-statistic
residual notation", and the high-level hierarchical summaries under "Direct
annotation-aware STBLR and MTBLR as hierarchical models".

**Appendix links:** `technical_summary_statistics.qmd`,
`technical_annotation_priors.qmd`, and `technical_sbayesrc_annotations.qmd`.

### `data_representations.qmd`

**Purpose:** Explain what data each workflow consumes, marker-order contracts,
and how individual-level data become summary statistics and sparse LD.

**Audience:** Analysts preparing inputs and developers connecting new data
sources.

**Detail level:** Operational plus conceptual; defer cross-trait derivations.

**Sections:**

1. Individual-level PLINK BED representation
2. Dense summary-statistic representation
3. Disk-backed CSR sparse-LD representation
4. Sufficient statistics: `wy`, `ww`, `yy`
5. Multi-trait cross-products, sample overlap, and `SSY`
6. Marker ordering and annotation alignment
7. Representation-specific limitations and validation checks

**Equations/concepts:** \(X^\top y\), diagonal \(X^\top X\), sparse LD,
`SSY`, overlap matrix \(N\), exact versus reconstructed cross-trait
quantities.

**Functions/workflows:** `bed_xtx_xty()`, `sparseLD_stream_CSR()`,
`sparseLD_read_CSR()`, `readLD_to_CSR()`, `stblr_bed_marker()`,
`stblr_csr()`, and the preparation steps in both principal workflow scripts.
This note defines representations and contracts; it defers the end-to-end
BED-to-CSR bridge to `sparse_ld_and_bed_workflows.qmd`.

**Draft reuse:** The overview and practical motivation from
`docs/sblr_theory.qmd`; annotation draft "Summary-statistic residual notation".

**Appendix links:** `technical_summary_statistics.qmd`,
`technical_multitrait_overlap.qmd`, and `technical_sparse_ld_csr.qmd`.

### `sparse_ld_and_bed_workflows.qmd`

**Purpose:** Explain the central bridge from one individual-level BED genotype
source to two fitting paths: a direct BED BLR path and a BED-derived
sufficient-statistics plus disk-backed CSR LD path.

**Audience:** Analysts choosing between BED-direct and CSR fitting, and
developers validating that both paths use aligned markers and comparable
statistics.

**Detail level:** Operational and architectural. Include the defining
sufficient-statistic equations and data-flow diagram, while moving CSR storage
and sparse-LD construction details to an appendix.

**Sections:**

1. One genotype source, two fitting paths
2. Selecting samples and markers from a BED-backed `Glist`
3. Simulating or supplying aligned phenotypes
4. Path A: direct individual-level fitting with `stblr_bed_marker()`
5. Path B, step 1: deriving `wy`, `ww`, and `yy` with `bed_xtx_xty()`
6. Path B, step 2: streaming sparse LD to disk with `sparseLD_stream_CSR()`
7. Inspecting and validating CSR files with `sparseLD_read_CSR()`
8. Path B, step 3: fitting summary-statistics BLR with `stblr_csr()`
9. Comparing BED-direct and CSR posterior outputs
10. Marker order, allele frequency, scaling, row selection, and LD-threshold
    contracts
11. Runtime, disk, reproducibility, and convergence cautions

**Equations/concepts:** \(X^\top y\), diagonal \(X^\top X\), \(y^\top y\),
sparse approximation to \(X^\top X\), direct versus sufficient-statistic
likelihood inputs, disk-backed CSR, and why aligned marker order is a
cross-layer invariant.

**Functions/workflows:** Use
`examples/workflows/sparse_ld_bed_workflow.R` as the primary worked source.
Mention `bed_xtx_xty()`, `sparseLD_stream_CSR()`, `sparseLD_read_CSR()`,
`stblr_csr()`, `stblr_bed_marker()`, and the qgg comparison calls only as an
external reference point.

**Source reuse:** Reuse the workflow's sequencing, comments, validation steps,
and CSR-versus-BED comparison logic. Extract only short code fragments into
the note and link to the complete script.

**Appendix links:** `technical_summary_statistics.qmd` and
`technical_sparse_ld_csr.qmd`.

### `single_trait_models.qmd`

**Purpose:** Explain single-trait BLR/ST-BLR models and the difference between
BED-direct and summary-statistics CSR fitting.

**Audience:** Analysts fitting one trait or running one ST-BLR per trait.

**Detail level:** Model definition, prior interpretation, workflow choices,
and output interpretation. Marker-update algebra moves to an appendix.

**Sections:**

1. Single-trait BLR model
2. Global BayesC-style prior
3. Individual-level BED workflow
4. Summary-statistics CSR workflow
5. Scheduled sparse updates and runtime controls
6. Posterior outputs and diagnostics
7. Transition from plain to annotation-informed priors

**Equations/concepts:** spike-and-slab prior, global \(\pi\), effect variance
\(v_b\), PIP, posterior marker effects, residual score.

**Functions/workflows:** `stblr_bed_marker()`, `stblr_csr()`, relevant parts of
`sblr()`, `bm`, `dm`, `vbs`, `vgs`, `ves`, `covb`, `covg`, `cove`, with the
operational comparison delegated to `sparse_ld_and_bed_workflows.qmd`.

**Draft reuse:** Annotation draft "Using marker-specific prior variances in
STBLR" and "Modified STBLR marker update" only as concise summaries; the full
derivations go to `technical_annotation_priors.qmd`.

**Appendix links:** `technical_summary_statistics.qmd`,
`technical_sparse_ld_csr.qmd`, and `technical_annotation_priors.qmd`.

### `multi_trait_models.qmd`

**Purpose:** Explain joint multi-trait modeling, covariance structures,
trait-sharing models, and the role of sample overlap.

**Audience:** Researchers using or extending MTBLR.

**Detail level:** Moderate-to-high. Include model equations and interpretation;
move overlap derivation and proposed annotation-aware MTBLR algebra to
appendices.

**Sections:**

1. Why model traits jointly
2. Effect vectors and covariance matrices
3. Trait-sharing model space
4. Residual covariance and cross-trait sufficient statistics
5. Full versus partial sample overlap
6. Annotation-informed MTBLR concepts
7. Current implementation status and proposed extensions
8. Interpretation: sharing is not causal direction

**Equations/concepts:** \(\mathbf b_i\), \(B\), \(E\), sharing models such as
`000` through `111`, model probabilities, cross-trait `SSY`, overlap \(N\).

**Functions/workflows:** `sblr()` and multi-trait components represented in
fit covariance outputs; clearly label annotation-aware MTBLR sharing priors as
research proposals unless a public wrapper exists.

**Draft reuse:** All of `docs/sblr_theory.qmd` at a summarized level; annotation
draft sections "Annotation-specific MTBLR prior covariance matrices",
"Annotation-specific MTBLR model priors", "Direct annotation-informed MTBLR",
"MTBLR model prior in the marker update", and "Annotation-informed MTBLR
effect covariance".

**Appendix links:** `technical_multitrait_overlap.qmd` and
`technical_annotation_priors.qmd`.

### `annotation_informed_models.qmd`

**Purpose:** Be the main conceptual guide to annotation-informed modeling,
clearly separating implemented prior families and empirical-Bayes strategies.

**Audience:** Analysts comparing annotation models and researchers designing
annotation priors.

**Detail level:** High-level equations, comparison tables, implementation
status, and practical cautions. Full derivations and C++ sketches move to
appendices.

**Sections:**

1. Why annotations inform priors
2. Annotation matrix, overlap, centering, and marker ordering
3. Plain/global prior baseline
4. Fixed marker-specific priors
5. Group priors
6. Direct learned annotation effects
7. SBayesRC-style overlapping-annotation mixture priors
8. Empirical-Bayes annotation workflows
9. Fixed versus learned versus empirical-Bayes comparison
10. Regularization, caps, shrinkage, and overfitting
11. Interpretation and causal limits

**Equations/concepts:** \(\operatorname{logit}(\pi_{it})\), log variance
multipliers, group-specific priors, probit stick-breaking, expected gamma,
overlapping annotation attribution, cross-fitting.

**Functions/workflows:** `mtsim_annotation()`, `summarize_annotation_signal()`,
`stblr_csr_prior_annot()`, `stblr_csr_learn_annot()`,
`stblr_csr_group_annot()`, `stblr_csr_sbayesrc_generic()`,
`make_sbayesrc_alpha_init()`, and the four SBayesRC diagnostic helpers.

**Draft reuse:** Annotation draft "Motivation", "Handling overlapping
annotations", "Deriving annotation-specific prior inclusion probabilities",
"Deriving annotation-specific prior effect variances", "Marker-specific prior
scales for overlapping annotations", "Empirical-Bayes priors from an
annotation-level STBLR or MTBLR", "Direct annotation modelling inside the
genomic STBLR or MTBLR sampler", "Centering, shrinkage, and caps", "Recommended
first direct model", and "Practical cautions".

**Appendix links:** `technical_annotation_priors.qmd` and
`technical_sbayesrc_annotations.qmd`.

### `workflows_and_outputs.qmd`

**Purpose:** Connect model choices to executable workflows and explain compact
posterior summaries without printing marker-scale matrices. This note indexes
the complete workflows; it does not repeat the detailed BED-to-CSR bridge from
`sparse_ld_and_bed_workflows.qmd`.

**Audience:** Users running package workflows.

**Detail level:** Practical and concise, with links to theory.

**Sections:**

1. Workflow decision table
2. Workflow scripts as executable architecture documents
3. BED-direct and summary-statistics CSR workflow index
4. Annotation-model comparison workflow
5. MCMC settings: demonstration versus real analysis
6. Fit object structure
7. Compact fit summaries
8. Top-marker summaries
9. Annotation-level summaries
10. SBayesRC prior-architecture diagnostics
11. Convergence, posterior predictive checks, and reproducibility

**Concepts:** marker ordering, runtime-heavy steps, full matrices versus
printed summaries, PIP totals, mean absolute effects, variance traces,
covariance/correlation summaries, expected gamma.

**Functions/workflows:** all public wrappers. Treat
`examples/workflows/sparse_ld_bed_workflow.R` and
`examples/workflows/annotation_based_models.R` as coequal primary worked
sources; also link to `examples/workflows/basic_sblr_summary_stats.R`.

**Draft reuse:** Annotation draft "Recommended modelling hierarchy", "Minimal
example workflow", "Recommended empirical-Bayes workflow", "Warm-started
regional MTBLR", and "Practical causal workflow". Keep warm-starting labeled as
experimental/proposed if not exposed by a stable wrapper.

**Appendix links:** all technical appendices as appropriate, especially
`technical_sparse_ld_csr.qmd` for the BED/CSR workflow and
`technical_annotation_priors.qmd` plus `technical_sbayesrc_annotations.qmd`
for the annotation workflow.

## Technical Appendices

### `technical_summary_statistics.qmd`

**Purpose:** Preserve detailed algebra connecting individual-level models to
summary-statistic residual calculations.

**Audience:** Method developers and readers validating sampler equations.

**Detail level:** Full derivation and numerical assumptions.

**Sections and content:**

1. Summary-statistic notation: `wy`, `ww`, LD, `r`
2. Identity \(X^\top Xb = wy-r\)
3. Genomic variance and covariance from residual scores
4. Set-specific and weighted genomic variance contributions
5. Assumptions and numerical caveats
6. R helper implementations and their package status

**Draft reuse:** Annotation draft "Summary-statistic residual notation",
"Set-specific genomic variance decomposition", "Combining effect size and
inclusion probability", and their R implementations.

**Functions/workflows and links:** Connect `bed_xtx_xty()`, `stblr_csr()`, and
both primary workflows back to `data_representations.qmd`,
`sparse_ld_and_bed_workflows.qmd`, and `single_trait_models.qmd`.

### `technical_multitrait_overlap.qmd`

**Purpose:** Preserve the complete `SSY` and sample-overlap derivation.

**Audience:** Readers implementing or auditing multi-trait summary-statistic
likelihoods.

**Detail level:** Full first-principles derivation, assumptions, and edge
cases.

**Sections and content:**

1. Cross-trait `SSY` definition
2. Phenotypic variance from summary statistics
3. Covariance reconstruction from \(R_y\)
4. Full-overlap and fractional-overlap \(N\)
5. Reconstruction \(\widehat{SSY}=N\odot\widehat{\operatorname{Cov}}(y)\)
6. Why full residual covariance estimation requires `SSY` and \(N\)
7. First-principles derivation and assumptions
8. Limitations, missingness, and sensitivity analysis

**Draft reuse:** Every section of `docs/sblr_theory.qmd`, including the
duplicated first-principles derivation. Consolidate repeated prose, but retain
all equations and assumptions.

**Functions/workflows and links:** Relate `SSY`, overlap matrices, and current
multi-trait interfaces back to `data_representations.qmd` and
`multi_trait_models.qmd`.

### `technical_sparse_ld_csr.qmd`

**Purpose:** Preserve the storage, construction, alignment, and computational
details behind disk-backed sparse LD used by CSR BLR workflows.

**Audience:** Method developers, maintainers, and analysts diagnosing
BED-to-CSR discrepancies or scaling behavior.

**Detail level:** Full storage and computational contract, including
implementation-sensitive caveats.

**Sections and content:**

1. From BED genotypes to sufficient statistics and sparse LD
2. Sparse approximation to \(X^\top X\): threshold and distance choices
3. CSR representation: row pointers, column indices, and values
4. Disk-backed prefix and file contract
5. One-based R marker identities versus internal CSR indexing
6. Marker ordering, sample rows, allele frequency, and scaling invariants
7. Streaming and blocking behavior, memory use, and disk tradeoffs
8. Reading and validating CSR with `sparseLD_read_CSR()`
9. How `stblr_csr()` consumes statistics and the CSR prefix
10. Scheduled sparse updates and performance controls
11. Numerical approximation, reproducibility, and BED-versus-CSR comparison

**Functions/workflow reuse:** Use
`examples/workflows/sparse_ld_bed_workflow.R` as the primary operational
source and align contracts with `bed_xtx_xty()`, `sparseLD_stream_CSR()`,
`sparseLD_read_CSR()`, `readLD_to_CSR()`, and `stblr_csr()` help. Preserve
implementation-sensitive details, but label them so they can be checked when
the CSR format changes.

**Main-note links:** `data_representations.qmd`,
`sparse_ld_and_bed_workflows.qmd`, `single_trait_models.qmd`, and
`workflows_and_outputs.qmd`.

### `technical_annotation_priors.qmd`

**Purpose:** Preserve detailed derivations for annotation summaries,
empirical-Bayes priors, direct annotation-aware STBLR/MTBLR, and proposed
implementation changes.

**Audience:** Method developers and researchers extending samplers.

**Detail level:** Full derivations, implementation sketches, uncertainty
considerations, and proposed extensions.

**Sections and content:**

1. Annotation-level posterior summaries and overlap
2. PIP enrichment and variance attribution
3. Empirical-Bayes inclusion and variance priors
4. Marker-specific prior scale construction
5. Full STBLR marker-update derivation with marker-specific priors
6. Annotation-specific MTBLR covariance and sharing-model priors
7. Annotation-level STBLR/MTBLR as a prior-learning stage
8. Direct hierarchical annotation models
9. Fixed versus jointly updated annotation coefficients
10. Proposed C++ changes and sparse annotation representation
11. Centering, shrinkage, caps, cross-fitting, and uncertainty
12. Warm-started regional MTBLR

**Draft reuse:** All detailed equations, R implementations, and C++ sketches
from the annotation draft except the SBayesRC-specific material and causal
interpretation, which have dedicated destinations.

**Functions/workflows and links:** Connect current annotation wrappers and
`examples/workflows/annotation_based_models.R` back to
`annotation_informed_models.qmd`, `multi_trait_models.qmd`, and
`workflows_and_outputs.qmd`.

### `technical_sbayesrc_annotations.qmd`

**Purpose:** Explain the implemented SBayesRC-style annotation prior and its
diagnostics in full mathematical detail.

**Audience:** Users interpreting SBayesRC fits and developers validating the
stick-breaking implementation.

**Detail level:** Full probability construction and diagnostic contracts, with
implementation-aligned examples.

**Sections and content:**

1. Mixture components and gamma multipliers
2. Annotation coefficient matrix `alpha`
3. Generalized probit stick-breaking probabilities
4. Marker predictors \(A\alpha\)
5. Annotation- and marker-level expected gamma
6. Initialization with `make_sbayesrc_alpha_init()`
7. Fit outputs: `alpha`, `sigmaSqAlpha`, `comp_prob`, `ncomp`
8. Diagnostic helper contracts and compact summaries
9. Interpretation and limits

**Functions:** `stblr_csr_sbayesrc_generic()`, `make_sbayesrc_alpha_init()`,
`sbayesrc_annotation_pi()`, `sbayesrc_annotation_gamma_mean()`,
`sbayesrc_marker_pi()`, and `sbayesrc_marker_gamma_mean()`.

**Draft reuse:** Use the annotation draft's direct annotation-model and
mixture-prior concepts where applicable, then align the note with the current
generalized probit stick-breaking implementation and workflow diagnostics.

**Main-note links:** `annotation_informed_models.qmd` and
`workflows_and_outputs.qmd`.

## Workflow Scripts as Documentation Sources

The two principal workflow scripts should be maintained as executable
architecture documents:

- `sparse_ld_bed_workflow.R` explains how a common BED genotype source supports
  both direct individual-level fitting and a summary-statistics CSR fitting
  path.
- `annotation_based_models.R` explains how the CSR path is extended with an
  aligned annotation matrix and compared across fixed, learned, group, and
  SBayesRC-style priors.

Conceptual notes should explain stable definitions, alternatives, assumptions,
and interpretation. Workflow scripts should show the current exact call order,
object handoffs, and integration points. Main notes may reuse short excerpts,
but should link to the complete scripts, mark runtime-heavy steps, avoid
machine-specific paths, and never imply that a demonstration chain is
adequate for real analysis. Scripts are worked integrations, not substitutes
for package help or automated tests.

### The BED-to-CSR Bridge

`sparse_ld_bed_workflow.R` is the central bridge between package layers:

1. A BED-backed `Glist` supplies individual-level genotypes and marker
   metadata.
2. `stblr_bed_marker()` fits the individual-level model directly from that
   source.
3. `bed_xtx_xty()` derives sufficient statistics from the same selected rows
   and markers.
4. `sparseLD_stream_CSR()` derives a thresholded, disk-backed sparse LD
   representation from the same BED source.
5. `sparseLD_read_CSR()` exposes the stored representation for validation.
6. `stblr_csr()` combines the sufficient statistics and CSR prefix to fit the
   summary-statistics model.
7. Comparing posterior effects and inclusion probabilities from both paths
   demonstrates which inputs and outputs are comparable and where sparse-LD
   approximation or scaling choices may matter.

The new notes should present this as a data-flow diagram and a set of explicit
contracts, not as two unrelated model examples.

## Existing Material Mapping

No substantive section below is discarded. Repetition may be consolidated
only after its equations, assumptions, cautions, and implementation details
are preserved at the listed destinations.

### `examples/workflows/sparse_ld_bed_workflow.R`

| Existing workflow section | Main-note destination | Detailed destination | Documentation role |
|---|---|---|---|
| Data setup, `Glist`, chromosome and marker selection | `data_representations.qmd`, `sparse_ld_and_bed_workflows.qmd` | `technical_sparse_ld_csr.qmd` | Establish the shared BED source and marker-order contract. |
| Trait simulation and scaling | `sparse_ld_and_bed_workflows.qmd`, `workflows_and_outputs.qmd` | `technical_summary_statistics.qmd` | Explain demonstration inputs and why scaling must match both paths. |
| qgg comparison | `workflows_and_outputs.qmd` | None required | Retain as an external comparison point, not a package architecture dependency. |
| BED-derived sufficient statistics with `bed_xtx_xty()` | `data_representations.qmd`, `sparse_ld_and_bed_workflows.qmd` | `technical_summary_statistics.qmd`, `technical_sparse_ld_csr.qmd` | Document the transition from individual data to summary-statistic inputs. |
| Stream sparse LD with `sparseLD_stream_CSR()` | `sparse_ld_and_bed_workflows.qmd` | `technical_sparse_ld_csr.qmd` | Document disk-backed construction, thresholding, and blocking. |
| Inspect CSR with `sparseLD_read_CSR()` | `sparse_ld_and_bed_workflows.qmd` | `technical_sparse_ld_csr.qmd` | Show validation of the stored representation. |
| Fit `stblr_csr()` | `single_trait_models.qmd`, `sparse_ld_and_bed_workflows.qmd` | `technical_summary_statistics.qmd`, `technical_sparse_ld_csr.qmd` | Explain the summary-statistics CSR path. |
| Fit `stblr_bed_marker()` | `single_trait_models.qmd`, `sparse_ld_and_bed_workflows.qmd` | None beyond model appendices | Explain the direct individual-level BED path. |
| Compare CSR and BED outputs | `sparse_ld_and_bed_workflows.qmd`, `workflows_and_outputs.qmd` | `technical_sparse_ld_csr.qmd` | Explain comparable outputs, approximation differences, and validation limits. |

### `examples/workflows/annotation_based_models.R`

| Existing workflow section | Main-note destination | Detailed destination | Documentation role |
|---|---|---|---|
| BED data setup, marker matching, and CSR prefix | `data_representations.qmd`, `sparse_ld_and_bed_workflows.qmd` | `technical_sparse_ld_csr.qmd` | Reuse the BED-to-CSR contract before adding annotations. |
| `mtsim_annotation()` and `summarize_annotation_signal()` | `annotation_informed_models.qmd`, `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` | Define simulated annotation signal and compact truth-aware diagnostics. |
| `bed_xtx_xty()` and `sparseLD_stream_CSR()` | `sparse_ld_and_bed_workflows.qmd` | `technical_summary_statistics.qmd`, `technical_sparse_ld_csr.qmd` | Link annotation models to the same CSR foundation. |
| Annotation matrix marker-order contract | `data_representations.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` | Make annotation alignment a required input invariant. |
| Fixed-prior model | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` | Demonstrate externally supplied marker-specific priors. |
| Learned annotation model | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` | Demonstrate direct annotation learning within the sampler. |
| Group annotation model | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` | Demonstrate mutually exclusive group priors. |
| SBayesRC-style model | `annotation_informed_models.qmd` | `technical_sbayesrc_annotations.qmd` | Demonstrate overlapping-annotation mixture priors. |
| Compact fit, top-marker, and annotation summaries | `workflows_and_outputs.qmd` | Relevant model appendices | Establish readable output patterns while retaining full matrices in fit objects. |
| Expected-gamma diagnostics | `annotation_informed_models.qmd`, `workflows_and_outputs.qmd` | `technical_sbayesrc_annotations.qmd` | Explain annotation- and marker-level prior architecture diagnostics. |

### `docs/sblr_theory.qmd`

| Existing section | Main-note destination | Detailed destination | Treatment |
|---|---|---|---|
| Theory: Constructing Cross-Trait Summary-Statistic Matrices; Overview | `data_representations.qmd`, `multi_trait_models.qmd` | `technical_multitrait_overlap.qmd` | Main notes summarize motivation; appendix preserves definitions and equations. |
| Phenotypic Variance From Summary Statistics | `data_representations.qmd` | `technical_multitrait_overlap.qmd` | Preserve variance approximation and centering assumption. |
| Phenotypic Covariance From Correlation Information | `multi_trait_models.qmd` | `technical_multitrait_overlap.qmd` | Preserve rescaling from correlations to covariance units. |
| The Sample-Overlap Matrix; Full-overlap assumption; Fractional-overlap assumption | `data_representations.qmd`, `multi_trait_models.qmd` | `technical_multitrait_overlap.qmd` | Preserve both overlap models and diagonal constraints. |
| Reconstructing the Cross-Trait SSY Matrix | `multi_trait_models.qmd` | `technical_multitrait_overlap.qmd` | Preserve \(\widehat{SSY}=N\odot\widehat{\operatorname{Cov}}(y)\). |
| Why SSY and Nmat Are Needed in Multivariate BLR | `multi_trait_models.qmd` | `technical_multitrait_overlap.qmd` | Preserve residual-covariance rationale and diagonal-only limitation. |
| Summary | `index.qmd` and cross-links | `technical_multitrait_overlap.qmd` | Convert to navigation/recap; no technical claims removed. |
| 5.1 Deriving SSY from first principles | brief link from `multi_trait_models.qmd` | `technical_multitrait_overlap.qmd` | Preserve full derivation, random-missingness and homogeneity assumptions. |

### `docs/annotation_informed_priors_stblr_mtblr.qmd`: posterior summaries and empirical-Bayes foundations

| Existing section group | Main-note destination | Detailed destination |
|---|---|---|
| Motivation | `index.qmd`, `model_overview.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Summary-statistic residual notation | `model_overview.qmd`, `data_representations.qmd` | `technical_summary_statistics.qmd` |
| Set-specific genomic variance decomposition; R implementation | `annotation_informed_models.qmd` | `technical_summary_statistics.qmd` |
| Inclusion-probability enrichment; R implementation | `annotation_informed_models.qmd`, `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |
| Combining effect size and inclusion probability; R implementation | `annotation_informed_models.qmd` | `technical_summary_statistics.qmd` |
| Handling overlapping annotations; linear, weighted, and fractional-response implementations | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Using STBLR or MTBLR instead of ridge regression | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Deriving annotation-specific prior inclusion probabilities; implementation | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Deriving annotation-specific prior effect variances; implementation | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Marker-specific prior scales for overlapping annotations; implementation | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |

### `docs/annotation_informed_priors_stblr_mtblr.qmd`: STBLR, MTBLR, and staged workflows

| Existing section group | Main-note destination | Detailed destination |
|---|---|---|
| Using marker-specific prior variances in STBLR | `single_trait_models.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Annotation-specific MTBLR prior covariance matrices; implementation | `multi_trait_models.qmd` | `technical_annotation_priors.qmd` |
| Annotation-specific MTBLR model priors; implementation | `multi_trait_models.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Prior scale derivation from annotation BLR | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Recommended modelling hierarchy, Stages 1-4 | `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |
| Warm-started regional MTBLR | `workflows_and_outputs.qmd` as proposed workflow | `technical_annotation_priors.qmd` |
| Practical cautions | `annotation_informed_models.qmd`, `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |
| Minimal example workflow | `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |
| Summary | `annotation_informed_models.qmd` | Cross-links to all relevant appendices |

### `docs/annotation_informed_priors_stblr_mtblr.qmd`: empirical-Bayes annotation-level BLR

| Existing section group | Main-note destination | Detailed destination |
|---|---|---|
| Empirical-Bayes priors from an annotation-level STBLR or MTBLR | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| First-stage genomic fit | `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |
| Annotation-level STBLR | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Annotation-level MTBLR | `multi_trait_models.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Deriving marker-specific inclusion priors | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Deriving marker-specific prior effect variances | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Deriving MTBLR marker-specific sharing-model priors | `multi_trait_models.qmd` | `technical_annotation_priors.qmd` |
| Deriving a global MTBLR prior covariance matrix from annotations | `multi_trait_models.qmd` | `technical_annotation_priors.qmd` |
| Recommended empirical-Bayes workflow and R sketch | `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |
| Interpretation | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Cross-fitting and overfitting | `annotation_informed_models.qmd`, `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |

### `docs/annotation_informed_priors_stblr_mtblr.qmd`: direct annotation models and implementation

| Existing section group | Main-note destination | Detailed destination |
|---|---|---|
| Direct annotation modelling inside the genomic STBLR or MTBLR sampler | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Direct annotation-informed STBLR | `single_trait_models.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Modified STBLR marker update and C++ sketch | brief equation in `single_trait_models.qmd` | `technical_annotation_priors.qmd` |
| Estimating annotation effects inside the sampler; fixed effects; updating effects | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Direct annotation-informed MTBLR | `multi_trait_models.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| MTBLR model prior in the marker update and C++ sketch | `multi_trait_models.qmd` | `technical_annotation_priors.qmd` |
| Annotation-informed MTBLR effect covariance | `multi_trait_models.qmd` | `technical_annotation_priors.qmd` |
| Direct annotation-aware STBLR and MTBLR as hierarchical models | `model_overview.qmd`, `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Practical implementation strategy, Steps 1-3 | `annotation_informed_models.qmd`, `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |
| Sparse annotation matrix implementation and C++ sketches | mention as scalability consideration | `technical_annotation_priors.qmd` |
| Centering, shrinkage, and caps | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Recommended first direct model | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |

### `docs/annotation_informed_priors_stblr_mtblr.qmd`: interpretation and causal hypotheses

| Existing section group | Main-note destination | Detailed destination |
|---|---|---|
| Causal interpretation and causal-hypothesis generation | `annotation_informed_models.qmd`, `multi_trait_models.qmd` | `technical_annotation_priors.qmd` |
| Genetic sharing versus causal direction | `multi_trait_models.qmd` | `technical_annotation_priors.qmd` |
| What can be inferred directly | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Connection to colocalization | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Connection to Mendelian randomization | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Annotation-specific causal hypotheses | `annotation_informed_models.qmd` | `technical_annotation_priors.qmd` |
| Practical causal workflow and manuscript statement | `workflows_and_outputs.qmd` | `technical_annotation_priors.qmd` |

## How Technical Depth Is Preserved

The main notes should use a three-layer pattern for every technical topic:

1. **Definition and decision:** state the model, why it is used, and which
   package interface implements it.
2. **Interpretation and caveats:** explain outputs, assumptions, and limits in
   the main note.
3. **Full derivation:** link directly to a technical appendix section that
   preserves all equations, derivation steps, R implementations, and
   implementation sketches.

Dense derivations should not be replaced by prose summaries. Instead, the main
note should show the defining equation and a short interpretation, followed by
a prominent link such as:

> For the full marker-update derivation, assumptions, and implementation
> sketch, see `technical_annotation_priors.qmd#marker-specific-stblr-update`.

Each migrated appendix section should carry a source comment or editorial note
identifying the original draft section, making completeness review possible.
Before retiring any source draft from active use, compare its heading list
against the mapping tables above and confirm every equation/code block is
either retained or explicitly identified as obsolete.

Operational depth from the workflow scripts should be preserved with the same
discipline. Each migrated workflow stage should retain its object handoff,
ordering/scaling assumptions, validation purpose, and relationship to the
next package layer. Main notes should explain the path and show compact
fragments; appendices should preserve storage and algorithmic contracts; the
full scripts should remain linked as executable references.

## Unresolved Questions

1. Should `docs/notes/` become a standalone Quarto website/book, or remain a
   linked collection of notes without a site configuration?
2. Should the existing rendered `docs/sblr_theory.html` and
   `docs/sblr_theory_files/` remain archived, or be removed only after a new
   rendered site exists?
3. Several functions named in `sblr_theory.qmd`, such as
   `build_SSY_full_overlap()`, `build_Nmat_full()`, and `build_Nmat_frac()`,
   are not currently public package wrappers. Should they be documented as
   conceptual helpers, internal/proposed functions, or future API work?
4. Which MTBLR components are considered stable public functionality versus
   research-stage internals? The notes must label this consistently.
5. Should empirical-Bayes set-summary helpers in `R/settest.R` become public,
   or remain technical examples?
6. Should `annotation_informed_models.qmd` describe only implemented wrappers
   in its main path, with all proposed MTBLR annotation models moved to a
   clearly marked "Research extensions" section?
7. What convergence diagnostics and minimum-chain guidance can be stated
   concretely for the current samplers?
8. Should examples use simulated package-internal data, or continue pointing
   to external qgg/PLINK data preparation?
9. The annotation draft contains exploratory C++ snippets. Should they remain
   as design documentation, or move later to a separate developer note?
10. Some draft text contains encoding artifacts. These should be corrected
    during migration without altering mathematical meaning.
11. Should `sparse_ld_bed_workflow.R` be made portable enough to run without
    local qgg/PLINK paths, or remain an explicitly environment-dependent
    integration workflow?
12. Which parts of the disk-backed CSR file contract are stable enough to
    document as user-facing guarantees, and which should be labeled as current
    implementation details?
13. What level of agreement between `stblr_bed_marker()` and `stblr_csr()` is
    expected after LD thresholding, scaling choices, and finite MCMC error?
14. Should the workflow notes define a small shared validation checklist for
    marker order, sample rows, allele frequencies, and annotation alignment?

## Recommended Writing Order

1. **Archive and inventory:** create `docs/notes/drafts/` copies and produce a
   heading/equation/code-block inventory before rewriting.
2. **`model_overview.qmd`:** establish terminology and the four-axis
   classification used everywhere else.
3. **`data_representations.qmd`:** define BED, summary-statistic, and CSR inputs
   before model-specific workflows.
4. **`sparse_ld_and_bed_workflows.qmd`:** document the central BED-direct and
   BED-to-CSR bridge from the current workflow script.
5. **`technical_sparse_ld_csr.qmd` and `technical_summary_statistics.qmd`:**
   preserve the storage, alignment, and sufficient-statistic contracts needed
   by that bridge.
6. **`technical_multitrait_overlap.qmd`:** migrate the focused theory draft
   nearly intact; this is the smallest complete migration and establishes the
   appendix pattern.
7. **`single_trait_models.qmd` and `multi_trait_models.qmd`:** describe the
   baseline models before annotation extensions.
8. **`technical_annotation_priors.qmd`:** migrate the long annotation draft by
   mapped section group, preserving equations and code.
9. **`technical_sbayesrc_annotations.qmd`:** align the current public helpers
   and wrapper with the generalized stick-breaking mathematics.
10. **`annotation_informed_models.qmd`:** write the readable comparison and
   decision guide after the technical appendices are stable.
11. **`workflows_and_outputs.qmd`:** connect both primary workflow scripts to
   current public interfaces and compact output summaries.
12. **`index.qmd`:** write last so its navigation reflects the completed note
    set.
13. **Completeness review:** compare both original draft heading lists and
    both workflow mapping tables, equations, code blocks, contracts, and
    cautions against the new notes before changing the status of the original
    drafts.
