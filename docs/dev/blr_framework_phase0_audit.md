# Unified BLR Framework: Phase 0 Repository Audit

This report audits the local `sblr` repository at commit `1b3e8d6`. It is a
repository-grounded description of current behavior, not an implementation
specification. Statements labelled **Observed** are supported by the cited
local files and functions. Statements labelled **Recommendation** are proposed
migration actions and do not describe code that already exists.

## 1. Executive summary

The proposed unified framework is compatible with the strongest parts of the
current repository, but it should be introduced incrementally around the
validated single-trait (ST) paths rather than by generalizing the legacy
multi-trait sampler. The best reusable infrastructure is:

- the named `stblr_raw` version 1 contract and canonical `.as_stblr_fit()`
  formatter in `R/sparse_ld_bed_helper.R`;
- marker/trait alignment and prior preparation in the public ST wrappers;
- the `CsrOperator`/`BlockEigenOperator` interface in
  `src/st_ld_operator.h`;
- packed-BED decoding and score/update helpers in `src/st_bed_decode.h` and
  `src/st_bed_bayesr_common.h`;
- chain-task seed helpers in `src/st_chain_utils.h` and the chain-local RNG
  designs in the current CSR BayesR/SBayesRC and BED BayesR/BayesRC backends;
- the likelihood-independent BayesRC stick-breaking utilities in
  `src/st_bayesrc_annotation_prior.h`; and
- the typed, standard-C++ sufficient-statistic preparation example in
  `src/mt_sufficient_statistics_core.h`.

The highest-risk legacy area is the public multi-trait `sblr()` route. Its
default native implementation in `src/mtblr.cpp` does not implement one
coherent covariance sampler: it performs multiple different `B` updates in a
single iteration, ignores `updateB` during per-set covariance draws, updates
only the diagonal of `E`, and uses only diagonal residual precision in the
active marker kernel. The alternative CSR multi-trait implementation in
`src/mt_cpg_omp_csr.cpp` has the more useful full-`YY`, shared-design residual
cross-product calculation, but its partial-pattern marker update does not use
the full residual score in the active subspace and its covariance prior under
partially active patterns needs a formal model decision. Neither path should be
declared canonical without new derivations and reduction tests.

There is also a reproducibility risk outside the multi-trait code. The active
scheduled CSR BayesC and two active scheduled individual/BED BayesC kernels
retain `static thread_local std::normal_distribution` objects. A normal
distribution may cache a variate, so reseeding only `std::mt19937` does not
reset all random state and job-to-thread assignment can affect results.

**Compatibility assessment.** The target separation into data operator, state
space, probability policy, marker-scale policy, covariance policy, MCMC
control, and result aggregation matches identifiable seams in the repository.
It is not yet a drop-in refactor because most marker kernels fuse those
policies and because the legacy multi-trait statistical contract is
unresolved.

**Recommended design-document revisions.** Keep the policy composition and
typed boundary goals, but:

1. make covariance selection and the meaning of a partially active
   multivariate Gaussian prior explicit prerequisites rather than implementation
   details;
2. require a full `Y'Y`/sample-overlap contract before a full residual
   covariance model can be selected;
3. treat the current `sblr()` algorithms as reference/experimental routes, not
   as the canonical starting implementation;
4. add an explicit scheduled-BayesC RNG remediation gate before claiming
   family-wide core-count invariance; and
5. defer a raw-schema version change until typed multi-trait orientations and
   pattern outputs are settled.

**Recommended first implementation phase.** Introduce an internal resolved R
specification and language-neutral typed C++ specification/result vocabulary
for one production reference route: unscheduled CSR BayesC with unit marker
scale, independent traits, and the existing scalar variance updates. The
resolved specification should be expanded back into the exact existing native
arguments, so the sampler, RNG trajectory, raw schema, public API, and
formatted output remain unchanged. Do not migrate multi-trait sampling in that
phase.

## 2. Repository baseline

### Local state

The initial commands were:

```text
git status --short
git branch --show-current
git log -5 --oneline
```

Observed baseline:

- branch: `master`;
- initial working tree: clean (`git status --short` produced no output);
- latest commits:
  - `1b3e8d6 FInal update before refactoring`
  - `adac7c6 Update sparse_ld_bed_helper.R`
  - `56be6dc Update sparse_ld_bed_helper.R`
  - `89c7e57 Sequnece 3`
  - `6f0763d updates sequence 3`
- no branch was created, and no commit or push was performed.

The three authoritative documents were read in full before implementation
evaluation:

- `docs/dev/blr_framework_implementation_plan.md`;
- `docs/dev/blr_model_capability_matrix.md`; and
- `docs/dev/blr_reduction_test_matrix.md`.

Relevant existing notes inspected include
`docs/dev/stblr_backend_naming.md`, `docs/dev/stblr_raw_schema.md`,
`docs/dev/stblr_annotation_backend_design.md`,
`docs/dev/individual_level_bayesrc_plan.md`,
`docs/dev/stblr_backend_computation_inventory.md`,
`docs/dev/stblr_block_eigen_operator.md`,
`docs/dev/archive/stblr_csr_multichain_design.md`, and the user-facing notes on
multi-trait models, multi-trait overlap, data representations, CSR, BED, and
annotation-informed models.

### Validation baseline

The environment is R 4.4.1. `pkgload` 1.4.1, `testthat` 3.3.2, and `devtools`
2.4.6 are installed. `testthat` reports that it was built under R 4.4.3.

The attempted compiled focused run was:

```text
Rscript -e "pkgload::load_all('.', compile=TRUE); testthat::set_max_fails(Inf); testthat::test_dir('tests/testthat', filter='individual-bayesrc|stblr-raw-schema|stblr-annotation|stblr-selection-s|backend-consistency|bayesr-bed-backend|bayesr-csr-backend', reporter='summary')"
```

It stopped before tests because R could not find Rtools. There was no existing
DLL or object file in `src/` to load. This is an environment limitation, not a
test failure.

Read-only fallback checks used `pkgload::load_all('.', compile=FALSE)`:

| Check | Test blocks | Passed expectations | Failures | Errors | Skips | Interpretation |
|---|---:|---:|---:|---:|---:|---|
| `test-stblr-raw-schema.R` | 7 | 81 | 0 | 0 | 0 | R-only schema/formatter baseline passes |
| `test-individual-bayesrc.R` | 12 | 27 | 0 | 0 | 8 | R-level checks pass; native tests skip without DLL |
| focused capability set | 161 | 465 | 0 | 79 | 16 | errors occur when tests call unavailable registered native symbols |
| complete `tests/testthat` directory | 238 | 725 | 0 | 105 | 16 | not a package pass; compile-free environmental baseline only |

The focused set included raw schema, annotation models and chains, BED
BayesRC, BayesR BED/CSR, `selection_s`, block-eigen, BED interface, and backend
consistency tests. A representative error is `object
'_sblr_stblr_cpg_omp_csr_bayesr' not found`, confirming the missing DLL cause.
No automated multi-trait BLR test file exists.

Additional read-only checks:

```text
Rscript -e 'parse every R/*.R file'
git diff --check
```

All 23 R source files parsed with zero errors. The initial `git diff --check`
passed. No package build/check was attempted after the compiler prerequisite
failed.

## 3. Public interface map

| Interface | Source/status | Families and representation | Alignment and initialization | Dispatch and formatting | Migration constraints / duplication |
|---|---|---|---|---|---|
| `sblr()` | `R/interface_mtblr.R:sblr`; exported | Experimental summary-statistic single/multi-trait BayesC only; dense/list `XX`, `Xy`, `yy` | Trait count is `length(Xy)`; marker count is mean vector length; names are taken from `names(yy)` and `names(XXvalues[[1]])`; defaults enumerate all `2^T` patterns and derive `B`, `E`, and scale matrices from `h2` | Direct `.Call` to `mtblr`, `mtblr_cpg`, `mtblr_cpg_arma`, `mtblr_cpg_omp`, or `mtblr_eigen`; manually formats 20 positional slots | Preserve the callable algorithms and 18-field BayesC fit while migrating. It duplicates setup, dispatch, and formatting and does not use `stblr_raw` |
| `stblr_csr()` | `R/sparse_ld_bed_helper.R:stblr_csr`; exported | ST/parallel-trait BayesC and BayesR on CSR LD | `.prepare_csr_stats()`-style checks align marker rows and trait columns; Glist/LD prefix metadata supply order; priors are initialized per trait | BayesC chooses scheduled or unscheduled native backend; BayesR delegates to `stblr_csr_bayesr()`; named raw output goes to `.as_stblr_fit()` | Strongest insertion point for a resolved spec: after input/alignment resolution and before native dispatch |
| `stblr_csr_bayesr()` | `R/sparse_ld_bed_helper.R:stblr_csr_bayesr`; exported | ST/parallel-trait BayesR on CSR LD | Validates `mixture_var`, global `pi`, Dirichlet `alpha`, initial components, residual state, chain seeds, and optional `selection_s` scale | `stblr_cpg_omp_csr_bayesr`; `.as_stblr_fit()` | BayesR setup is partly duplicated in BED and block-eigen wrappers; preserve component order and `component_0` semantics |
| `stblr_bed()` | `R/sparse_ld_bed_helper.R:stblr_bed`; exported | ST/parallel-trait BayesC, BayesR, BayesRC on packed PLINK BED | `.make_bed_marker_data()` aligns phenotype row names to `Glist$ids`, resolves `chr`/`cls`, concatenates `Glist$rsids[[chr]][cls]`, and aligns allele frequencies | Separate model branches call packed-BED backends, then `.as_stblr_fit()` | Public default remains BayesC. Common BED preparation should remain shared; model branch arguments can be represented by a resolved spec |
| `stblr_bed_marker()` | `R/sparse_ld_bed_helper.R`; exported lower-level helper | BayesC packed-BED marker interface | Accepts already selected files/columns/rows/AF and prepared priors | Chooses sparse, scheduled, or scheduled-chain BED backend | Important compatibility route and test fixture; not the initial unified public API |
| `stblr_csr_annot()` | `R/stblr-csr-annot.R:stblr_csr_annot`; exported facade | Fixed prior, learned annotation, group annotation, or SBayesRC on CSR | Checks annotation model and forwards model-specific arguments | Dispatches to the four explicit wrappers and standardizes annotation aliases | Good R-level probability/scale-policy resolver seam; do not collapse distinct posterior models into one generic coefficient object |
| `stblr_csr_prior_annot()` | `R/stblr-csr-prior-annot.R`; exported | BayesC with fixed marker inclusion probabilities and/or fixed variance multipliers | `.stblr_prepare_annotation_matrix()` and `.stblr_make_prior_from_annotations()` align/construct marker priors | `stblr_cpg_omp_csr_prior`; canonical formatter | Fixed probability and scale are independent policies despite sharing one backend |
| `stblr_csr_learn_annot()` | `R/stblr-csr-learn-annot.R`; exported | BayesC with learned logit inclusion and/or learned log-variance annotation effects | Shared annotation preprocessing; initializes `eta_pi`, `eta_vb` and MH controls | `stblr_cpg_omp_csr_annot`; canonical formatter plus aliases | The learned probability and learned scale updates are model-specific MH policies |
| `stblr_csr_group_annot()` | `R/stblr-csr-group-annot.R`; exported | BayesC with categorical group inclusion probabilities and group variance multipliers | Converts group labels to 0-based indices and initializes group-by-trait parameters | `stblr_cpg_omp_csr_group_annot`; canonical formatter | One categorical layer exists; it is not an implementation of arbitrary overlapping hierarchy |
| `stblr_csr_sbayesrc_generic()` | `R/stblr-csr-sbayesrc.R`; exported | SBayesRC mixture with annotation-dependent probit stick probabilities on CSR | Shared annotation preprocessing adds intercept by default; initializes `alpha`, `sigmaSqAlpha`, gamma, chain state, and optional `selection_s` | `stblr_cpg_omp_csr_sbayesrc`; canonical formatter | Source of truth for stick-breaking and annotation output semantics |
| block-eigen wrappers | `R/sparse_ld_bed_helper.R:.stblr_csr_bayesc_block_eigen`, `.stblr_csr_bayesr_block_eigen`; `R/stblr-csr-sbayesrc.R:.stblr_csr_sbayesrc_block_eigen`; internal | BayesC, BayesR, SBayesRC with a BED-derived block-eigen approximation | Validate contiguous block starts, selected BED markers, eigen filter and transformed score | Native `_block_eigen` variants; canonical formatter | They already demonstrate one sampler with two operator implementations; remain internal during first framework phase |
| `finemap_stblr_csr()` | `R/finemap-stblr-csr.R`; exported | Local CSR BayesC re-fits around selected loci | Aligns fit marker/trait names with stats/Glist, writes subset CSR, prepares one-chain priors | Calls `stblr_cpg_omp_csr`; canonical raw formatter | Preserve locus/result contracts; use a resolved spec only after base CSR equivalence is protected |
| `mtsim()` | `R/mtsim.R:mtsim`; exported | Multi-trait individual-level simulation | Uses a common genotype matrix and individual set; aligns supplied matrices and trait covariance targets | Pure R simulation | Reusable known-truth fixture; not evidence that a native joint sampler is correct |
| `mtsim_annotation()` / `summarize_annotation_signal()` | `R/mtsim.R`; exports documented in `R/mtsim-annotation.R` | Multi-trait simulation with overlapping annotations and marker priors | Marker rows follow simulated/supplied W; annotation sets and enriched traits are explicit | Pure R | Reusable fixture for future policy tests; current fitted workflows are chiefly parallel ST |

**Resolved-spec insertion point.** For the ST APIs, resolution belongs after
all marker/trait/annotation/data alignment and before construction of native
argument lists. Initially it should be an internal immutable list/class whose
legacy adapter reproduces the exact existing arguments. For `sblr()`, a
separate compatibility adapter is needed later because its data and output
contracts are positional and statistically unresolved.

## 4. Native backend map

All active ST entries below return named `stblr_raw` version 1 and are formatted
by `.as_stblr_fit()`. Legacy multi-trait entries return 20 positional slots and
are formatted manually by `sblr()` unless noted.

| Backend / source | Exported native function | Data/operator | Model and covariance | Chains, OpenMP, RNG | Callers/tests/status |
|---|---|---|---|---|---|
| `src/mtblr.cpp` | `mtblr` | dense/list summary `XX` | joint pattern BayesC; full latent `B`, effectively diagonal-updated `E`; conflicting covariance updates | one chain; serial; `mt19937`, with old R/Armadillo RNG helpers in same file | public `sblr(algorithm="default")`; no automated tests; experimental/legacy production route |
| `src/mtblr.cpp` | `mtblr_hybrid` | dense/list summary | alternative masked/hybrid pattern sampler | one chain; serial; mixed legacy helpers | generated native callable, not selected by `sblr()`; experimental |
| `src/mtblr.cpp` | `mtblr_eigen` | rotated/eigen summary | pattern BayesC; diagonal covariance traces; alternative marker update | one chain; serial; marker draws use `mvrnormARMA()`/Armadillo global RNG despite a supplied `mt19937` seed | public `sblr(algorithm="eigen")`; no automated tests; experimental and migration-required |
| `src/mt_cpg.cpp` | `mtblr_cpg` | dense/list summary | legacy CPG multi-trait patterns | one chain; serial; `mt19937` | public experimental algorithm; no automated tests; reference/legacy |
| `src/mt_cpg_arma.cpp` | `mtblr_cpg_arma` | Armadillo dense/list summary | alternative CPG pattern implementation | one chain; serial; `mt19937` | public experimental algorithm; no automated tests; reference/duplicated |
| `src/mt_cpg_omp.cpp` | `mtblr_cpg_omp` | sparse-list summary with graph-coloured independent sets | latent pattern BayesC | one chain; parallel within colours; per-iteration/thread seeds depend on thread number | public experimental algorithm; no tests; core-count-dependent by construction |
| `src/mt_cpg_omp_csr.cpp` | `mtblr_cpg_omp_csr` | one shared disk CSR LD | joint pattern BayesC with full `YY`, full IW `E`, full `B` attempt | one serial chain despite name; `mt19937`; no chain aggregation | called directly by exploratory examples, not `sblr()`; no automated tests; strongest full-`YY` reference but statistically unresolved |
| `src/st_cpg_omp_csr.cpp` | `stblr_cpg_omp_csr` | `CsrOperator` | binary BayesC, scalar trait-specific variances; optional LD-swap and `selection_s` | trait-chain tasks, static OpenMP, chain-local `mt19937`, explicit chain seeds | public CSR/finemap; extensive tests; production reference |
| `src/st_cpg_omp_csr.cpp` | `stblr_cpg_omp_csr_block_eigen` | `BlockEigenOperator` | same BayesC kernel, approximate likelihood operator | same chain infrastructure | internal block-eigen tests; experimental operator behind production kernel |
| `src/st_cpg_omp_csr_scheduled.cpp` | `stblr_cpg_omp_csr_scheduled` | CSR | scheduled binary BayesC | trait-chain tasks; static OpenMP; `mt19937` but persistent thread-local normal/uniform distributions | public scheduled CSR; tested functionally; production with reproducibility remediation required |
| `src/st_cpg_omp_csr_bayesr.cpp` | `stblr_cpg_omp_csr_bayesr` | CSR operator helpers | BayesR components, global Dirichlet pi, scalar variance, optional LD-swap/`selection_s` | trait-chain static OpenMP; chain-local engine/distribution state | public; reduction/component/RNG tests; production |
| `src/st_cpg_omp_csr_bayesr.cpp` | `stblr_cpg_omp_csr_bayesr_block_eigen` | block-eigen | same BayesR policy | same chain design | internal tests; experimental operator variant |
| `src/st_sbayesrc_omp_csr.cpp` | `stblr_cpg_omp_csr_sbayesrc` | CSR | BayesRC components; probit stick annotation probability; scalar component scale | trait-chain static OpenMP; chain-local engine; shared header calls R normal CDF/quantile | public annotation route; extensive tests; production, mixed core/binding |
| `src/st_sbayesrc_omp_csr.cpp` | `stblr_cpg_omp_csr_sbayesrc_block_eigen` | block-eigen | same SBayesRC policy | same | internal tests; experimental operator variant |
| `src/st_cpg_omp_csr_prior.cpp` | `stblr_cpg_omp_csr_prior` | CSR | binary BayesC with fixed marker probability/scale | trait-chain static OpenMP; chain-local `mt19937` | public annotation route; tests; production |
| `src/st_cpg_omp_csr_annot.cpp` | `stblr_cpg_omp_csr_annot` | CSR | binary BayesC with learned logit probability and/or learned log scale | trait-chain static OpenMP; chain-local `mt19937`, local MH distributions | public annotation route; chain/LD-swap tests; experimental statistical extension but active |
| `src/st_cpg_omp_csr_group.cpp` | `stblr_cpg_omp_csr_group_annot` | CSR | group Beta inclusion and group inverse-chi-square scale | trait-parallel OpenMP; chain wrapper; `mt19937` | public annotation route; tests; active experimental extension |
| `src/st_cpg_omp_individual.cpp` | `stblr_cpg_omp_bed_marker_sparse` | packed BED, individual residual | binary BayesC with sparse marker scheduling | trait parallelism; `mt19937`; mixed decode/binding | lower-level BED route; BED tests; retained implementation/reference |
| `src/st_cpg_omp_individual_scheduled.cpp` | `stblr_cpg_omp_bed_marker_scheduled` | packed BED | scheduled binary BayesC | trait OpenMP; persistent thread-local distribution state | public/lower-level BED; tests; production with RNG remediation required |
| `src/st_cpg_omp_individual_scheduled_chains.cpp` | `stblr_cpg_omp_bed_marker_scheduled_chains` | packed BED | scheduled binary BayesC | trait-chain static OpenMP; `mt19937`; persistent thread-local distribution state | public `stblr_bed()` main BayesC route; tests; production with RNG remediation required |
| `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp` | `stblr_cpg_omp_bed_marker_scheduled_chains_bayesr` | shared packed-BED helpers | BayesR components with optional adaptive scheduling | trait-chain static OpenMP; chain-local engine and distributions | public BED BayesR; exact reduction/RNG tests; production |
| `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp` | `stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc` | shared packed-BED helpers | full-sweep BayesRC probit stick components | trait-chain static OpenMP; chain-local engine and distributions | public BED BayesRC; exact BayesR reduction/RNG/alignment tests; production |

Sampler/helper headers are part of the inventory even though they are not
separate exported backends:

- `src/st_bed_bayesr_common.h`: packed-BED decode/standardization, residual
  scores, variance/CPO helpers, and BayesR/BayesRC chain result structures;
- `src/st_bayesrc_annotation_prior.h`: stick probabilities, latent probit and
  annotation-coefficient/sigma updates shared by CSR and BED;
- `src/st_csr_common.h`: CSR loading, residual and `selection_s` utilities;
- `src/st_ld_operator.h`: common CSR/block-eigen operator interface;
- `src/st_block_eigen.{h,cpp}`: block eigen construction/application; and
- `src/cpg_samplers.{h,cpp}` and `src/distributions.{h,cpp}`: legacy
  multivariate samplers and Wishart/random utilities.

`src/mt_sufficient_statistics.cpp` and
`src/mt_sufficient_statistics_core.h` prepare BED-based sufficient statistics;
they are data helpers, not BLR samplers.

## 5. Current statistical model map

### Single-trait/parallel-trait models

The active ST backends fit traits independently, optionally in parallel. For
trait `t`, the core likelihood is represented either by individual residuals
`e_t = y_t - X b_t` (BED) or score residuals
`r_t = X'y_t - X'X b_t` (CSR/block-eigen). There is no sampled cross-trait
covariance in these paths; `covb`, `covg`, and `cove` are diagonal matrices
assembled from trait-specific summaries.

Current prior families are:

- BayesC: `d_j ~ Bernoulli(pi)` and
  `b_j | d_j=1 ~ N(0, vb*q_j)`;
- BayesR: component `k` with `gamma_0=0` and
  `b_j | k>0 ~ N(0, vb*gamma_k*q_j)`;
- SBayesRC/BayesRC: the same component scale, with marker-specific component
  probabilities generated by probit stick-breaking from annotation
  coefficients;
- fixed-prior annotation BayesC: marker-specific `pi_j` and/or `q_j`;
- learned annotation BayesC: centered logit probability and centered log-scale
  predictors updated by random-walk MH; and
- group BayesC: group-specific inclusion probabilities and group variance
  multipliers.

### Legacy joint multi-trait model

In the default `mtblr()` path, marker `j` has a latent vector
`beta_j in R^T`, a binary trait pattern `d_j`, and an effective effect
`b_j = D_j beta_j`, where `D_j=diag(d_j)`. Pattern probabilities are global
`pi_k`, and the supplied/default patterns are categorical. The active native
code treats `beta_j ~ N(0,B)` and updates marker score residuals trait by
trait. The residual input is a covariance matrix `E`, but the active marker
kernel and active residual update reduce its statistical role to diagonal
precision/variance.

`sblr()` currently stops unless `method == "bayesC"`
(`R/interface_mtblr.R:48`); the subsequent BayesR numeric mapping is
unreachable. Therefore user documentation or examples suggesting a public
multi-trait BayesR route are not supported by this wrapper.

## 6. Multivariate BLR audit

### Public setup and state space

**Observed.** `R/interface_mtblr.R:sblr()` determines `T` from `length(Xy)` and
when `models` is `NULL` creates every binary pattern with `expand.grid`
(`R/interface_mtblr.R:55-58`). There is no maximum-pattern or maximum-trait
guard. `models="restrictive"` keeps only all-null and all-active. Initial
pattern mass is `(1-pi, pi/(K-1),...)`. Trait and marker dimension/name
validation is substantially weaker than in the ST wrappers.

The documented and implemented `yy` contracts conflict. The public
documentation describes a trait cross-product matrix, but `sblr()` passes
`yy` to native default/CPG backends whose generated signatures take
`std::vector<double>`, initializes variances with `diag(yy/(n-1), nt)`, and,
when names are absent, creates trait names with `seq_len(length(yy))`
(`R/interface_mtblr.R:67,80-88`). For an unnamed `T x T` matrix,
`length(yy)=T^2`, so assigning those names to `T x T` covariance matrices is
dimensionally inconsistent. The separate `mtblr_cpg_omp_csr()` export accepts
a full matrix, but `sblr()` does not dispatch to it. Therefore the public
wrapper cannot currently be treated as a validated full-`Y'Y` interface.

`sblr()` has no seed argument. It samples an internal seed from R's global RNG
at `R/interface_mtblr.R:113`, so a caller cannot directly supply a native seed.

### Default marker update

The active update is `sampleBetaCPG_Mt_latent()` at
`src/mtblr.cpp:2978-3175`, called at `src/mtblr.cpp:4227`. For each pattern
`k`, it computes, in current notation,

```text
rhs_t = E^{-1}_{tt} (r_tj + w_tj b_tj)   for active t, otherwise 0
C_k   = B^{-1} + diag(d_k * w_j * diag(E^{-1}))
log w_k = log(pi_k) - 1/2 log|C_k| + 1/2 rhs' C_k^{-1} rhs.
```

It then draws a full latent `beta_j ~ N(C_k^{-1}rhs, C_k^{-1})` and stores
`b_j=D_k beta_j`.

Consequences:

- off-diagonal residual precision is not used in pattern probabilities or
  effect draws (`src/mtblr.cpp:2999-3002,3027-3031`);
- full `B^{-1}` is used, so the effect prior couples latent traits even when
  some effective effects are masked;
- null/partially inactive latent coordinates are still drawn and enter the
  latent covariance update; and
- the alternative `sampleBetaCPG_Mt_latent_fullR()` at
  `src/mtblr.cpp:3178` attempts a full-precision construction but is not called
  by the default sampler.

### Marker-effect covariance

Three incompatible updates occur in one default iteration:

1. `sampleBset()` draws trait diagonals by scaled inverse chi-square and builds
   shrunken empirical correlations from active shared effects
   (`src/mtblr.cpp:475-594`).
2. Immediately afterward, `sampleB_latent()` overwrites `B` with an
   inverse-Wishart draw based on all latent marker vectors
   (`src/mtblr.cpp:148-187`).
3. After marker updates, `sampleB()` overwrites `B` again with the heuristic
   diagonal/shrunken-correlation rule (`src/mtblr.cpp:596-709,4260-4268`).

The first two calls occur once per user `set`, use all latent markers in the
second call, and are not guarded by `updateB`
(`src/mtblr.cpp:4185-4196`). Thus `updateB=FALSE` does not fix `B`; it only
skips the final heuristic update. The first draw is overwritten but still
consumes random numbers. Reported `covb` is accumulated from the final
heuristic `B`, while marker updates used the preceding inverse-Wishart `B`.
This is a statistical-contract conflict, not merely duplicated code.

### Residual covariance and likelihood precision

The active `sampleE()` at `src/mtblr.cpp:882-923` calculates a separate SSE
for each trait and updates only `E(t,t)` by scaled inverse chi-square. Existing
off-diagonal entries remain whatever was initialized, are accumulated into
`cove`, and are not learned. Regardless, the active marker kernel uses only
`diag(E^{-1})`. `sampleE_full()` (heuristic), the overload using a fixed
correlation matrix, and `sampleE_exact_sameX()` are present but inactive in the
default route.

`sampleE_exact_sameX()` at `src/mtblr.cpp:1207-1228` contains the correct
shared-design residual cross-product identity:

```text
S_e = Y'Y - (X'Y)'B - B'(X'Y) + B'(X'X)B,
E | B ~ IW(nue + N, S0 + S_e).
```

It is useful reference code but does not establish that the active default
sampler implements that model.

### CSR multi-trait alternative

`mtblr_cpg_omp_csr()` requires equal `n` and identical `ww` across traits and
one shared LD prefix (`src/mt_cpg_omp_csr.cpp:1096-1174`). Its residual
covariance update reconstructs a full cross-product from full `YY`, `WY`,
current effects, and shared LD, then uses an inverse-Wishart draw. This is the
closest current implementation to a full-overlap shared-design likelihood.

However, `sampleBetaCPG_Mt_arma_ld()`
(`src/mt_cpg_omp_csr.cpp:759-958`) restricts the likelihood RHS and matrices to
the active traits. With full residual precision, the linear term for an active
effect generally depends on residual scores from all traits; inactive-effect
coordinates do not imply that the corresponding observed residual score
vanishes. The implementation therefore does not consistently use full
off-diagonal residual precision for partial patterns.

It also selects active entries from the full prior precision `B^{-1}` rather
than inverting the active covariance submatrix. These are different priors:
one is a conditional precision and the other is a marginal active-effect
covariance. The correct choice depends on whether the intended model is a
latent full vector with a mask or a pattern-specific active subvector. The
repository currently contains both ideas without one authoritative contract.

The CSR `B` inverse-Wishart update accumulates zero/masked full vectors with a
degree-of-freedom rule based on active markers. Its conjugacy under arbitrary
partial patterns is not established in code or tests. This requires a model
derivation, not an implementation-only decision.

### Genetic and residual covariance summaries

`computeG()` in `src/mtblr.cpp:1147-1203` symmetrizes terms of the form
`b_t'(wy_s-r_s)/sqrt(n_t n_s)`. Under one shared design and compatible scaling,
this is related to `B'X'XB/n`; under trait-specific LD or samples it is not an
exact cross-trait genetic covariance without cross-design statistics.

The default sampler saves states when `it > nburn && it % nthin == 0`, whereas
current ST code generally uses a post-burn offset convention. Covariance means
and mean pattern probabilities are divided by `nit`, not by the actual saved
sample count (`src/mtblr.cpp:4371-4382`). Traces retain `nit+nburn` rows.
These semantics are insufficiently specified for direct migration.

### One-trait reduction

A one-trait fixed-variance marker update can algebraically resemble BayesC,
but an exact reduction to current ST-BLR is not presently plausible for the
full public default because:

- `updateB=FALSE` does not prevent the per-set `B` draws;
- latent effects for null markers enter `sampleB_latent()`;
- pi pseudo-count/update and saved-sample conventions differ; and
- output and chain semantics differ.

No one-trait reduction test exists. The legacy multivariate route should be
retained as reference until such a contract is selected and tested.

## 7. Data and likelihood operator map

| Representation | Current operations | Separation/coupling | Important assumptions | Candidate boundary |
|---|---|---|---|---|
| Dense individual level | Mostly simulation and BED sufficient-statistic preparation; no clean canonical dense joint sampler | legacy marker kernels directly manipulate matrices/lists | shared individuals/design required for full covariance | typed `score(j,t)`, `diag(j,t)`, `apply_delta(j, delta_vector)`, `residual_crossprod()` |
| Packed BED | `src/st_bed_decode.h`, `src/st_bed_bayesr_common.h`: decode selected columns/rows, standardize from AF, dot residual, update residual, reconstruct `Xb` | decoding and numeric operations are reusable, but the common header also converts Rcpp inputs and uses `NA_REAL` | exact concatenated marker order; PLINK two-bit coding; selected rows; AF/scaling contract; missing genotypes use zero on standardized scale | retain language-neutral byte reader/marker view; move Rcpp list/vector conversion to binding adapter |
| Dense sufficient statistics | legacy `sblr()` passes dense/list `XX`; marker kernels own traversal; `mt_sufficient_statistics_core.h` creates typed BED statistics | preparation core is clean; legacy sampler is tightly coupled | full `Y'Y` needed for full residual covariance; common design/scaling | first canonical multivariate reference should use exact shared-design `X'X`, `X'Y`, `Y'Y` |
| CSR LD | `CsrOperator::diag`, `apply_offdiag`, and `rebuild` in `src/st_ld_operator.h`; residual update subtracts diagonal and stored neighbours | best existing sampler/operator separation | disk values are correlations converted with `sqrt(xx_i xx_j)`; metadata does not itself identify marker IDs; R alignment is authoritative | retain as canonical ST operator; extend only with explicitly typed multi-trait shared/trait-specific views |
| Block eigen | `BlockEigenOperator` supplies the same interface using block factors; R transforms inputs according to hard/ridge filtering | already separated from BayesC/BayesR/SBayesRC policies | contiguous blocks, BED-derived selected order, approximation/filter contract | retain as optional approximate operator; not first multivariate correctness reference |
| Legacy MT CSR | private `LDCSR` reader/update logic in `src/mt_cpg_omp_csr.cpp` | duplicates ST CSR reader/residual propagation and is fused to pattern code | equal sample sizes, equal marker diagonals, one shared LD | retain as reference; do not make its private operator canonical unchanged |

For CSR, `r = X'y-X'Xb`; a marker effect change updates the marker diagonal and
stored LD neighbours. For BED, `e=y-Xb`; the marker score is `x'e+x'x b_j`,
and an effect change applies `e <- e-x delta`. Block eigen represents an
approximate symmetric operator `A` and reconstructs `r=wy-Ab`.

LD-swap requires neighbour/friend enumeration, residual SSE evaluation, and a
model-specific prior ratio. It should be an optional operator capability plus
a state/scale-policy acceptance calculation, not a mandatory method on every
operator.

**Recommendation.** The first canonical small-`T` multivariate BayesC
evaluation should use an exact dense shared-design sufficient-statistic
operator (`X'X`, `X'Y`, full `Y'Y`) derived from a common individual sample.
It makes the full residual cross-product testable without BED/CSR
approximations. `mt_sufficient_statistics_core.h` is the cleanest existing
typed preparation example, while the residual identity in
`sampleE_exact_sameX()` and the CSR MT implementation are reference formulas,
not ready-made canonical samplers. Selection of dense sufficient statistics
versus dense individual data remains a maintainer decision.

## 8. State, probability, and scale-policy map

### State and probability policies

| Policy | Probability/state location | Counts/update | Returned posterior | Shareable versus specific |
|---|---|---|---|---|
| Binary BayesC | duplicated `sampleBetaC*` functions across CSR, scheduled CSR, BED, prior/group/learned annotation files | active/null count; Beta/global pi where enabled | `dm=P(d=1)`, pi traces/final/mean | stable log-weight/categorical utilities shareable; scheduling and prior ratio model-specific |
| BayesR components | `src/st_cpg_omp_csr_bayesr.cpp` and BED BayesR file | component counts and Dirichlet global pi | marker-by-component `comp_prob`, `dm`, `dm_component_mean`, `ncomp` | component identity/count aggregation shareable; adaptive scheduling model-specific |
| BayesRC stick components | `src/st_bayesrc_annotation_prior.h` plus CSR/BED marker kernels | latent probit, alpha normal updates, sigma inverse-chi-square | component probabilities, alpha/sigma, marker priors | stick construction/update is already shared; likelihood weight remains operator-specific |
| Trait patterns | legacy MT files, especially `sampleBetaCPG_Mt_latent` and `sampleBetaCPG_Mt_arma_ld` | categorical pattern counts, Dirichlet pi | positional `dm`, current `d`, global `pi/pim` | pattern enumeration/categorical sampling shareable; active-subspace covariance definition unresolved |
| Group inclusion | `src/st_cpg_omp_csr_group.cpp` | group active/null counts; group Beta updates | group pi and included counts | categorical indexing shareable; group posterior is specific |
| Fixed marker pi | `src/st_cpg_omp_csr_prior.cpp` | no probability update when fixed | resolved marker prior namespace | clean fixed probability policy |
| Learned annotation pi | `src/st_cpg_omp_csr_annot.cpp` | centered-logit RW-MH | `eta_pi` mean/chain output | model-specific MH policy |

The proposed `state space` / `state-probability policy` / `marker-scale policy`
split is compatible conceptually. It is not yet a mechanical extraction:
current kernels calculate likelihood weights, state prior weights, scale,
draws, residual changes, counts, and diagnostics in one function. The first
shared utilities should be value-only operations (stable normalization,
categorical draw, state labels, count accumulation), not virtual calls inside
every marker update.

### Scale mechanisms

| Mechanism | Parameterization/update | Identifiability/normalization | Output/tests |
|---|---|---|---|
| Unit/global | active variance `vb` | scalar inverse-chi-square update | `vbs`, `vb`; extensive BayesC tests |
| Component | `vb * gamma_k` for BayesR/BayesRC | fixed ordered gamma, exact null at zero | `mixture_var`, component identities/reductions |
| Fixed marker multiplier | `vb*q_j` | caller supplies positive finite q; vb update uses `b_j^2/q_j` | prior namespace; fixed-prior tests |
| `selection_s` | `q_j=h_j^(S+1)`, `h_j=2p_j(1-p_j)` | fixed or bounded RW-MH `S in [-3,2]`; no mean normalization because functional MAF scale plus vb defines model | traces/acceptance; CSR BayesC/R/SBayesRC tests |
| Group multiplier | `vb*theta_g`; scaled inverse-chi-square from active `b^2/vb` | optional marker-count-weighted arithmetic mean normalization to one at `src/st_cpg_omp_csr_group.cpp:341-355`; vb is not counter-rescaled in that function | group means/chains; group tests |
| Learned annotation multiplier | `q_i=exp(center(A eta_vb)_i)`, clipped to bounds; RW-MH | centering gives geometric mean one before clipping; clipping may alter exact normalization | `eta_vb`; learned annotation/LD-swap tests |
| Trait-specific scale | all ST jobs own trait-specific `vb`, optional S, annotation/group state | independent traits, not a cross-trait hierarchy | trait columns and chains |

The current group arithmetic normalization is not the scale-preserving
geometric normalization proposed for a composable hierarchy. It should be a
reference behavior protected by tests, not silently reused as a general
hierarchy rule.

Scale also enters LD-swap acceptance. The ordinary CSR BayesC swap uses the
marker-specific normal prior density under `vb*q_j`
(`src/st_cpg_omp_csr.cpp:765-809`); fixed, learned-annotation, and group
backends implement corresponding scale- and probability-aware prior ratios in
their own files. BayesR/SBayesRC additionally condition the active prior on the
current non-null component multiplier. This duplication is evidence that the
scale policy must expose a log-prior contribution to optional moves, not only
the variance used for an effect draw.

**Smallest reusable scale-policy boundary.** A value policy needs to provide
`q(marker, trait, state/component)` and accumulate the sufficient statistic
required by its own update, plus an explicit normalization/identifiability
rule. Concrete initial policies are unit, fixed marker, and MAF/`selection_s`.
A categorical layer can provide one multiplier per level; multiple hierarchy
layers compose multiplicatively only after a common normalization convention
is approved. Component composition is `gamma_state * q_marker`. Arbitrary
overlapping annotation scales are future work and are not implied by the
existing categorical group backend.

## 9. Covariance-policy map

| Exact implementation | Prior/update | RNG/caller | Assessment |
|---|---|---|---|
| scalar ST `vb`/`ve` helpers across active ST kernels | scaled inverse-chi-square using active-effect SS or residual SSE and explicit degrees of freedom | chain-local `mt19937`; all production ST callers | canonical validated scalar reference |
| `src/distributions.cpp:rwishart/rinvertedwishart` and equivalents in legacy MT files | Bartlett Wishart / inverse-Wishart | `mt19937` in modern overloads | reusable reference after Armadillo-only/core separation and parameterization unification |
| `src/mtblr.cpp:sampleB_latent` | `B ~ IW(nub+m, S0 + sum beta beta')` | `mt19937`; active default but overwritten | exact conditional only for the explicitly defined full latent-beta prior; not the reported covariance update |
| `src/mtblr.cpp:sampleBset`, `sampleB` | inverse-chi-square diagonals; empirical shared-marker correlation shrunk by `n_shared/(n_shared+20)`, clamped and eigen-floored | `mt19937`; active default | heuristic, not a Gibbs draw from a stated joint covariance prior |
| `src/mtblr.cpp:sampleE` | independent scaled inverse-chi-square diagonals | `mt19937`; active default | exact diagonal variance updates only; off-diagonals are unchanged inputs |
| `src/mtblr.cpp:sampleE_full` | exact-ish diagonals plus residual-score empirical correlations, shrink/clamp/eigen floor | inactive | heuristic reference only |
| `src/mtblr.cpp:sampleE_exact_sameX` | full inverse-Wishart from exact shared-design residual cross-product | inactive | mathematically useful reference for full overlap |
| `src/mt_cpg_omp_csr.cpp` full `E` update | inverse-Wishart from full `YY` and shared LD residual cross-product, followed by SPD correction | `mt19937`; direct native caller only | candidate reference; exact only under its shared-design input contract and before target-altering corrections |
| `src/mt_cpg_omp_csr.cpp` full `B` update | inverse-Wishart from masked effects with active-count degrees of freedom | `mt19937`; direct native caller | conjugacy unresolved for partial patterns |
| positive-definite corrections | ridge, correlation clamp, eigenvalue floor scattered through `mtblr.cpp` and CSR MT | deterministic | numerical safeguards may alter a posterior transition; must be explicit diagnostics/policy, not hidden cleanup |

There is no factor or low-rank covariance sampler in the active repository.
Correlation summaries are derived in R with `cov2cor` when covariance
diagonals are positive.

**Recommendation.** First evaluate a small-`T`, full-overlap covariance model
with a full inverse-Wishart residual covariance and either (a) a formally
defined latent full-effect covariance with masking, or (b) a pattern-specific
active-subvector covariance. Do not mix these. Include fixed covariance and
diagonal covariance modes as reductions. The scale-matrix convention
(`S0` versus `nu*S0`) and degrees of freedom must be named in the resolved
specification.

**Maintainer decision required.** The repository does not determine which of
the two active-effect prior interpretations is scientifically intended, nor
whether full inverse-Wishart residual covariance should be the canonical
default. Those choices cannot be made by code extraction alone.

## 10. Chain, RNG, and parallelism audit

### Current canonical ST infrastructure

`src/st_chain_utils.h` defines deterministic task mapping and seeds:

```text
trait-only: seed + 1000003*(trait+1) + 9176
trait-chain: seed + 1000003*(trait+1) + 9176*(chain+1)
explicit chain base: chain_seed[chain] + 1000003*(trait+1)
```

Current unscheduled CSR backends run trait-chain tasks with OpenMP
`schedule(static)`, construct one `std::mt19937` per task, aggregate plain C++
chain results after the parallel region, and optionally materialize compact
chains. Burn-in/thinning and `nsamples` are recorded in the raw metadata. This
is the best canonical infrastructure.

Marker order is part of each RNG contract. Unscheduled ST kernels use a fixed
marker traversal; scheduled BayesC/BayesR may skip or wake null markers and
perform periodic full sweeps; BED BayesRC deliberately performs full sweeps.
The default legacy MT sampler sorts markers by the maximum trait-specific
`x2` score (`src/mtblr.cpp:4153-4167`) and then filters that order through the
user sets. Any adapter must preserve these existing orders until a separately
tested model defines a different schedule.

BED BayesR/BayesRC also own uniform and normal distribution state inside each
chain runner. Their fixed-seed and one-core/two-core tests protect the recent
chain-local distribution fix. BED backends currently reject explicit
`chain_seeds` at the public wrapper.

### Migration-required occurrences

The requested kernel search covered `R::rnorm`, `R::rchisq`, `R::runif`,
`GetRNGstate`, `PutRNGstate`, `arma::randn`, `arma::randu`, `Rcpp::Rcout`,
`Rprintf`, static distributions, `thread_local`, and `std::mt19937`.

- No `GetRNGstate`, `PutRNGstate`, `R::runif`, `Rprintf`, or `arma::randu`
  occurrence was found in `src`.
- `src/mtblr.cpp:14-41` contains old overloads using `arma::randn`,
  `R::rchisq`, and `R::rnorm`; some old covariance functions call them.
  `mtblr_eigen` reaches `mvrnormARMA()` through
  `sampleBetaCMt_eigen()` (`src/mtblr.cpp:3800,3897,4975`), so its native seed
  does not control all marker draws.
- `src/mt_cpg_omp.cpp:548-590` seeds marker engines with iteration and OpenMP
  thread number. Changing core count or static assignment changes which
  engine a marker receives.
- Active persistent distribution state occurs at:
  - `src/st_cpg_omp_csr_scheduled.cpp:85-86`;
  - `src/st_cpg_omp_individual_scheduled.cpp:311-312`; and
  - `src/st_cpg_omp_individual_scheduled_chains.cpp:537-538`.
  These are statistical-core occurrences. A cached normal variate can survive
  a chain/call and migrate with a worker thread, so they are potentially
  non-deterministic and require chain-local ownership.
- `src/st_bayesrc_annotation_prior.h:13,20` uses `R::pnorm`/`R::qnorm` from
  trait-chain workers. These are R-specific math calls inside the statistical
  core and violate the intended language-neutral/thread boundary even though
  they are not RNG calls.
- `Rcpp::Rcout` occurs in many exported kernels. The inspected production ST
  paths generally print setup/diagnostics outside worker loops, but output is
  still mixed into the core translation units and should move to bindings or
  returned diagnostics.

### Coverage gaps

BED BayesRC has strong repeated-call, intervening-call, reversed core-order,
different-seed, and one/two-core tests. CSR SBayesRC and annotation chains have
explicit seed tests. There is no corresponding automated coverage for any
joint multi-trait backend, for `mtblr_eigen`, or for all scheduled BayesC
families. Therefore current tests do not establish family-wide core-count or
chain-order invariance.

## 11. Language-neutral C++ boundary audit

### Classification

| Classification | Files/functions | Reason |
|---|---|---|
| Suitable core now | `src/mt_sufficient_statistics_core.h`; most value algorithms in `src/st_chain_utils.h`; typed parts of `src/st_block_eigen.*`; standard-RNG Wishart logic in `src/distributions.*` | standard C++/Armadillo values, explicit ownership, exceptions instead of R object construction |
| Thin binding already | generated `src/RcppExports.cpp`; exported argument validation portions of `src/mt_sufficient_statistics.cpp` and `src/sparse_ld_bed_core.cpp` | R object conversion/registration is appropriate here |
| Mixed core/binding requiring separation | all active ST backend `.cpp` files; `src/st_bed_bayesr_common.h`; `src/st_csr_common.h`; `src/st_bayesrc_annotation_prior.h`; legacy MT files | statistical loops coexist with `Rcpp::List`, `NumericMatrix`, `Nullable`, `Rcout`, `NA_REAL`, or R math/RNG |
| Appropriately R-specific | final named raw-list assembly, R attributes/names, nullable argument decoding, error translation | should remain in a thin Rcpp adapter |

Across `src`, `Rcpp::List`, `NumericVector`, `NumericMatrix`, and `Nullable`
appear in 15-17 translation units, including most statistical backends.
`NA_REAL` appears in 11 files, including selection/optional output assembly.
These are not all unsafe: most R objects are read before OpenMP and constructed
after it. They nevertheless prevent a direct Python binding to the core.

### Ownership, orientation, and storage findings

- Canonical ST marker matrices are `m x nt`; traces are iteration x trait;
  component probabilities are a trait list of `m x K` matrices.
- Legacy MT native results are trait-major nested vectors that R transposes or
  coerces positionally. This orientation must not leak into a new core result.
- Rcpp/Armadillo inputs are commonly copied into Armadillo or standard vectors
  before workers. BED/CSR disk readers own standard buffers; these formats are
  language-neutral.
- Packed BED is disk-backed by PLINK files. CSR uses binary arrays plus text
  metadata. Block eigen is reconstructed from BED and R-provided block/filter
  metadata. Core statistical workflows do not require R serialization.
- R-backed `NumericMatrix y` and R list arguments exist at binding boundaries;
  worker code must continue to use preconverted/owned buffers and must not
  allocate or mutate SEXP-backed objects.

### Proposed typed adapter boundary

The following are boundary concepts, not Phase 0 code:

```text
BlrDataSpec
  representation, marker_count, trait_count, marker_ids/order token,
  trait_ids, sample/overlap contract, scaling contract, operator handle/spec

BlrModelSpec
  state-space tag and states, probability-policy tagged payload,
  scale-policy tagged payload, covariance-policy tagged payload,
  fixed/update flags and validated dimensions

McmcControl
  iterations, burn-in, thinning, chain count, task seeds,
  core count, scheduling/rebuild controls

ChainResult
  marker means/probabilities/current state, traces, covariance summaries,
  diagnostics, and tagged optional model results

BlrResult
  aggregate marker/trace/covariance results, chain summaries,
  metadata, and tagged optional model results

ModelSpecificResult
  component, pattern, annotation, group, selection, LD-swap payloads
```

Use enums/tagged variants and owned standard/Armadillo values. Rcpp should only
convert R objects into these values, invoke a typed function, and convert the
typed result into the existing `stblr_raw` representation.

### Work now versus Python-deferred work

Do now: define orientations/ownership, remove R objects from new typed types,
standardize exception/error translation, make seeds/distribution ownership
explicit, and isolate normal CDF/quantile behind a language-neutral math
utility. Protect exact existing output through adapters.

Defer until an actual Python binding: Python package layout, NumPy ownership
and zero-copy policy, Python exceptions, wheels/toolchains, Python-facing disk
object classes, and a Python-specific serialization format. Do not distort the
R API or introduce a second result schema solely for hypothetical bindings.

## 12. Raw schema and formatter map

### Current schemas

| Family | Native return | Dimensions/optional namespaces | Formatter |
|---|---|---|---|
| Legacy joint MT | positional 20-slot nested vector list | marker slots trait-major; traces include burn-in; covariance `T x T`; pattern pi; no schema metadata/chains/diagnostics contract | manual code in `R/interface_mtblr.R:252-294` |
| ST BayesC | named `stblr_raw` v1 | marker `m x nt`, trace `n_trace x nt`, variance `nt x nt`, two-column pi; optional chains/selection/LD-swap | `.validate_stblr_raw()` then `.as_stblr_fit()` |
| ST BayesR | named v1 + `component` | trait list of `m x K` probabilities; `component_0`; optional component counts/chains | same formatter, which derives `dm` and component mean |
| SBayesRC/BED BayesRC | named v1 + `component` and `annotation` | `gamma_0.00`; alpha trait list `P x (K-1)`; sigma `(K-1) x nt`; marker priors | same formatter plus annotation aliases/summaries |
| Fixed/group/learned annotation BayesC | named v1 + `prior`, `group`, or `annotation` | backend-specific matrices and compact-chain payloads | same formatter plus standardized public aliases |

`.is_stblr_raw()` accepts only class `stblr_raw`, version 1
(`R/sparse_ld_bed_helper.R:541-546`). `.stblr_raw_structure_problems()` checks
meta-implied marker, trace, covariance, pi, component, annotation, group, and
chain dimensions. Production wrappers reject positional fallback. Older
positional formatter functions remain only for direct compatibility tests.

The central raw-to-fit mapping is explicit at
`R/sparse_ld_bed_helper.R:802-967`:

- `marker.{bm,dm,wy,r,b,state}` -> `bm,dm,wy,r,b,d`;
- `trace.{vbs,vgs,ves,vle,vld,pis}` -> like-named public traces;
- `variance.{covb,covg,cove,vb,vg,ve}` -> like-named matrices;
- `pi.final/pi.mean` -> `pi/pim` (and BayesR aliases);
- `component.prob` -> `comp_prob`, with `dm` recomputed from the named null;
- annotation alpha/sigma/final marker priors -> `alpha`, `sigmaSqAlpha`, and
  annotation fields; and
- diagnostics and compact chains -> public diagnostic/chain structures.

### Historical/optional differences

- `Rcpp::List chains = R_NilValue` remains in
  `src/st_cpg_omp_csr.cpp:2095`,
  `src/st_cpg_omp_csr_bayesr.cpp:2176`, and
  `src/st_sbayesrc_omp_csr.cpp:2432`. Newer annotation/group/prior and BED
  BayesRC paths use `Rcpp::RObject` to preserve actual R `NULL`. Raw `NULL`
  representation is therefore not uniform across all backends.
- `.as_stblr_fit()` assigns `out$chains <- NULL` in the no-chain branch
  (`R/sparse_ld_bed_helper.R:1219-1222`), which removes the list element in R.
  This differs from the repository's present-but-`NULL` formatted-field rule,
  although `$chains` still evaluates to `NULL`. `ld_swap` fields are repaired
  through a dedicated ensure helper; chains are not.
- Optional model fields are naturally model-specific; differences in null
  component names, annotation coefficients, group outputs, and pattern
  probabilities are substantive. Positional slot reuse and trait-major legacy
  orientation are historical accidents.

### Version recommendation

`stblr_raw_v1` can support another independent-trait policy by adding an
optional namespace without changing canonical dimensions. A canonical joint
multivariate result needs explicit pattern posterior dimensions, full
covariance traces/means/final states, overlap metadata, and possibly
trait-pair diagnostics. Because legacy MT has no stable named contract, a
versioned extension is feasible. Whether it is called `stblr_raw_v2` or
`blr_raw_v1` should wait until those dimensions are approved. Typed C++
results can be converted into v1 for unchanged ST routes and a future v2 for
joint MT; no schema change is required in Phase 1.

## 13. Test-coverage map

| Reduction/capability from `blr_reduction_test_matrix.md` | Current coverage | Gap before migration |
|---|---|---|
| Deterministic seed/core/chain order | strong for BED BayesRC; seed tests for CSR SBayesRC, annotation chains, `selection_s` | no joint MT coverage; scheduled BayesC persistent distribution state is not covered across preceding calls/core order |
| Chain aggregation/`keep_chains` | CSR BayesC/R, annotation families, BED BayesRC; chain mean identities and single-chain degenerate summaries | actual raw `NULL` not uniform; not all BED backends retain compact chains |
| BayesC reductions (RED-01..04) | scalar/ST behavior is extensively tested within backend families | no canonical MT `T=1`, diagonal covariance, fixed covariance, or all-active reduction |
| BayesR reductions (RED-05) | exact CSR/BED BayesR backend tests and component identities | add explicit one-non-null-component-to-BayesC canonical reduction if not already isolated |
| BayesRC fixed-prior reduction (RED-06) | exact BED BayesRC -> fixed-pi BED BayesR at strict tolerance; CSR helper/reduction coverage | preserve as canonical policy reduction |
| Annotation alignment/enrichment (ANN) | public BED BayesRC shuffled IDs, `chr`/`cls`, missing/duplicate IDs, preprocessing, enrichment; CSR annotation interface tests | reusable fixtures are strong; add hierarchy reductions only when hierarchy exists |
| Fixed/group/learned scale | prior/group/learned annotation tests; group enrichment and LD-swap; `selection_s` tests | no multiple-layer hierarchy or component x hierarchy test |
| Block eigen (OP) | construction/filter/smoke and wrapper tests | no complete exact-operator equality suite for every policy; approximation error contract needs explicit tests |
| LD-swap | CSR BayesC/R/annotation/SBayesRC tests and diagnostics | not available for BED or canonical MT; policy/operator split needs acceptance-ratio reductions |
| BED/CSR consistency | backend field/identity and low-level BED score/pairwise checks | no exact end-to-end likelihood-operator equivalence across dense/BED/CSR under one fixture |
| Multivariate/covariance (COV) | no automated `sblr()`/`mtblr*` tests; exploratory scripts only | all MT reductions, full precision behavior, covariance recovery, SPD, seed/core invariance, and output schema are missing |
| Raw schema/formatting/public fields | `test-stblr-raw-schema.R`, backend inventory/consistency and public-interface tests | joint MT is positional and outside schema; present-but-NULL `chains` deserves a contract test |
| Future hierarchy/factor/Python | none, appropriately | tests should be added with implementation, not pre-assert invented behavior |

There are 238 `test_that` blocks in the current suite and zero test files that
call `sblr()` or a joint `mtblr*` sampler. Exploratory fixtures in
`examples/Evaluation_mtblr_revised.R`, `examples/mtsim_sparse_ld.R`, and
`examples/mtsim_prior_sparse_ld.R` are useful data sources but are not
regression protection.

Canonical reusable fixtures include the tiny BED/CSR marker-order fixtures,
the BED BayesRC fixed-prior reduction, annotation-ID shuffle/subset fixtures,
component identity checks, `selection_s` fixed/sampled fixtures, and raw schema
malformation tests. New framework tests should compare stable formatted fields
and named raw namespaces, not legacy positional slot numbers.

### Documentation and example inventory

The principal user-facing model documents are `README.md` and the QMD/rendered
HTML pairs under `docs/notes/`: `model_overview`, `single_trait_models`,
`multi_trait_models`, `data_representations`,
`technical_multitrait_overlap`, `technical_summary_statistics`,
`technical_sparse_ld_csr`, `annotation_informed_models`,
`technical_annotation_priors`, `technical_sbayesrc_annotations`,
`sparse_ld_and_bed_workflows`, `stblr_interface_and_output_map`, and
`workflows_and_outputs`. Developer architecture is primarily documented by
the three unified-framework documents and the backend naming, computation
inventory, raw schema, block-eigen, annotation, individual BayesRC, and chain
design notes cited in section 2.

Executable/user examples are concentrated in `examples/workflows/01` through
`06`, the compatibility workflow scripts in the same directory,
`examples/Evaluation_mtblr_revised.R`, `examples/mtsim_sparse_ld.R`, and
`examples/mtsim_prior_sparse_ld.R`. The numbered workflows chiefly exercise
the validated ST interfaces. The multi-trait scripts are exploratory and do
not establish a public joint-sampler contract.

Observed documentation inconsistencies:

- `man/sblr.Rd` calls `yy` a trait cross-product matrix, while the public
  default/CPG native signatures accept a vector and the wrapper's unnamed
  matrix name handling is inconsistent, as detailed in section 6.
- `docs/notes/model_overview.qmd` lists only BayesC/BayesR for `stblr_bed()`;
  the local public interface also supports BayesRC.
- `docs/dev/stblr_backend_naming.md` and the leading backend map in
  `docs/dev/stblr_backend_computation_inventory.md` do not fully inventory the
  now-public BED BayesRC route, although later computation-inventory text and
  `docs/dev/stblr_raw_schema.md` discuss it.
- `docs/notes/multi_trait_models.qmd` and
  `docs/notes/technical_multitrait_overlap.qmd` appropriately call current MT
  algorithms experimental and distinguish the CSR full-`YY` export from other
  contracts. Those caveats should be retained and strengthened with the exact
  covariance findings from this audit.
- Source QMD and rendered HTML coexist. Later documentation work must render
  the selected source pages so their checked-in HTML does not silently remain
  stale; no rendering was appropriate in this no-modification audit.
- Draft and `old_workflows_and_outputs` documents should remain explicitly
  non-authoritative. They contain historical ideas and should not be used to
  infer current supported methods or output fields.

## 14. Duplication and coupling analysis

### Critical

1. **Conflicting active covariance transitions in public default MT.** Files:
   `R/interface_mtblr.R:sblr`, `src/mtblr.cpp:mtblr`, `sampleBset`,
   `sampleB_latent`, `sampleB`, `sampleE`, and
   `sampleBetaCPG_Mt_latent`. Consequence: marker transitions, reported `B`,
   and `updateB` do not describe one posterior; residual off-diagonals are not
   learned or used. Migration action: freeze as legacy reference, select and
   derive one covariance/pattern model before extraction. Protection: full
   MT reductions, conditional-moment checks, covariance known-truth tests, and
   fixed-update-flag tests.

2. **No automated joint MT contract.** Files: all `mt*.cpp`,
   `R/interface_mtblr.R`; absence under `tests/testthat`. Consequence: any
   refactor can silently change mathematics, RNG, orientation, or outputs.
   Migration action: add characterization tests before routing production.
   Protection: RED-01..04, COV, RNG/core, and positional compatibility tests.

### High

1. **Persistent random-distribution state in active scheduled BayesC.** Exact
   locations are listed in section 10. Consequence: repeated-seed and
   core-count results can depend on prior calls/thread assignment. Migration
   action: make distributions chain-local without changing draw order, in a
   separately scoped correction. Protection: both core orders and intervening
   calls for scheduled CSR and BED BayesC.

2. **Duplicated/fused marker kernels.** Files: CSR BayesC/R/annotation/group,
   BED BayesC/R/RC, legacy MT. Consequence: fixes and identities drift between
   representations/policies. Migration action: extract value utilities and
   operator calls one validated family at a time. Protection: exact
   backend/policy reductions.

3. **Legacy MT RNG is not one engine contract.** Files:
   `R/interface_mtblr.R:113`, `src/mt_cpg_omp.cpp:548-590`,
   `src/mtblr.cpp:mvrnormARMA` and `mtblr_eigen`. Consequence: users cannot
   directly set the native seed; core count and Armadillo global RNG can change
   results. Migration action: explicit MCMC seed spec and chain-local engine,
   after characterization. Protection: repeated/core/algorithm-order tests.

4. **MT CSR duplicates the operator and has unresolved full-precision partial
   patterns.** File: `src/mt_cpg_omp_csr.cpp`. Consequence: apparent full
   covariance support does not imply a coherent partial-pattern likelihood.
   Migration action: derive active-subspace contract, then use a shared typed
   operator. Protection: brute-force tiny Gaussian conditional comparisons.

### Medium

1. **Raw and formatted optional-null inconsistency.** Files cited in section
   12. Consequence: schema consumers can observe empty list versus actual
   `NULL`, and `chains` can be absent rather than present-but-NULL. Migration
   action: normalize only under an explicit schema/formatter correction.
   Protection: names plus `is.null` tests for every family.

2. **R-specific math and object assembly mixed with cores.** Files:
   `src/st_bayesrc_annotation_prior.h` and active backend translation units.
   Consequence: unsafe/awkward Python reuse and harder OpenMP review.
   Migration action: typed results and language-neutral math adapter.
   Protection: exact BayesRC probability/reduction tests.

3. **Repeated R setup/metadata code.** Files: `R/sparse_ld_bed_helper.R` and
   all annotation wrappers. Consequence: defaults, seed rules, and scale names
   can diverge. Migration action: resolved immutable spec plus legacy argument
   adapters. Protection: snapshot/structural resolved-spec tests and unchanged
   public results.

4. **Group versus learned-scale normalization differs.** Files:
   `src/st_cpg_omp_csr_group.cpp` and
   `src/st_cpg_omp_csr_annot.cpp`. Consequence: a naive hierarchy composition
   would change effective global variance. Migration action: require an
   approved scale normalization contract. Protection: one-layer exact
   reductions and scale-product identities.

### Low

1. **Stale inventories/user notes.** `docs/dev/stblr_backend_naming.md` and the
   leading table in `docs/dev/stblr_backend_computation_inventory.md` omit or
   incompletely list public BED BayesRC; `docs/notes/model_overview.qmd` says
   `stblr_bed()` supports BayesC/BayesR but not BayesRC. Consequence: readers
   can select the wrong capability assumptions. Migration action: update in a
   documentation phase, not this audit. Protection: documentation link/check
   review.

2. **Large commented/alternative code in production files.** Most pronounced
   in `src/mtblr.cpp`, `src/mt_cpg_omp_csr.cpp`,
   `src/st_sbayesrc_omp_csr.cpp`, and scheduled individual files. Consequence:
   active call paths are hard to audit and old R RNG code appears beside modern
   code. Migration action: retain temporarily, retire only after reference
   tests and history are available. Protection: call graph and reduction suite.

## 15. Proposed migration map

| Classification | Current code | Migration treatment |
|---|---|---|
| Retain unchanged | public function signatures/defaults; marker/trait/annotation alignment; `stblr_raw_v1`; `.as_stblr_fit()`; validated ST samplers; BED decoding and CSR disk formats | build adapters around them first |
| Extract into shared core | `st_chain_utils.h`; stable categorical/log-weight helpers; component identities/counts; typed operator value methods; standard-RNG Wishart utility after parameterization decision | extract incrementally with exact tests |
| Wrap through adapter | `stblr_csr()` unscheduled BayesC first, then BayesR/SBayesRC, BED, annotation models | resolved spec expands to identical legacy native arguments and existing raw formatter |
| Retain temporarily as reference | all `mtblr*` implementations; `sampleE_exact_sameX`; MT CSR full-`YY` formulas; old positional formatters | characterize before any retirement or canonical reuse |
| Deprecate later | unreachable BayesR mapping in `sblr()`, unused alternative covariance/marker functions, duplicated private MT CSR reader, positional MT output | only after a supported replacement and compatibility policy |
| Unresolved | active-effect covariance under partial patterns; canonical residual covariance prior; partial overlap; maximum pattern space; raw schema version; first MT operator; legacy algorithm migration | maintainer decisions listed below |

Migration order should be ST CSR BayesC specification boundary, other validated
ST probability/scale policies, BED adapters, then a new separately validated
canonical MT BayesC. It should not be a package-wide rewrite.

## 16. Recommended Phase 1 scope

One bounded Phase 1 task:

> Introduce an internal resolved BLR specification and typed C++
> specification/result boundary for unscheduled CSR BayesC with unit marker
> scale and independent trait jobs, while routing the existing public call
> through an argument-preserving legacy adapter and leaving sampler mathematics
> and `stblr_raw_v1` unchanged.

Concrete scope:

1. Add an internal R resolver that converts the already aligned
   `stblr_csr(method="bayesc", scheduled=FALSE)` inputs into explicit data,
   state, fixed/global probability, unit-scale, scalar-variance, and MCMC
   sub-specifications.
2. Validate the resolved object and expand it back into exactly the current
   `stblr_cpg_omp_csr()` arguments. Do not expose a new public class yet.
3. Add language-neutral C++ POD/tagged types for `BlrDataSpec`,
   `BlrModelSpec`, `McmcControl`, `ChainResult`, and `BlrResult`, plus a thin
   Rcpp conversion/round-trip adapter. Do not put `Rcpp` or `SEXP` in those
   types.
4. Keep the existing native sampler as the execution target. Typed result
   conversion should demonstrate the future boundary but must convert to the
   existing named raw schema.
5. Add exact fixed-seed before/after comparisons for every formatted field,
   raw namespaces/dimensions, chain aggregation, explicit chain seeds,
   one/two-core equality, `selection_s=0` fixed reduction, LD-swap-off behavior,
   and malformed resolved-spec errors.

Candidate files for Phase 1 (names may follow repository conventions): a new
`R/blr-model-spec.R`, new `src/blr_spec.h`, `src/blr_result.h`, and a thin
`src/blr_rcpp_adapter.cpp`; narrowly scoped changes to
`R/sparse_ld_bed_helper.R`; and a new focused test file. Generated Rcpp files
would be regenerated only if a new exported adapter is actually required.

Files that must remain untouched in Phase 1: `src/mtblr.cpp`, all other
`src/mt_*.cpp`, all BED samplers, BayesR/SBayesRC/annotation/group samplers,
`src/st_bayesrc_annotation_prior.h`, block-eigen mathematics, public roxygen
signatures/defaults, and `stblr_raw_v1` field semantics. The scheduled BayesC
RNG issue should be a separate correction with its own trajectory tests, not
folded into the specification boundary.

## 17. Decision points requiring maintainer approval

1. **Active-effect covariance model.** Is a partially active marker a masked
   latent full Gaussian vector or a Gaussian on the active subvector? This
   determines both marker probabilities and covariance conjugacy.
2. **Canonical covariance prior.** Select diagonal scaled inverse-chi-square,
   full inverse-Wishart, fixed covariance, or another explicitly derived prior;
   define degrees of freedom and scale parameterization.
3. **Residual covariance contract.** Decide whether the first joint model
   requires full shared-sample `Y'Y`, permits diagonal `E`, and how partial/no
   overlap is represented. `Nmat` alone is not sufficient.
4. **Default multi-trait sharing modes.** Decide whether all `2^T` patterns are
   allowed, what maximum is enforced, and which restricted defaults are
   supported.
5. **Raw schema version/name.** Decide between extending `stblr_raw_v1`,
   `stblr_raw_v2`, or `blr_raw_v1` only after pattern/covariance orientations
   are fixed.
6. **Initial canonical multivariate operator.** Choose exact dense
   shared-design sufficient statistics versus dense individual data. This
   audit recommends the former for transparent reductions.
7. **Legacy `sblr()` migration policy.** Decide whether each `default`, `cpg`,
   `cpg_arma`, `cpg_omp`, and `eigen` algorithm is characterized and preserved,
   explicitly marked experimental, or deprecated after a new route exists.
8. **Hierarchy normalization.** Select the identifiability rule for one and
   multiple categorical variance layers before treating current group
   multipliers as a hierarchy reduction.

## 18. Risk register

| Risk type | Current risk | Likelihood/impact | Mitigation gate |
|---|---|---|---|
| Statistical | conflicting MT `B` updates, diagonal-only active `E`, unresolved partial-pattern covariance | high / critical | formal model contract plus RED/COV tests before canonical MT |
| Reproducibility | persistent scheduled-BayesC distribution caches; MT thread/Armadillo RNG; no public MT seed | high / high | chain-local RNG ownership and repeated/core/order tests |
| Compatibility | public positional `sblr()` output versus named ST raw; defaults and optional fields | medium / high | adapter-first migration and exact formatted-field comparisons |
| Performance | virtual/policy abstraction in marker loop, pattern explosion, copies at R boundary | medium / medium-high | value/tag dispatch outside hot loop, benchmark after correctness, explicit pattern cap |
| Build system | Rcpp/OpenMP/Rtools dependence; current audit environment cannot compile | medium / medium | CI/compiler matrix and generated-file discipline; do not infer code failure from missing toolchain |
| Language binding | Rcpp/R math/`NA_REAL` mixed into core; orientation and lifetime assumptions | high / medium-high | owned typed core values, thin adapters, Python work deferred until real binding phase |
| Schema | joint MT lacks named raw contract; actual-null/optional fields vary | medium / high | approve dimensions/version, typed result conversion, schema validation tests |
| Documentation | capability tables lag BED BayesRC and legacy MT caveats | medium / medium | dedicated documentation synchronization after architecture decisions |

## 19. Readiness assessment

Phase 0 is complete as a repository-grounded audit. The current repository has
enough validated ST infrastructure to begin a bounded specification/adapter
phase without changing sampler mathematics. It is not ready to declare a
canonical multivariate sampler: covariance and active-pattern semantics require
maintainer decisions and new tests first. Those are Phase 1/ later design
inputs, not blockers to completing this audit.

PHASE 0 AUDIT COMPLETE — READY TO DESIGN PHASE 1
