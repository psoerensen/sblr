# selection_s implementation status

The `s` prefix identifies summary-statistics models; it does not request MAF
scaling. `selection_s` and `estimate_selection_s` independently control the
MAF-dependent variance policy for any route that explicitly supports it.
Without that option BayesC uses `unit` and BayesR/BayesRC use `component`;
with it they use `maf_s` and `component_maf_s`, respectively.

The package argument is `selection_s`. It corresponds to the BayesS-style
MAF architecture parameter `S`, applied as marker-specific prior variance
scaling in the CSR samplers that currently support it.

## Scale convention

CSR effects are on the standardized-genotype scale. The MAF-dependent prior
scale is

```text
q_j(S) = h_j^(S + 1)
h_j = 2 p_j (1 - p_j)
```

where `p_j` is the minor allele frequency. The `S + 1` exponent is used
because CSR marker effects are standardized-genotype-scale effects rather than
allele-scale effects.

## Backend support

| Public route and model | Fixed `selection_s` | Sampled `selection_s` |
|---|---:|---:|
| `stblr_csr(sbayesc)` | yes | yes |
| `stblr_csr(sbayesr)` | yes | yes |
| `stblr_block_eigen(sbayesc)` | yes | yes |
| `stblr_block_eigen(sbayesr)` | yes | yes |
| `stblr_block_eigen(sbayesrc)` | yes | yes |
| `stblr_csr_annot(sbayesrc, annotation_probit_stick)` | yes | yes |
| `stblr_csr_annot(sbayesc, fixed_marker/group/learned_logistic)` | no | no |
| `stblr_bed(bayesc/bayesr/bayesrc)` | no | no |
| `mtblr_bed(bayesr)` | yes | no |
| `mtblr_csr(sbayesr)` | yes | no |
| `mtblr_block_eigen(sbayesr)` | yes | no |
| MT BayesC routes | no | no |

## User-facing modes

```r
selection_s = NULL
estimate_selection_s = FALSE
```

fits the ordinary model.

```r
selection_s = -0.5
estimate_selection_s = FALSE
```

fits a fixed global `selection_s` model. The same supplied value is used for
all traits and chains.

```r
selection_s = NULL
estimate_selection_s = TRUE
```

samples `selection_s` by Metropolis-Hastings. Internally, sampled `S` is
trait-specific and chain-specific; the formatted fit reports posterior
summaries across retained iterations and chains.

Supplying both a numeric `selection_s` and `estimate_selection_s = TRUE` is
invalid.

## Fixed-S prior forms

For CSR BayesC:

```text
b_j | d_j = 1, vb, S ~ N(0, vb * q_j(S))
```

For CSR BayesR:

```text
b_j | component_j = m, vb, S ~ N(0, vb * gamma_m * q_j(S))
```

For CSR SBayesRC:

```text
b_j | component_j = m, vb, S ~ N(0, vb * gamma_m * q_j(S))
```

For the SBayesRC backend, annotations affect component probabilities. The
BayesS-style MAF-scaling argument affects marker-specific effect-size prior
variance.

## Sampled-S defaults

```r
selection_s_init = 0
selection_s_prior = c(-3, 2)
selection_s_proposal_sd = 0.35
```

## Sampled-S MH kernels

For BayesC:

```text
log p(S_t | b_t, d_t, vb_t)
  = log p(S_t)
    - 0.5 sum_{j: d_jt = 1} [
        log q_j(S_t) + b_jt^2 / (vb_t q_j(S_t))
      ]
```

For BayesR and SBayesRC:

```text
log p(S_t | b_t, gamma_t, vb_t)
  = log p(S_t)
    - 0.5 sum_{j: gamma_jt > 0} [
        log q_j(S_t) + b_jt^2 / (vb_t gamma_jt q_j(S_t))
      ]
```

Posterior summaries are averaged or summarized across chains in the returned
fit object.

## Component naming

BayesR component probabilities use `component_0` for the null component, and
`dm = 1 - P(component_0)`.

SBayesRC component probabilities are named by gamma values. The null component
column is `gamma_0.00`, and `dm = 1 - P(gamma_0.00)`.

## Interpretation

`selection_s = -1` gives `q_j(S) = 1`, so it reproduces the ordinary unscaled
model under the same seed for fixed-S fits.

Valid compact examples:

```r
fitC <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesc",
  estimate_selection_s = TRUE
)

fitR <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "sbayesr",
  estimate_selection_s = TRUE
)

fitS <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = annotations,
  annotation_model = "sbayesrc",
  estimate_selection_s = TRUE
)
```

Unsupported examples include sampled `selection_s` with
`annotation_model = "prior"`, `"learned"`, or `"group"`, and all BED backends.

## MTBLR Phase 19 status

MT SBayesR supports one fixed scalar `selection_s` shared by the joint model.
The default is zero and `selection_s = -1` gives exactly unit marker scale.
`estimate_selection_s = TRUE` fails explicitly because a validated joint-MT
Metropolis-Hastings update is not yet part of the retained prior contract.
## Phase 20 annotation models

Fixed `selection_s` is supported independently for MT BayesRC/SBayesRC on all
three operators. `NULL` means no MAF scale; `-1` explicitly requests the
unit-scale reduction. Sampled MT S remains unsupported. MAF-named annotations
plus selection-S are allowed with explicit overlap metadata and one advisory.

Phase 21 diagnoses S only when the sampler genuinely updates it. A fixed finite
S is `not_updated`; `NULL` is `not_applicable`; acceptance rate remains a
sampler diagnostic rather than an R-hat/ESS quantity. Sampled MT S is still
rejected.

The sampled ST trace is task-private and recorded at every post-burn completed
iteration. It coexists with component-probability, annotation-parameter, and
explicit selected-marker capture without diagnostic thinning or RNG use.
