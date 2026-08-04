# SBayesRC reference crosswalk

## Sources and scope

This crosswalk records the scalar contracts inspected during the Study 06
stabilization. The scientific reference is Zheng et al., *Leveraging functional
genomic annotations and genome coverage to improve polygenic prediction of
complex traits within and between ancestries*, **Nature Genetics** 56,
767--777 (2024), DOI: 10.1038/s41588-024-01704-y. The manuscript repository was
`zhilizheng/SBayesRC` `main` at
`2881f9a6f57374b80b4d19a09dbf9387939a2f46` (package 0.2.6). The public GCTB
repository was `jianzeng/GCTB` `master` at
`cc7fa7d765c83a89c6375946cf77fe50ba1a317e` (version 2.5.5).

Primary implementation files included manuscript `src/SBayesRC.cpp`,
`src/AnnoProb.cpp`, `src/BlockLDeig.cpp`, and `src/dist.cpp`; GCTB
`scr/model.hpp`, `scr/model.cpp`, `scr/data.cpp`, `scr/eigen.cpp`, and
`scr/stat.cpp`; and the corresponding `sblr` scalar BED, CSR, and retained
block-eigen sources. Classification terms below are:

- **same**: same model and parameterization;
- **equivalent**: a mathematically equivalent coordinate or scale change;
- **different-valid**: a deliberate different model or prior;
- **heuristic**: a robustness action not implied by the posterior;
- **defect**: implementation does not satisfy its stated contract;
- **unresolved**: evidence does not establish equivalence or validity.

## Contract matrix

| Contract | Original SBayesRC | Manuscript implementation | GCTB | `sblr` BED BayesRC | `sblr` exact CSR | `sblr` hard-sparse CSR | `sblr` retained block eigen |
|---|---|---|---|---|---|---|---|
| Data level | GWAS summary | Summary | Summary | Individual | Summary | Summary | Summary |
| Likelihood | Block projected LD likelihood | Block eigen likelihood | Eigen/block summary likelihood | Supplied selected genotypes | Full `X'X` summary likelihood | Supplied sparse `R` likelihood | Projected retained-factor likelihood |
| LD representation | Block eigen, retained positive mass | Block eigen files | Eigen LD | Packed BED | Complete disk CSR | Threshold/window disk CSR | BED-built block factors |
| Residual coordinates | Retained eigen coordinates | Block eigen coordinates | Block eigen coordinates | Sample coordinates | Marker scores | Marker scores | Concatenated retained coordinates |
| Residual variance | Eigen likelihood | Configurable implementation update | Block-specific procedures and current robustness logic | One global individual-level variance | One global summary variance | One global sparse-operator variance | One global projected variance |
| Global vs block residual variance | Block likelihood | Implementation-specific | Block-specific | Global | Global | Global | Global (**different-valid**) |
| Effect-variance parameterization | Genetic-variance/component formulation | Includes derived/robust options | Current robust derived variance available | Sampled common variance | Sampled common variance | Sampled common variance | Sampled common variance (**different-valid**) |
| Mixture scaling | Component multipliers | `gamma` multipliers | `gamma` multipliers | `gamma * vb` | `gamma * vb` | `gamma * vb` | `gamma * vb` (**same within `sblr`**) |
| Annotation stick direction | Sequential continuation probit | Component `k` uses preceding continuation sticks | Same continuation direction | Same | Same | Same | Same (**same**) |
| Alpha initialization | Marginal mixture mapped to intercepts | Reverse-stick intercept mapping | Reverse-stick intercept mapping | Same helper | Same helper | Same helper | Same helper (**same**) |
| Alpha prior | Intercept treated separately; Gaussian annotation effects | Flat intercept, Gaussian non-intercepts | Flat intercept, Gaussian non-intercepts | Flat intercept by default | Same | Same | Same (**same**, but separation qualification below) |
| `sigmaSqAlpha` prior | Variance prior on non-intercepts | Inverse-chi-square form | df 4, scale 1 | User `a,b`; Study 06 used 2,2 | Same | Same | Same (**different-valid** from current GCTB) |
| Empty-stick handling | Not a posterior definition in the paper | Empty eligible set sets non-intercepts to zero and intercept to -10 | Same family of fallback | Leaves previous value | Same | Same | Same (manuscript/GCTB action is a **heuristic**, not copied) |
| Continuous-annotation scaling | Standardized annotations | Preprocessed annotation matrix | Preprocessed annotation matrix | Explicit R preprocessing controls | Same | Same | Same (**same for Study 06**) |
| Residual rebuild | Projected likelihood state | Implementation rebuilds/checks | Rebuilds plus robustness checks | Maintained sample residual | Optional marker rebuild | Optional marker rebuild | Periodic plus final reduced rebuild |
| Operator approximation | Retained block spectrum and omitted cross-block LD | Explicit eigen approximation | Explicit eigen approximation | None for supplied samples | None when CSR is complete | Threshold/window omission | Retained positive spectral mass and block omission |
| Numerical fallback | Not a model term | Resampling and removal controls exist | Effect blow-up/problem-SNP handling exists | No residual-scale fallback | No residual-scale fallback | No residual-scale fallback | No residual-scale fallback |
| Problematic-SNP handling | Quality-control workflow | Can remove/resample flagged states | Documented removal/fallback tools | None silently | None silently | None silently | None silently |

## Annotation update comparison

The manuscript `AnnoProb::sampleFromFC` and GCTB
`ApproxBayesRC::AnnoEffects::sampleFromFC_Gibbs` use Albert--Chib latent
variables followed by coefficient-wise Gaussian updates. `sblr` uses the same
stick direction, eligible set, latent-variable model, flat-intercept convention,
and within-sweep residual updating. GCTB shuffles non-intercept coefficient
order; fixed order in `sblr` remains a valid Gibbs order.

Before this audit, `sblr` sampled extreme truncated normals by clipping normal
probabilities at `1e-12`. That can return a latent value on the wrong side of
zero when the predictor is sufficiently extreme. This was an implementation
defect, not a reparameterization. It is replaced by an exact normal/exponential
rejection kernel.

The manuscript empty-set assignment (`intercept=-10`, other coefficients zero)
is a robustness heuristic. It is not a draw from the stated flat-intercept
posterior and is therefore not copied. Study 06 exposes a different, more
fundamental case: later sticks have nonempty eligible sets but all eligible
markers continue, creating complete separation. With a flat intercept the
conditional posterior has no finite stabilizing tail. Blocking or reordering
the same Gibbs update cannot repair that model contract.

## Operator conclusion

The evidence supports four distinct roles, with qualifications:

1. Packed BED is the individual-level reference for the supplied selected
   genotypes.
2. Complete CSR is the summary-statistics reference when complete aligned LD is
   actually supplied. A CSR filename alone is not evidence of exactness.
3. Hard-sparse CSR is an explicit approximation. PSD is necessary for operator
   integrity but does not establish corrected-score or quadratic fidelity. The
   Study 05 operator was positive definite and still materially unfaithful; the
   Study 06 hard-sparse operator additionally failed sparse Cholesky and is
   therefore not a valid covariance/LD likelihood operator as constructed.
4. Retained block eigen is the canonical scalable scalar summary SBayesRC route
   in `sblr`, and is closest to the paper/GCTB likelihood. It remains distinct
   from GCTB because `sblr` uses one global projected residual variance and
   omits cross-block LD.
