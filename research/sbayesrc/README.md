# SBayesRC research archive

SBayesRC uses annotations to modify marker-specific mixture membership. In the
investigated settings it improved causal-variant prioritisation and prediction.
Marker-level posterior quantities were more stable than raw annotation
coefficients, while reliable estimation and joint exploration of unrestricted
raw alpha remained difficult under LD.

Joint Gibbs, EM/MCEM, higher-information and sparsity controls, prior and
intercept changes, implementation comparisons, and the continuous-alpha HMC
prototype did not provide a satisfactory general solution. Same-posterior
sampler development ended at PMA-R3, and further unrestricted continuous-alpha
development is paused. This boundary does not mean SBayesRC is scientifically
invalid or unsuccessful; it separates useful downstream genomic inference from
precision and mixing claims about raw alpha.

- [`continuous_alpha_hmc/`](continuous_alpha_hmc/) contains the qualified
  research prototype formerly at `research/annotation_alpha/`.
- [`studies/`](studies/) contains historical supporting Study 06 investigations
  moved from `sblrbench`; they are not accepted benchmark capsules.
- [`tools/`](tools/) contains research and development scripts, not maintained
  package utilities or package tests.

The supported package implementation and public API remain defined by package
code, tests, and current developer contracts outside this research tree.

