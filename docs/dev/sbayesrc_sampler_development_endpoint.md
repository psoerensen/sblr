# SBayesRC sampler-development endpoint

## Purpose

This page is the navigation record for the continuous-alpha SBayesRC sampler
work completed before the `sblr` 0.2.0 baseline. It distinguishes supported
production behavior from posterior-correct development references. Detailed
derivations and experiment records remain in the linked documents.

## Production status

The supported scalar BayesRC/SBayesRC implementation remains the ordinary
marker/allocation and annotation-hierarchy Gibbs sampler. Retained block-eigen
SBayesR/SBayesRC uses `residual_policy = "gctb_block"` and
`block_ve_mode = "allMixVe"` by default; the block residual mathematics and
defaults are specified in
[`block_eigen_gctb_residual_contract.md`](block_eigen_gctb_residual_contract.md).
The established LD-swap controls remain supported localization diagnostics and
moves where their documented route contract applies. They are not global
learned-alpha transitions.

No PX, conditional particle-Gibbs, particle-marginal, coupling-tempering,
partial-exchange, pairwise-allocation, or coordinated alpha/allocation method
is a supported production sampler option. None is exported, enabled by default,
or part of the ordinary sampler RNG path.

## Methods evaluated

| Method | Posterior target preserved? | Purpose | Scale tested | Result | Production status | Reference document |
|---|---:|---|---|---|---|---|
| Ordinary Gibbs | Yes | Production BayesRC/SBayesRC transition | Tiny fixtures through Study 06 designs | Scientific reference; learned-alpha global mixing remains difficult | Supported default | [`bayesrc_probit_intercept_prior.md`](bayesrc_probit_intercept_prior.md) |
| Repeated hierarchy updates | Yes | Equilibrate alpha conditional on allocations | 1,500-marker BED and block schedules | Improved some alpha quantities; did not resolve allocation/annotation mixing | Development control only | [`study06_allocation_hierarchy_kernel_composition.md`](study06_allocation_hierarchy_kernel_composition.md) |
| Repeated allocation sweeps | Yes | Equilibrate allocations conditional on marker priors | 1,500-marker BED and block schedules | Unfavorable; did not resolve regime separation | Development control only | [`study06_allocation_hierarchy_kernel_composition.md`](study06_allocation_hierarchy_kernel_composition.md) |
| Coupling tempering and exchange audits | Yes for the validated complete-state construction | Move complete allocation/annotation regimes between coupling levels | Tiny oracle and short 1,500-marker BED screen | Exact target passed; all Study 06 exchanges rejected; partial-overlap audit retained | Development reference only | [`study06_bed_coupling_tempering_screen.md`](study06_bed_coupling_tempering_screen.md), [`study06_partial_exchange_feasibility.md`](study06_partial_exchange_feasibility.md) |
| Exact pair update | Yes | Coordinate two marker component/effect states | Tiny exact cases and Study 06 development runs | Too sparse to resolve occupancy mixing | Not promoted | [`bayesrc_pairwise_allocation_update.md`](bayesrc_pairwise_allocation_update.md) |
| Coordinated subset transition | Yes for the validated finite-subset construction | Move alpha with compatible allocation/effect subsets | Tiny oracle through marker-count scaling | Exponential enumeration and extensive acceptance loss prevent scaling | Development reference only | [`bayesrc_coordinated_alpha_allocation_result.md`](bayesrc_coordinated_alpha_allocation_result.md) |
| PX/sandwich | Yes | Improve latent-probit/alpha movement without changing the model | Tiny oracle and 1,500-marker retained-block diagnostic | Exact but insufficient for joint alpha/occupancy convergence | Development reference only | [`sbayesrc_block_px_transition.md`](sbayesrc_block_px_transition.md) |
| Conditional particle Gibbs | Yes | Refresh allocation/effect paths within retained LD blocks | Tiny oracle, 100-marker and 500-marker blocks | Exact and diverse, but block-local moves did not solve global sparsity movement | Development reference only | [`sbayesrc_block_particle_transition.md`](sbayesrc_block_particle_transition.md) |
| Particle-marginal alpha | Yes for the validated selected-path reference | Marginalize compatible block allocation/effect states during global alpha moves | Tiny oracle, 1,500-marker screen, representative 76-block feasibility | PMA-R3: exact global reference transition, computationally impractical at qualification scale | Development reference only | [`sbayesrc_particle_marginal_alpha_result.md`](sbayesrc_particle_marginal_alpha_result.md) |

## Endpoint

No further unrestricted continuous-alpha same-posterior sampler engineering is
currently recommended. PMA-R3 established an exact global reference transition
but production-scale computation was not justified. The ordinary sampler and
validated likelihood contracts remain the package baseline; the negative and
reference results above are retained as scientific provenance, not advertised
features.

Future annotation-selection and annotation-PIP work is a different posterior
model and must not be confused with a sampler fix for standard SBayesRC. Its
design belongs to a separate task.

## Separate MCEM inference line

The PMA-R3 endpoint applies to full-joint continuous-alpha sampler
development. Phase 5A subsequently introduced the separate SBayesRC-EM
inference line, initially qualified on CSR in
[`sbayesrc_mcem_phase5a.md`](sbayesrc_mcem_phase5a.md) and completed for the
retained block backend in
[`sbayesrc_em_phase5b.md`](sbayesrc_em_phase5b.md). MCEM estimates an
observed-data annotation MAP by alternating fixed-prior genomic Monte Carlo
blocks with a soft-probit M-step. It is not a replacement sampler for the
joint alpha posterior and does not change ordinary SBayesRC behavior.

Phase 5C adds the separate shared-annotation-selection line SBayesRC-S-EM.
Its internal implementation uses MCEM for genomic latent states,
discrete/continuous MAP annotation fitting, and a responsibility-conditioned
Laplace model distribution for `annotation_pip_eb`. CSR remains the reference
backend and retained block-eigen remains the scalable backend. See
[`sbayesrc_s_em_phase5c.md`](sbayesrc_s_em_phase5c.md). This inference line
does not replace or alter joint SBayesRC-S posterior sampling.
