# Developer-document authority map

This page records which developer documents govern the current repository and
which are evidence, plans, or historical research. Documentation is not proof
that a capability exists: current R/C++ source, executable dispatch, schemas,
generated registration where relevant, and tests are the highest authority for
implemented behavior.

## Authority hierarchy

Use evidence in this order:

1. **Current implementation truth:** source, executable dispatch, schemas,
   generated registration, and tests.
2. **Current development contracts:** maintained implementation, schema,
   convergence, backend, and test-ownership documents listed below.
3. **Approved future annotation-prior architecture:** the audit, capability
   matrix, and phased implementation plan. Proposed entries are not current
   functionality.
4. **Canonical statistical theory:** [`../methods/`](../methods/index.qmd).
5. **Qualification records:** evidence for a particular current feature; they
   do not override general contracts.
6. **Research, experimental, historical, and decision records:** provenance,
   including negative results; they never override current code, tests,
   contracts, or the approved annotation-prior direction.
7. **Practical guidance:** `docs/notes/`, which is not architecture authority.

Status terms used here are `CURRENT_CONTRACT`, `CURRENT_PLAN`,
`CURRENT_QUALIFICATION`, `RESEARCH_EXPERIMENTAL`, `HISTORICAL`, `SUPERSEDED`,
and `DECISION_RECORD`.

## Current contracts

These are the normal starting set before changing the BLR framework.

| Document | Authoritative for |
|---|---|
| [`README.md`](README.md) | developer-document authority and status |
| [`blr_development_guide.md`](blr_development_guide.md) | change discipline and documentation ownership |
| [`blr_architecture.md`](blr_architecture.md) | current high-level backend/raw/formatter architecture |
| [`blr_backend_inventory.md`](blr_backend_inventory.md) | current backend and operator ownership |
| [`blr_model_contracts.md`](blr_model_contracts.md) | current model semantics and distinctions |
| [`blr_output_schema.md`](blr_output_schema.md) | named raw schemas and formatted-output ownership |
| [`blr_convergence_contract.md`](blr_convergence_contract.md) | convergence modes, traces, diagnostics, and warnings |
| [`blr_phase2_provider_operator_checkpoint.md`](blr_phase2_provider_operator_checkpoint.md) | implemented shared marker/resource/provider adapters and native-view qualification boundary |
| [`blr_phase3_execution_checkpoint.md`](blr_phase3_execution_checkpoint.md) | activated logical-task, seed-v1, retention-v1, convergence-capture, and sampler-worker boundary |
| [`blr_phase4a_cheng_mt_bayesc_checkpoint.md`](blr_phase4a_cheng_mt_bayesc_checkpoint.md) | qualification-only corrected two-trait common-sample BED Cheng MT-BayesC$\Pi$ slice with fixed full $V_e$ |
| [`blr_phase4b_sampled_residual_covariance_checkpoint.md`](blr_phase4b_sampled_residual_covariance_checkpoint.md) | qualification-only extension with one authoritative sampled full inverse-Wishart $V_e$ transition |
| [`blr_phase5a_general_t_cheng_mt_checkpoint.md`](blr_phase5a_general_t_cheng_mt_checkpoint.md) | qualification-only general-$T$ complete-pattern Cheng MT-BayesC$\Pi$ implementation, checked $M2^T$ plus aligned packed-BED fit-memory boundary, validated compact diagnostics, and $T=2$--$4$ evidence |
| [`blr_test_ownership.md`](blr_test_ownership.md) | permanent tests and validation tiers |
| [`blr_block_eigen_contract.md`](blr_block_eigen_contract.md) | current scalar retained/dense block-eigen route contract |
| [`stblr_low_rank_operator_design.md`](stblr_low_rank_operator_design.md) | retained low-rank scalar operator mathematics |
| [`block_eigen_gctb_residual_contract.md`](block_eigen_gctb_residual_contract.md) | current scalar GCTB-compatible block residual policy |

The final three are feature-specific contracts. Read them when changing scalar
block-eigen construction, residual policy, or operator semantics.

## Approved annotation-prior restructuring

All three documents have status `CURRENT_PLAN`. They are authoritative for the
approved direction, while their current-capability statements remain subject
to executable evidence.

| Document | Role |
|---|---|
| [`annotation_prior_architecture_audit.md`](annotation_prior_architecture_audit.md) | audited current observations and approved P/Q/H/$V_b$ conclusions |
| [`annotation_prior_architecture_matrix.md`](annotation_prior_architecture_matrix.md) | audited current-versus-proposed capability map |
| [`annotation_prior_architecture_implementation_plan.md`](annotation_prior_architecture_implementation_plan.md) | gated future reorganization sequence |

In particular, proposed structured prior specifications, future
`fit$architecture`, informative theta means, external-q composition, MT-LV,
BED-LV, and arbitrary P/Q composition are not current implementation claims.

## Current qualification and reference records

These records provide evidence for implemented features. They are not general
architecture owners.

| Document | Evidence for |
|---|---|
| [`annotation_log_variance_implementation.md`](annotation_log_variance_implementation.md) | qualified version-1 ST CSR/block-eigen LV implementation and reduction gates |
| [`bayesrc_probit_intercept_prior.md`](bayesrc_probit_intercept_prior.md) | current proper probit-stick intercept decision and qualification |
| [`block_eigen_gctb_residual_implementation_result.md`](block_eigen_gctb_residual_implementation_result.md) | implementation qualification for the block residual contract |
| [`blr_prior_calibration_audit.md`](blr_prior_calibration_audit.md) | selected current calibration contract, with clearly marked historical defect audit |
| [`sbayesrc_alpha_reference_validation.md`](sbayesrc_alpha_reference_validation.md) | standard SBayesRC alpha conditional and permanent oracle guards |
| [`sbayesrc_reference_crosswalk.md`](sbayesrc_reference_crosswalk.md) | pinned external/scalar SBayesRC crosswalk |
| [`sbayesrc_sampler_development_endpoint.md`](sbayesrc_sampler_development_endpoint.md) | supported scalar SBayesRC endpoint and non-promoted sampler inventory |
| [`stblr_low_rank_gctb_crosswalk.md`](stblr_low_rank_gctb_crosswalk.md) | pinned retained-factor/GCTB scale crosswalk |
| [`stblr_low_rank_performance.md`](stblr_low_rank_performance.md) | retained-operator performance and memory qualification |

## Research and experimental records

The following 27 Markdown files have status `RESEARCH_EXPERIMENTAL`. They
record model research, sampler experiments, internal inference lines,
diagnostics, or rejected approaches. They are evidence only, even where an
internal implementation or exact transition was validated.

- `bayesrc_annotation_mixing_review.md`,
  `bayesrc_coordinated_alpha_allocation_result.md`,
  `bayesrc_coordinated_alpha_allocation_transition.md`, and
  `bayesrc_pairwise_allocation_update.md`;
- all five `sbayesrc_block_*` records;
- `sbayesrc_particle_marginal_alpha_design.md` and
  `sbayesrc_particle_marginal_alpha_result.md`;
- `sbayesrc_mcem_phase5a.md`, `sbayesrc_em_phase5b.md`, and
  `sbayesrc_mcem_phase5d_study07_diagnosis.md`;
- `sbayesrc_s_mathematical_specification.md`,
  `sbayesrc_s_cpp_implementation.md`, `sbayesrc_s_em_phase5c.md`,
  `sbayesrc_s_information_flow.md`, and `sbayesrc_s_joint_mixing.md`;
- every `study06_*.md` record: allocation/hierarchy composition and audit,
  BED coupling-tempering, large-B0 and compact-trace work, execution-unblock,
  partial-exchange feasibility, and scalar stabilization.

No PX, particle, MCEM/EM, annotation-selection, coupling-tempering, exchange,
pairwise, or coordinated experimental transition becomes a supported model or
general contract merely because its record or internal code exists.

## Historical, superseded, and decision records

| Status | Documents | Meaning |
|---|---|---|
| `HISTORICAL` | [`blr_cleanup_manifest.md`](blr_cleanup_manifest.md), [`history.md`](history.md), [`package_cleanup_after_study06.md`](package_cleanup_after_study06.md), [`history/phase22_stabilization_report.md`](history/phase22_stabilization_report.md) | cleanup and phase provenance anchored to older commits |
| `SUPERSEDED` | [`blr_model_capability_matrix.md`](blr_model_capability_matrix.md) | replaced for annotation/prior truth by the audited matrix; retained only for scoped historical context |
| `DECISION_RECORD` | `study06_allocation_hierarchy_kernel_decision.json`, `study06_alpha_hierarchy_decision.json`, `study06_bed_coupling_tempering_decision.json`, `study06_partial_exchange_decision.json` | machine/human-readable Study06 evidence; not current architecture authority |

The old cleanup manifest recorded the then-existing
`blr_block_eigen_contract.md` as transferred/retired. The file was subsequently
reintroduced and revised with the retained low-rank implementation. Its current
contract identifiers are used by source and tests, so its present status is
`CURRENT_CONTRACT`, not the historical disposition in that manifest.

## Phase 0/1 mandatory reading

Before annotation-prior Phase 0/1, read:

1. repository-root `AGENTS.md` and this authority map;
2. `blr_development_guide.md`, `blr_architecture.md`,
   `blr_backend_inventory.md`, `blr_model_contracts.md`,
   `blr_output_schema.md`, `blr_convergence_contract.md`, and
   `blr_test_ownership.md`;
3. the three approved annotation-prior architecture documents;
4. [`../methods/model_theory.qmd`](../methods/model_theory.qmd) and
   [`../methods/annotation_priors.qmd`](../methods/annotation_priors.qmd).

Also read `annotation_log_variance_implementation.md` when changing learned Q,
and the three block-eigen contracts/qualification records when the phase
touches retained operators. Research records should be consulted only for the
specific scientific question they document.

## Inventory checkpoint

At this checkpoint the classification covers all 60 developer records: 59
files directly under `docs/dev/` (55 Markdown and four decision JSON files)
and one file under `docs/dev/history/`. Counts are:

| Status | Count |
|---|---:|
| `CURRENT_CONTRACT` | 12 |
| `CURRENT_PLAN` | 3 |
| `CURRENT_QUALIFICATION` | 9 |
| `RESEARCH_EXPERIMENTAL` | 27 |
| `HISTORICAL` | 4 |
| `SUPERSEDED` | 1 |
| `DECISION_RECORD` | 4 |
