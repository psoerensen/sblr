# BLR convergence contract

## Modes and controls

Modes are `auto`, `none`, `core`, and `extended`. `auto` computes only the
five core trait quantities when at least two chains are available; it does not
allocate extended traces. `none` captures no formal traces. `core` explicitly
requests `vbs`, `vgs`, `ves`, `vle`, and `vld`. `extended` adds applicable
low-dimensional groups and explicitly selected markers.

Common controls are `warn`, `rhat_threshold`, `ess_per_chain_threshold`,
`mcse_mean_over_sd_threshold`, and `keep_traces`. Extended controls are
`extended_groups`, `selected_markers`, `selected_marker_quantities`,
`full_probability_states`, `max_trace_gb`, and `allow_large_traces`.
`extended_groups` selects covariance, probability, sampled `maf_effect_s`, and
annotation groups. Selected marker IDs or one-based indices must be explicit;
`NULL` never means all markers and no all-marker shortcut exists.

## Trace ownership

Every physical trace is logical-chain-private, post-burn, and captured at the
completed-iteration posterior-accumulation checkpoint. Capture occurs every
post-burn iteration without diagnostic thinning, before task state is
discarded. It is observational: it does not consume RNG, rerun a chain, depend
on `nthin`, or depend on `keep_chains`. `keep_traces` controls only whether the
formal array remains in the fit.

The `blr_convergence_trace_bundle` schema remains version 1 with numeric values
ordered iteration × chain × scalar quantity plus stable identity descriptors.

## Mathematics

One shared scalar engine implements Blom rank normalization with deterministic
average ranks for ties, rank-normalized split R-hat, folded R-hat, their
maximum, FFT autocovariances, Geyer initial positive/monotone sequences, bulk
ESS, tail ESS, mean ESS, posterior SD, MCSE of the posterior mean, and relative
mean MCSE. Discrete state traces are not jittered. No ESS for posterior SD or
quantile/median MCSE is claimed.

## Extended quantities

MT covariance diagnostics use strict R column-major lower-triangle order and
never duplicate diagonal core quantities. Probability diagnostics consume
sampled low-dimensional simplexes or marginals, not marker-level posterior
probability matrices. Annotation diagnostics use sampled coefficients and
variances. Fixed parameters are `not_updated`; absent quantities are
`not_applicable`; diagonal-policy off-diagonal residual covariances are
`structural_zero`.

Selected marker quantities are current effective `b`, binary activity `d`,
and ordered component code where the model defines components. They are direct
indexed chain-state captures, not posterior means or latent effects.

## Memory, summaries, and warnings

The resolved plan accounts separately for core, extended, selected-state,
workspace, descriptor, and retained-trace memory before numerical execution.
Requests above `max_trace_gb` fail before sampling unless
`allow_large_traces = TRUE`, which produces one explicit memory warning and
does not thin or truncate the request.

The summary has one scalar row per diagnostic key and group overviews for core,
covariance, probability, `maf_effect_s`, annotations, and selected markers.
Advisories aggregate R-hat, ESS-per-chain, tail ESS, and relative-MCSE flags at
most once per fit. Status is `ok`, `warning`, `partial`, `unavailable`, or
`not_requested`—never a claim that a model has definitively converged.

## Phase 0 task, retention, and seed target

The preceding sections describe current convergence capture. It is unthinned,
post-burn, observational, and remains unchanged until production migration.
The exact target contracts and current-route evidence are in Sections 23--26
of [the unified framework design](blr_unified_framework_design.md).

Target retention-contract version 1 indexes post-burn transitions by
$u=1,\ldots,n_{\mathrm{sampling}}$ and retains exactly those satisfying

$$
u\bmod n_{\mathrm{thin}}=0.
$$

Thus retained count is
$\lfloor n_{\mathrm{sampling}}/n_{\mathrm{thin}}\rfloor$. The exact retained
indices are stored in the resolved specification. Current native summaries
generally retain $u=1,1+n_{\mathrm{thin}},\ldots$, while formal convergence
captures every post-burn transition. Phase 1 labels and preserves those legacy
contracts; it does not silently change current trajectories or draw counts.

Logical task IDs are `chain` for `single_trait` and `joint_multitrait`, and
`trait_id × chain` for `independent_traits`. Seed-contract version 1 uses the
exact FNV-1a/SplitMix64 algorithm and reference vectors in the unified design
and standalone fixtures. Current arithmetic seed rules remain explicit legacy
versions during Phase 1. Their eventual replacement is a deliberate trajectory
migration requiring deterministic tests and release notes; it is not described
as a posterior-target change.

Raw v2 convergence arrays preserve `draw × chain` axes plus fixed quantity
axes. Diagnostics consume no RNG. Serial and parallel execution resolve the
same logical task-seed table, and worker diagnostics separately report
requested cores, configured workers, actual sampler team size, and task worker
IDs.

## Phase 1 implementation status

Phase 1 records, but does not replace, the current ST contracts. The resolved
specification uses `seed_contract_version = 0` with final legacy task seeds and
`retention_contract_version = 0` with exact indices
$1,1+n_{\mathrm{thin}},\ldots$. Scientific traces that currently include
burn-in remain available under the explicitly named
`derived$legacy_iteration_quantities`; formal convergence capture remains the
existing unthinned post-burn sequence. Activating unified seed-contract version
1 or target retention-contract version 1 remains a Phase 3 trajectory
migration.

The Phase 1 compute contract represents `memory_limit_bytes` as `NULL`, a
finite nonnegative scalar (including zero), or positive `Inf`. It does not
alter native scheduling or retention. Original public call spellings are
validated exactly before abbreviated iteration, chain, thinning, seed, or core
controls can be partially matched by R.
