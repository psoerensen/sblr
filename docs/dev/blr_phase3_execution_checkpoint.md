# Phase 3 logical-task execution checkpoint

## Status and scope

**Status:** `READY FOR INDEPENDENT VERIFICATION`

This checkpoint records the bounded Phase 3 implementation of shared logical
tasks, seed-contract version 1, retention-contract version 1, observational
convergence capture, and sampler-local worker diagnostics. It does not change
a posterior target, marker transition, likelihood operator, variance update,
or within-chain update order. Current MT samplers remain legacy.

## Source trace

| Route | Pre-Phase-3 topology | Seed and retention evidence | Phase 3 decision |
|---|---|---|---|
| ordinary CSR BayesC/BayesR | static trait-major trait × chain tasks | task-private `std::mt19937`; legacy first-post-burn retention | activate version 1 |
| packed-BED BayesC/BayesR | static chain/trait work items with task-private state | existing per-job RNG and legacy local retention | activate with canonical trait × chain identity mapping |
| block-eigen BayesC/BayesR | existing CSR scalar task loop with representation adapter | same scalar task ownership; representation-specific residual policy | activate version 1 |
| fixed-marker and learned-logistic BayesC | chains assembled outside trait-parallel sampler regions | chain-local seeds and sampler-local trait workers | activate version 1 |
| scheduled CSR, log-variance, group, BayesRC/SBayesRC | route-specific scheduling or additional policies | not qualified under a common version-1 comparison in this phase | retain version 0 |
| current MT routes | one joint chain per task, with legacy covariance behavior | inventory only | retain legacy; no dispatch change |

The canonical ordinary CSR BayesC path is the typed
`CsrBayesCExecutionInput`/`run_csr_bayesc()` engine. Phase 3 controls therefore
enter that engine as plain validated values rather than being added only to an
obsolete or bypassed loop. BED, CSR BayesR, block-eigen, and annotation routes
reuse their existing task loops and mutable task-owned state.

## Logical tasks and seeds

Canonical task order is:

- `single_trait`: increasing zero-based chain index;
- `independent_traits`: declared stable trait-ID order, then increasing chain
  index;
- `joint_multitrait`: increasing chain index.

R resolves the final task table before native dispatch. Seed-contract version
1 uses the frozen UTF-8 FNV-1a, ordered SplitMix64, and uint32 xor-fold
algorithm. Values above signed integer range remain exact R doubles until the
native parser validates and converts them to `std::uint32_t`. User seed zero is
valid. Worker identity, core count, and operator representation never enter
seed derivation. An explicit final task-seed table bypasses derivation only
after identity, shape, name, integer, and uint32-range validation.

The eight frozen R/native reference values in the Phase 0 fixtures pass
unchanged. Stable trait IDs, rather than trait positions, identify
independent-trait tasks.

## Retention and convergence

For post-burn transition $u=1,\ldots,n_{\mathrm{sampling}}$, version 1 retains
exactly

$$
u\bmod n_{\mathrm{thin}}=0.
$$

Resolved indices are validated before execution and used for posterior
accumulation. The final state is still the state after the last transition.
Version 0 remains available only on explicitly declared legacy routes and for
controlled qualification.

Convergence capture remains every completed post-burn iteration, unthinned,
task-private, observational, and zero-RNG. `keep_traces` affects only whether
the formal bundle remains in the formatted fit. A one-chain trace can be
retained for inspection but cannot provide between-chain convergence evidence.

## Scheduler and worker diagnostics

Existing static OpenMP task loops are retained. Every task owns its RNG,
effects, states, residuals, workspaces, accumulators, convergence capture, and
diagnostic slots. R/Rcpp objects are parsed before and constructed after worker
regions.

Native diagnostics distinguish:

- requested cores;
- configured workers;
- actual team size observed inside the sampler region;
- zero-based worker ID for every canonical logical task;
- scheduler version and logical task order;
- runtime OpenMP availability and maximum;
- zero diagnostic RNG draws.

The actual team and worker identifiers are written to task-owned slots, so the
capture is race-free and does not influence scheduling or scientific output.

## Qualification evidence

| Check | Result |
|---|---|
| frozen seed vectors | all eight R and native uint32 values exact |
| retention examples | `(1,1) -> 1`, `(5,2) -> 2,4`, `(4,4) -> 4`, `(3,4) -> none` |
| scheduler neutrality | canonical CSR raw scientific output bitwise identical under equal final seed and equal retained indices |
| serial/parallel | identical posterior summaries and resolved task seeds for the same task table |
| diagnostics | actual CSR sampler team size and per-task worker IDs observed; zero RNG draws |
| convergence | retained one-chain trace has every post-burn iteration while posterior draws use indices 2 and 4 |
| earlier contracts | Phase 0 fixture, Phase 1 schema, and Phase 2 provider/operator suites pass |

Timing fields and the new observational worker record are excluded from the
scientific equality comparison. No frozen version-0 trajectory is rewritten
as a version-1 expectation.

## Activation and legacy boundary

Version 1 is active for newly created eligible ST fits on:

- ordinary CSR BayesC and BayesR;
- packed-BED BayesC and BayesR;
- block-eigen BayesC and BayesR without log-variance policy;
- fixed-marker BayesC;
- learned-logistic BayesC.

Scheduled CSR, log-variance, group, BayesRC/SBayesRC, and all current MT routes
remain explicitly version 0. Joint-MT task-plan construction is structural
only and never invokes the current MT covariance sampler. Multi-chain raw-v2
promotion remains gated because Phase 3 does not fabricate unavailable native
final-chain state.

Eligible raw-v2 objects record exact task IDs, final seeds, retained indices,
convergence indices, contract versions, and worker diagnostics. Legacy
full-iteration scientific traces remain explicitly identified and are not
silently relabelled as canonical retained draws.

## Deferred work

- migrate additional ST policies only after route-specific scheduler
  neutrality qualification;
- build-time Git provenance injection;
- corrected Cheng joint-MT sampling and full residual-covariance policy;
- heterogeneous-provider likelihood accumulation.

These items do not weaken or broaden the Phase 3 checkpoint.
