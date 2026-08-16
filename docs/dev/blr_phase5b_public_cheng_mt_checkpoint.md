# Phase 5B public Cheng MT-BayesC checkpoint

Status: implemented for focused verification.

## Public capability

`mtblr_bed()` now resolves directly to the corrected general-$T$ Cheng
MT-BayesC$\Pi$ implementation qualified in Phases 4a, 4b, and 5A. The public
boundary supports:

- one common-sample packed-BED provider;
- complete finite phenotype matrices with stable sample and trait IDs;
- feasible $T\geq2$ complete binary pattern spaces with the first trait
  changing fastest;
- joint categorical activity-pattern updates, null collapse, and conditional
  completion;
- sampled Dirichlet activity-pattern probabilities;
- one authoritative inverse-Wishart $V_b$;
- fixed full SPD $V_e$ or sampled full inverse-Wishart $V_e$;
- Phase 3 task seeds, retention, static chain scheduling, convergence capture,
  and worker diagnostics;
- the Phase 5A pre-provider memory gate and compact transition diagnostics.

The public covariance-prior arguments use explicit
`degrees_of_freedom` and `scale` terminology. For dimension $T$, an
inverse-Wishart prior is proper when $\nu_0>T-1$; finite prior mean additionally
requires $\nu_0>T+1$. A proper prior without a finite prior mean remains valid.
Phenotypes must be centred or residualized by the caller when required by the
analysis; the sampler does not estimate fixed effects or subtract their
degrees of freedom.

The shared public execution spellings `nit`, `nburn`, `nthin`, `seed`,
`nchains`, `ncores`, `chain_seeds`, and `keep_chains` retain their established
meanings. Optional `chain_seeds` are chain base seeds and still pass through
seed-contract version 1; they are not copied directly into the native RNG.
`convergence = "core"` optionally retains the zero-RNG completed-iteration
capture through `keep_traces`, while `convergence = "none"` requires
`keep_traces = FALSE`. The exact incremental allocation gate is controlled by
`memory_limit_bytes`.

## Execution path

The supported path is:

```text
mtblr_bed()
  -> exact public argument validation
  -> Phase 5A BED preparation and memory preflight
  -> global marker map and common-sample provider
  -> blr_resolved_spec v1
  -> corrected general-T native sampler
  -> validated blr_raw v2
  -> shared formatted-fit finalization
```

The formatted object is a view of the raw-v2 scientific fields. The complete
validated raw object remains available through `attr(fit, "blr_raw")`; no
second result API or legacy result formatter is used.

## Output

Raw output retains named non-dropped draw, chain, marker, trait, observation,
covariance, and activity-pattern axes. It includes realised and latent effect
draws, joint and traitwise activity, markerwise pattern probabilities,
traitwise PIPs, all-active pleiotropic probabilities, sampled
$\boldsymbol\Pi$, authoritative $V_b$, truthful fixed or sampled $V_e$,
predictions, final states, retained and convergence indices, exact task seeds,
workers, compact pattern occupancy, and pattern-change counts.

The formatted fit exposes the corresponding explicit fields and the stable
unambiguous common aliases. Fixed $V_e$ has no fabricated residual-covariance
draws or posterior mean. Generic `pi`, `pis`, and `pim` fields are not created.

## Legacy replacement

The previous BED covariance-hybrid R execution body is retained temporarily
only as an explicitly named internal historical implementation because the
same source file still owns preparation helpers used by the corrected route.
No public or unified dispatch calls it. Its native entry points remain
registered temporarily for historical internal tests, but are not selected by
any supported public function.

The former internal positional `sblr()` hybrid route was renamed as historical
internal code and no longer provides a unified dispatch spelling. Corrected MT
CSR and block-eigen likelihoods are unavailable: `mtblr_csr()` and
`mtblr_block_eigen()` are no longer exported, and their internal boundary
spellings fail before legacy execution. Old hybrid fits cannot be converted to
corrected covariance output and must be rerun.

This retained unreachable source is a narrow cleanup follow-up. Removing the
large shared legacy translation unit in Phase 5B would create unrelated churn
across historical fixtures and generated interfaces without changing public
scientific behavior.

## Scientific invariants

Phase 5B changes no native transition. The completed-iteration order remains:

1. complete marker sweep;
2. Dirichlet $\boldsymbol\Pi$ update;
3. inverse-Wishart $V_b$ update;
4. optional inverse-Wishart $V_e$ update;
5. unthinned convergence capture;
6. retained capture.

Fixed $V_e$ remains a zero-residual-covariance-RNG path. Every marker is
visited in every sweep. Parallelism is across complete chains only.

## Unsupported and deferred

Phase 5B does not provide corrected CSR or block-eigen MT likelihoods,
heterogeneous or overlap-aware summary providers, missing phenotypes, fixed
effects, MT-BayesR, MT-BayesRC, annotation-informed MT models, regional
covariance, covariance templates, factor models, restricted large-$T$
patterns, or adaptive marker scheduling. Marker scheduling is potential Phase
5C work, not part of this checkpoint.

## Focused qualification

Permanent tests compare public and internal fixed-$V_e$ $T=2$ and sampled-$V_e$
$T=3$ runs under identical inputs and task seeds. They also cover raw/formatted
consistency, two-chain $T=4$ serial/parallel identity, exact covariance-axis
validation, and the unreachable CSR/block-eigen boundaries. Existing Phase
4/5 scientific, memory, provider, raw-schema, and execution tests remain the
scientific qualification foundation.
