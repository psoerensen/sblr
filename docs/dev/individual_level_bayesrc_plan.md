# Individual-Level BayesRC Plan

This developer note tracks the staged internal work toward a possible
individual-level BayesRC backend. It does not define or expose a public model.

## Sequence 1 Status

- likelihood-independent probit stick-breaking annotation-prior utilities were
  extracted into `src/st_bayesrc_annotation_prior.h`
- the existing CSR SBayesRC sampler remains the only consumer
- no individual-level BED BayesRC backend was added
- there are no public R API changes and no new statistical model
- the existing SBayesRC component ordering, annotation updates, RNG draw order,
  and formatted output are preserved

## Sequence 2 Status

- added the internal native full-sweep BED BayesRC backend
  `stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc`
- the backend uses the packed-BED BayesR likelihood with the shared probit
  stick-breaking annotation prior
- the public `stblr_bed()` API is not extended in this sequence
- every MCMC iteration visits every selected marker; adaptive null scheduling
  is intentionally absent
- annotation rows are supplied in the final selected BED marker order defined
  by the concatenated `cls` marker map; the low-level caller is responsible for
  applying that same index to `A`
- `selection_s`, LD-swap, annotation-dependent effect variances, and
  multivariate covariance sampling remain deferred

## Sequence 2A Status

- extracted shared packed-BED BayesR-family utilities into
  `src/st_bed_bayesr_common.h`
- ordinary BED BayesR and internal BED BayesRC now include the shared header
- removed the macro-renamed inclusion of the complete BayesR `.cpp` file from
  the BayesRC translation unit
- no public API changes were made
- ordinary BayesR adaptive scheduling remains unchanged
- internal BayesRC remains structurally full-sweep with adaptive skipping off
- raw schemas and formatted outputs remain unchanged

## Sequence 2B Status

- added fixed-prior reduction coverage against ordinary fixed-pi BED BayesR
- added learned annotation-prior and directional enrichment smoke coverage
- added multiple-trait, multiple-chain, retained-chain, repeated-seed, and
  OpenMP thread-count reproducibility checks
- added component probability, `dm`, `dm_component_mean`, `ncomp`, and `pis`
  identity checks
- added annotation output dimension and final marker-prior checks
- added variance trace, CPO, thinning, optional `wy`/residual, native input,
  full-sweep diagnostic, and raw schema validation
- retained regression coverage for ordinary BED BayesR and CSR SBayesRC through
  their existing focused test groups
- annotation rows remain an internal caller responsibility: `A` must already
  follow the final selected BED marker order represented by `cls`; a future
  public wrapper must perform marker-ID alignment before calling native code
- the backend remains internal; public `stblr_bed(method = "bayesrc")` routing
  is not implemented
