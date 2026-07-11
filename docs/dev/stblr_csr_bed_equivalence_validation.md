# CSR / individual-level (BED) equivalence validation

This note records a validation exercise checking whether the CSR
summary-statistic backend (`stblr_csr()`) and the individual-level BED
backend (`stblr_bed()`) implement the same underlying per-marker BayesC
update when run on matched in-sample data. It does not change sampler math,
public arguments, or default behavior. `src/` was not modified and nothing
here was added to the automated test suite (the checks below run real MCMC
chains and are too slow for `testthat`).

## The invariant

The individual-level BayesC marker score is

```text
score_j = x_j' e + xx_j * b_j
```

(`st_cpg_omp_individual.cpp:402`, where `e` is the current residual and
`xx_j` is the marker's `X'X` diagonal), and the CSR marker score is

```text
score_j = r_j + ww_j * b_j
```

(`st_cpg_omp_csr.cpp:70`, where `r_j` is the CSR residual-adjusted
sufficient statistic and `ww_j` is the marker's diagonal from the sparse LD
input). These are algebraically the same quantity when:

- (a) the LD reference used to build the CSR sparse-LD input *is* the GWAS
  sample (in-sample), so `ww` and the off-diagonal LD entries are exactly
  `X'X` rather than an external-reference approximation;
- (b) the sparse LD is the *full* `X'X` (`r2_threshold = 0`, no distance or
  windowing filters), so no entries are dropped; and
- (c) both backends are seeded and swept in a way that makes their sampling
  paths comparable.

Under (a) and (b), the two backends' posterior-mean marker effects (`bm`)
should agree closely. A real mismatch here would point at a bug in the
score/diagonal/operator wiring shared between the two backends, not at
approximation error, and should not be papered over with a looser tolerance.

## Method

A standalone script built a tiny in-sample fixture (250 individuals, 20
markers, 4 causal, no missing genotypes), wrote it as a real `.bed` file,
and derived both backends' inputs from the *same* underlying data:

```r
stats    <- make_summary_stats(Glist, y, scale = TRUE)
Glist_ld <- make_sparse_ld(
  Glist, rows = seq_len(n),        # LD reference == GWAS sample
  r2_threshold = 0,                # full X'X, no thresholding
  max_distance_bp = 0,
  max_distance_variants = 0,       # no windowing
  allow_full_ld = TRUE
)
```

`method = "bayesc"`, `full_sweep_every = 1` (disables BED's scheduled
candidate-skipping so both backends visit every marker every iteration),
and default priors (`h2 = 0.3`, `adjE = 0.9`, `nub = nue = 4`,
`pi_prior_mean = 0.001`, `pi_prior_strength = 5e5`) were used on both sides.

Three diagnostics were run, in order of increasing reliance on RNG:

1. **Diagnostic 0 (deterministic, no RNG).** `stblr_bed(..., updateB =
   FALSE, return_wy = TRUE)` keeps `b` fixed at 0 for the whole run, so the
   backend's internal `wy` (the `x_j'y` term feeding the score) never moves
   from its initial value. Comparing that directly against `stats$wy` from
   `make_summary_stats()` isolates the shared sufficient-statistic wiring
   from all sampling behavior.
2. **Diagnostic A (matched seed).** Both backends run with the same `seed`
   and `full_sweep_every = 1`, compared after 1 iteration and after 50
   iterations with no burn-in.
3. **Diagnostic B (multi-seed ensemble).** 20 independent seeds per
   backend, `nit = 800`, `nburn = 200`. Per-marker posterior means are
   compared against each other's Monte Carlo standard error
   (`combined_se = sqrt(se_bed^2 + se_csr^2)`, `z = diff / combined_se`).

## Findings

**Diagnostic 0 — exact match.** `fit_bed$wy` vs `stats$wy`: max absolute
difference was `0`. The deterministic sufficient statistic that feeds the
score is bit-identical between the two backends under in-sample, full-LD
conditions. This is the actual operator/diagonal wiring the invariant is
about, and it checks out.

**Diagnostic A — diverges from iteration 1, but not because of (1).** With
a matched `seed`, `bm` differs by up to ~0.13 after a single iteration
(tolerance was 1e-4) and the divergence persists at 50 iterations. Given
Diagnostic 0, this is not a scoring bug. Reading the sampler source
explains it directly: the two backends derive their internal RNG seed from
the public `seed` argument differently.

- CSR, single trait/chain: `seed + 1000003 * (trait + 1)`
  (`stblr_trait_seed()` in `src/st_chain_utils.h:7-11`, used from
  `src/st_cpg_omp_csr.cpp:1277`).
- BED scheduled-chain BayesC: `seed + 1000003 * (trait + 1) + 9176 * (chain
  + 1)` (`src/st_cpg_omp_individual_scheduled_chains.cpp:837-839`).

For a single chain these differ by a constant `9176 * (chain + 1)` offset,
so a "matched" `seed` argument does not produce a matched `mt19937` stream
between the two backends. This is a seed-derivation mismatch, not a
score/diagonal/operator bug, and it falls squarely under "impossible
because of draw-order/update-order differences" rather than something to
work around by hand-deriving a compensating seed.

**Diagnostic B — agrees within Monte Carlo error.** Across 20 seeds per
backend, `max|bed_mean - csr_mean| = 0.0028`, `max|z| = 2.08`, and 0/20
markers exceeded `|z| > 3`. This is consistent with the two backends
sampling from the same posterior.

## Conclusion

The invariant holds. No bug was found in the CSR / individual-level
score, diagonal, or operator wiring. The two backends' sufficient
statistics are exactly equal under in-sample, full-LD conditions
(Diagnostic 0), and their posterior-mean marker effects agree within Monte
Carlo error across independent seeds (Diagnostic B). The apparent
matched-seed divergence (Diagnostic A) is fully explained by the two
backends' RNG seed derivation formulas, not by the update math.

## Known related asymmetry (not fixed here)

The BED scheduled-chain seed formula includes a `9176 * (chain + 1)` term
that has no counterpart in CSR's single-chain seed formula
(`stblr_trait_seed()`, used when CSR runs with `nchains == 1`; CSR's
multi-chain path, `stblr_chain_seed()`, does include an analogous
`9176 * (chain + 1)` term). This means CSR and BED never produce comparable
RNG streams for the "same" `seed`, even in the single-chain case. This may
be intentional (to decorrelate chains across backends), but it is worth a
deliberate decision rather than an implicit side effect. No change was made
as part of this validation; see the native-code rules in `AGENTS.md`
(no sampler/RNG changes without an explicit request and explanation).

## Reproducing this check

The validation script is not part of the repository or the test suite (per
`AGENTS.md`, long-running MCMC should not be added to `testthat`). To
reproduce, write a script that:

1. builds a small in-sample `.bed` fixture and phenotype (see
   `tests/testthat/test-stblr-bed-interface.R` for the BED-packing
   convention used by other fixtures in this repo);
2. calls `make_summary_stats()` and `make_sparse_ld(r2_threshold = 0,
   max_distance_bp = 0, max_distance_variants = 0, allow_full_ld = TRUE,
   rows = <all in-sample rows>)` on the same `Glist`;
3. runs the three diagnostics above with `stblr_bed()` and `stblr_csr()`.

Keep `nit`/`nseeds` small enough to run in seconds to minutes; the ensemble
diagnostic does not need long chains to detect a systematic wiring bug,
only enough replicates to bound Monte Carlo error.
