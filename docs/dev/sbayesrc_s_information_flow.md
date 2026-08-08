# SBayesRC / SBayesRC-S Phase 4D information flow

## Scope and starting conclusion

Phase 4D starts from commit `8c8c2e88986eb0fc4d10b733ebf406006a4e89`.
Phase 4C identified `COUPLING-D1`: continuous annotation alpha and genomic
allocation feedback was the first failing coupling layer, before dynamic
annotation selection or its hyperparameters. Phase 4D is diagnostic only. It
changes no posterior, transition probability, random draw, standard-SBayesRC
default, or public output contract.

## Three distinct probability layers

For marker `i`, component `k`, and the state present when the single-marker
Gibbs update visits that marker:

1. `prior component probabilities`, pi_ik, come from the current annotation
   hierarchy;
2. `RB component probabilities`, p_ik proportional to pi_ik times the marker
   Bayes factor, are the normalized conditional probabilities used by the
   component draw;
3. the `hard allocation` is the realized categorical draw from p_i.

The aggregate active summaries are correspondingly named
`prior_active_trace`, `rb_active_trace`, and `hard_active_trace`. The RB rows
are sequential-Gibbs conditionals: each is evaluated against the state at the
moment its marker is visited. They are not simultaneous genome-wide marginal
posterior probabilities.

Existing `dm` and `comp_prob` remain hard allocation-frequency summaries.
The internal development binding adds `rb_dm` and `rb_comp_prob` only when the
Phase-4D flag is explicitly enabled.

## Soft continuation and annotation information

For stick `j`, the diagnostic uses

```text
soft eligible = sum_{k >= j} p_ik
soft success  = sum_{k >  j} p_ik
```

alongside the existing hard indicators. It streams per-iteration expected
component counts, active counts, stick eligible/success/rate summaries, and
annotation cross-moments

```text
E_aj = sum_i A_ia e_ij
S_aj = sum_i A_ia s_ij
Q_aj = sum_i A_ia^2 e_ij.
```

These are information diagnostics, not fractional probit observations or a
replacement sufficient statistic for Albert--Chib updating. Conditional
entropy, KL from the annotation prior to the marker conditional, and total
variation measure how much the current marker likelihood changes its prior.

The supplied Jian Zeng/GCTB R reference constructs the per-component
`probDelta` vector before drawing hard membership. This provides conceptual
precedent for retaining Rao--Blackwell-type component summaries. It does not
establish that active production GCTB uses those probabilities to update
alpha, and no GCTB transition is copied here.

## Instrumentation validation

The deterministic oracle passed probability normalization, non-negativity,
`rb_dm = 1 - p0`, hard and soft stick algebra, monotone eligibility, annotation
cross-moments, and finite entropy/KL/TV. A same-seed internal genomic
regression compared diagnostics off and on. Marker effects, hard allocations,
alpha, variance traces, delta, `pi_A`, `tau2`, and all other RNG-dependent
scientific state were identical. The only added object was
`chain$information_flow`. Thus the scientific kernel and RNG consumption are
unchanged.

## 160-marker hard-versus-soft audit

Four chains used seeds 20271101--20271104, 1,800 retained iterations, 400
burn-in iterations, and the committed 160-marker fixture. Ranges below are
ranges of chain means.

| configuration | prior active range | RB active range | hard active range | RB / hard R-hat | RB / hard minimum ESS | max alpha R-hat |
|---|---:|---:|---:|---:|---:|---:|
| B0 frozen allocations | 0.498 | 0.456 | 0 | 1.001 / 1.000 | 371 / 1800 | 1.002 |
| B1 fixed all included | 21.004 | 20.884 | 21.214 | 1.077 / 1.078 | 9.15 / 9.35 | 1.134 |
| B2 fixed informative subset | 37.485 | 37.113 | 37.558 | 1.151 / 1.151 | 8.42 / 8.62 | 1.136 |
| B3 full hierarchy | 36.730 | 36.566 | 36.586 | 1.165 / 1.162 | 8.33 / 8.36 | 1.192 |

B1--B3 soft stick continuation ranges, R-hats, and ESS were likewise nearly
identical to their hard counterparts. None of the 12 preregistered active/stick
comparisons satisfied at least two of: 50% range reduction, R-hat reduction to
at most 1.05, or twofold ESS improvement. Annotation cross-moment separation
also persisted: median hard/soft chain ranges were 4.79/4.73 (B1), 4.93/4.99
(B2), and 6.59/6.51 (B3); upper-tail and maximum ranges were similarly close.

At marker level, RB summaries agreed closely with hard posterior averages
within each pooled fit (mean absolute PIP differences about 0.0033--0.0035),
but did not reduce between-chain separation. For B1, median hard/RB SNP-PIP
chain ranges were 0.165/0.163; for B2, 0.248/0.237; and for B3, 0.296/0.294.
The corresponding 90th percentiles were 0.228/0.226, 0.355/0.347, and
0.342/0.340.

Mean conditional entropy was 0.763, 0.846, and 0.902 for B1--B3. Mean
conditional-to-prior KL was only 0.034, 0.038, and 0.039. The RB layer closely
tracks the chain-specific annotation prior (`cor(prior, soft)` about 0.999),
while soft and hard active traces correlate about 0.983--0.985.

The primary interpretation is `INFO-D2`: chains occupy genuinely different
conditional genomic-information regions. Hard categorical realization is not
the main explanation for the chain separation.

## Standard SBayesRC bridge

The same internal diagnostic path with the selection policy disabled follows
the standard continuous-alpha hierarchy. Its prior/RB/hard active ranges were
31.945/31.878/32.089, RB/hard active R-hats were 1.247/1.242, and minimum ESS
was 6.61/6.88. Median hard/RB SNP-PIP chain ranges were 0.483/0.478. Thus
standard SBayesRC shows the same pattern: RB conditional information remains
separated when continuous alpha/allocation chains occupy different regions.
No standard model source, RNG path, or public behavior changed.

## Information-scale ladder

The marker-count ladder used newly simulated deterministic fixtures rather
than marker replication. Seeds were derived as `20271200 + M`; annotation
prevalence, mixture components, sample-size policy, and contrast logic were
held fixed.

| M | prior active range | RB active range | hard active range | RB / hard R-hat | max alpha R-hat | mean KL | seconds |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 160 | 17.855 | 18.084 | 18.284 | 1.028 / 1.028 | 1.136 | 0.0465 | 7.59 |
| 500 | 77.115 | 76.520 | 76.670 | 1.112 / 1.110 | 1.211 | 0.0200 | 22.31 |
| 2000 | 840.421 | 839.461 | 839.667 | 2.041 / 2.042 | 2.050 | 0.00615 | 85.50 |

The increasing marker-count sequence did not stabilize either layer. RB and
hard chain separation stayed nearly identical while alpha/occupancy
convergence worsened strongly. This is not evidence that the original
160-marker fixture alone created the problem.

At fixed M=160, the preregistered annotation-contrast multipliers 0, 0.5, 1,
and 2 gave RB/hard active ranges 38.43/38.76, 23.98/24.18, 18.08/18.28, and
25.27/25.43. Their maximum alpha R-hats were 1.112, 1.078, 1.136, and 1.077.
The non-monotone response is architecture-specific, but no level showed a
hard-versus-soft divergence. Under the prespecified classification the scale
result is `SCALE-D3`: the continuous alpha/allocation geometry persists as
information increases.

## Limited full-hierarchy confirmation

The follow-up scenario was fixed before its result: use M=2000 because it
maximizes aggregate annotation information without changing the baseline
contrast or inference priors. Full-hierarchy prior/RB/hard active ranges were
580.16/579.36/580.53; RB/hard R-hats were 1.473/1.474 and alpha R-hat was
1.435. The baseline C5 conclusion therefore remains visible when selection
and hyperparameters are restored.

## Decision and next task

**SBS4D-R2:** RB/soft genomic information remains materially chain-separated;
hard allocation noise is not the main explanation.

**SCALE-D3:** the problem persists as marker information increases.

Overall Phase 4 remains **SBS4-R3**. No production sampler is qualified and
Study 07 remains blocked. The evidence supports **NEXT-B**: if development
continues, return to a genuinely global posterior-preserving continuous-alpha
/ allocation transition. Phase 4D does not implement that transition and
does not implement an outer/inner soft-allocation algorithm.
