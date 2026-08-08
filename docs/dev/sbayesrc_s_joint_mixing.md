# SBayesRC-S Phase 4C joint genomic mixing

## Scope and frozen target

Phase 4C starts from commit `ccce30621de449da837f83e1717517ce52d8333f`,
where the proper-intercept R hierarchy, its C++ implementation, and the
internal CSR genomic target were validated, but the 160-marker multichain
screen ended at `SBS4B2-R2`. This audit changes neither the SBayesRC-S target
nor standard SBayesRC. It adds internal-only hyperparameter-freeze controls
and compact selection traces used by the development qualification binding.

The preregistered final gates were maximum annotation-PIP range 0.10, key
R-hat 1.05, and minimum per-chain ESS 100. All genomic comparisons used the
same 160-marker fixture, scientific inputs, and four dispersed annotation
initial states.

## 4C-A coupling ladder

| rung | max PIP range | max alpha R-hat | min alpha ESS | M R-hat | min M ESS | active R-hat | min active ESS | max component R-hat |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| C0 frozen allocations | 0.039 | 1.002 | 87.5 | 1.002 | 89.8 | n/a | n/a | n/a |
| C1 fixed delta, all included | 0 | 1.098 | 6.27 | 1.000 | 1800 | 1.086 | 7.03 | 1.086 |
| C1 fixed delta, informative subset | 0 | 1.169 | 5.83 | 1.000 | 1800 | 1.113 | 6.02 | 1.228 |
| C2 dynamic delta, fixed pi_A/tau2 | 0.348 | 1.138 | 4.69 | 1.024 | 18.0 | 1.017 | 5.16 | 1.017 |
| C3 dynamic pi_A, fixed tau2 | 0.587 | 1.030 | 6.65 | 1.227 | 6.41 | 1.005 | 6.67 | 1.036 |
| C4 dynamic tau2, fixed pi_A | 0.394 | 1.120 | 4.90 | 1.037 | 6.60 | 1.127 | 4.71 | 1.127 |
| C5 full hierarchy | 0.515 | 1.192 | 4.70 | 1.171 | 4.31 | 1.162 | 8.36 | 1.162 |

C0 is the positive control: conditional on known allocations, chains agree and
selection states switch repeatedly. The first material deterioration is C1,
before delta, pi_A, or tau2 feedback is permitted. Phase 4C therefore assigns
`COUPLING-D1`: continuous alpha/allocation coupling. Dynamic selection and the
hyperparameters amplify or redistribute the problem (`COUPLING-D2/D3/D4/D5`)
but do not create it.

State geometry supports this classification. In the fixed-informative-delta
C1 chains, mean prior expected active counts were 87.0--100.7 and realized
active counts were 86.6--100.6; within-chain correlations were 0.982--0.988.
In C5, the corresponding chain means were 93.7--112.8 and 93.8--112.8, with
correlations 0.983--0.991. Delta switches changed expected activity by about
3.88--4.39 markers and realized activity by 5.62--5.91 markers, similar to
ordinary iteration-to-iteration movement. The separated modes are therefore
carried by continuous hierarchy/occupancy feedback, not an absorbing delta
state or unusually discontinuous delta toggle.

## Chain-length screen

| retained-length multiplier | max PIP range | max alpha R-hat | min alpha ESS | M R-hat | min M ESS | active R-hat | min active ESS | seconds |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1x | 0.515 | 1.192 | 4.70 | 1.171 | 4.31 | 1.162 | 8.36 | 14.4 |
| 2x | 0.321 | 1.183 | 7.10 | 1.063 | 12.6 | 1.192 | 10.0 | 30.1 |
| 4x | 0.171 | 1.055 | 13.0 | 1.014 | 23.9 | 1.052 | 13.5 | 59.4 |

ESS growth is far below the length multiplier for the limiting alpha and
active-count summaries. Four times the work still misses the PIP, R-hat, and
ESS gates. This is structural local inefficiency, not evidence for
`SBS4C-R1`.

## 4C-B invariant schedule screen

Every schedule is a composition of existing exact Gibbs kernels and therefore
preserves the same target.

| H:G | max PIP range | max alpha R-hat | min alpha ESS | M R-hat | min M ESS | active R-hat | min active ESS | max component R-hat | seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1:1 | 0.515 | 1.192 | 4.70 | 1.171 | 4.31 | 1.162 | 8.36 | 1.162 | 14.4 |
| 2:1 | 0.255 | 1.034 | 6.93 | 1.028 | 15.5 | 1.007 | 7.01 | 1.013 | 23.0 |
| 5:1 | 0.163 | 1.038 | 16.5 | 1.005 | 22.4 | 1.025 | 17.5 | 1.064 | 47.7 |
| 1:2 | 0.326 | 1.062 | 6.38 | 1.017 | 8.62 | 1.034 | 7.16 | 1.053 | 18.5 |
| 1:5 | 0.489 | 1.080 | 8.29 | 1.017 | 9.19 | 1.018 | 9.33 | 1.074 | 13.8 |

No schedule passes. H=5 improves several R-hats but retains a 0.163 PIP
range, component R-hat 1.064, and limiting ESS below 23 while taking 3.3 times
the baseline runtime. Repeated genomic sweeps do not resolve the state
separation. The decision is `SCHEDULE-R2`.

## 4C-C exact-kernel assessment

### Existing collapsed annotation block: KERNEL-R2

The current update draws `(delta_j, alpha_j1:K)` from its exact conditional
given latent z and the remaining hierarchy, then performs the exact blocked
Gaussian redraw of all selected coefficients. Its posterior correctness is
protected by the Phase-1 finite-state oracle and the proper-intercept refresh.
Frequent switches in C0--C5 show that it is active, but C1 proves that no
selection-only improvement can remove the first failing coupling layer.

### Global delta block conditional on z: KERNEL-R2

For tiny J, enumerate the `2^J` models and let every row of the transition
matrix equal the exact conditional model-probability vector `p(delta | z)`.
Then

```text
p(delta) K(delta, delta') = p(delta) p(delta')
                         = p(delta') K(delta', delta),
```

so detailed balance and stationarity are exact. The permanent tiny guard has
stationarity and detailed-balance error below `1e-14`; 25,000 sampled states
match every exact model probability within 0.01. The underlying Phase-1 oracle
had PIP error 0.00643 and total-variation distance 0.00666. This block is not
integrated into the genomic backend: the existing sampler already redraws the
continuous coefficients jointly, C1 fails with delta fixed, and enumeration
scales as `2^J`. It is exact but not a relevant practical remedy.

### Annotation plus allocation subset: KERNEL-R4

For a proposal from state `x=(delta,alpha,c,beta)` to `x'`, an exact MH ratio
would require

```text
min(1, target(x') q(x | x') / (target(x) q(x' | x))).
```

If only allocation subset `S` is refreshed, every unchanged marker outside S
still contributes `p(c_i | alpha') / p(c_i | alpha)` when alpha moves. Thus a
useful global alpha change retains the genome-wide mismatch unless S is
essentially global. No scalable selected-state proposal with computable
forward and reverse density was established here, so no heuristic was coded.

### Global compatible allocation/effect refresh: KERNEL-R3/R4 boundary

The existing standard-SBayesRC development references already validate the
two relevant exact constructions: finite coordinated subsets lost acceptance
and required exponential enumeration, while selected-path particle-marginal
alpha updates were exact but computationally impractical (`PMA-R3`). The C1
all-included bridge shows that SBayesRC-S contains the same continuous
alpha/allocation geometry even before selection is dynamic. Extending the
particle-marginal target to random delta and the proper intercept would require
a fresh selected-path derivation, so it is not claimed exact for SBayesRC-S in
this task (`KERNEL-R4`). Repeating a known production-scale-impractical global
particle construction is not justified by the 160-marker screen.

No new coordinated genomic kernel is promoted. This is a principled stopping
point: the only exact blocks validated specifically for SBayesRC-S do not
address the diagnosed layer, while the relevant global continuous-alpha move
has no demonstrated practical implementation.

## Decision

There is no selected schedule or kernel for 4C-D requalification. Therefore
the correlated-annotation and moderate-J genomic screens, which are gated on
the primary fixture passing, were not run. Empty-stick support remains proper;
the tested schedules recorded finite entry/exit diagnostics and introduced no
absorbing state.

**Phase 4C: SBS4C-R5 — joint genomic mixing remains unresolved.**

**Overall Phase 4: SBS4-R3 — valid genomic target, practical production
sampler unresolved.** Study 07 remains blocked. A future task should proceed
only with a newly justified scalable exact global continuous-alpha/allocation
strategy; it should not add more local sweep tuning or change the validated
SBayesRC-S model.
