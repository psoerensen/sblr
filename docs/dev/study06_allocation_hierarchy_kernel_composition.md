# Study 06 allocation--hierarchy kernel-composition audit

## Status and scope

This package-development audit tests whether repeating the existing scalar
BayesRC/SBayesRC Gibbs blocks improves communication between annotation-prior
and component-occupancy regimes. It does not change the posterior, package
defaults, Study 06, or the role of any LD operator. The formal Study 06
qualification remains failed and the final benchmark remains unauthorized.

The audit started from `sblr` commit
`8d0ad4c811c5543b62d88e5c264bf27f49aea4ed` (version 0.2.0) and the clean,
read-only `sblrbench` evidence checkout at
`de31f62e182d8540488d4135df4c58f052a515d9`. The immutable informative v2
truth has specification hash
`241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56` and
truth hash
`169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.
Marker order, causal set, effects, phenotype, summary statistics, BED data,
blocks, eigen inputs, and annotations were checked by the existing Study 06
semantic checkpoint before fitting.

## Existing transition and invariant target

Write the state as non-annotation state `x` and hierarchy state `h`. Kernel A
is the existing complete non-annotation transition: the ordinary marker
allocation/effect sweep and residual bookkeeping followed by the existing
effect-, residual-, and other BayesR variance updates. Throughout A, `alpha`,
`sigmaSqAlpha`, and the marker mixture probabilities are fixed. Kernel H is
the existing hierarchy transition: construct stick outcomes and eligible
sets, sample latent probit variables, update the proper-prior intercept and
non-intercept coefficients, sample `sigmaSqAlpha`, and reconstruct normalized
marker component probabilities. Throughout H, allocations and the remaining
model state are fixed.

For target density `pi(x,h)`, a complete A transition is invariant for
`pi(x | h)` and hence for `pi`; a complete H transition is invariant for
`pi(h | x)` and hence for `pi`. Therefore powers `A^a` and `H^h`, and their
fixed state-independent composition, leave the same joint posterior invariant.
No extra transition is a deterministic recomputation: every A and H repeat
resamples all random quantities belonging to that block.

The historical routes differ only in statement order. BED used A then H;
CSR/block eigen placed H between its marker transition and variance updates.
The implementation retains an exact legacy branch for `(a,h)=(1,1)`, including
RNG order. Non-default diagnostic schedules use complete A repeats followed by
complete H repeats and retain one state after both blocks. This is a different
valid scan, not a different model.

## Implementation contract

Two dot-prefixed diagnostic controls were added to the scalar BED and scalar
CSR/block-eigen R entry points:

```text
.diagnostic_allocation_updates_per_cycle = 1L
.diagnostic_annotation_updates_per_cycle = 1L
```

They are positive integers and are passed in the already resolved annotation
control matrix, avoiding another native argument. The native shared annotation
control owner records both counts. Values of one take the exact legacy path;
the controls are not documented public tuning parameters and no default or
public result field changed. A disabled/default run and an explicit 1/1 run
are bitwise identical under the same seed. Component-trace collection also
leaves draws unchanged. Diagnostics do not consume sampler RNG.

Tiny BED and retained-low-rank tests covered S1, H5, H20, A5, and A20. They
verified finite annotation state, finite normalized component probabilities,
valid component indices and trace dimensions, retained-operator residual
integrity, deterministic logical-chain results, and rejection of invalid
schedule counts.

## Registered experiment

All fits used four chains, fit seed 701020, chain seeds 701121, 701222, 701323,
and 701424, 9,000 recorded outer iterations, burn-in 3,000, and 6,000 retained
states. Complete 1,500-marker component traces were retained locally. The
retained block-eigen fits used `representation="low_rank"` and
`eigen_prop=0.995`; all 100 modes were retained in each of 15 100-marker
blocks. Consequently this audit does not test meaningful eigenvalue
truncation.

| Route | Schedule | A/state | H/state | Total A | Total H | Runtime (s) |
|---|---:|---:|---:|---:|---:|---:|
| BED | S1 | 1 | 1 | 36,000 | 36,000 | 433.0 |
| BED | H5 | 1 | 5 | 36,000 | 180,000 | 696.4 |
| BED | H20 | 1 | 20 | 36,000 | 720,000 | 1,755.6 |
| BED | A5 | 5 | 1 | 180,000 | 36,000 | 2,074.8 |
| BED | A20 | 20 | 1 | 720,000 | 36,000 | 4,165.3 |
| Block eigen | S1 | 1 | 1 | 36,000 | 36,000 | 678.2 |
| Block eigen | H5 | 1 | 5 | 36,000 | 180,000 | 886.4 |
| Block eigen | H20 | 1 | 20 | 36,000 | 720,000 | 1,673.2 |
| Block eigen | A5 | 5 | 1 | 180,000 | 36,000 | 901.9 |
| Block eigen | A20 | 20 | 1 | 720,000 | 36,000 | 1,716.4 |

## Convergence

The unchanged contract is maximum R-hat <= 1.01, minimum bulk and tail ESS >=
400, and maximum relative MCSE <= 0.05. The registered table contains 35 raw
hierarchy, occupancy, and variance quantities per fit.

| Route | Schedule | Passed | Failed | max R-hat | min bulk ESS | min tail ESS | max rel. MCSE |
|---|---:|---:|---:|---:|---:|---:|---:|
| BED | S1 | 4 | 31 | 1.084 | 34.6 | 33.0 | 0.177 |
| BED | H5 | 12 | 23 | 1.023 | 278.0 | 180.6 | 0.078 |
| BED | H20 | 21 | 14 | 1.011 | 320.2 | 183.2 | 0.091 |
| BED | A5 | 4 | 31 | 1.062 | 51.0 | 67.4 | 0.142 |
| BED | A20 | 3 | 32 | 1.173 | 17.3 | 17.2 | 0.340 |
| Block eigen | S1 | 1 | 34 | 1.141 | 23.9 | 39.0 | 0.206 |
| Block eigen | H5 | 3 | 32 | 1.092 | 32.8 | 50.5 | 0.194 |
| Block eigen | H20 | 6 | 29 | 1.043 | 103.3 | 137.7 | 0.105 |
| Block eigen | A5 | 5 | 30 | 1.169 | 17.5 | 13.2 | 0.257 |
| Block eigen | A20 | 0 | 35 | 1.234 | 13.5 | 14.1 | 0.341 |

BED H20 brought all 12 raw alpha coefficients inside the contract (maximum
R-hat 1.0085, minimum bulk ESS 560, minimum tail ESS 624, maximum relative
MCSE 0.045). It did not bring occupancy or the implied expected active count
inside the contract: their minimum bulk ESS values were approximately 327 and
346, and effect-variance tail ESS was approximately 194. Block-eigen H20 also
improved the hierarchy but retained deficient alpha mixing (minimum bulk/tail
ESS 222/218 and maximum R-hat 1.0206). Derived stick-prior summaries improved
from 0/15 passing to 12/15 for BED H20 and from 0/15 to 3/15 for block H20.
Thus hierarchy repetition can equilibrate `alpha` conditional on a local
allocation regime without making the joint regime exploration reliable.

The stick-specific maxima/minima below make the late-stick behavior explicit.
Each cell is `maximum R-hat / minimum bulk ESS`; sticks 0, 1, and 2 are the
sequential component-0, component-1, and component-2 continuation sticks.

| Fit | alpha stick 0 | alpha stick 1 | alpha stick 2 | sigmaSqAlpha stick 0 | sigmaSqAlpha stick 1 | sigmaSqAlpha stick 2 |
|---|---:|---:|---:|---:|---:|---:|
| BED S1 | 1.084/34.6 | 1.074/57.1 | 1.076/50.6 | 1.006/746.8 | 1.038/256.7 | 1.044/113.7 |
| BED H5 | 1.023/278.0 | 1.020/339.1 | 1.011/335.2 | 1.003/2094.8 | 1.005/921.6 | 1.007/761.1 |
| BED H20 | 1.009/560.0 | 1.003/843.3 | 1.005/611.7 | 1.001/4671.7 | 1.003/2410.9 | 1.002/1684.7 |
| BED A5 | 1.062/51.0 | 1.033/107.5 | 1.047/79.5 | 1.006/534.0 | 1.013/265.2 | 1.030/99.8 |
| BED A20 | 1.071/46.2 | 1.173/17.3 | 1.061/72.5 | 1.012/320.0 | 1.159/17.5 | 1.020/149.2 |
| Block S1 | 1.141/23.9 | 1.060/78.2 | 1.021/166.0 | 1.017/243.6 | 1.010/505.9 | 1.005/351.1 |
| Block H5 | 1.049/94.0 | 1.040/129.4 | 1.015/245.4 | 1.015/409.7 | 1.011/263.1 | 1.005/852.3 |
| Block H20 | 1.021/222.2 | 1.018/257.3 | 1.014/364.1 | 1.007/857.5 | 1.004/1031.6 | 1.007/510.7 |
| Block A5 | 1.169/17.5 | 1.084/39.1 | 1.017/146.1 | 1.047/69.1 | 1.036/97.3 | 1.006/581.5 |
| Block A20 | 1.096/39.7 | 1.234/13.5 | 1.091/35.0 | 1.003/339.7 | 1.158/17.1 | 1.016/352.4 |

Allocation repetition did not help. A5 and especially A20 reduced chain
agreement and ESS on both routes. This rules out the simple explanation that
ordinary allocation sweeps merely need to run longer under each frozen set of
marker priors.

## Occupancy regimes and movement

Regimes were registered before inspection as 0--49, 50--74, 75--99,
100--149, 150--199, and 200+ active markers. The retained chain mean ranges and
lag-50 active-count autocorrelation ranges were:

| Route | Schedule | Chain mean active-count range | lag-50 ACF range |
|---|---:|---:|---:|
| BED | S1 | 53.9--59.4 | 0.319--0.499 |
| BED | H5 | 54.3--57.9 | 0.117--0.311 |
| BED | H20 | 55.0--56.3 | 0.147--0.281 |
| BED | A5 | 54.3--58.1 | 0.282--0.592 |
| BED | A20 | 54.3--65.8 | 0.274--0.644 |
| Block eigen | S1 | 122.5--142.3 | 0.545--0.691 |
| Block eigen | H5 | 127.3--181.4 | 0.422--0.628 |
| Block eigen | H20 | 128.6--136.9 | 0.390--0.567 |
| Block eigen | A5 | 135.4--151.8 | 0.589--0.774 |
| Block eigen | A20 | 104.2--142.8 | 0.511--0.612 |

H20 materially improved overlap, especially for BED, but did not meet the
joint convergence contract. Local marker movement was not absent: BED fits
changed about 59--62 allocations per recorded state and block fits about
150--164; null entries/exits averaged about 27 for BED and 67--72 for block.
The failure is slow communication between aggregate regimes despite frequent
local changes. Fixed-band proportions, transition matrices, dwell lengths,
and equivalent expected-active summaries are retained in the ignored local
evidence directory.

## Runtime-adjusted efficiency

For BED alpha, median bulk ESS/s increased from about 0.249 (S1) to 0.604
(H5) and 0.617 (H20). Expected-active bulk ESS/s was highest at H5 (about
0.705 versus 0.532 for S1) and fell to about 0.197 for H20. For block eigen,
H20 improved raw hierarchy and occupancy diagnostics but had approximately
neutral alpha efficiency (median bulk ESS/s about 0.227 versus 0.238 for S1);
expected-active bulk ESS/s changed from about 0.098 to 0.110. A schedules were
both statistically and computationally inefficient. Raw ESS alone therefore
does not justify a new schedule default.

## Scientific summaries and route comparison

Annotation-informed causal ranking was preserved across schedules despite the
joint convergence failures. BED PIP AUPRC ranged 0.592--0.596 and AUROC
0.848--0.853; block-eigen AUPRC ranged 0.541--0.553 and AUROC 0.830--0.836.
At 10 markers every fit recovered 5.85% of causal markers with precision 1.00.
At 25, BED schedules except A20 had recall/precision 0.140/0.96 (A20:
0.135/0.92); at 50 the ranges were 0.240--0.246 and 0.82--0.84. BED recall at
100 was 0.409 with precision 0.70. Block-eigen recall/precision was
0.222--0.228/0.76--0.78 at 50 and 0.374--0.386/0.64--0.66 at 100.

At Bayesian FDR 5%, BED selected 19--20 markers and all were true discoveries;
block eigen selected 23--24 with 22--23 true discoveries. At FDR 10%, BED
selected 23 with 22 true discoveries; block eigen selected 32--33 with 27--28.
Validation genetic-value correlations were 0.9178--0.9183 for BED and
0.8955--0.8958 for block eigen; phenotype-prediction correlations were
0.6545--0.6550 and 0.6316--0.6320. These stable ranking results do not turn a
non-converged hierarchy into a qualified posterior analysis.

No schedule passed, so posterior agreement among converged schedules cannot be
claimed. BED and block eigen continued to occupy different active-count and
heritability regimes. The previously documented BED/block variance-calibration
offset remains separate from annotation-hierarchy mixing and was not modified.

## Decision and next task

The primary decision is **K4: schedule recomposition does not solve the
problem**. Additional H transitions are diagnostically useful and improve
conditional hierarchy equilibration and regime overlap, but even H20 fails the
joint contract. Additional A transitions consistently fail and often worsen
chain-location disagreement. This is not merely a raw-ESS gain: H schedules
change regime overlap, but the improvement is insufficient and is not
consistently favorable per second.

The next package task should review a larger coordinated posterior-preserving
transition, beginning with a derivation and exact small-state validation of a
partially collapsed or block allocation/hierarchy move. Tempering, split--merge,
mode jumping, or interweaving remain later design candidates, not conclusions
of this audit. Package defaults must remain unchanged pending that method
review and a separate qualification-length assessment.

## Validation

The focused BED/block scheduling tests, independent annotation-hierarchy tests,
and component-trace tests passed. The source-tree testthat run completed with
0 failures and 0 errors (one pre-existing covariance warning and one opt-in
reproducibility skip). A clean isolated source install succeeded. Packaged
`R CMD check --no-manual` completed with 0 errors, 0 warnings, and one installed
size NOTE; its installed-package tests reported 4,050 passes, one known warning,
and two skips (the opt-in reproducibility test and a source-tree-only assertion).
The 24-page Quarto documentation render succeeded. `git diff --check` passed.

## Local evidence and limitations

The reproducible harness is `research/sbayesrc/tools/study06_kernel_composition_audit.R`. Raw fit
objects, complete component traces, stick-prior draws, convergence tables,
regime summaries and transitions, allocation-change summaries, representative
marker traces, power metrics, residual drift, and runtime tables are under the
ignored `results/local/study06_kernel_composition_audit/` directory. Raw chains
are intentionally untracked.

This audit uses one informative Study 06 truth, fixed schedules, and the
registered four chains. It does not test combined A/H schedules, another truth,
substantial eigen truncation, or any new transition. It does not rerun formal
qualification and does not authorize the final benchmark.
