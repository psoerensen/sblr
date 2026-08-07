# Conditional particle-Gibbs block transition

## Target sequence

For a retained block factor `Q`, transformed score `w`, fixed block residual
variance `Ve`, fixed annotation probabilities, and effects outside the block,
the block conditional is proportional to

```text
exp{-||w - Q beta||^2 / (2 Ve)}
prod_j pi[j, c[j]] p(beta[j] | c[j], vb, gamma).
```

The development kernel uses the valid prefix sequence obtained by setting all
not-yet-visited marker effects to zero. Adding marker `j` changes the log
likelihood by

```text
(beta[j] q[j]' r - 0.5 beta[j]^2 q[j]'q[j]) / Ve.
```

The component/effect proposal is the exactly normalized version of this
incremental factor. Its incremental importance weight is therefore only that
normalizing constant. Conditional SMC retains the current complete block path
as particle one, resamples only the other ancestors, and selects one terminal
path by normalized final weight. This is the standard particle-Gibbs extended
target, so its marginal transition preserves the existing block conditional.

No block likelihood, prior, alpha state, or residual-variance rule changes.
The first implementation is an R development reference, not a production
sampler option. It deliberately omits ancestor sampling until ordinary
conditional SMC path diversity has been measured.

## Gate

The reference must match exact allocation enumeration on tiny blocks. It is
then screened at 100 and 500 markers for particle ESS, ancestry collapse,
changed allocations, active-count jumps, and runtime. Production integration
is permitted only if fixed-alpha block paths remain diverse at useful particle
counts and costs.
