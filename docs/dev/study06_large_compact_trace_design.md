# Study 06 large compact allocation traces

## Contract

Extended convergence control now accepts `aggregate_component_states = TRUE`
for BayesR/BayesRC kernels. The option retains integer histories for component
counts, realized active count, and the eligible/continue/stop counts of every
sequential stick. It is implemented for packed BED and the shared scalar
summary operator owner used by CSR and retained block eigen. The multivariate
owners use the same count semantics.

For component labels `0, ..., K - 1`, retained state `d`, and stick index
`j = 0, ..., K - 2`:

```text
component_count[k] = sum(d == k)
realized_active_count = sum(d > 0)
eligible[j] = sum(d >= j)
continue[j] = sum(d > j)
stop[j] = eligible[j] - continue[j]
```

Thus the first stick is eligible for every fitted marker, its continuation
count is the realized active count, and its stop count is the null-component
count. These are sampled occupancies, not expected counts from marker prior
probabilities.

## Public result

When retained traces are requested, the fit exposes:

```text
component_count_trace
realized_active_count_trace
stick_eligible_count_trace
stick_continue_count_trace
stick_stop_count_trace
```

Each is an integer array with dimensions `draw x chain x quantity`. Quantity
metadata remain attached. Scalar fits require `keep_chains = TRUE` when this
option is enabled so chain boundaries cannot be silently lost.

The native capture reads the current in-memory component vector after each
retained state. It neither draws random numbers nor changes marker update
order, residual state, allocation probabilities, or sampler transitions.

## Storage

For Study 06 large feasibility (`m = 37,991`, four chains, 12,000 package
retained draws, four components), a full double component array would require
14,588,544,000 bytes per fit (13.59 GiB; 14.59 GB). The compact integer arrays
contain 14 integers per draw/chain and require 2,688,000 bytes (2.56 MiB) per
fit, a 5,427-fold reduction. Six fits require about 15.38 MiB for these arrays,
excluding ordinary scalar and the unchanged 300-marker selected histories.

## Validation

The focused tests compare every compact count with an independent R tabulation
of full selected-marker component histories at every draw and chain. They
cover four-component BayesR and BayesRC, packed BED, retained block eigen,
multiple chains, burn-in, simultaneous selected-marker capture, and explicit
chain-boundary enforcement. Paired tracing-disabled/enabled fits require exact
equality (`tolerance = 0`) for existing RNG-sensitive scientific outputs.

