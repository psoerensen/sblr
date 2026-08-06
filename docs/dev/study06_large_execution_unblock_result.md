# Study 06 large execution-unblock result

## Status

The compact-retention blocker is removed and independently validated. The B0
residual-scale blocker remains a mathematical operator-contract issue, so the
frozen six-fit experiment remains blocked before scientific execution.

No sampler transition, probability, prior, initialization, RNG order, marker
order, eigen policy, or posterior quantity was changed. In particular, no
clipping and no 0.995 fallback was introduced.

## Compact histories

The public extended-diagnostic API now returns component count, realized active
count, and stick eligible/continue/stop histories using about 2.56 MiB per
large fit rather than a 13.59 GiB full component array. Exact native-versus-R
oracles and tracing on/off comparisons pass for the required BED BayesR,
BayesRC, block-eigen SBayesR, and shared SBayesRC owner.

## B0

The exact frozen state was reproduced. The block-factor calculation agrees
with an independent oracle, but omitting cross-block quadratic contributions
changes the valid direct SSE of 13,105.26 to an invalid block SSE of -4,223.39.
This is classified B0-C4, not a numerical-tolerance case.

## Consequence

Only non-inferential diagnostic states were run. No registered full fit was
launched. The next package task must explicitly choose and validate the
block-eigen residual-variance likelihood contract; only then can the frozen
experiment be resumed without disguising a model change as an implementation
fix.
