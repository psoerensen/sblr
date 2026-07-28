# BLR convergence contract

Modes are `auto`, `none`, `core`, and `extended`. `auto` computes only the
five core trait quantities when at least two chains are available. Core rows
are `vbs[trait]`, `vgs[trait]`, `ves[trait]`, `vle[trait]`, and `vld[trait]`.

Extended diagnostics add applicable covariance, probability, sampled
`selection_s`, annotation, and explicitly selected marker quantities. Marker
selection must be an explicit unique vector of IDs or one-based indices; there
is no all-marker shortcut. Trace capture is post-burn, every iteration,
unthinned, chain private, RNG neutral, and independent of `keep_chains`,
`keep_traces`, worker count, and posterior thinning.

All scalar traces use the single rank-normalized split/folded R-hat, Geyer ESS,
tail ESS, posterior-SD, MCSE-mean, and relative-MCSE implementation. Diagnostic
memory is resolved and guarded before sampling. Warnings are aggregated and
advisory; diagnostics assess chain mixing rather than prove model correctness.
