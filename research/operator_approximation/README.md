# Operator-approximation research handoff

Hard entrywise CSR sparsification can underestimate the accumulated effect-
state quadratic $b^\mathsf{T}Cb$ even when every omitted correlation is weak.
The SBayesRV Study 02-coordinate `q500_h050` attempt exposed this boundary in
ordinary SBayesR with 1,006 active markers:

- sparse $b^\mathsf{T}Cb$: 9072.37;
- minimum quadratic required for nonnegative SSE: 11674.59;
- missing quadratic: at least 2602.22.

The failure belongs to the approximation of the LD operator, not to gsim,
SBayesRV, or the SBayesRV theta transition. Future operator work may compare a
wider or lower-threshold CSR representation, retained block eigen, or an exact
matrix-free residual correction. Those alternatives require their own
prespecified qualification; none is implemented or tested here.
