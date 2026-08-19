# SBayesRV research boundary

For marker $j$ assigned to a non-null BayesR component $k$, SBayesRV studies a
relative variance modifier $q_j$ derived from annotations or other signals:

$$
b_j \mid c_j=k>0, q_j \sim N(0, v_b \gamma_k q_j).
$$

The null component remains a point mass at zero. The research question is how
to specify, identify, and validate $q_j$ without conflating marker-level
variance modulation with SBayesRC's annotation-dependent mixture membership.

This document sets naming and scope only. It specifies no public interface,
sampler transition, implementation plan, or benchmark claim. Historical
Study 12 used the name SBayesR-LV; that wording is retained inside the archived
files only for traceability.

