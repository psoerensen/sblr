# Study 12 — log-variance annotation evaluation

**Status: `paused_pending_redesign`.** Phase 8B preserves scientific
traceability but does not maintain Study 12 as an executable benchmark.

Study 12 remains separate from Study 06. Aspect 1 compares SBayesR,
SBayesR-LV, and SBayesRC on the immutable
`canonical_annotation_simulation_v1`; Aspect 2 is the reciprocal SBayesR-LV
truth comparison on the same genomic backbone.

The scientific contracts, design documents, reports, historical seeds,
thresholds, decisions, and local evidence identity are preserved. Evidence
remains under ignored `results/local/12_logvar_annotation_evaluation/` and is
not promoted. [`analysis.R`](analysis.R) validates status and traceability only
and writes nothing. Fit, truth-generation, and summarization gates fail
immediately because rerun orchestration is deliberately deferred to the future
redesign.

The frozen backbone retains its truthful historical source identity: the
protected legacy `human_independent` qgdata fixture at
`6cca5819e711d326cfb2614d7e9d9f34942612cd`. It is not relabelled as the
current `human_1000g_eur` panel or the current qgdata HEAD. The historical
executable implementation can be recovered from Git checkpoint
`7a5cbf947af9b3441538fc64d7966c693e1433f8` if a later redesign needs it; Git
history, not a dormant compatibility workflow, is the archive.
