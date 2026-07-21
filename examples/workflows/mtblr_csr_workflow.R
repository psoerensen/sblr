## Example 1: one BED/Glist construction and a multi-trait stats object.
## stats <- make_summary_stats(Glist, y = cbind(trait1, trait2))
## Glist <- make_sparse_ld(Glist, out_prefix = "shared_ld")
## fit <- mtblr_csr(stats, Glist = Glist, nit = 20, nburn = 10, seed = 1)

## Example 2: one raw LD correlation reference with trait-specific ww.
## stats_list <- list(study1 = stats1, study2 = stats2)
## fit <- mtblr_csr(stats_list, ld_prefix = "shared_ld",
##                  ld_metadata = shared_descriptor,
##                  nit = 20, nburn = 10, seed = 1)

## Example 3: trait/study-specific LD resources and provenance.
## fit <- mtblr_csr(
##   list(EUR = eur_stats, AFR = afr_stats),
##   Glist = list(eur_glist, afr_glist),
##   trait_metadata = data.frame(
##     trait_id = c("EUR", "AFR"), study_id = c("study_eur", "study_afr"),
##     ancestry = c("European", "African"), population = c("EUR", "AFR"),
##     ld_reference = c("eur_panel", "afr_panel")),
##   sample_overlap = "not_modeled", nit = 20, nburn = 10, seed = 1)
