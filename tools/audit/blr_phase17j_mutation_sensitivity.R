root <- normalizePath(".", winslash = "/")
r <- paste(readLines(file.path(root, "R", "mtblr-csr.R"), warn = FALSE), collapse = "\n")
cpp <- paste(readLines(file.path(root, "src", "mtblr.cpp"), warn = FALSE), collapse = "\n")
core <- paste(readLines(file.path(root, "src", "blr_mt_default_core_impl.h"), warn = FALSE), collapse = "\n")
checks <- c(
  marker_id_validation = grepl("Marker IDs must be unique", r, fixed = TRUE),
  ld_order_validation = grepl("LD resources must already use identical marker IDs and order", r, fixed = TRUE),
  allele_swap = grepl("Swapped effect/other alleles", r, fixed = TRUE),
  unknown_scale = grepl("standardized_genotype", r, fixed = TRUE),
  offdiag_yy = grepl("Nonzero off-diagonal yy", r, fixed = TRUE),
  overlap_policy = grepl("sample_overlap must be exactly", r, fixed = TRUE),
  named_raw_validation = grepl(".validate_mtblr_raw(raw)", r, fixed = TRUE),
  no_research_route = !grepl("mtblr_eigen|mtblr_cpg_omp_csr", r),
  one_gibbs_loop = lengths(regmatches(core, gregexpr("for ( int it =", core, fixed=TRUE))) == 1L,
  model_names_preserved = grepl("apply(models, 1L, paste", r, fixed = TRUE),
  one_raw_wrapper = lengths(regmatches(cpp, gregexpr("mtblr_csr_raw_internal(", cpp, fixed=TRUE))) == 1L
)
print(checks)
if (!all(checks)) stop("Phase 17J mutation-sensitivity invariant was not detected.")
