read_text <- function(path) paste(readLines(path, warn=FALSE), collapse="\n")
core <- read_text("src/blr_mt_default_core_impl.h")
access <- read_text("src/blr_mt_ld_access.h")
adapter <- read_text("src/mtblr.cpp")
legacy <- read_text("src/blr_mt_default_legacy_adapter.h")
rroute <- read_text("R/interface_mtblr.R")

checks <- list(
  trait_zero = !grepl("trait_ld\\[0\\].*apply_offdiag", access),
  diagonal_omitted = grepl("residual\\[static_cast<std::size_t>\\(marker\\)\\] -=", access),
  duplicate_gibbs = lengths(regmatches(core, gregexpr("for \\( int it =", core))) == 1L,
  unequal_shared = grepl("shared LD requires identical trait diagonals", adapter, fixed=TRUE),
  traversal_order = grepl("position = start; position < end; ++position",
    read_text("src/blr_sparse_ld_csr.h"), fixed=TRUE),
  public_csr = !grepl("mtblr_csr_internal", rroute, fixed=TRUE),
  duplicate_adapter = lengths(regmatches(legacy,
    gregexpr("MtDefaultLegacyResult result\\(20\\)", legacy))) == 1L
)
for (name in names(checks)) {
  cat(sprintf("%-24s detected=%s\n", name, checks[[name]]))
}
stopifnot(all(unlist(checks)))
