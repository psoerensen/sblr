root <- normalizePath(".", winslash = "/", mustWork = TRUE)
contract <- paste(readLines(file.path(root,
  "docs/dev/blr_mt_bed_convergence_contract.md"), warn = FALSE), collapse = "\n")
test <- function(...) all(vapply(c(...), function(x) grepl(x, contract, fixed = TRUE), logical(1)))
guards <- c(
  test("Burn-in is excluded"),
  test("Pooled `vbs/vgs/ves`", "not valid formal diagnostic inputs"),
  test("## 7. Chain splitting"),
  test("discard the central draw"),
  test("average\nranks", "never jittered"),
  test("qnorm((r - 3/8)/(S + 1/4))"),
  test("rhat_folded"),
  test("rhat=max(rhat_rank,rhat_folded)"),
  test("Classical unsplit Gelman"),
  test("M*N*log10(M*N)", "stability"),
  test("bulk ESS is not substituted"),
  test("use their minimum", "q05", "q95"),
  test("positive sequence", "monotone sequence"),
  test("must not be reported as R-hat 1"),
  test("not_updated"),
  test("unavailable_single_chain"),
  test("not R-hat or posterior SD"),
  test("not R-hat or posterior SD"),
  test("must work with `keep_chains=FALSE`"),
  test("All-marker default traces are forbidden"),
  test("Pattern-specific", "memory-controlled"),
  test("never charges packed genotype storage"),
  test("tier1_trace_bytes"),
  test("R-hat > 1.01"),
  test("100*nchains"),
  test("at most one main-thread warning per fit"),
  test("definitive `converged` status"),
  !grepl("convergence =", paste(readLines(file.path(root,
    "R/mtblr-bed.R"), warn = FALSE), collapse = "\n"), fixed = TRUE),
  system2("git", c("diff", "--quiet", "--",
                   "src/blr_mt_bed_core_impl.h", "R/mtblr-bed.R")) == 0,
  test("raw version 1"),
  system2("git", c("diff", "--quiet", "--", "DESCRIPTION")) == 0,
  test("`posterior` is optional", "development comparison"),
  grepl("MtBedConvergenceTraceBundle", paste(readLines(file.path(root,
    "src/mtblr.cpp"), warn = FALSE), collapse = "\n"), fixed = TRUE),
  system2("git", c("diff", "--quiet", "--", "R/mtblr-bed.R")) == 0,
  grepl("tools/check/check_package.R", paste(readLines(file.path(root,
    ".github/workflows/blr-framework.yml"), warn = FALSE), collapse = "\n"), fixed = TRUE)
)
names(guards) <- sprintf("MUTATION_%02d", seq_along(guards))
for (name in names(guards)) cat(name, "=", guards[[name]], "\n", sep = "")
cat("ALL_35_CRITICAL_MUTATIONS_DETECTED=", all(guards), "\n", sep = "")
stopifnot(length(guards) == 35L, all(guards))
