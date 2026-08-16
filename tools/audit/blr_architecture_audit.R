root <- normalizePath(if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else ".",
                      winslash = "/", mustWork = TRUE)
old <- setwd(root); on.exit(setwd(old), add = TRUE)
suppressPackageStartupMessages(pkgload::load_all(".", compile = FALSE, quiet = TRUE))

checks <- c()
record <- function(name, value) {
  checks[[name]] <<- isTRUE(value)
  cat(sprintf("%-48s %s\n", name, if (isTRUE(value)) "PASS" else "FAIL"))
}
text <- function(paths) paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")

exports <- getNamespaceExports("sblr")
fitters <- c("stblr_csr", "stblr_csr_annot", "stblr_block_eigen", "stblr_bed",
             "mtblr_bed", "mtblr_csr", "mtblr_block_eigen")
record("seven canonical fitting exports", all(fitters %in% exports))
record("obsolete fitting exports absent", !any(c("sblr", "stblr_bed_marker",
  "check_stblr_convergence", "stblr_csr_bayesr", "stblr_csr_prior_annot",
  "stblr_csr_group_annot", "stblr_csr_learn_annot") %in% exports))
matrix <- sblr:::.blr_model_capability_matrix()
record("six canonical model semantics", setequal(unique(matrix$model),
  c("bayesc", "sbayesc", "bayesr", "sbayesr", "bayesrc", "sbayesrc")))
record("three canonical operators", setequal(unique(matrix$operator),
  c("packed_bed", "csr", "block_eigen")))
record("maf_effect_s independent", !grepl("sbayes.*maf_effect_s = 0",
  text(list.files("R", full.names = TRUE)), ignore.case = TRUE))

phase_tests <- list.files("tests/testthat", "phase[0-9].*[.]R$", ignore.case = TRUE)
phase_audits <- list.files("tools/audit", "phase[0-9].*[.]R$", ignore.case = TRUE)
record("no phase-numbered active tests", !length(phase_tests))
record("no phase-numbered active audits", !length(phase_audits))
record("permanent test owners", all(file.exists(file.path("tests/testthat", c(
  "test-blr-unified-public-contract.R", "test-blr-model-semantics.R",
  "test-blr-operator-reductions.R", "test-blr-unified-convergence.R",
  "test-blr-selected-marker-diagnostics.R")))))

conv <- text(c("R/mtblr-convergence.R", "R/blr-extended-convergence.R"))
record("one rank-normalized R-hat owner", length(gregexpr(
  "[.]blr_convergence_rhat_basic <- function", conv, fixed = TRUE)[[1L]]) == 1L)
record("one ESS owner", length(gregexpr(
  "[.]blr_convergence_ess <- function", conv, fixed = TRUE)[[1L]]) == 1L)
record("one MCSE-mean owner", length(gregexpr(
  "mcse_mean <-", conv, fixed = TRUE)[[1L]]) >= 1L)
record("one diagnostic plan implementation", length(grep(
  "^[.]blr_convergence_controls <- function", readLines("R/blr-unified.R"))) == 1L)

native_files <- setdiff(list.files("src", "[.](cpp|h)$", full.names = TRUE),
                        "src/RcppExports.cpp")
native <- text(native_files)
native_code <- gsub("(?m)^[[:space:]]*//.*$", "", native, perl = TRUE)
record("no unsolicited native stdout", !grepl("std::cout|std::cerr", native))
record("no R console calls in OpenMP source", !grepl(
  "#pragma omp[\\s\\S]{0,500}(Rcpp::Rcout|Rprintf|R_FlushConsole)", native_code,
  perl = TRUE))
record("direct selected-marker capture", grepl("selected_marker", native, fixed = TRUE))
record("no all-marker diagnostic shortcut", grepl(
  "all-marker shortcuts are not supported", text(list.files("R", full.names = TRUE)),
  fixed = TRUE))

workflow <- text(list.files(".github/workflows", full.names = TRUE))
record("CI uses permanent owners", !grepl("phase[0-9]", workflow, ignore.case = TRUE))
record("CI retains package check", grepl("check_package[.]R|R CMD check", workflow))
record("authoritative architecture docs", all(file.exists(file.path("docs/dev", c(
  "blr_architecture.md", "blr_model_contracts.md", "blr_output_schema.md",
  "blr_convergence_contract.md", "blr_backend_inventory.md",
  "blr_test_ownership.md", "blr_development_guide.md")))))

passed <- unlist(checks, use.names = TRUE)
if (!all(passed)) stop("BLR architecture audit failed: ",
  paste(names(passed)[!passed], collapse = ", "), call. = FALSE)
cat(sprintf("architecture_audit=%d/%d PASS\n", sum(passed), length(passed)))
