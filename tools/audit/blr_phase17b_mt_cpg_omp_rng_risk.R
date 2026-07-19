root <- normalizePath(if (file.exists("DESCRIPTION")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE)
source(file.path(root, "tests/testthat/fixtures/blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))

run_risk_case <- function(root, threads) {
  setwd(root)
  pkgload::load_all(".", compile = FALSE, quiet = TRUE)
  source("tests/testthat/fixtures/blr_phase17b_mt_default/blr-phase17b-mt-default-reference.R")
  x <- phase17b_mt_native_arguments(phase17b_mt_config(1L))
  m <- length(x$wy[[1]])
  x$XXindices <- lapply(seq_len(m), function(i) i - 1L)
  x$XXvalues <- lapply(x$wy, function(unused)
    lapply(seq_len(m), function(i) x$ww[[1]][i]))
  names(x)[names(x) == "b"] <- "b_init"
  do.call(sblr:::mtblr_cpg_omp, x)
}

if (!requireNamespace("callr", quietly = TRUE)) stop("callr is required")
one <- callr::r(run_risk_case, list(root, 1L), env = c(OMP_NUM_THREADS = "1"))
two <- callr::r(run_risk_case, list(root, 2L), env = c(OMP_NUM_THREADS = "2"))
source(file.path(root, "tests/testthat/helper-source-architecture.R"))
difference <- reference_first_difference(one, two)
print(difference)
if (is.null(difference)) stop("Expected worker-sensitive difference was not observed")
