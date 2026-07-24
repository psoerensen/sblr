root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pkgload::load_all(root, compile = FALSE, quiet = TRUE)
source(file.path(root, "tests/testthat/helper-mtblr-block-eigen-fixtures.R"),
       local = TRUE)
case <- phase17m_public_case()
rejects <- function(...) inherits(try(phase17m_call(case, ...), silent = TRUE),
                                  "try-error")

external <- case$stats
external$source <- "external"
rows <- case$stats
rows$rows <- rev(rows$rows)
columns <- case$stats
columns$cls[[1L]] <- c(2L, 1L, 3L, 4L)
frequencies <- case$stats
frequencies$af[[1L]][1L] <- .2
source_r <- paste(readLines(file.path(root, "R/mtblr-block-eigen.R"),
                            warn = FALSE), collapse = "\n")
source_cpp <- paste(readLines(file.path(root, "src/mtblr.cpp"),
                              warn = FALSE), collapse = "\n")
namespace <- paste(readLines(file.path(root, "NAMESPACE"), warn = FALSE),
                   collapse = "\n")

checks <- c(
  EXTERNAL_STATS = rejects(stats = external),
  BED_ROWS = rejects(stats = rows),
  MARKER_COLUMNS = rejects(stats = columns),
  ALLELE_FREQUENCIES = rejects(stats = frequencies),
  TRANSFORMED_WY = grepl("result.transformed_wy", source_cpp, fixed = TRUE),
  WW_NOT_RUNTIME_DIAGONAL =
    !grepl("mtblr_block_eigen_raw_internal(\n    st$wy, st$ww", source_r,
           fixed = TRUE),
  PUBLIC_ONE_BASED_BLOCKS =
    rejects(block_start = c(0L, 3L)),
  SHARED_INCOMPATIBILITY =
    rejects(operator_sharing = "shared",
            eigen_filter = c("hard_truncate", "ridge_fixed")),
  ONE_RAW_CONVERTER =
    length(gregexpr("Rcpp::List mtblr_legacy_to_raw(", source_cpp,
                    fixed = TRUE)[[1L]]) == 1L,
  CSR_RAW_USES_CONVERTER =
    grepl('legacy, models, "mt_csr_bayesc"', source_cpp, fixed = TRUE),
  ACTUAL_BUILD_DIAGNOSTICS =
    grepl("&owner_diagnostics[owner]", source_cpp, fixed = TRUE),
  SINGLE_ADAPTER_EXECUTION =
    length(gregexpr("MtBlockEigenAdapterResult adapter=run_mt_block_eigen_adapter(",
                    source_cpp, fixed = TRUE)[[1L]]) == 1L,
  NO_RESEARCH_ROUTE = !grepl("mtblr_eigen(", source_r, fixed = TRUE),
  PUBLIC_EXPORT = grepl("export(mtblr_block_eigen)", namespace, fixed = TRUE),
  SAMPLE_OVERLAP = rejects(sample_overlap = "modeled"),
  DIAGONAL_RESIDUAL =
    rejects(ve = matrix(c(1, .1, .1, 1), 2L)),
  ONE_GENERAL_FORMATTER =
    length(gregexpr(".as_mtblr_fit <- function(",
      paste(readLines(file.path(root, "R/mtblr-csr.R"), warn = FALSE),
            collapse = "\n"), fixed = TRUE)[[1L]]) == 1L)

for (name in names(checks)) cat(name, "=", checks[[name]], "\n", sep = "")
if (!all(checks)) stop("Undetected Phase 17M critical mutation: ",
                       paste(names(checks)[!checks], collapse = ", "))
cat("ALL_PHASE17M_CRITICAL_MUTATIONS_DETECTED=TRUE\n")
