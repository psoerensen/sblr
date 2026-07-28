root <- normalizePath(if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else ".",
                      winslash = "/", mustWork = TRUE)
old <- setwd(root); on.exit(setwd(old), add = TRUE)
r <- paste(readLines("R/RcppExports.R", warn = FALSE), collapse = "\n")
cpp <- paste(readLines("src/RcppExports.cpp", warn = FALSE), collapse = "\n")
calls <- unique(regmatches(r, gregexpr("_sblr_[A-Za-z0-9_]+", r))[[1L]])
registered <- regmatches(cpp, gregexpr('"_sblr_[A-Za-z0-9_]+"', cpp))[[1L]]
registered <- gsub('"', "", registered, fixed = TRUE)
if (anyDuplicated(registered)) stop("duplicate native registration", call. = FALSE)
if (!setequal(calls, registered)) stop("R/native generated wrapper mismatch", call. = FALSE)
exports <- sub("^export[(]([^)]*)[)]$", "\\1",
  grep("^export[(]", readLines("NAMESPACE", warn = FALSE), value = TRUE))
canonical <- c("stblr_csr", "stblr_csr_annot", "stblr_block_eigen", "stblr_bed",
               "mtblr_csr", "mtblr_block_eigen", "mtblr_bed")
if (!all(canonical %in% exports)) stop("canonical fitter export missing", call. = FALSE)
if (any(c("sblr", "stblr_bed_marker", "check_stblr_convergence") %in% exports))
  stop("obsolete public export", call. = FALSE)
cat(sprintf("generated_interfaces=PASS wrappers=%d registrations=%d exports=%d\n",
            length(calls), length(registered), length(exports)))
