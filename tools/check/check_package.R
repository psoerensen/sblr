args <- commandArgs(trailingOnly = TRUE)
source_root <- normalizePath(if (length(args)) args[[1L]] else ".",
                             winslash = "/", mustWork = TRUE)
check_root <- tempfile("sblr-check-")
dir.create(check_root, recursive = TRUE)
on.exit(unlink(check_root, recursive = TRUE, force = TRUE), add = TRUE)
r <- file.path(R.home("bin"), "R")
run <- function(arguments, label) {
  status <- system2(r, arguments, stdout = "", stderr = "")
  if (!identical(status, 0L)) stop(label, " failed with status ", status,
                                    call. = FALSE)
}
old <- setwd(check_root); on.exit(setwd(old), add = TRUE)
message("R CMD build ", source_root)
run(c("CMD", "build", shQuote(source_root), "--no-build-vignettes"),
    "R CMD build")
tarball <- list.files(check_root, pattern = "[.]tar[.]gz$", full.names = TRUE)
if (length(tarball) != 1L) stop("Expected exactly one source tarball.")
message("R CMD check ", tarball)
run(c("CMD", "check", "--no-manual", "--no-build-vignettes",
      shQuote(tarball)), "R CMD check")
log <- readLines(file.path(check_root, "sblr.Rcheck", "00check.log"), warn = FALSE)
status <- grep("^Status:", log, value = TRUE)
cat(paste(log, collapse = "\n"), "\n")
if (!length(status) || grepl("ERROR|WARNING", status)) {
  stop("Package check reported ERROR or WARNING: ", paste(status, collapse = "; "))
}
