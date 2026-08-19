args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else
  file.path("research", "sbayesrv", "study02_coordinate", "analysis.R")
script <- normalizePath(script, winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script), "..", "..", ".."),
                           winslash = "/", mustWork = TRUE)
setwd(repo_root)

if (!requireNamespace("pkgload", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("qgg", quietly = TRUE)) {
  stop("This research run requires the existing pkgload, digest, and qgg packages.",
       call. = FALSE)
}
bench_root <- normalizePath(file.path(dirname(repo_root), "sblrbench"),
                            winslash = "/", mustWork = TRUE)
suppressMessages(pkgload::load_all(repo_root, compile = FALSE, quiet = TRUE))
suppressMessages(pkgload::load_all(bench_root, compile = FALSE, quiet = TRUE))
source(file.path(repo_root, "research", "sbayesrv", "prototype.R"),
       local = .GlobalEnv)
source(file.path(repo_root, "research", "sbayesrv", "study02_coordinate",
                 "helpers.R"), local = .GlobalEnv)

status <- tryCatch(
  run_study02_coordinate(repo_root),
  study02_coordinate_error = function(e) {
    output <- file.path(repo_root, "research", "sbayesrv",
                        "study02_coordinate", "output")
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    completed_signal <- identical(e$label, "BLOCKED_Q500_H050") &&
      file.exists(file.path(output, "metrics.csv")) &&
      all(c("q50_h030", "q500_h030") %in%
            utils::read.csv(file.path(output, "metrics.csv"),
                            stringsAsFactors = FALSE)$stage)
    checkpoint_status <- if (completed_signal) {
      "SBayesRV_PREDICTION_SIGNAL_OBSERVED_AT_STUDY02_COORDINATE"
    } else {
      e$label
    }
    limitation <- if (completed_signal) {
      "ORDINARY_SBayesR_SPARSE_LD_BOUNDARY_AT_Q500_H050"
    } else {
      NA_character_
    }
    utils::write.csv(data.frame(
      status = checkpoint_status,
      limitation = limitation,
      execution_error = conditionMessage(e),
      stringsAsFactors = FALSE),
      file.path(output, "status.csv"), row.names = FALSE)
    writeLines(c(
      "# Study 02-coordinate SBayesRV research run",
      "",
      paste0("Status: `", checkpoint_status, "`"),
      if (!is.na(limitation)) paste0("Limitation: `", limitation, "`") else "",
      "",
      conditionMessage(e),
      "",
      "Completed stages remain in compact output tables. A failed stage has no",
      "held-out evaluation. LD was reused from the exact cache, not rebuilt.",
      "This one-replicate diagnostic is not benchmark evidence."),
      file.path(output, "summary.md"))
    message(conditionMessage(e))
    checkpoint_status
  })
cat(status, "\n")
