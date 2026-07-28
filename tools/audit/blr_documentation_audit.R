root <- normalizePath(if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else ".",
                      winslash = "/", mustWork = TRUE)
old <- setwd(root); on.exit(setwd(old), add = TRUE)
docs <- c(list.files("docs/dev", "[.]md$", recursive = TRUE, full.names = TRUE),
          list.files("docs/notes", "[.](md|qmd)$", recursive = TRUE,
                     full.names = TRUE))
active <- docs[!grepl("docs/dev/history", docs, fixed = TRUE) &
               !grepl("blr_cleanup_manifest", docs, fixed = TRUE)]
txt <- paste(unlist(lapply(active, readLines, warn = FALSE)), collapse = "\n")
obsolete <- c("stblr_bed_marker", "check_stblr_convergence",
              "stblr_csr_bayesr", "blr_framework_phase")
bad <- obsolete[vapply(obsolete, grepl, logical(1), x = txt, fixed = TRUE)]
if (length(bad)) stop("obsolete documentation terms: ", paste(bad, collapse = ", "),
                      call. = FALSE)
required <- c("bayesc", "sbayesc", "bayesr", "sbayesr", "bayesrc",
              "sbayesrc", "selection_s", "extended")
if (!all(vapply(required, grepl, logical(1), x = txt, fixed = TRUE)))
  stop("canonical documentation terminology incomplete", call. = FALSE)
link_records <- unlist(lapply(docs, function(file) {
  x <- paste(readLines(file, warn = FALSE), collapse = "\n")
  hit <- regmatches(x, gregexpr("(?<=\\]\\()[^)]+(?=\\))", x,
                                perl = TRUE))[[1L]]
  hit <- hit[hit != "" & !grepl("^(https?:|#|mailto:)", hit)]
  if (!length(hit)) return(character())
  setNames(file.path(dirname(file), hit), rep(file, length(hit)))
}))
missing <- link_records[!file.exists(link_records)]
if (length(missing)) stop("broken documentation links: ",
  paste(paste0(names(missing), " -> ", missing), collapse = ", "), call. = FALSE)
cat(sprintf("documentation_audit=PASS active_documents=%d candidate_links=%d\n",
            length(active), length(link_records)))
