pkgload::load_all(".", compile=FALSE, quiet=TRUE)
test_lines <- readLines("tests/testthat/test-blr-framework-phase17i.R", warn=FALSE)
stop_line <- grep("^test_that", test_lines)[1L] - 1L
eval(parse(text=test_lines[seq_len(stop_line)]), envir=.GlobalEnv)

bench <- function(label, case, repetitions=5L) {
  invisible(do.call(sblr:::mtblr, case$dense))
  invisible(do.call(sblr:::mtblr_csr_internal, case$csr))
  dense <- replicate(repetitions,
    system.time(do.call(sblr:::mtblr, case$dense))[["elapsed"]])
  csr <- replicate(repetitions,
    system.time(do.call(sblr:::mtblr_csr_internal, case$csr))[["elapsed"]])
  data.frame(label=label, traits=length(case$dense$wy),
    markers=length(case$dense$wy[[1]]), nit=case$dense$nit,
    nburn=case$dense$nburn, nthin=case$dense$nthin,
    dense_mean=mean(dense), dense_median=median(dense), dense_min=min(dense),
    dense_max=max(dense), csr_mean=mean(csr), csr_median=median(csr),
    csr_min=min(csr), csr_max=max(csr), completed_fit_rss=NA_real_)
}
out <- rbind(
  bench("shared_all_updates", phase17i_case(updates=TRUE)),
  bench("trait_specific", phase17i_case(trait_specific=TRUE, updates=TRUE)),
  bench("independent_patterns", phase17i_case(trait_specific=TRUE,
    independent=TRUE, updates=FALSE))
)
print(out, row.names=FALSE)
cat("Tiny timings are regression signals; completed-fit RSS is not peak RSS.\n")
cat("No MCMC-time I/O occurs. Moderate scaling is reported by the memory audit.\n")
