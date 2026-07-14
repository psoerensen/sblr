pkgload::load_all(".", compile=FALSE)
source("tests/testthat/fixtures/blr-phase5a-bayesr-reference.R")
configs<-phase5a_bayesr_configs[c("one_chain","two_chains_one_core","two_chains_two_cores")]
cat("Phase 5B CSR BayesR post-migration baseline\n",R.version.string,"\n")
for(nm in names(configs)) {
  cfg<-configs[[nm]]; invisible(phase5a_bayesr_run(cfg,FALSE))
  elapsed<-replicate(3,system.time(phase5a_bayesr_run(cfg,FALSE))[["elapsed"]])
  cat(nm,"markers=4 traits=",cfg$traits,"components=",length(cfg$mixture),
      "chains=",cfg$nchains,"cores=",cfg$ncores,"times=",paste(elapsed,collapse=","),
      "mean=",mean(elapsed),"median=",median(elapsed),"range=",paste(range(elapsed),collapse=":"),"\n")
}
cat("RSS: whole-process sampling not available in this compact in-process baseline; exact fixtures and repeated timings are primary.\n")
