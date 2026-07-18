# Phase 13D migration-closure baseline. Completed-fit RSS is not peak RSS.
pkgload::load_all(".",compile=FALSE,quiet=TRUE)
source(file.path("tests","testthat","fixtures","blr-phase13a-bed-bayesr-reference.R"))
rss_mb <- function() if (requireNamespace("ps",quietly=TRUE))
 as.numeric(ps::ps_memory_info()[["rss"]])/1024^2 else NA_real_
bench <- function(label,cores,chains,full=10L,base=50L,progress=0,reps=5L) {
 phase13a_capture(cores,chains,1300L,full,base,progress_every=progress)
 times <- vapply(seq_len(reps),function(i) system.time(phase13a_capture(
  cores,chains,1300L+i,full,base,progress_every=progress))[["elapsed"]],numeric(1))
 data.frame(label,samples=6L,markers=2L,components=4L,scales="0,.01,.1,1",
  iterations=8L,burnin=2L,thinning=1L,chains,cores,full_sweep_every=full,
  null_skip_base=base,progress_every=progress,times=paste(times,collapse=","),
  mean=mean(times),median=median(times),minimum=min(times),maximum=max(times),
  range=diff(range(times)),completed_fit_rss_mb=rss_mb(),
  memory_method="completed-fit RSS; not peak")
}
results <- do.call(rbind,list(bench("dense_1x1",1L,1L,0L,1L),
 bench("dense_2x1",1L,2L,0L,1L),bench("dense_2x2",2L,2L,0L,1L),
 bench("aggressive_skip",2L,2L,25L,100L),bench("conservative_skip",2L,2L,5L,2L),
 bench("progress_enabled",1L,1L,10L,50L,1)))
cat("Phase 13D packed-BED BayesR migration-closure baseline\n")
cat("Tiny BED size:",file.info(phase11a_fixture()$Glist$bedfiles)$size,"bytes\n")
cat("First call is warm-up; OS page cache affects later reads.\n")
cat("Moderate 2,000-marker/200-sample and peak-RSS workloads remain opt-in resource runs.\n")
cat("R:",R.version.string,"\n"); print(results,row.names=FALSE)
