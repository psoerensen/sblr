pkgload::load_all(".",compile=FALSE,quiet=TRUE)
library(testthat)
source(file.path("tests","testthat","fixtures","blr-phase14a-bed-bayesrc-reference.R"))
rss_mb <- function() if(requireNamespace("ps",quietly=TRUE)) as.numeric(ps::ps_memory_info()[["rss"]])/1024^2 else NA_real_
bench <- function(label,cores,chains,update=TRUE,multi=TRUE,reps=5L){
 phase14a_capture(cores,chains,1400L,update,multi)
 times<-vapply(seq_len(reps),function(i) system.time(phase14a_capture(cores,chains,1400L+i,update,multi))[["elapsed"]],numeric(1))
 data.frame(label,samples=6L,markers=2L,annotations=if(multi)3L else 1L,components=4L,
  scales="0,.01,.1,1",iterations=8L,burnin=2L,thinning=1L,chains,cores,
  updateAlpha=update,times=paste(times,collapse=","),mean=mean(times),median=median(times),
  minimum=min(times),maximum=max(times),range=diff(range(times)),completed_fit_rss_mb=rss_mb(),
  memory_method="completed-fit RSS; not peak")
}
results<-do.call(rbind,list(bench("intercept_fixed_1x1",1,1,FALSE,FALSE),
 bench("multi_annotation_1x1",1,1,TRUE,TRUE),bench("multi_annotation_2x1",1,2,TRUE,TRUE),
 bench("multi_annotation_2x2",2,2,TRUE,TRUE)))
cat("Phase 14A packed-BED BayesRC audit baseline\n")
cat("Warm-up precedes timing; page cache affects repeated reads.\n")
cat("Tiny workload only; moderate/large and sampled peak RSS are opt-in.\n")
cat("R:",R.version.string,"\n")
print(results,row.names=FALSE)
