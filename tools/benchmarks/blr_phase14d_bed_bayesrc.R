pkgload::load_all(".",compile=FALSE,quiet=TRUE)
source(file.path("tests","testthat","fixtures","blr-phase14a-bed-bayesrc-reference.R"))
rss_mb<-function() if(requireNamespace("ps",quietly=TRUE)) as.numeric(ps::ps_memory_info()[["rss"]])/1024^2 else NA_real_
bench<-function(label,cores,chains,update=TRUE,multi=TRUE,reps=5L){
 phase14a_capture(cores,chains,1400L,update,multi)
 times<-vapply(seq_len(reps),function(i) system.time(phase14a_capture(cores,chains,1400L+i,update,multi))[["elapsed"]],numeric(1))
 data.frame(label,samples=6L,markers=2L,annotations=if(multi)3L else 1L,components=4L,
  scales="0,.01,.1,1",iterations=8L,burnin=2L,thinning=1L,chains,cores,
  updateAlpha=update,times=paste(times,collapse=","),mean=mean(times),median=median(times),
  minimum=min(times),maximum=max(times),range=diff(range(times)),completed_fit_rss_mb=rss_mb(),
  memory_method="completed-fit RSS; not sampled peak")
}
results<-do.call(rbind,list(bench("intercept_fixed_1x1",1,1,FALSE,FALSE),
 bench("multi_annotation_1x1",1,1,TRUE,TRUE),bench("multi_annotation_2x1",1,2,TRUE,TRUE),
 bench("multi_annotation_2x2",2,2,TRUE,TRUE)))
cat("Phase 14D packed-BED BayesRC migrated baseline\n")
cat("Warm-up precedes five repetitions; page cache affects repeated BED reads.\n")
cat("Completed-fit RSS is not sampled peak RSS; moderate/large workloads remain opt-in.\n")
cat("R:",R.version.string,"\n")
print(results,row.names=FALSE)
