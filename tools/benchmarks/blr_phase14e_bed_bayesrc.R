pkgload::load_all(".",compile=FALSE,quiet=TRUE)
source(file.path("tests","testthat","fixtures","blr-phase14a-bed-bayesrc-reference.R"))
rss_mb<-function() if(requireNamespace("ps",quietly=TRUE)) as.numeric(ps::ps_memory_info()[["rss"]])/1024^2 else NA_real_
bench<-function(label,cores,chains,update=TRUE,multi=TRUE,reps=5L){
 phase14a_capture(cores,chains,1500L,update,multi)
 times<-vapply(seq_len(reps),function(i) system.time(phase14a_capture(cores,chains,1500L+i,update,multi))[["elapsed"]],numeric(1))
 data.frame(label,samples=6L,markers=2L,annotations=if(multi)3L else 1L,
  annotation_density=1,components=4L,scales="0,.01,.1,1",iterations=8L,burnin=2L,
  thinning=1L,chains,cores,updateAlpha=update,times=paste(times,collapse=","),
  mean=mean(times),median=median(times),minimum=min(times),maximum=max(times),
  range=diff(range(times)),completed_fit_rss_mb=rss_mb(),bed_size_bytes=NA_real_,
  annotation_size_bytes=as.numeric(object.size(if(multi)matrix(0,2,3) else matrix(1,2,1))),
  memory_method="completed-fit RSS; not peak")
}
results<-do.call(rbind,list(bench("intercept_fixed_1x1",1,1,FALSE,FALSE),
 bench("multi_annotation_1x1",1,1,TRUE,TRUE),bench("multi_annotation_2x1",1,2,TRUE,TRUE),
 bench("multi_annotation_2x2",2,2,TRUE,TRUE)))
cat("Canonical Phase 14E packed-BED BayesRC baseline\n")
cat("Warm-up precedes five repetitions; page cache affects repeated BED reads.\n")
cat("completed-fit RSS is not peak RSS\n")
cat("sampled peak RSS is opt-in unless successfully measured\n")
cat("tiny timings are regression signals, not performance claims\n")
cat("Moderate 2,000-marker/200-sample and larger workloads remain opt-in.\n")
cat("R:",R.version.string,"\n")
cat("sblr:",as.character(utils::packageVersion("sblr"))," Rcpp:",as.character(utils::packageVersion("Rcpp")),"\n")
print(results,row.names=FALSE)
