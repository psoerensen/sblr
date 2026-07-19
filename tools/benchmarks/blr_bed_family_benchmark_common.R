bed_family_benchmark_metadata<-function(model,samples,markers,iterations,burnin,
 thinning,chains,cores,times,completed_fit_rss_mb,bed_size_bytes,controls){
 data.frame(model,samples,markers,iterations,burnin,thinning,chains,cores,
  times=paste(times,collapse=","),mean=mean(times),median=median(times),
  minimum=min(times),maximum=max(times),range=diff(range(times)),
  completed_fit_rss_mb,bed_size_bytes,controls,
  r_version=R.version.string,stringsAsFactors=FALSE)
}
bed_family_benchmark_notice<-function(){
 cat("no cross-model speed ranking\nworkloads are model-specific\n")
 cat("completed-fit RSS is not peak RSS\npage-cache effects apply\n")
 cat("tiny timings are regression signals\n")
}
