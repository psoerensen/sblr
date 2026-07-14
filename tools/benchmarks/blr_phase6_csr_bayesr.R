pkgload::load_all(".",compile=FALSE)
source("tests/testthat/fixtures/blr-phase5a-bayesr-reference.R")
make_case<-function(m=2000L,traits=2L){
 p<-tempfile("phase6_bayesr_"); sblr:::.stblr_write_uint64_file(paste0(p,".row_ptr.u64.bin"),rep(0,m+1L));file.create(paste0(p,".col_idx.u32.0based.bin"));file.create(paste0(p,".values.f32.bin"));writeLines(c(paste0("n_variants=",m),"nnz=0"),paste0(p,".meta.txt"));ids<-paste0("m",seq_len(m));tn<-paste0("T",seq_len(traits));wy<-lapply(seq_len(traits),function(t)stats::setNames(20*sin(seq_len(m)/17+t),ids));ww<-replicate(traits,stats::setNames(rep(499,m),ids),simplify=FALSE);names(wy)<-names(ww)<-tn;list(prefix=p,stats=list(yy=stats::setNames(rep(499,traits),tn),wy=wy,ww=ww,n=500L,m=m),Glist=list(rsidsLD=list(ids),rsids=list(ids),maf=list(rep(.25,m))))
}
run<-function(x,chains,cores,keep=FALSE)stblr_csr_bayesr(x$stats,x$Glist,x$prefix,nit=80L,nburn=20L,nthin=1L,nchains=chains,ncores=cores,keep_chains=keep,seed=31L,updateE=FALSE,updateLDswap=FALSE)
bench<-function(label,x,chains,cores,keep=FALSE){invisible(run(x,chains,cores,keep));z<-replicate(5,system.time(run(x,chains,cores,keep))[["elapsed"]]);rss<-if(requireNamespace("ps",quietly=TRUE))unname(ps::ps_memory_info()[["rss"]])/1024^2 else NA_real_;cat(label,"times",paste(z,collapse=","),"mean",mean(z),"median",median(z),"min",min(z),"max",max(z),"IQR",IQR(z),"rss_mb_sample",rss,"\n")}
cat("R",R.version.string,"Rcpp",as.character(packageVersion("Rcpp")),"\n")
tiny<-phase5a_bayesr_configs$one_chain;invisible(phase5a_bayesr_run(tiny,FALSE));cat("tiny_elapsed",system.time(phase5a_bayesr_run(tiny,FALSE))[["elapsed"]],"\n")
x<-make_case();bench("moderate_1c1core",x,1L,1L);bench("moderate_2c1core",x,2L,1L);bench("moderate_2c2core",x,2L,2L);bench("moderate_2c2core_chains",x,2L,2L,TRUE)
cat("RSS is sampled whole-process resident memory after each configuration, not instantaneous sampler-only peak; interval sampling was unavailable.\n")
