pkgload::load_all(".",compile=FALSE)
source(file.path("tests","testthat","fixtures","blr-phase7a-sbayesrc-reference.R"))

bench<-function(label,cfg,reps=3L) {
  invisible(phase7a_sbayesrc_run(cfg,FALSE))
  rss<-numeric(reps); elapsed<-numeric(reps)
  for(i in seq_len(reps)) {
    elapsed[i]<-system.time(invisible(phase7a_sbayesrc_run(cfg,FALSE)))[["elapsed"]]
    rss[i]<-if(requireNamespace("ps",quietly=TRUE)) as.numeric(ps::ps_memory_info()[["rss"]])/1024^2 else NA_real_
  }
  data.frame(configuration=label,elapsed=paste(elapsed,collapse=","),mean=mean(elapsed),median=median(elapsed),minimum=min(elapsed),maximum=max(elapsed),rss_mib=max(rss,na.rm=TRUE),chains=cfg$nchains,cores=cfg$ncores,updateAlpha=cfg$update,alpha_update_every=cfg$every)
}

configs<-phase7a_sbayesrc_configs[c("fixed_one_chain","fixed_two_chains","learned_two_cores","learned_explicit_keep")]
results<-do.call(rbind,Map(bench,names(configs),configs))
print(results,row.names=FALSE)
cat("fixture: 4 markers, 1 trait, 3 annotations, 3 components, 10 iterations, 8 retained\n")

moderate_run<-function(updateAlpha,nchains,ncores) {
 m<-500L; ids<-paste0("m",seq_len(m)); prefix<-tempfile("phase7a_moderate_")
 sblr:::.stblr_write_uint64_file(paste0(prefix,".row_ptr.u64.bin"),rep(0,m+1L)); file.create(paste0(prefix,".col_idx.u32.0based.bin")); file.create(paste0(prefix,".values.f32.bin"))
 writeLines(c("format=sparse_ld_csr","storage=streamed_upper_triangle","n_bed=NA","n_used=NA","n_samples_used=NA",paste0("n_variants=",m),"nnz=0","triangle=upper","diagonal=implicit_1",paste0("row_ptr_file=",prefix,".row_ptr.u64.bin"),paste0("col_idx_file=",prefix,".col_idx.u32.0based.bin"),paste0("values_file=",prefix,".values.f32.bin"),"row_ptr_type=uint64","col_idx_type=uint32","values_type=float32","index_base=0","value=r"),paste0(prefix,".meta.txt"))
 stats<-list(wy=list(T1=stats::setNames(8*sin(seq_len(m)/17),ids)),ww=list(T1=stats::setNames(rep(200,m),ids)),yy=c(T1=200),n=200L,m=m)
 A<-cbind(intercept=1,coding=as.numeric(seq_len(m)%%5==0),qtl=as.numeric(seq_len(m)%%11==0),continuous=scale(seq_len(m))[,1]); rownames(A)<-ids
 sblr::stblr_csr_sbayesrc_generic(stats,prefix,A,gamma=c(0,.01,.1,1),add_intercept=FALSE,standardize_annotations=FALSE,updateAlpha=updateAlpha,alpha_update_every=5L,updateE=FALSE,nit=40L,nburn=10L,nchains=nchains,ncores=ncores,seed=811L)
}
for(z in list(c(FALSE,1L,1L),c(TRUE,2L,1L),c(TRUE,2L,2L))) {
 invisible(moderate_run(z[[1]],z[[2]],z[[3]])); times<-replicate(3,system.time(invisible(moderate_run(z[[1]],z[[2]],z[[3]])))[["elapsed"]])
 cat(sprintf("moderate updateAlpha=%s chains=%d cores=%d times=%s mean=%.3f median=%.3f\n",z[[1]],z[[2]],z[[3]],paste(times,collapse=","),mean(times),median(times)))
}
cat("moderate fixture: 500 markers, 1 trait, 4 annotations, 4 components, 40 iterations, 30 retained\n")
cat("environment:",R.version.string,"; sblr",as.character(utils::packageVersion("sblr")),"\n")
cat("memory: whole-process RSS sampled after each run; not sampler-only peak; no interval peak sampler available\n")
