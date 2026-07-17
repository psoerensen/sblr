# Phase 10A pre-migration scheduled ordinary-CSR benchmark (BayesC only).
suppressPackageStartupMessages(pkgload::load_all(".",compile=FALSE))

make_prefix <- function(m) {
  p <- tempfile("phase10a_bench_")
  sblr:::.stblr_write_uint64_file(paste0(p,".row_ptr.u64.bin"),rep(0,m+1L))
  file.create(paste0(p,".col_idx.u32.0based.bin")); file.create(paste0(p,".values.f32.bin"))
  writeLines(c("format=sparse_ld_csr","storage=streamed_upper_triangle","n_bed=NA",
    "n_used=NA","n_samples_used=NA",paste0("n_variants=",m),"nnz=0",
    "triangle=upper","diagonal=implicit_1",paste0("row_ptr_file=",p,".row_ptr.u64.bin"),
    paste0("col_idx_file=",p,".col_idx.u32.0based.bin"),
    paste0("values_file=",p,".values.f32.bin"),"row_ptr_type=uint64",
    "col_idx_type=uint32","values_type=float32","index_base=0","value=r"),
    paste0(p,".meta.txt")); p
}
rss_mb <- function() {
  if (requireNamespace("ps",quietly=TRUE)) {
    return(as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])/(1024^2))
  }
  NA_real_
}
run_case <- function(name,m,nchains,ncores,full,base,max,threshold,lifetime,wakeup,reps=5L) {
  ids <- paste0("m",seq_len(m)); prefix <- make_prefix(m)
  z <- seq_len(m); wy <- sin(z/17)*1.5 + cos(z/31)*.5
  stats <- list(wy=list(trait1=setNames(wy,ids)),ww=list(trait1=setNames(rep(100,m),ids)),
    yy=setNames(100,"trait1"),n=100L,m=m,marker_names=ids,trait_names="trait1")
  fit <- function() sblr::stblr_csr(stats=stats,ld_prefix=prefix,scheduled=TRUE,
    pi_init=.02,pi_prior_mean=.02,pi_prior_strength=50,updateB=FALSE,updateE=FALSE,
    updatePi=FALSE,nit=30L,nburn=10L,nthin=2L,seed=501L,nchains=nchains,
    ncores=ncores,chain_seeds=if(nchains==2L)c(601L,602L) else NULL,
    full_sweep_every=full,null_skip_base=base,null_skip_max=max,
    candidate_threshold=threshold,candidate_lifetime=lifetime,
    skip_nulls_burnin_only=FALSE,wakeup_ld_neighbors=wakeup,
    wakeup_diff_threshold=.01,wakeup_max_neighbors=5L,updateLDswap=FALSE)
  invisible(fit())
  before <- rss_mb(); elapsed <- numeric(reps)
  for(i in seq_len(reps)) elapsed[i] <- system.time(invisible(fit()))[["elapsed"]]
  after <- rss_mb()
  data.frame(model="scheduled CSR BayesC",case=name,markers=m,traits=1L,
    iterations=30L,retained_samples=10L,chains=nchains,cores=ncores,
    full_sweep_every=full,null_skip_base=base,null_skip_max=max,
    candidate_threshold=threshold,candidate_lifetime=lifetime,
    wakeup_ld_neighbors=wakeup,attempted_updates=NA_integer_,skipped_updates=NA_integer_,
    full_sweeps=NA_integer_,times=paste(elapsed,collapse=","),mean=mean(elapsed),
    median=median(elapsed),minimum=min(elapsed),maximum=max(elapsed),range=diff(range(elapsed)),
    completed_fit_rss_before_mb=before,completed_fit_rss_after_mb=after,
    memory_method="ps process RSS after completed fits; not peak RSS")
}
cases <- list(
  c("tiny_dense",60,1,1,1,1,1,0,0,FALSE),
  c("moderate_dense_1c1",2000,1,1,1,1,1,0,0,FALSE),
  c("moderate_dense_2c1",2000,2,1,1,1,1,0,0,FALSE),
  c("moderate_dense_2c2",2000,2,2,1,1,1,0,0,FALSE),
  c("moderate_aggressive",2000,1,1,10,20,200,.02,10,FALSE),
  c("moderate_conservative",2000,1,1,5,3,20,.002,20,FALSE),
  c("moderate_wakeup",2000,1,1,10,20,200,.02,10,TRUE)
)
results <- do.call(rbind,lapply(cases,function(x) run_case(x[[1]],as.integer(x[[2]]),
  as.integer(x[[3]]),as.integer(x[[4]]),as.integer(x[[5]]),as.integer(x[[6]]),
  as.integer(x[[7]]),as.numeric(x[[8]]),as.integer(x[[9]]),as.logical(x[[10]]))))
cat("Phase 10A scheduled CSR pre-migration baseline\n")
cat("R:",R.version.string,"\n")
cat("Compiler:",R.version$compiler,"\n")
cat("OpenMP task scheduling: static; production scheduler diagnostics are not returned\n")
print(results,row.names=FALSE)
