phase5a_bayesr_starting_commit <- "b307db3"

phase5a_bayesr_configs <- list(
  one_chain = list(traits=1L, nchains=1L, ncores=1L, keep=FALSE, seeds=NULL,
                   mixture=c(0, .01, .1, 1), update_pi=TRUE),
  two_chains_one_core = list(traits=1L, nchains=2L, ncores=1L, keep=FALSE, seeds=NULL,
                             mixture=c(0, .01, .1, 1), update_pi=TRUE),
  two_chains_two_cores = list(traits=1L, nchains=2L, ncores=2L, keep=FALSE, seeds=NULL,
                              mixture=c(0, .01, .1, 1), update_pi=TRUE),
  multiple_traits = list(traits=2L, nchains=2L, ncores=2L, keep=FALSE, seeds=NULL,
                         mixture=c(0, .01, .1, 1), update_pi=TRUE),
  explicit_seeds_keep = list(traits=1L, nchains=2L, ncores=1L, keep=TRUE,
                             seeds=c(401L, 402L), mixture=c(0, .01, .1, 1), update_pi=TRUE),
  explicit_scales_fixed_pi = list(traits=1L, nchains=1L, ncores=1L, keep=FALSE,
                                  seeds=NULL, mixture=c(0, .025, .25), update_pi=FALSE)
)

phase5a_bayesr_write_csr <- function() {
  prefix <- tempfile("blr_phase5a_bayesr_")
  sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), c(0,1,2,2,2))
  sblr:::.stblr_write_uint32_file(paste0(prefix, ".col_idx.u32.0based.bin"), c(1L,2L))
  writeBin(c(.35,-.20), paste0(prefix, ".values.f32.bin"), size=4, endian="little")
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle", "n_bed=NA",
    "n_used=NA", "n_samples_used=NA", "n_variants=4", "nnz=2",
    "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=",prefix,".row_ptr.u64.bin"),
    paste0("col_idx_file=",prefix,".col_idx.u32.0based.bin"),
    paste0("values_file=",prefix,".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"), paste0(prefix,".meta.txt"))
  prefix
}

phase5a_bayesr_inputs <- function(traits) {
  ids <- paste0("m",1:4); tn <- paste0("T",seq_len(traits))
  scores <- list(c(42,-27,18,5), c(-16,35,9,-22))[seq_len(traits)]
  wy <- lapply(scores, stats::setNames, nm=ids)
  ww <- replicate(traits, stats::setNames(rep(99,4),ids), simplify=FALSE)
  names(wy) <- names(ww) <- tn
  list(stats=list(yy=stats::setNames(rep(99,traits),tn), ww=ww, wy=wy,
                  n=100L, m=4L),
       Glist=list(rsidsLD=list(ids),rsids=list(ids),maf=list(rep(.25,4))))
}

phase5a_bayesr_run <- function(config, raw=FALSE) {
  prefix <- phase5a_bayesr_write_csr(); x <- phase5a_bayesr_inputs(config$traits)
  pi <- c(.70, rep(.30/(length(config$mixture)-1L),length(config$mixture)-1L))
  alpha <- pi*20
  args <- list(stats=x$stats,Glist=x$Glist,ld_prefix=prefix,mixture_var=config$mixture,
    pi=pi,alpha=alpha,h2=.4,adjE=.9,nit=14L,nburn=4L,nthin=1L,ncores=config$ncores,
    seed=31L,nchains=config$nchains,keep_chains=config$keep,chain_seeds=config$seeds,
    updateB=TRUE,updateE=FALSE,updatePi=config$update_pi,updateLDswap=FALSE,
    method="sbayesr",convergence="none")
  fit <- do.call(stblr_csr,args)
  if (!raw) return(fit)
  nt <- config$traits; m <- 4L; vy <- rep(1,nt); active <- sum(pi[-1])
  B <- diag((vy*.4)/(m*active),nt); E <- diag(vy*.6,nt)
  ssb <- diag(((4-2)/4)*(vy*.4)/(m*active),nt)
  sse <- diag(((4-2)/4)*(vy*.6),nt)
  sblr:::stblr_cpg_omp_csr_bayesr(x$stats$wy,x$stats$ww,x$stats$yy,
    replicate(nt,rep(0,m),simplify=FALSE),replicate(nt,rep(0,m),simplify=FALSE),FALSE,
    x$stats$wy,FALSE,FALSE,prefix,B,E,
    split(ssb,rep(seq_len(nt),each=nt)),split(sse,rep(seq_len(nt),each=nt)),
    pi,config$mixture,alpha,4,4,TRUE,FALSE,config$update_pi,.9,rep(100L,nt),
    14L,4L,1L,config$ncores,31L,config$nchains,config$keep,
    if(is.null(config$seeds)) integer() else config$seeds,
    updateE_start=0L,updateE_every=1L,updateLDswap=FALSE)
}

phase5a_bayesr_normalize <- function(x) {
  if (!is.list(x)) return(x)
  nms <- names(x)
  for (i in seq_along(x)) {
    nm <- if (is.null(nms)) NA_character_ else nms[[i]]
    if (!is.na(nm) && nm %in% c("seconds_mean","seconds_max")) x[[nm]][] <- 0
    else if (!is.na(nm) && identical(nm,"ld_prefix")) x[[nm]] <- "<fixture>"
    else x[i] <- list(phase5a_bayesr_normalize(x[[i]]))
  }
  x
}
phase5a_bayesr_pick <- function(x,current,historical=current) if(!is.null(x[[current]])) x[[current]] else x[[historical]]
phase5a_bayesr_fit_science <- function(x) list(
 bm=x$bm,dm=x$dm,wy=x$wy,r=x$r,b=x$b,d=x$d,vbs=x$vbs,vgs=x$vgs,
 ves=x$ves,vle=x$vle,vld=x$vld,component=x$component,
 pi_final=as.numeric(phase5a_bayesr_pick(x,"pi_final","final_pi")),
 pi_mean=as.numeric(phase5a_bayesr_pick(x,"pi_mean","mean_pi")),
 cov_b_mean=phase5a_bayesr_pick(x,"cov_b_mean","covb"),
 cov_g_mean=phase5a_bayesr_pick(x,"cov_g_mean","covg"),
 cov_e_mean=phase5a_bayesr_pick(x,"cov_e_mean","cove"),
 cov_b_final=phase5a_bayesr_pick(x,"cov_b_final","vb"),
 cov_g_final=phase5a_bayesr_pick(x,"cov_g_final","vg"),
 cov_e_final=phase5a_bayesr_pick(x,"cov_e_final","ve"),
 component_probabilities=phase5a_bayesr_pick(x,"component_probabilities","comp_prob"),
 dm_component_mean=x$dm_component_mean,ncomp=x$ncomp,mixture_var=x$mixture_var,
 updateE_diagnostics=x$updateE_diagnostics,rb=x$rb,rg=x$rg,re=x$re)
phase5a_bayesr_drop_null <- function(x) {
 if(!is.list(x)) return(x)
 x <- x[!vapply(x,is.null,logical(1))]
 lapply(x,phase5a_bayesr_drop_null)
}
phase5a_bayesr_raw_science <- function(x)
 phase5a_bayesr_drop_null(x[c("marker","trace","variance","pi","component")])

phase5a_bayesr_metadata <- function(name,config) list(
  starting_commit=phase5a_bayesr_starting_commit,R_version=R.version.string,
  compiler="Rtools44 GCC 13.2.0 C++17",fixture="4-marker ordinary CSR, two off-diagonal entries",
  marker_count=4L,trait_count=config$traits,component_count=length(config$mixture),
  component_scales=config$mixture,seed=31L,chain_seeds=config$seeds,nit=14L,nburn=4L,
  nthin=1L,nchains=config$nchains,ncores=config$ncores,keep_chains=config$keep,
  update_component_probabilities=config$update_pi,ld_swap=FALSE,
  schema_name="stblr_raw",schema_version=1L,configuration=name)
