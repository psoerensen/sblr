phase7a_sbayesrc_starting_commit <- "f6785f7"

phase7a_sbayesrc_configs <- list(
  fixed_one_chain=list(traits=1L,nchains=1L,ncores=1L,keep=FALSE,seeds=NULL,update=FALSE,every=10L,gamma=c(0,.1,1)),
  fixed_two_chains=list(traits=1L,nchains=2L,ncores=1L,keep=FALSE,seeds=NULL,update=FALSE,every=10L,gamma=c(0,.1,1)),
  learned_two_cores=list(traits=1L,nchains=2L,ncores=2L,keep=FALSE,seeds=NULL,update=TRUE,every=1L,gamma=c(0,.1,1)),
  learned_explicit_keep=list(traits=1L,nchains=2L,ncores=1L,keep=TRUE,seeds=c(701L,702L),update=TRUE,every=2L,gamma=c(0,.1,1)),
  multiple_traits=list(traits=2L,nchains=2L,ncores=2L,keep=FALSE,seeds=c(711L,712L),update=TRUE,every=1L,gamma=c(0,.1,1)),
  explicit_scales=list(traits=1L,nchains=1L,ncores=1L,keep=FALSE,seeds=NULL,update=FALSE,every=3L,gamma=c(0,.025,.25,1))
)

phase7a_sbayesrc_write_csr <- function() {
  prefix<-tempfile("blr_phase7a_sbayesrc_")
  sblr:::.stblr_write_uint64_file(paste0(prefix,".row_ptr.u64.bin"),c(0,1,2,2,2))
  sblr:::.stblr_write_uint32_file(paste0(prefix,".col_idx.u32.0based.bin"),c(1L,2L))
  writeBin(c(.25,-.15),paste0(prefix,".values.f32.bin"),size=4,endian="little")
  writeLines(c("format=sparse_ld_csr","storage=streamed_upper_triangle","n_bed=NA","n_used=NA","n_samples_used=NA","n_variants=4","nnz=2","triangle=upper","diagonal=implicit_1",paste0("row_ptr_file=",prefix,".row_ptr.u64.bin"),paste0("col_idx_file=",prefix,".col_idx.u32.0based.bin"),paste0("values_file=",prefix,".values.f32.bin"),"row_ptr_type=uint64","col_idx_type=uint32","values_type=float32","index_base=0","value=r"),paste0(prefix,".meta.txt"))
  prefix
}

phase7a_sbayesrc_inputs <- function(traits=1L) {
  ids<-paste0("m",1:4); tn<-paste0("T",seq_len(traits)); scores<-list(c(30,-20,12,5),c(-14,25,8,-18))[seq_len(traits)]
  wy<-lapply(scores,stats::setNames,nm=ids); ww<-replicate(traits,stats::setNames(rep(80,4),ids),simplify=FALSE); names(wy)<-names(ww)<-tn
  A<-matrix(c(1,0,1, 1,1,0, 1,0,0, 1,1,1),4,3,byrow=TRUE,dimnames=list(ids,c("intercept","coding","qtl")))
  list(stats=list(wy=wy,ww=ww,yy=stats::setNames(rep(80,traits),tn),n=80L,m=4L),A=A)
}

phase7a_sbayesrc_run <- function(config,raw=FALSE) {
  x<-phase7a_sbayesrc_inputs(config$traits); captured<-NULL
  native<-function(...) { captured<<-do.call(sblr:::stblr_cpg_omp_csr_sbayesrc,list(...)); captured }
  fit<-sblr:::.stblr_csr_sbayesrc_generic_impl(stats=x$stats,ld_prefix=phase7a_sbayesrc_write_csr(),A=x$A,gamma=config$gamma,pi_init=.35,pi_prior_mean=.35,pi_prior_strength=2,add_intercept=FALSE,standardize_annotations=FALSE,updateAlpha=config$update,alpha_update_every=config$every,updateB=TRUE,updateE=FALSE,nit=10L,nburn=2L,nthin=1L,ncores=config$ncores,seed=73L,nchains=config$nchains,chain_seeds=config$seeds,keep_chains=config$keep,updateLDswap=FALSE,.native_fun=native)
  if(raw) captured else fit
}

phase7a_sbayesrc_normalize <- function(x) {
  if(!is.list(x)) return(x)
  for(nm in names(x)) if(nm %in% c("seconds_mean","seconds_max")) x[[nm]][]<-0 else if(nm=="ld_prefix") x[[nm]]<-"<fixture>" else x[[nm]]<-phase7a_sbayesrc_normalize(x[[nm]])
  x
}

phase7a_sbayesrc_metadata <- function(name,cfg) {
 init<-sblr::make_sbayesrc_alpha_init(phase7a_sbayesrc_inputs(1L)$A,gamma=cfg$gamma,pi_init=.35)
 list(starting_commit=phase7a_sbayesrc_starting_commit,R_version=R.version.string,compiler="Rtools44 GCC 13.2 C++17",fixture="4-marker CSR, 3 annotations including explicit intercept",marker_count=4L,trait_count=cfg$traits,annotation_count=3L,annotation_names=c("intercept","coding","qtl"),intercept_convention="explicit first column; intercept_flat",component_count=length(cfg$gamma),component_scales=cfg$gamma,null_component=0L,initial_component_probabilities=init$component_prob_init,alpha_dimensions=dim(init$alpha_init),alpha_initial_values=init$alpha_init,alpha_prior=c(a=2,b=2),updateAlpha=cfg$update,alpha_update_every=cfg$every,seed=73L,chain_seeds=cfg$seeds,nit=10L,nburn=2L,nthin=1L,nchains=cfg$nchains,ncores=cfg$ncores,keep_chains=cfg$keep,ld_swap=FALSE,schema_name="stblr_raw",schema_version=1L,configuration=name)
}
