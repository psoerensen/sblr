phase9a_starting_commit <- "4e7d6b4"

phase9a_prefix <- function(m=4L) {
 p<-tempfile("phase9a_annot_"); sblr:::.stblr_write_uint64_file(paste0(p,".row_ptr.u64.bin"),rep(0,m+1L)); file.create(paste0(p,".col_idx.u32.0based.bin")); file.create(paste0(p,".values.f32.bin"))
 writeLines(c("format=sparse_ld_csr","storage=streamed_upper_triangle","n_bed=NA","n_used=NA","n_samples_used=NA",paste0("n_variants=",m),"nnz=0","triangle=upper","diagonal=implicit_1",paste0("row_ptr_file=",p,".row_ptr.u64.bin"),paste0("col_idx_file=",p,".col_idx.u32.0based.bin"),paste0("values_file=",p,".values.f32.bin"),"row_ptr_type=uint64","col_idx_type=uint32","values_type=float32","index_base=0","value=r"),paste0(p,".meta.txt")); p
}
phase9a_inputs <- function(nt=1L) {
 ids<-paste0("m",1:4); tn<-paste0("trait",seq_len(nt)); wy<-lapply(seq_len(nt),function(t) stats::setNames(c(2,-1,.5,1)*(1+.1*(t-1)),ids)); ww<-lapply(seq_len(nt),function(t) stats::setNames(rep(50+5*(t-1),4),ids)); names(wy)<-names(ww)<-tn
 A<-cbind(intercept=1,coding=c(1,0,1,0),qtl=c(0,1,1,0)); rownames(A)<-ids
 list(stats=list(wy=wy,ww=ww,yy=stats::setNames(50+5*(seq_len(nt)-1),tn),n=rep(50L,nt),m=4L,marker_names=ids,trait_names=tn),A=A,group=stats::setNames(c("coding","background","coding","background"),ids))
}
phase9a_configs <- list(
 prior=list(fixed_one=list(nt=1L,nchains=1L,ncores=1L,keep=FALSE,seeds=NULL),fixed_chains=list(nt=1L,nchains=2L,ncores=2L,keep=TRUE,seeds=c(11L,12L)),fixed_explicit=list(nt=1L,nchains=2L,ncores=1L,keep=FALSE,seeds=c(21L,22L))),
 group=list(group_one=list(nt=1L,nchains=1L,ncores=1L,keep=FALSE,seeds=NULL,normalize=TRUE),group_chains=list(nt=1L,nchains=2L,ncores=2L,keep=TRUE,seeds=c(31L,32L),normalize=FALSE),group_explicit=list(nt=1L,nchains=2L,ncores=1L,keep=FALSE,seeds=c(41L,42L),normalize=TRUE)),
 annotation=list(annot_fixed=list(nt=1L,nchains=1L,ncores=1L,keep=FALSE,seeds=NULL,learn_pi=FALSE,learn_vb=FALSE,every=2L),annot_learned=list(nt=1L,nchains=2L,ncores=2L,keep=TRUE,seeds=c(51L,52L),learn_pi=TRUE,learn_vb=TRUE,every=1L),annot_explicit=list(nt=1L,nchains=2L,ncores=1L,keep=FALSE,seeds=c(61L,62L),learn_pi=TRUE,learn_vb=TRUE,every=2L)))

phase9a_run <- function(backend,cfg,raw=FALSE) {
 z<-phase9a_inputs(cfg$nt); common<-list(stats=z$stats,ld_prefix=phase9a_prefix(),h2=.5,pi_init=.3,pi_prior_mean=.3,pi_prior_strength=2,updateB=FALSE,updateE=FALSE,updatePi=FALSE,nit=6L,nburn=2L,nthin=1L,ncores=cfg$ncores,seed=909L,nchains=cfg$nchains,keep_chains=cfg$keep,chain_seeds=cfg$seeds,updateLDswap=FALSE,convergence="none")
 if(raw){ trace(".as_stblr_fit",where=asNamespace("sblr"),tracer=quote(assign(".phase9a_raw_capture",raw,envir=.GlobalEnv)),print=FALSE); on.exit(untrace(".as_stblr_fit",where=asNamespace("sblr")),add=TRUE) }
 fit<-switch(backend,
  prior=do.call(sblr::stblr_csr_annot,c(common,list(annotation_model="prior",annotations=list(A=z$A,fixed_pi_marker=rep(list(c(.2,.35,.5,.65)),cfg$nt),fixed_vb_multiplier=rep(list(c(.7,1,1.4,2)),cfg$nt),use_pi_marker=TRUE,use_vb_multiplier=TRUE)))),
  group=do.call(sblr::stblr_csr_annot,c(common,list(annotation_model="group",annotations=z$group,group_names=c("coding","background"),group_pi_init=c(.35,.2),group_vb_multiplier_init=c(1.4,.7),updateGroupVb=TRUE,normalize_group_vb=cfg$normalize))),
  annotation=do.call(sblr::stblr_csr_annot,c(common,list(annotation_model="learned",annotations=z$A,add_intercept=FALSE,standardize_annotations=FALSE,learn_pi_annot=cfg$learn_pi,learn_vb_annot=cfg$learn_vb,eta_pi_init=matrix(0,ncol(z$A),cfg$nt),eta_vb_init=matrix(0,ncol(z$A),cfg$nt),rw_sd_eta_pi=.02,rw_sd_eta_vb=.02,annot_update_every=cfg$every,pi_min=.01,pi_max=.9,vb_multiplier_min=.1,vb_multiplier_max=10))))
 if(raw){ out<-get(".phase9a_raw_capture",envir=.GlobalEnv); rm(".phase9a_raw_capture",envir=.GlobalEnv); out } else fit
}
phase9a_normalize <- function(x){
 if(!is.list(x)) return(x)
 nms<-names(x)
 for(i in seq_along(x)){
  nm<-if(is.null(nms)) NA_character_ else nms[[i]]
  if(!is.na(nm)&&nm%in%c("seconds_mean","seconds_max")) x[[nm]][]<-0
  else if(!is.na(nm)&&identical(nm,"ld_prefix")) x[[nm]]<-"<fixture>"
  else x[i]<-list(phase9a_normalize(x[[i]]))
 }
 x
}
phase9a_pick <- function(x,current,historical=current) if(!is.null(x[[current]])) x[[current]] else x[[historical]]
phase9a_drop_null <- function(x){
 if(!is.list(x)) return(x)
 x<-x[!vapply(x,is.null,logical(1))]
 lapply(x,phase9a_drop_null)
}
phase9a_raw_science <- function(x) phase9a_drop_null(
 x[c("marker","trace","variance","pi","prior","group","annotation")])
phase9a_fit_science <- function(x) {
 extras <- c("prior","group","group_pi","group_vb_multiplier","group_nincluded",
  "group_size","annotation","annotation_prior","annotation_summary","annotation_pi",
  "annotation_variance","eta_pi","eta_vb","annotation_effects")
 c(list(bm=x$bm,dm=x$dm,wy=x$wy,r=x$r,b=x$b,d=x$d,vbs=x$vbs,vgs=x$vgs,
  ves=x$ves,vle=x$vle,vld=x$vld,pi_trace=phase9a_pick(x,"pi_trace","pis"),
  pi_final=phase9a_pick(x,"pi_final","pi"),pi_mean=phase9a_pick(x,"pi_mean","pim"),
  cov_b_mean=phase9a_pick(x,"cov_b_mean","covb"),cov_g_mean=phase9a_pick(x,"cov_g_mean","covg"),
  cov_e_mean=phase9a_pick(x,"cov_e_mean","cove"),cov_b_final=phase9a_pick(x,"cov_b_final","vb"),
  cov_g_final=phase9a_pick(x,"cov_g_final","vg"),cov_e_final=phase9a_pick(x,"cov_e_final","ve")),
  phase9a_drop_null(x[intersect(extras,names(x))]))
}
phase9a_metadata <- function(backend,name,cfg) list(starting_commit=phase9a_starting_commit,backend=backend,R_version=R.version.string,compiler="Rtools44 GCC 13.2 C++17",marker_count=4L,trait_count=cfg$nt,annotation_count=3L,group_count=if(backend=="group")2L else NA_integer_,marker_order=paste0("m",1:4),annotation_order=c("intercept","coding","qtl"),group_order=if(backend=="group")c("coding","background") else NULL,intercept_convention="explicit input column; no implicit intercept",seed=909L,chain_seeds=cfg$seeds,nit=6L,nburn=2L,nthin=1L,nchains=cfg$nchains,ncores=cfg$ncores,keep_chains=cfg$keep,LD_swap=FALSE,schema="stblr_raw_v1",configuration=name)
