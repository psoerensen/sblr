source("tests/testthat/helper-mtblr-csr-fixtures.R",local=TRUE)
source("tests/testthat/helper-mtblr-block-eigen-fixtures.R",local=TRUE)
pkgload::load_all(".",compile=FALSE,quiet=TRUE)
run <- function(nt,shared,filters,blocks) {
  x <- phase17l_case(nt=nt,shared=shared,filters=filters,blocks=blocks)
  block_time <- system.time(do.call(sblr:::mtblr_block_eigen_internal,x$block))[["elapsed"]]
  dense_time <- system.time(do.call(sblr:::mtblr,x$dense))[["elapsed"]]
  sizes <- unlist(lapply(x$inspections,`[[`,"block_size"))
  data.frame(traits=nt,markers=length(x$transformed[[1]]),reference_samples=8L,
    blocks=length(sizes),sharing=if(shared) "shared" else "trait_specific",
    filters=paste(filters,collapse=","),iterations=x$block$nit,burnin=x$block$nburn,
    thinning=x$block$nthin,build_and_mcmc_seconds=block_time,
    dense_mcmc_seconds=dense_time,packed_bytes=sum(sizes*(sizes+1)/2)*4,
    dense_estimated_bytes=nt*length(x$transformed[[1]])^2*8)
}
print(rbind(
  run(2L,TRUE,rep("hard_truncate",2),rep(list(c(0L,2L)),2)),
  run(3L,FALSE,c("hard_truncate","ridge_fixed","ridge_lw"),
      list(c(0L,2L),c(0L,1L,3L),c(0L,3L)))),row.names=FALSE)
set.seed(17L); m <- 500L; nt <- 2L; nref <- 40L
dosage <- matrix(sample(c(0:2,NA),m*nref,replace=TRUE,prob=c(.3,.4,.28,.02)),m,nref)
descriptor <- phase17l_descriptor(dosage,seq.int(0L,m-1L,by=25L),"ridge_fixed",eta=.25,af=rep(.35,m))
wy <- list(seq(-1,1,length.out=m),seq(.8,-.8,length.out=m))
inspection <- phase17l_inspect(descriptor,wy)
inspections <- lapply(seq_len(nt),function(t){z<-inspection;z$transformed_wy<-inspection$transformed_wy[t,,drop=FALSE];z})
op <- phase17l_dense_operator(inspections); transformed <- lapply(inspections,function(x)as.numeric(x$transformed_wy[1,]))
models <- list(c(0L,0L),c(1L,0L),c(0L,1L),c(1L,1L)); mat_list <- function(x)split(x,rep(seq_len(ncol(x)),each=nrow(x)))
common <- list(yy=c(50,52),b=rep(list(rep(0,m)),nt),sets=list(0:(m-1L)),B=diag(.15,nt),E=diag(.8,nt),
 ssb_prior=mat_list(diag(.05,nt)),sse_prior=mat_list(diag(.3,nt)),models=models,pi=c(.7,.1,.1,.1),nub=4,nue=4,
 updateB=FALSE,updateE=FALSE,updatePi=FALSE,n=c(41L,53L),nit=2L,nburn=1L,nthin=1L,seed=17012L,method=4L)
block_args <- c(list(wy=wy),common[c("yy","b")],list(operator_descriptors=list(descriptor)),common[setdiff(names(common),c("yy","b"))])
dense_args <- c(list(wy=transformed,ww=op$ww),common[c("yy","b")],list(XXvalues=op$XXvalues,XXindices=op$XXindices),common[setdiff(names(common),c("yy","b"))])
cat(sprintf("moderate traits=%d markers=%d reference_samples=%d blocks=%d build_mcmc=%.3f dense=%.3f packed_bytes=%d dense_bytes=%d\n",
 nt,m,nref,length(descriptor$block_start),system.time(do.call(sblr:::mtblr_block_eigen_internal,block_args))[["elapsed"]],
 system.time(do.call(sblr:::mtblr,dense_args))[["elapsed"]],sum(25L*26L/2L)*4L*(m/25L),nt*m*m*8L))
cat("Regression signal only; completed-fit RSS is not peak RSS.\n")
