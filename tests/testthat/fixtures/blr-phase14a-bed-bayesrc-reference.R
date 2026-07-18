phase14a_fixture_dir <- if (file.exists(file.path("tests","testthat","fixtures","blr-phase11a-bed-reference.R")))
 file.path("tests","testthat","fixtures") else file.path(testthat::test_path(),"fixtures")
source(file.path(phase14a_fixture_dir, "blr-phase11a-bed-reference.R"))

phase14a_capture <- function(ncores=1L,nchains=1L,seed=141L,updateAlpha=TRUE,
                             multiple_annotations=TRUE) {
 x <- phase11a_fixture(); ids <- x$Glist$rsids[[1L]]
 annotation <- if (multiple_annotations)
  cbind(intercept=1,enriched=c(0,1),continuous=c(-0.5,0.5)) else
  matrix(1,2L,1L,dimnames=list(ids,"intercept"))
 rownames(annotation) <- ids
 ns <- asNamespace("sblr"); assign(".phase14a_raw",NULL,envir=.GlobalEnv)
 suppressMessages(trace(".as_stblr_fit",tracer=quote(assign(".phase14a_raw",raw,envir=.GlobalEnv)),where=ns,print=FALSE))
 on.exit(suppressMessages(try(untrace(".as_stblr_fit",where=ns),silent=TRUE)))
 fit <- sblr::stblr_bed(y=x$y,Glist=x$Glist,method="bayesrc",annotation=annotation,
  add_intercept=FALSE,standardize_annotations=FALSE,updateAlpha=updateAlpha,
  mixture_var=c(0,.01,.1,1),pi=c(.95,.03,.015,.005),nit=8L,nburn=2L,nthin=1L,
  seed=seed,ncores=ncores,nchains=nchains,updateB=FALSE,updateE=FALSE,
  rebuild_every=2L,read_block_size=2L,return_wy=TRUE,return_r=TRUE)
 list(raw=get(".phase14a_raw",envir=.GlobalEnv),fit=fit)
}

phase14a_normalize <- phase11a_normalize
