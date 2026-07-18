# Manual maintenance tool; ordinary tests never call this file.
pkgload::load_all(".",compile=FALSE,quiet=TRUE)
library(testthat)
source(file.path("tests","testthat","fixtures","blr-phase14a-bed-bayesrc-reference.R"))
out <- file.path("tests","testthat","fixtures","blr_phase14a_bed_bayesrc")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
cfg <- list(one_chain_one_core=c(1,1,141,1,0),two_chains_one_core=c(1,2,143,1,1),two_chains_two_cores=c(2,2,143,1,1))
for(nm in names(cfg)) {
 z<-cfg[[nm]]
 observed<-phase14a_normalize(phase14a_capture(z[1],z[2],z[3],as.logical(z[4]),as.logical(z[5])))
 reference<-list(raw=observed$raw,fit=observed$fit,metadata=list(starting_commit="c446d1c",R=R.version.string,
  reference_mode="normalized raw and formatted",cores=z[1],chains=z[2],seed=z[3],updateAlpha=as.logical(z[4]),
  multiple_annotations=as.logical(z[5]),components=c(0,.01,.1,1),null_component=0L))
 saveRDS(reference,file.path(out,paste0(nm,".rds")),version=3)
}
