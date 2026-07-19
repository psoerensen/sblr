source("tools/benchmarks/blr_bed_family_benchmark_common.R")
bed_family_benchmark_notice()
scripts<-c(BayesC="tools/benchmarks/blr_phase11d_bed_scheduled_bayesc.R",
 BayesR="tools/benchmarks/blr_phase13e_bed_bayesr.R",
 BayesRC="tools/benchmarks/blr_phase14e_bed_bayesrc.R")
catalog<-data.frame(model=names(scripts),script=unname(scripts),
 controls=c("adaptive binary scheduler","ordered global mixture scheduler",
 "annotation probit-stick full sweep"),stringsAsFactors=FALSE)
cat("Phase 15B canonical packed-BED common reporting convention\n")
print(catalog,row.names=FALSE)
cat("Each canonical benchmark follows in its own R process for within-model comparison.\n")
for(i in seq_len(nrow(catalog))){
 cat("\n=== ",catalog$model[i]," canonical benchmark ===\n",sep="")
 status<-system2(file.path(R.home("bin"),"Rscript"),catalog$script[i])
 if(!identical(status,0L))stop(catalog$model[i]," benchmark failed with status ",status)
}
