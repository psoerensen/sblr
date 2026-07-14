pkgload::load_all(".",compile=FALSE)
source(file.path("tests","testthat","fixtures","blr-phase7a-sbayesrc-reference.R"))
out<-file.path("tests","testthat","fixtures","blr_phase7a_sbayesrc"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(nm in names(phase7a_sbayesrc_configs)) {
  cfg<-phase7a_sbayesrc_configs[[nm]]
  saveRDS(list(metadata=phase7a_sbayesrc_metadata(nm,cfg),raw=phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg,TRUE)),fit=phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg,FALSE))),file.path(out,paste0(nm,".rds")),version=3)
}
