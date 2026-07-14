# Manual maintenance tool. Never run automatically from tests.
pkgload::load_all(".",compile=FALSE)
source(file.path("tests","testthat","fixtures","blr-phase9a-annotation-reference.R"))
for(backend in names(phase9a_configs)){
 dir.create(file.path("tests","testthat","fixtures",paste0("blr_phase9a_",backend)),recursive=TRUE,showWarnings=FALSE)
 for(nm in names(phase9a_configs[[backend]])){
  cfg<-phase9a_configs[[backend]][[nm]]
  saveRDS(list(metadata=phase9a_metadata(backend,nm,cfg),raw=phase9a_normalize(phase9a_run(backend,cfg,TRUE)),fit=phase9a_normalize(phase9a_run(backend,cfg,FALSE))),file.path("tests","testthat","fixtures",paste0("blr_phase9a_",backend),paste0(nm,".rds")),version=3)
 }
}
