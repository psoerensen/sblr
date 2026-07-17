# Manual maintenance tool. Never sourced by ordinary tests.
source(file.path("tests","testthat","fixtures","blr-phase10b-scheduled-reference.R"))
dir.create(file.path("tests","testthat","fixtures","blr_phase10b_scheduled_csr"),
           recursive=TRUE,showWarnings=FALSE)
for (nm in names(phase10b_configs)) {
  result <- callr::r(function(root,name) {
    setwd(root); pkgload::load_all(".",compile=FALSE,quiet=TRUE)
    source(file.path("tests","testthat","fixtures","blr-phase10b-scheduled-reference.R"))
    phase10b_run(phase10b_configs[[name]])
  },list(root=normalizePath(".",winslash="/"),name=nm),show=TRUE)
  saveRDS(c(list(metadata=phase10b_metadata(nm,phase10b_configs[[nm]])),result),
    file.path("tests","testthat","fixtures","blr_phase10b_scheduled_csr",paste0(nm,".rds")),version=3)
}
