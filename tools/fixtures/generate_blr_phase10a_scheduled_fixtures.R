# Manual maintenance tool. Never sourced by ordinary tests.
source(file.path("tests","testthat","fixtures","blr-phase10a-scheduled-reference.R"))
dir.create(file.path("tests","testthat","fixtures","blr_phase10a_scheduled"),
           recursive=TRUE,showWarnings=FALSE)
for (nm in names(phase10a_configs)) {
  result <- callr::r(function(root,name) {
    setwd(root)
    pkgload::load_all(".", compile=FALSE, quiet=TRUE)
    source(file.path("tests","testthat","fixtures","blr-phase10a-scheduled-reference.R"))
    phase10a_run(phase10a_configs[[name]])
  }, list(root=normalizePath(".",winslash="/"),name=nm), show=TRUE)
  saveRDS(c(list(metadata=phase10a_metadata(nm,phase10a_configs[[nm]])),result),
          file.path("tests","testthat","fixtures","blr_phase10a_scheduled",paste0(nm,".rds")),
          version=3)
}
