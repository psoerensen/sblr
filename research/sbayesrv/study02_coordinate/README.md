# SBayesRV at the Study 02 coordinate

This directory preserves one descriptive SBayesRV replicate using the exact
validated Study 02 data and operator coordinate. It reuses qgdata
`human_independent`, the historical 37,991-marker QC panel, fixed 3,500/1,500
split, training-subset allele frequencies, and a content-identified cached CSR
operator. The code does not rebuild LD or replace accepted Study 02 evidence.

The three frozen stages and fitting controls are documented in `design.md`.
The compact results already obtained are tracked in `findings.md`. Generated
rerun output belongs only below ignored `output/`; fit objects, genotype
matrices, marker-by-draw arrays, and duplicate LD files are not retained.

The analysis is intentionally not launched by the root SBayesRV entry point.
An explicit rerun, if scientifically authorized later, is:

```powershell
Rscript research/sbayesrv/study02_coordinate/analysis.R
```

This is one research replicate, not benchmark evidence or a production
qualification.
