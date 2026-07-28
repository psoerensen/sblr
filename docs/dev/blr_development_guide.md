# BLR development guide

Changes must preserve posterior targets, operation order, RNG draw order,
logical-task seed mapping, operator preparation, and convergence mathematics
unless a documented correctness issue explicitly requires otherwise.

Use canonical model and operator identifiers. Add no compatibility aliases.
Public behavior belongs in portable tests; source structure belongs in one of
the four permanent audits. New scientific fixtures require a documented owner,
tolerance, and regeneration policy in `tests/testthat/fixtures/README.md`.

Before handoff run focused owners, the neutrality matrix when numerical code or
formatting changes, generated-interface and architecture audits, documentation
audit, workflow parsing, compilation/roxygen/Rd validation, the fast filter,
the full source suite, and the built-package check. Remove verified generated
objects, DLLs, tarballs, check directories, and temporary libraries.
