# ISDSA simulated acceptance fixture

`isdsa_jpark26_fixture.rds` is a simulated 20-subject, 4,000-row fixture
created for the ISDSA practice analysis. It contains two parameter-clone
subgroups of ten subjects, five modeled variables, and the generating group,
subgroup, and individual drift matrices.

- Simulation seed: `31072026`
- Generating runtime: R 4.4.2, OpenMx 2.21.13
- Fixture SHA-256:
  `5ADD83906470012F9DCE2B0F89098683B02370A42E698EDD018C9A8B24412985`
- Original generator: `QuartoFile/ISDSA_JPark26.qmd` in the ISDSA practice
  project (generator SHA-256
  `9F12641712852375941C288D30FA2448A32B5513A4C0AB88493BFCDBD0CB36D7`)

The fixture is repository-only validation data. The complete `validation/`
directory is excluded from the CRAN source archive and ordinary package
checks. Its embedded generator path is repository-relative; the generator
hash above preserves provenance without publishing a workstation path.
