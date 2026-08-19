# Publishing ctsgimme 0.0.11

The GitHub repository is named `ctgimme`; the R package is named `ctsgimme`.
The repository can be renamed to `ctsgimme` later for consistency, but all URLs
in this release must use the actual repository name at publication time.

## Local release checks

From the repository root:

```powershell
Rscript --vanilla -e "devtools::document(); devtools::test(stop_on_failure = TRUE)"
New-Item -ItemType Directory -Force dist | Out-Null
Push-Location dist
R CMD build ..
R CMD check --no-manual ctsgimme_0.0.11.tar.gz
Pop-Location
```

Confirm that `dist/ctsgimme.Rcheck/00check.log` ends with `Status: OK` and
install the exact archive into a temporary library before publishing it.

The final Windows validation for this repository used R 4.4.2, passed all 158
test expectations, completed `R CMD check --no-manual` with `Status: OK`, and
successfully installed and loaded the exact source archive in a fresh library.
The release archive is `dist/ctsgimme_0.0.11.tar.gz` (SHA-256
`29ABB94D6D703D25FBAC33C2947776A51C0ECA6B165063D53F537816FB3B6A08`).

## Git publication

Create and review the release commit, merge it locally, and tag the checked
commit before publishing:

```sh
git switch release/0.0.11
git add -A
git commit -m "Release ctsgimme 0.0.11"
git switch main
git merge --ff-only release/0.0.11
git tag -a v0.0.11 -m "ctsgimme 0.0.11"
```

Then publish the already-reviewed branch and tags:

```sh
git push origin main
git push origin pre-package-0.0.6 v0.0.11
```

Create a GitHub release for `v0.0.11` and attach:

- `dist/ctsgimme_0.0.11.tar.gz`;
- `dist/ctsgimme_0.0.11_multisubject_figures.pdf`, the validated subgroup
  parameter and transition figures (SHA-256
  `3C806C35E33210FD384E9B046B9DAD491F6D08E2383876AB7F72EFB1352C07EE`).

Suggested release title: `ctsgimme 0.0.11`

Suggested release summary:

> ctsgimme 0.0.11 replaces chained subgroup trajectories with one
> shared-parameter multisubject likelihood per subgroup. Each member retains
> its own observation times, likelihood, and filter initialization, while
> subgroup parameters are jointly estimated once. `time.intervals` continues
> to generate discrete-time transition plots, and each subgroup writes one
> fitted `Subgroup_<g>Model.RDS`. The obsolete chaining selector and
> concatenation controls have been removed; see NEWS.md for migration details.

## GitHub repository metadata

Suggested About description:

> Continuous-time group, subgroup, and individual network estimation for
> intensive longitudinal data.

Suggested topics:

`r-package`, `continuous-time`, `time-series`, `state-space-models`,
`vector-autoregression`, `intensive-longitudinal-data`, `gimme`, `openmx`

If the GitHub repository is renamed to `ctsgimme`, update `DESCRIPTION`, README
badges and install commands, `RELEASE.md`, and the local remote before tagging:

```sh
git remote set-url origin https://github.com/JPark93/ctsgimme.git
```
