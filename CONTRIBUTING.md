# Contributing to ctsgimme

Contributions and reproducible bug reports are welcome. Please open an issue
before beginning a large API or statistical-method change so its intended
scope can be discussed first.

## Development workflow

1. Fork the repository and create a focused branch from `main`.
2. Make changes in the package root; do not commit generated analysis output,
   `*.Rcheck` directories, or source archives.
3. Regenerate documentation with `devtools::document()` after editing roxygen
   comments.
4. Run the package tests and source-package check:

```r
devtools::test(stop_on_failure = TRUE)
```

```sh
R CMD build .
R CMD check --no-manual ctsgimme_*.tar.gz
```

5. Describe the behavioral change and validation performed in the pull
   request. Add or update tests for any user-visible change.

Please avoid including identifiable participant data in issues, fixtures, or
pull requests.
