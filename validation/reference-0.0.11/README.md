# Frozen multisubject reference

These CSV files are the frozen numerical results from the validated
`ctsgimme` 0.0.11 multisubject implementation. They are used only to detect
regressions in the renamed `ctgimme` implementation; the historical chained
estimator is not a reference target.

- `fit_summary.csv` SHA-256:
  `36A24EB4EDAE4AB8460C48BA2678A8EAD0EA6C597369C000B8FA402F3F51E3EF`
- `parameter_estimates.csv` SHA-256:
  `F1A7EF1BC628C7E1A3C874A96880E2BC76F26359FD2DCC10BCC26504064DB31F`
- `subject_order_invariance.csv` SHA-256:
  `29602E001B766270E476475434F66B23C59636C77521300BB821DF472609DDC1`

The CSV hashes are calculated after normalizing line endings to LF so the
gate is invariant to Git checkout settings and operating system.

The acceptance test compares parameter keys and values, likelihoods,
architecture, subject-order invariance, and truth-recovery metrics. It does
not compare serialized RDS bytes or cross-platform PNG hashes.
