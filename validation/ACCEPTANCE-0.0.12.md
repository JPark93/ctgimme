# ctgimme 0.0.12 installed-archive acceptance

This record covers the exact checked source archive:

- Archive: `ctgimme_0.0.12.tar.gz`
- Size: 51,021 bytes
- SHA-256:
  `485E698808387C8FB15E30EC9B464878BCF88B82BEF710ED2D756CE955AF1ABA`
- Check runtime: Windows 11, R 4.6.1
- Check result: 0 errors, 0 warnings, 2 local notes (new submission and no
  local HTML Tidy executable)
- Tests: 190 passed, 0 failed/warned/skipped
- Main example: 0.68 seconds

## Truth-controlled multisubject parity

The archive was installed into a fresh library under R 4.6.1/OpenMx 2.22.11
and evaluated against the frozen R 4.4.2/OpenMx 2.21.13 results. The simulated
fixture and reference inputs were checksum-verified before fitting.

- Four of four fits converged with status 0.
- Each subgroup had ten shared free parameters and ten independently
  initialized subject likelihood blocks.
- Maximum estimate difference from 0.0.11: `5.0675e-9`.
- Maximum standard-error difference: `2.9558e-6`.
- Maximum -2LL difference: `4.0018e-11`.
- Maximum forward/reverse subject-order difference: `2.8323e-10`.
- Subgroup 1 MAE/RMSE: `0.03584/0.04344`.
- Subgroup 2 MAE/RMSE: `0.01889/0.02413`.
- Two fitted subgroup RDS files and eight parameter/transition PNGs passed
  semantic and file-signature checks.

## Complete public ISDSA workflow

The installed public `ctgimme()` workflow was then run with two workers,
`subgroup.model = TRUE`, `time.intervals = c(0.25, 1, 2)`, and
`verbose = FALSE`.

- Elapsed time: 433.28 seconds.
- Standard output, package messages, and warnings: none.
- Membership, group support, and both subgroup supports: exact.
- All 20 individual support masks matched the retained installed 0.0.10
  practice run exactly.
- Joint subgroup fits: both status 0, each with ten shared parameters and ten
  subject likelihood blocks.
- Individual fits: 20/20 converged; 19/20 supports were exact. Subject 5
  retained one additional path, and no true path was missing for any subject.
- All requested membership, structure, individual-model, subgroup-model, RDS,
  parameter-plot, and delta-t transition-plot artifacts were present.

The fixture is one stochastic realization whose within-subgroup subjects are
parameter clones. These results establish release regression parity and
finite-sample recovery for this practice case; they do not establish
frequentist unbiasedness, coverage, or general structure-recovery rates.
