# ctgimme 0.1.0 installed-archive acceptance

## Accepted archive

- Archive: `ctgimme_0.1.0.tar.gz`
- Size: 52,096 bytes
- SHA-256:
  `C23BC88796F2207FC16F4F7E0CD0A39C00A78D534BE9EDAB1F8AB330921C0005`
- Package version: 0.1.0
- Exports: `ctgimme`, `ctgimme_demo`

The archive was installed into a fresh library. The validation process loaded
that installed copy and did not load package source from the repository.

## ISDSA multisubject regression

The truth-controlled 20-subject ISDSA clone fixture was evaluated with two
subgroups of ten subjects in forward and reversed subject order. All four
fits converged with status 0. Each model had ten shared free parameters, ten
independently initialized subject likelihood blocks, and a multigroup fit
equal to the sum of the child likelihoods.

Against the frozen 0.0.11 multisubject reference:

- maximum estimate difference: `4.4408920985e-16`;
- maximum standard-error difference: `4.85722573274e-17`;
- maximum `-2LL` difference: `4.36557456851e-11`; and
- maximum forward/reverse estimate difference: `5.17531573152e-9`.

Subgroup 1 had mean signed error `-0.00686929`, MAE `0.03583639`, RMSE
`0.04344323`, and maximum absolute error `0.08407568`. Subgroup 2 had mean
signed error `-0.00033432`, MAE `0.01889391`, RMSE `0.02412866`, and maximum
absolute error `0.05143360`.

All quick-demo, local-time, architecture, truth-recovery, and subgroup PNG/RDS
artifact gates passed. Elapsed time was 63.92 seconds.

## Parallel-worker regression

A separate installed-package smoke test created four distinct PSOCK worker
processes. The same public toy analysis was then run at `cores = 1` and
`cores = 4`. Membership and group structure were identical, and the maximum
difference among fitted individual drift estimates was zero.

The package examples use `cores = 1`; real process-spawning package tests use
at most two workers. The four-worker test was run outside `R CMD check`.

## Environment

- Windows 11 x64 (build 26200)
- R 4.4.2
- OpenMx 2.21.13

This exact-support fixture is a regression test. It does not establish
frequentist unbiasedness or replace the slower public search/orchestration
acceptance recorded for 0.0.12.
