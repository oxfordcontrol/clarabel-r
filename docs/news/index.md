# Changelog

## clarabel 0.11.9003

- Added persistent solver API for warm starts:
  [`clarabel_solver()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md),
  [`solver_solve()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_solve.md),
  [`solver_update()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_update.md),
  [`solver_is_update_allowed()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_is_update_allowed.md)
- Robustified Rust interface: type coercion, error handling, CscMatrix
  validation, and regex ordering fixes
- Switched all R-side error messaging to
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
  with markup
- Upgraded `savvy` crate from 0.8.13 to 0.9.2
- Added examples to
  [`clarabel_control()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_control.md)
  and solver functions
- Added vignette section on updating problem data (warm starts)

## clarabel 0.11.1

CRAN release: 2025-09-24

- Synced up to v0.11.1 of `Clarabel.rs`
- Fixed `flang` issues shown on CRAN for fedora-clang-rdevel

## clarabel 0.10.1

CRAN release: 2025-04-17

- Minor fixes to fix CRAN issues

## clarabel 0.10.0

CRAN release: 2025-02-18

- Synced up to v0.10.0 of `Clarabel.rs`
- Fixed up to use new `r-src` that detects R being used correctly (thank
  you, Ivan Krylov!)

## clarabel 0.9.0.1

CRAN release: 2024-09-03

- Explicit system requirements for Rust and Cargo added.

## clarabel 0.9.0

CRAN release: 2024-06-22

- Synced up to version 0.9.0 of Clarabel.rs
- Added all applicable tests from Clarabel.rs
- Updated documentation of cone specification
- Updated vignette
- Switched to savvy from rextendr

## clarabel 0.5.1

CRAN release: 2023-06-24

- Clarabel now supports semidefinite programs (syncing up to version
  0.5.1 of Clarabel.rs)
- Added tests
- Updated documentation of cone specification

## clarabel 0.4.1

CRAN release: 2023-04-25

- First public release
