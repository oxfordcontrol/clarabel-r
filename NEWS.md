# clarabel 0.11.2

- Added persistent solver API for warm starts: `clarabel_solver()`,
  `solver_solve()`, `solver_update()`, `solver_is_update_allowed()`
- Robustified Rust interface: type coercion, error handling,
  CscMatrix validation, and regex ordering fixes
- Switched all R-side error messaging to `cli::cli_abort()` with
  markup
- Upgraded `savvy` crate from 0.8.13 to 0.9.2
- Added examples to `clarabel_control()` and solver functions
- Added vignette section on updating problem data (warm starts)

# clarabel 0.11.1

- Synced up to v0.11.1 of `Clarabel.rs`
- Fixed `flang` issues shown on CRAN for fedora-clang-rdevel

# clarabel 0.10.1

- Minor fixes to fix CRAN issues

# clarabel 0.10.0

- Synced up to v0.10.0 of `Clarabel.rs`
- Fixed up to use new `r-src` that detects R being used correctly (thank you, Ivan Krylov!)

# clarabel 0.9.0.1

- Explicit system requirements for Rust and Cargo added.

# clarabel 0.9.0

- Synced up to version 0.9.0 of Clarabel.rs
- Added all applicable tests from Clarabel.rs
- Updated documentation of cone specification
- Updated vignette
- Switched to savvy from rextendr

# clarabel 0.5.1

- Clarabel now supports semidefinite programs (syncing up to version 0.5.1 of Clarabel.rs)
- Added tests
- Updated documentation of cone specification

# clarabel 0.4.1

- First public release
