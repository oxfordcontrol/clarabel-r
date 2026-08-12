# clarabel 0.11.3

- `clarabel_control()` defaults now match `Clarabel.rs` exactly. Several
  had fallen behind upstream:
  - Chordal decomposition is now **enabled by default** for semidefinite
    programs: `chordal_decomposition_enable`,
    `chordal_decomposition_compact` and
    `chordal_decomposition_complete_dual` are `TRUE`, and
    `chordal_decomposition_merge_method` is `"clique_graph"`. On a
    semidefinite program with a decomposable sparsity pattern this is
    dramatically faster — on a banded example with a 120x120 semidefinite
    cone, 0.007 seconds against 422 seconds. Two consequences: the dual
    variables reported for semidefinite constraints are not in general the
    same as before, and data updates are blocked when a decomposition
    actually takes place. Linear, quadratic and second-order cone problems
    are unaffected and remain updatable. Pass
    `chordal_decomposition_enable = FALSE` for the previous behavior.
  - `reduced_tol_infeas_abs` is now `5e-12`, not `5e-5`, matching upstream
    since `Clarabel.rs` v0.10.0. Smaller means *less* margin required, so
    the previous value was stricter than upstream; the change makes
    `AlmostPrimalInfeasible` and `AlmostDualInfeasible` marginally easier
    to reach on borderline problems.
  - `direct_solve_method` is now `"auto"`, not `"qdldl"`, matching upstream
    since v0.11.0. This build links only the QDLDL solver, so `"auto"`
    resolves to the same code path and results are unchanged.
- Removed two `direct_solve_method` options that never worked:
  `"cholmod"` (never existed in `Clarabel.rs`; it is a `Clarabel.jl`
  option) and `"mkl"` (requires a Pardiso feature this build does not
  enable). The accepted values are `"auto"` and `"qdldl"`.
- Corrected the documented defaults for `tol_ktratio` (`1e-6`, was
  documented as `1e-7`) and `iterative_refinement_reltol` (`1e-13`, was
  documented as `1e-12`). The code was already correct; only the help page
  was wrong.

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
- Reduced CRAN test suite to LP and QP only; full suite runs locally
  via `tinytest::at_home()`
- Fixed macOS CRAN linker warning by forwarding `MACOSX_DEPLOYMENT_TARGET`
  to Rust `cc` crate during `cargo build`

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
