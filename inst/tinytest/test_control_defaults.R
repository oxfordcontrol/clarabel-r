if (at_home()) {

## clarabel_control() hard-codes its defaults in the function signature so that
## they appear in the help page. That copy can silently fall out of step with
## the Rust solver: reduced_tol_infeas_abs was changed upstream in Clarabel.rs
## v0.10.0 and this package went two releases without noticing. These tests
## compare the R defaults against DefaultSettings::<f64>::default() field by
## field, so the next upstream change fails here instead of drifting.

rust <- clarabel:::rust_default_settings()
ctrl <- clarabel_control()

## Every setting the Rust solver knows about should be offered by
## clarabel_control(), and vice versa.
expect_equal(sort(names(ctrl)), sort(names(rust)))

## direct_solve_method is the one legitimate departure: the accepted values
## depend on which Cargo features this package is built with. "faer" needs
## faer-sparse and "mkl"/"panua" need a Pardiso feature, none of which are
## enabled here, so only "auto" and "qdldl" are selectable. The *default* must
## still agree with upstream.
expect_equal(ctrl$direct_solve_method, rust$direct_solve_method)
expect_equal(formals(clarabel_control)$direct_solve_method,
             quote(c("auto", "qdldl")))

## Everything else must match upstream exactly. Integer-vs-double storage
## differs by design (R's integer settings are passed as integers), so compare
## values rather than types.
for (nm in names(rust)) {
  expect_equal(unname(unlist(ctrl[[nm]])), unname(unlist(rust[[nm]])),
               info = sprintf("clarabel_control() default for '%s' disagrees with Clarabel.rs", nm))
}

## Spot-check the specific values this release corrected, so a future edit that
## silently reverts one of them is caught by name rather than only in bulk.
expect_equal(ctrl$reduced_tol_infeas_abs, 5e-12)
expect_equal(ctrl$reduced_tol_gap_abs, 5e-5)
expect_equal(ctrl$tol_ktratio, 1e-6)
expect_equal(ctrl$iterative_refinement_reltol, 1e-13)
expect_true(ctrl$chordal_decomposition_enable)
expect_equal(ctrl$chordal_decomposition_merge_method, "clique_graph")

## Invalid solver methods must be rejected by match.arg() rather than reaching
## the Rust side, where they would fail at solver construction.
expect_error(clarabel_control(direct_solve_method = "cholmod"))
expect_error(clarabel_control(direct_solve_method = "mkl"))

} ## end at_home()
