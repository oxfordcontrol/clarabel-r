# Create a persistent Clarabel solver object

Creates a persistent solver that can be reused across multiple solves
with updated problem data (warm starts). This avoids the overhead of
reallocating the solver's internal data structures when only the problem
data changes but the sparsity pattern stays the same.

## Usage

``` r
clarabel_solver(
  A,
  b,
  q,
  P = NULL,
  cones,
  control = list(),
  strict_cone_order = TRUE
)
```

## Arguments

- A:

  a matrix of constraint coefficients.

- b:

  a numeric vector giving the primal constraints

- q:

  a numeric vector giving the primal objective

- P:

  a symmetric positive semidefinite matrix, default `NULL`

- cones:

  a named list giving the cone sizes, see “Cone Parameters” below for
  specification

- control:

  a list giving specific control parameters to use in place of default
  values, with an empty list indicating the default control parameters.
  Specified parameters should be correctly named and typed to avoid Rust
  system panics as no sanitization is done for efficiency reasons

- strict_cone_order:

  a logical flag, default `TRUE` for forcing order of cones described
  below. If `FALSE` cones can be specified in any order and even
  repeated and directly passed to the solver without type and length
  checks

## Value

a `ClarabelSolver` environment object with methods
[`solve()`](https://rdrr.io/r/base/solve.html),
`update_data(Px, Ax, q, b)`, and `is_update_allowed()`

## Details

For data updates to work, the solver settings must have
`presolve_enable = FALSE`, `chordal_decomposition_enable = FALSE`, and
`input_sparse_dropzeros = FALSE`. Use
[`solver_is_update_allowed()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_is_update_allowed.md)
to check after construction.

## See also

[`solver_solve()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_solve.md),
[`solver_update()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_update.md),
[`solver_is_update_allowed()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_is_update_allowed.md),
[`clarabel()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
P <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = c(2, 1), dims = c(2, 2))
A <- matrix(c(1, 0, 0, 1), nrow = 2)
b <- c(1, 1)
q <- c(-2, -3)
cones <- list(l = 2L)
ctrl <- clarabel_control(presolve_enable = FALSE, verbose = FALSE)
s <- clarabel_solver(A, b, q, P, cones, control = ctrl)
sol1 <- solver_solve(s)
solver_update(s, q = c(-4, -1))
sol2 <- solver_solve(s)
} # }
```
