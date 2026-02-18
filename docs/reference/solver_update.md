# Update problem data on a persistent Clarabel solver

Update one or more of P (objective), q (linear objective), A
(constraints), b (constraint RHS) on an existing solver. The sparsity
pattern of P and A must remain the same as the original problem; only
the nonzero values can change.

## Usage

``` r
solver_update(solver, P = NULL, q = NULL, A = NULL, b = NULL)
```

## Arguments

- solver:

  a `ClarabelSolver` object created by
  [`clarabel_solver()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md)

- P:

  new upper-triangular P matrix (same sparsity), or `NULL` to leave
  unchanged

- q:

  new linear objective vector, or `NULL` to leave unchanged

- A:

  new constraint matrix (same sparsity), or `NULL` to leave unchanged

- b:

  new constraint RHS vector, or `NULL` to leave unchanged

## Value

invisible `NULL`

## See also

[`clarabel_solver()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md),
[`solver_solve()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_solve.md)

## Examples

``` r
if (FALSE) { # \dontrun{
solver_update(s, q = c(-4, -1))   # update linear objective only
solver_update(s, b = c(2, 2))     # update constraint RHS only
sol2 <- solver_solve(s)           # re-solve with updated data
} # }
```
