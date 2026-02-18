# Solve using a persistent Clarabel solver

Solve using a persistent Clarabel solver

## Usage

``` r
solver_solve(solver)
```

## Arguments

- solver:

  a `ClarabelSolver` object created by
  [`clarabel_solver()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md)

## Value

the same named list as
[`clarabel()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md):
solution vectors `x`, `z`, `s` and solver information

## See also

[`clarabel_solver()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md),
[`solver_update()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_update.md)

## Examples

``` r
if (FALSE) { # \dontrun{
s <- clarabel_solver(A, b, q, P, cones,
                     control = clarabel_control(presolve_enable = FALSE,
                                                verbose = FALSE))
sol <- solver_solve(s)
sol$status
} # }
```
