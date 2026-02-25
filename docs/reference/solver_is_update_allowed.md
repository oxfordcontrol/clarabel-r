# Check if data updates are allowed on a persistent solver

Returns `FALSE` if presolve, chordal decomposition, or
`input_sparse_dropzeros` is enabled, which prevents data updates.

## Usage

``` r
solver_is_update_allowed(solver)
```

## Arguments

- solver:

  a `ClarabelSolver` object created by
  [`clarabel_solver()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md)

## Value

logical scalar

## See also

[`clarabel_solver()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md),
[`solver_update()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_update.md)

## Examples

``` r
if (FALSE) { # \dontrun{
solver_is_update_allowed(s)  # TRUE if presolve and chordal decomp are off
} # }
```
