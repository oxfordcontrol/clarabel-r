# Check if data updates are allowed on a persistent solver

Returns `FALSE` if presolve, chordal decomposition or
`input_sparse_dropzeros` has actually restructured this problem, which
prevents data updates. Enabling one of those settings is not by itself
enough to prevent updates; the transformation has to have taken place.

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
solver_is_update_allowed(s)  # FALSE if the problem was restructured
} # }
```
