# Return the solver status description as a named character vector

Return the solver status description as a named character vector

## Usage

``` r
solver_status_descriptions()
```

## Value

a named list of solver status descriptions, in order of status codes
returned by the solver

## Examples

``` r
solver_status_descriptions()[2] ## for solved problem
#>                               Solved 
#> "Solver terminated with a solution." 
solver_status_descriptions()[8] ## for max iterations limit reached
#>                                                                MaxIterations 
#> "Iteration limit reached before solution or infeasibility certificate found" 
```
