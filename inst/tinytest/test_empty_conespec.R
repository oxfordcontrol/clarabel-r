if (at_home()) {

##                                                                                                                                      ## Minimal reproducer: R clarabel rejects empty cone spec.
##                                                                                                                                      ## Problem: minimize -x1 - x2  (no constraints)           
## This is unbounded. Python clarabel returns DualInfeasible.                                                                           ## R clarabel crashes with:                                                                                                             ##   "sanitize_cone_spec: no cone parameters specified"
##
## The Rust solver itself handles this fine (as shown by Python).
## The bug is in the R wrapper's `sanitize_cone_spec()` function,
## which rejects an empty list before passing it to the Rust solver.
##
## Expected behavior: solve successfully, return status = DualInfeasible
## (integer code for dual infeasibility, i.e. unbounded primal).
##

n <- 2L

P <- NULL                        # no quadratic term (LP)
q <- c(-1, -1)                   # minimize -x1 - x2

## Zero constraints: A is 0-by-n, b is empty
A <- Matrix::sparseMatrix(
  i = integer(0), j = integer(0), x = numeric(0),
  dims = c(0L, n)
)
b <- numeric(0)

cones <- list()                  # empty cone spec — triggers the bug

## This used to crash.
##   Error: sanitize_cone_spec: no cone parameters specified
sol <- clarabel(A = A, b = b, q = q, P = P, cones = cones)
expect_equal(status_codes[[sol$status]], status_codes[["DualInfeasible"]])

} ## end at_home()
