eq_constrained_A1 <- function() {
  ## A =
  ##[ 0. 1.  1.;
  ##  0. 1. -1.]
  Matrix::sparseMatrix(
    dims = c(2, 3),  ## m x n
    p = c(0, 0, 2, 4),      ##colptr
    i = c(0, 1, 0, 1),      ##rowval
    x = c(1., 1., 1., -1.), ##nzva;
    index1 = FALSE
  )
}

eq_constrained_A2 <- function() {
  ## A = [
  ## 0    1.0   1.0;
  ## 0    1.0  -1.0;
  ##1.0   2.0  -1.0l
  ##2.0  -1.0   3.0l
  ##]
  Matrix::sparseMatrix(
    dims = c(4, 3),                                   ## m x n
    p = c(0, 2, 6, 10),                               ##colptr
    i = c(2, 3, 0, 1, 2, 3, 0, 1, 2, 3),              ##rowval
    x = c(1., 2., 1., 1., 2., -1., 1., -1., -1., 3.), ##nzva;
    index1 = FALSE
    )
}

## test_eq_constrained_feasible
P <- diag(3)
q <- numeric(3)
A <- eq_constrained_A1() ## <- two constraints
b <- c(2., 0.)
cones <- list(z = 2)
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones)
refsol <- c(0., 1., 1.)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
expect_true(l2_dist(solution$x, refsol) <= 1e-6)


## test_eq_constrained_primal_infeasible
P <- diag(3)
q <- numeric(3)
A <- eq_constrained_A2() ## <- 4 constraints, 3 vars
b <- rep(1.0, 4)
cones <- list(z = 4)
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones)
expect_equal(status_codes[[solution$status]], status_codes[["PrimalInfeasible"]])


## test_eq_constrained_dual_infeasible
P <- diag(3); P[1, 1] <- 0;
q <- rep(1.0, 3)
A <- eq_constrained_A1() ## <- two constraints
b <- c(2., 0.)
cones <- list(z = 2L)
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones)
expect_equal(status_codes[[solution$status]], status_codes[["DualInfeasible"]])
