if (at_home()) {

basic_expcone_data <- function() {
  ## produces data for the following exponential cone problem
  ## max  x
  ## s.t. y * exp(x / y) <= z
  ##      y == 1, z == exp(5)
  
  P <- matrix(0, nrow = 3, ncol = 3)
  q <- c(-1., 0., 0.)
  A1 <- -Matrix::Diagonal(3)
  b1 <- rep(0., 3)
  A2 <- Matrix::sparseMatrix(
    dims = c(2, 3),    ## m x n
    p = c(0, 0, 1, 2), ##colptr
    i = c(0, 1),       ##rowval
    x = c(1., 1.),     ##nzval
    index1 = FALSE
    )
  b2 <- c(1.0, exp(5.0))
  A <- rbind(A1, A2)
  b <- c(b1, b2)

  cones <- list(ep = 1L, z = 2L)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

## test_expcone_feasible
## solve the following exponential cone problem
## max  x
## s.t. y * exp(x / y) <= z
##      y == 1, z == exp(5)
##
## This is just the default problem data above

problem_data <- basic_expcone_data()
refsol <- c(5.0, 1., exp(5.0))
refobj <- -5.0
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
expect_true(abs(solution$info$cost_primal - refobj) <= 1e-6)
expect_true(l2_dist(solution$x, refsol) <= 1e-6)

## test_expcone_primal_infeasible
## solve the following exponential cone problem
## max  x
## s.t. y * exp(x / y) <= z
##      y == 1, z == -1
##
## Same as default, but last element of b is different
problem_data <- basic_expcone_data()
problem_data$b[5] = -1.
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["PrimalInfeasible"]])

## test_expcone_dual_infeasible
## solve the following exponential cone problem
## max  x
## s.t. y * exp(x / y) <= z
##
## Same as default, but no equality constraint
P <- NULL
q <- c(-1., 0., 0.)
A <- -diag(3)
b <- numeric(3)
cones <- list(ep = 1L)
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones)
expect_equal(status_codes[[solution$status]], status_codes[["DualInfeasible"]])

} ## end at_home()
