
#source("inst/tinytest/preamble.R")

basic_lp_data <- function() {
  P <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(3, 3))
  #P <- NULL
  I1 <- diag(3)
  I2 <- -diag(3)
  A <- rbind(I1, I2)
  A <- A * 2
  q <- c(3., -2., 1.)
  b <- rep(1.0, 6)
  cones <- list(l1 = 3L, l2 = 3L)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

## test_lp_feasible
problem_data <- basic_lp_data();
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
refsol <- c(-0.5, 0.5, -0.5)
expect_true(l2_dist(solution$x, refsol) <= 1e-8)
refobj <- -3.0
expect_true(abs(solution$obj_val - refobj) <= 1e-8)
expect_true(abs(solution$obj_val_dual- refobj) <= 1e-8)


## test_lp_primal_infeasible
problem_data <- basic_lp_data();
problem_data$b[1] <- -1; problem_data$b[4] <- -1;
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["PrimalInfeasible"]])
expect_true(is.nan(solution$obj_val))
expect_true(is.nan(solution$obj_val_dual))

## test_lp_dual_infeasible
problem_data <- basic_lp_data();
problem_data$A[4, 1] <- 1. ##swap lower bound on first variable to redundant upper bound
problem_data$q <- c(1., 0., 0.)
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["DualInfeasible"]])
expect_true(is.nan(solution$obj_val))
expect_true(is.nan(solution$obj_val_dual))

## test_lp_dual_infeasible_ill_cond() {
problem_data <- basic_lp_data();
problem_data$A[1, 1] <- .Machine$double.eps
problem_data$A[4, 1] <- 0
problem_data$q <- c(1., 0., 0.)
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["DualInfeasible"]])
expect_true(is.nan(solution$obj_val))
expect_true(is.nan(solution$obj_val_dual))
