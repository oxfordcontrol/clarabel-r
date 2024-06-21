
## Simple QP test (old remnant)

P <- Matrix::Matrix(2 * c(3, 0, 0, 2), nrow = 2, ncol = 2, sparse = TRUE)
P <- as(P, "symmetricMatrix")  # P needs to be a symmetric matrix
q <- c(-1, -4)
A <- Matrix::Matrix(c(1, 1, 0, -1, 0, -2, 0, 1, 0, -1), ncol = 2, sparse = TRUE)
b <- c(0, 1, 1, 1, 1)
cones <- list(z = 1L, l = 4L)  ## 1 equality and 4 inequalities, in order
sol <- clarabel(A = A, b = b, q = q, P = P, cones = cones)
expect_equal(sol$status, 2L)
expect_equal(sol$x,
             c(0.428571428198843, 0.214285714099422),
             tolerance = 1e-7)

basic_qp_data <- function() {
  P <- matrix(c(4, 1, 1, 2), byrow = TRUE, ncol = 2)
  A <- matrix(c(1., 1, 1, 0, 0, 1), byrow = TRUE, ncol = 2)
  A <- rbind(-A, A)
  q <- c(1., 1.)
  b <- c(-1., 0., 0., 1., 0.7, 0.7)
  cones = list(l1 = 3L, l2 = 3L)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

basic_qp_data_dual_inf <- function() {
  P <- matrix(1, nrow = 2, ncol = 2)
  A <- matrix(c(1., 1, 1, 0), byrow = TRUE, ncol = 2)
  q <- c(1., -1.)
  b <- c(1., 1.)
  cones = list(l = 2L)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

## test_qp_univariate NEEDS FIXING because of univariate matrices!! The S3 methods for
## coercion don't work.
## P <- as.matrix(1.0)
## q <- 0.0
## A <- P
## b <- 1.0
## cones <- list(l = 1L)
## solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones)
## expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
## expect_true(abs(solution$x) < 1e-6)
## expect_true(abs(solution$obj_val) < 1e-6)
## expect_true(abs(solution$obj_val_dual) < 1e-6)

## test_qp_feasible
problem_data <- basic_qp_data();
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
refsol <- c(0.3, 0.7)
expect_true(l2_dist(solution$x, refsol) <= 1e-6)
refobj <- 1.8800000298331538;
expect_true(abs(solution$obj_val - refobj) <= 1e-6)
expect_true(abs(solution$obj_val_dual- refobj) <= 1e-6)

## test_qp_primal_infeasible
problem_data <- basic_qp_data();
problem_data$b[1] <- -1; problem_data$b[4] <- -1;
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["PrimalInfeasible"]])
expect_true(is.nan(solution$obj_val))
expect_true(is.nan(solution$obj_val_dual))

## test_qp_dual_infeasible
problem_data <- basic_qp_data_dual_inf()
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["DualInfeasible"]])
expect_true(is.nan(solution$obj_val))
expect_true(is.nan(solution$obj_val_dual))

## test_qp_dual_infeasible_ill_cond
problem_data <- basic_qp_data_dual_inf()
problem_data$A <- Matrix::sparseMatrix(i = rep(0, 2), p = c(0, 1, 2), x = rep(1.0, 2), index1 = FALSE, dims = c(1, 2))
problem_data$cones <- list(l = 1L)
problem_data$b <- 1.0
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["DualInfeasible"]])
expect_true(is.nan(solution$obj_val))
expect_true(is.nan(solution$obj_val_dual))

