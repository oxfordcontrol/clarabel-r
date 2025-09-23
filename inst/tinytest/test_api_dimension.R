
api_dim_check_data <- function() {
  P <- matrix(0.0, nrow = 4, ncol = 4)
  q <- numeric(4)
  A <- matrix(0.0, nrow = 6, ncol = 4)
  b <- numeric(6)
  cones = list(z = 1L, l = 5L)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

## api_dim_check_working()
## This example should work because dimensions are
## all compatible.  All following checks vary one
## of these sizes to test dimension checks
problem_data <- api_dim_check_data()
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])

## api_dim_check_bad_P() {
## Check that bad P dimension is detected
problem_data$P <- matrix(0.0, nrow = 3, ncol = 3)
expect_error(solution <- do.call(clarabel, problem_data))

## api_dim_check_bad_A_rows()
## Check that bad A rows are detected
problem_data <- api_dim_check_data()
problem_data$A <- matrix(0.0, nrow = 5, ncol = 4)
expect_error(solution <- do.call(clarabel, problem_data))

## api_dim_check_bad_A_cols()
## Check that bad A columns are detected
problem_data <- api_dim_check_data()
problem_data$A <- matrix(0.0, nrow = 6, ncol = 3)
expect_error(solution <- do.call(clarabel, problem_data))

## api_dim_check_P_not_square()
## Check that P not square is detected
problem_data <- api_dim_check_data()
problem_data$P <- matrix(0.0, nrow = 5, ncol = 4)
expect_error(solution <- do.call(clarabel, problem_data))

## Check that repeated specification of cone types with strict_cone_order TRUE fails
problem_data <- api_dim_check_data()
problem_data$cones <- list(z = 1L, l = 2L, l = 4L)
expect_error(solution <- do.call(clarabel, problem_data))

## Check that the cone dimensions mismatch with matrices is detected
problem_data <- api_dim_check_data()
problem_data$cones <- list(z = 1L, l = 6L)
expect_error(solution <- do.call(clarabel, problem_data))

## Check that repeated specification of cone types succeeds with strict_cone_order FALSE
problem_data <- api_dim_check_data()
problem_data$cones <- list(z = 1L, l = 2L, l = 3L)
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])

## api_dim_check_bad_cones()
## Check that bad dimensions for cones fails
problem_data <- api_dim_check_data()
problem_data$cones <- list(z = 1L, l = 2L, l = 4L)
expect_error(solution <- do.call(clarabel, problem_data))

