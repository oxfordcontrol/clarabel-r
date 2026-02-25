if (at_home()) {

## This is not the complete suite since we don't expose the solver object!
## TBD
presolve_test_data <- function() {
  n <- 3
  P <- diag(n)
  A <- 2 * rbind(P, -diag(n))
  q <- c(3., -2., 1.)
  b <- rep(1.0, 2 * n)
  cones <- list(l1 = 3L, l2 = 3L)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

## test_presolve_single_unbounded
problem_data <- presolve_test_data()
problem_data$b[4] <- 1e30
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
expect_equal(solution$z[4], 0.)

## test_presolve_completely_redundant_cone
problem_data <- presolve_test_data()
problem_data$b[1:3] <- 1e30
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
expect_equal(solution$z[1:3], numeric(3))
refsol <- c(-0.5, 2., -0.5)
expect_true(l2_dist(solution$x, refsol) <= 1e-6)

##test_presolve_every_constraint_redundant
problem_data <- presolve_test_data()
problem_data$b <- rep(1e30, length(problem_data$b))
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
expect_true(l2_dist(solution$x, -problem_data$q) <= 1e-6)

} ## end at_home()

