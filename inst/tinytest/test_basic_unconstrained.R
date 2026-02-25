if (at_home()) {

## test_unconstrained_feasible
P <- diag(3)
q <- c(1.0, 2.0, -3.0)
A <- matrix(0, nrow = 0, ncol = 3) ## <- no constraints
b <- numeric(0)
cones <- list()
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
refsol <- -q
expect_true(l2_dist(solution$x, refsol) <= 1e-6)

## test_unconstrained_dual_infeasible
P <- matrix(0, nrow = 3, ncol = 3)
q <- c(1.0, 0, 0)
A <- matrix(0, nrow = 0, ncol = 3) ## <- no constraints
b <- numeric(0)
cones <- list()
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
expect_equal(status_codes[[solution$status]], status_codes[["DualInfeasible"]])

} ## end at_home()
