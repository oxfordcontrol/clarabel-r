if (at_home()) {

## test_mixed_conic_feasible
## solves a problem with a mix of symmetric and asymmetric
## cones.   This exercises the barrier methods and unit
## initializations of the symmetric cones

n <- 3
P <- diag(3)
q <- c(1., 1., 1.)
## put a 3 dimensional vector into the composition of multiple
## cones, all with b = 0 on the RHS
cones <- list(z = 3L,
              l = 3L,
              q = 3L,
              p = 0.5,
              ep = 1L
              )
A <- rbind(diag(3), diag(3))
A <- rbind(A, A, diag(3))
b <- numeric(5 * n)
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
expect_true(abs(solution$info$cost_primal - 0.) <= 1e-8)

} ## end at_home()

