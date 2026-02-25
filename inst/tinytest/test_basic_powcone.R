
#source("inst/tinytest/preamble.R")

if (at_home()) {

## test_powcone
## solve the following power cone problem
## max  x1^0.6 y^0.4 + x2^0.1
## s.t. x1, y, x2 >= 0
##      x1 + 2y  + 3x2 == 3
## which is equivalent to
## max z1 + z2
## s.t. (x1, y, z1) in K_pow(0.6)
##      (x2, 1, z2) in K_pow(0.1)
##      x1 + 2y + 3x2 == 3

## x = (x1, y, z1, x2, y2, z2)
n <- 6;
P <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
q <- c(0., 0., -1., 0., 0., -1.)

## (x1, y, z1) in K_pow(0.6)
## (x2, y2, z2) in K_pow(0.1)
A1 <- -diag(n)
b1 <- numeric(n)
cones1 <- list(p1 = 0.6, p2 = 0.1)

## x1 + 2y + 3x2 == 3
## y2 == 1
A2 <- Matrix::sparseMatrix(i = c(rep(1, 3), 2), j = c(1, 2, 4, 5), x = c(1.0, 2.0, 3.0, 1.0),
                           dims = c(2, 6))
b2 <- c(3., 1.)
cones2 <- list(z = 2L)
A <- rbind(A1, A2)
b <- c(b1, b2)
cones <- c(cones1, cones2)
solution <- clarabel(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
refobj <- -1.8458;
expect_true(abs(solution$info$cost_primal - refobj) <= 1e-3)

} ## end at_home()
