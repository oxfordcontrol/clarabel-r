## Simple QP test

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


