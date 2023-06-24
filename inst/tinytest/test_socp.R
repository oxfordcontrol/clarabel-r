## Basic SOCP

P <- Matrix::Matrix(2 * c(0, 0, 0, 1), nrow = 2, ncol = 2, sparse = TRUE)
P <- as(P, "symmetricMatrix") # P needs to be a symmetric matrix
q <- c(0, 0)
A <- Matrix::Matrix(c(0, -2.0, 0, 0, 0, 1.0), nrow = 3, ncol = 2, sparse = TRUE)
b <- c(1, -2, -2)
cones <- list(q = 3L)
sol <- clarabel(A = A, b = b, q = q, P = P, cones = cones)
expect_equal(sol$status, 2L)
expect_equal(sol$x,
             c(1, -1),
             tolerance = 1e-7)
