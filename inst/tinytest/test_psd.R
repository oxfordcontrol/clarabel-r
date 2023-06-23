## Semidefinite cone example

## Basic example

##  problem will be 3x3, so upper triangle
##  of problem data has 6 entries

P <- as(Matrix::Diagonal(6), "symmetricMatrix")
A <- diag(6)
c <- numeric(6)
b <- c(-3., 1., 4., 1., 2., 5.)
cones <- list(s = 3L)
sol <- clarabel(P = P, A = A, b = b, q = c, cone = cones)
expect_equal(sol$status, 2L)
expect_equal(sol$x,
             c(-3.0729833267361095, 0.3696004167288786, -0.022226685581313674,
               0.31441213129613066, -0.026739700851545107, -0.016084530571308823),
             tolerance = 1e-7)

## Example from scs

#' Return an vectorization of symmetric matrix using the upper triangular part,
#' still in column order.
#' @param S a symmetric matrix
#' @return vector of values
vec <- function(S) {
  n <- nrow(S)
  sqrt2 <- sqrt(2.0)
  upper_tri <- upper.tri(S, diag = FALSE)
  S[upper_tri] <- S[upper_tri] * sqrt2
  S[upper.tri(S, diag = TRUE)]
}

#' Return the symmetric matrix from the [vec] vectorization
#' @param v a vector
#' @return a symmetric matrix
mat <- function(v) {
  n <- (sqrt(8 * length(v) + 1) - 1) / 2
  sqrt2 <- sqrt(2.0)
  S <- matrix(0, n, n)
  upper_tri <- upper.tri(S, diag = TRUE)
  S[upper_tri] <- v / sqrt2
  S <- S + t(S)
  diag(S) <- diag(S) / sqrt(2)
  S
}

q <- c(1, -1, 1) # objective: x_1 - x2 + x_3
A11 <- matrix(c(-7, -11, -11, 3), nrow = 2)
A12 <- matrix(c(7, -18, -18, 8), nrow = 2)
A13 <- matrix(c(-2, -8, -8, 1), nrow = 2)

A21 <- matrix(c(-21, -11, 0, -11, 10, 8, 0, 8, 5), nrow = 3)
A22 <- matrix(c(0, 10, 16, 10, -10, -10, 16, -10, 3), nrow = 3)
A23 <- matrix(c(-5, 2, -17, 2, -6, 8, -17, 8, 6), nrow = 3)

B1 <- matrix(c(33, -9, -9, 26), nrow = 2)
B2 <- matrix(c(14, 9, 40, 9, 91, 10, 40, 10, 15), nrow = 3)

A <- rbind(
  cbind(vec(A11), vec(A12), vec(A13)), # first psd constraint
  cbind(vec(A21), vec(A22), vec(A23))  # second psd constraint
)
b <- c(vec(B1), vec(B2)) # stack both psd constraints
cones <- list(s = c(2, 3)) # cone dimensions
sol <- clarabel(A = A, b = b, q = q, cone = cones)

expect_equal(sol$status, 2L)
expect_equal(sol$x,
             c(-0.367751526880478, 1.89833324679172, -0.887460228128519),
             tolerance = 1e-7)

