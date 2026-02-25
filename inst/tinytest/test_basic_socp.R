#source("inst/tinytest/preamble.R")

if (at_home()) {

basic_socp_data <- function() {
  ## P matrix data taken from corresponding Julia unit test.
  ## These nzvals form a 3x3 positive definite matrix

  P <- Matrix::sparseMatrix(
    p = c(0, 3, 6, 9),                ## colptr
    i = c(0, 1, 2, 0, 1, 2, 0, 1, 2), ## rowval
    x = c(                            ## nzval
        1.4652521089139698,
        0.6137176286085666,
        -1.1527861771130112,
        0.6137176286085666,
        2.219109946678485,
        -1.4400420548730628,
        -1.1527861771130112,
        -1.4400420548730628,
        1.6014483534926371
        ),
    index1 = FALSE,
    dims = c(3, 3)
  )
  ## A = [2I;-2I;I]
  A <- rbind(2*diag(3), -2*diag(3), diag(3))
  q <- c(0.1, -2.0, 1.0)
  b <- c(1., 1., 1., 1., 1., 1., 0., 0., 0.)
  cones <- list(l1 = 3L, l2 = 3L, q = 3L)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

## test_socp_feasible
problem_data <- basic_socp_data()
refsol <- c(-0.5, 0.435603, -0.245459)
refobj <- -8.4590e-01
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["Solved"]])
expect_true(l2_dist(solution$x, refsol) <= 1e-4)
expect_true(abs(solution$obj_val - refobj) <= 1e-4)
expect_true(abs(solution$obj_val_dual - refobj) <= 1e-4)

## test_socp_infeasible
problem_data <- basic_socp_data()
## make the cone constraint unsatisfiable
problem_data$b[7] <- -10.
solution <- do.call(clarabel, problem_data)
expect_equal(status_codes[[solution$status]], status_codes[["PrimalInfeasible"]])
expect_true(is.nan(solution$obj_val))
expect_true(is.nan(solution$obj_val_dual))


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

} ## end at_home()
