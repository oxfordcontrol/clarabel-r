SQRT_2 <- sqrt(2.0)

sdp_chordal_data <- function() {
  P <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(8, 8))
  q <- c(-1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0)
  A <- Matrix::sparseMatrix(
    dims = c(28, 8),
    p = c(0, 1, 4, 5, 8, 9, 10, 13, 16),
    i = c(24, 7, 10, 22, 8, 12, 15, 25, 9, 13, 18, 21, 26, 0, 23, 27),
    x = c(
            -1.0,
            -SQRT_2,
            -1.0,
            -1.0,
            -SQRT_2,
            -SQRT_2,
            -1.0,
            -1.0,
            -SQRT_2,
            -SQRT_2,
            -SQRT_2,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0
      ),
    index1 = FALSE
  )
  b <- c(
    0.0,
    3.0,
    2. * SQRT_2,
    2.0,
    SQRT_2,
    SQRT_2,
    3.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  )
  cones <- list(l = 1L, s = 6L, p1 = 1.0 / 3.0, p2 = 0.5)
  list(P = P, q = q, A = A, b = b, cones = cones, strict_cone_order = FALSE)
}

##test_sdp_chordal() {
problem_data <- sdp_chordal_data()
problem_data$control <- clarabel_control(
  chordal_decomposition_enable = TRUE,
  chordal_decomposition_compact = TRUE,
  chordal_decomposition_complete_dual = TRUE,
  chordal_decomposition_merge_method = "clique_graph",
  max_iter = 50L)

for (compact in c(FALSE, TRUE)) {
  for (complete_dual in c(FALSE, TRUE)) {
    for (merge_method in c("clique_graph", "parent_child", "none")) {
      problem_data$control <- clarabel_control(
        chordal_decomposition_enable = TRUE,
        chordal_decomposition_compact = compact,
        chordal_decomposition_complete_dual = complete_dual,
        chordal_decomposition_merge_method = merge_method,
        max_iter = 50L)
      solution <- do.call(clarabel, problem_data)
      print(expect_equal(status_codes[[solution$status]], status_codes[["Solved"]]))
    }
  }
}

