## Benchmark: chordal decomposition on / off for a banded (chordally sparse) SDP.
##
## Problem (a linear matrix inequality):
##
##   minimize    q'x
##   subject to  F0 + sum_i x_i Fi  >=  0   (positive semidefinite)
##
## All Fi share a banded sparsity pattern of bandwidth w, so the aggregate
## sparsity pattern is banded and therefore chordal, with cliques of size w+1.
## Chordal decomposition should replace one n x n PSD cone with ~n small ones.
##
## In Clarabel form:  A x + s = b,  s in PSD cone of side n
##   A = -[svec(F1) ... svec(Fm)],  b = svec(F0)

suppressMessages(library(Matrix))

SQRT2 <- sqrt(2)

## Clarabel svec: upper triangle, column-major, off-diagonal scaled by sqrt(2).
svec <- function(M) {
  n <- nrow(M)
  idx <- which(upper.tri(M, diag = TRUE), arr.ind = TRUE)
  idx <- idx[order(idx[, "col"], idx[, "row"]), , drop = FALSE]
  scale <- ifelse(idx[, "row"] == idx[, "col"], 1, SQRT2)
  M[idx] * scale
}

## Random symmetric matrix supported on a band of half-bandwidth w.
band_sym <- function(n, w, seed) {
  set.seed(seed)
  M <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in i:min(i + w, n)) {
      v <- rnorm(1)
      M[i, j] <- v
      M[j, i] <- v
    }
  }
  M
}

make_problem <- function(n, w, m, seed = 1L) {
  Fs <- lapply(seq_len(m), function(k) band_sym(n, w, seed * 1000L + k))
  ## F0 = alpha*I + band, with alpha large enough that x = 0 is strictly feasible
  F0 <- band_sym(n, w, seed * 1000L + m + 1L)
  F0 <- F0 + diag(n) * (abs(min(eigen(F0, only.values = TRUE)$values)) + 1)
  A <- -do.call(cbind, lapply(Fs, function(Fi) svec(Fi)))
  b <- svec(F0)
  ## objective: push against the constraint
  q <- vapply(Fs, function(Fi) sum(diag(Fi)), numeric(1))
  list(A = Matrix::Matrix(A, sparse = TRUE), b = b, q = q,
       cones = list(s = as.integer(n)))
}

run_one <- function(prob, chordal) {
  ctrl <- clarabel::clarabel_control(
    verbose = FALSE,
    chordal_decomposition_enable = chordal,
    chordal_decomposition_merge_method = if (chordal) "clique_graph" else "none",
    chordal_decomposition_compact = chordal,
    chordal_decomposition_complete_dual = chordal
  )
  t <- system.time(
    sol <- clarabel::clarabel(A = prob$A, b = prob$b, q = prob$q,
                              cones = prob$cones, control = ctrl)
  )
  list(elapsed = unname(t[["elapsed"]]), sol = sol)
}

cat(sprintf("clarabel %s\n\n", packageVersion("clarabel")))
cat(sprintf("%5s %4s %6s %10s %10s %7s %6s %6s %12s %12s\n",
            "n", "w", "psdlen", "t_off(s)", "t_on(s)", "speedup",
            "it_off", "it_on", "obj_reldiff", "z_maxdiff"))

for (n in c(20L, 40L, 60L, 80L, 120L)) {
  prob <- make_problem(n, w = 2L, m = 5L, seed = 7L)
  off <- run_one(prob, FALSE)
  on  <- run_one(prob, TRUE)
  objd <- abs(off$sol$obj_val - on$sol$obj_val) / max(1, abs(off$sol$obj_val))
  zd <- max(abs(off$sol$z - on$sol$z))
  cat(sprintf("%5d %4d %6d %10.3f %10.3f %7.2f %6d %6d %12.2e %12.2e\n",
              n, 2L, length(prob$b), off$elapsed, on$elapsed,
              off$elapsed / on$elapsed,
              off$sol$iterations, on$sol$iterations, objd, zd))
}
