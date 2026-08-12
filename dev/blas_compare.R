## Solve the banded SDP family under whichever BLAS R is currently linked to,
## saving timings AND solutions so accuracy can be compared across BLAS choices.
## Usage: Rscript blas_compare.R <output.rds>
suppressMessages(library(clarabel))

out_path <- commandArgs(trailingOnly = TRUE)[1]

vec <- function(S) {
  u <- upper.tri(S, diag = FALSE)
  S[u] <- S[u] * sqrt(2)
  S[upper.tri(S, diag = TRUE)]
}

banded_sdp <- function(n, w = 2, m = 5, seed = 7) {
  band_sym <- function(s) {
    set.seed(s)
    M <- matrix(0, n, n)
    for (i in seq_len(n)) for (j in i:min(i + w, n)) M[i, j] <- M[j, i] <- rnorm(1)
    M
  }
  Fs <- lapply(seq_len(m), function(k) band_sym(1000 * seed + k))
  F0 <- band_sym(1000 * seed + m + 1)
  F0 <- F0 + diag(n) * (abs(min(eigen(F0, only.values = TRUE)$values)) + 1)
  list(A = -do.call(cbind, lapply(Fs, vec)),
       b = vec(F0),
       q = vapply(Fs, function(Fi) sum(diag(Fi)), numeric(1)),
       cones = list(s = as.integer(n)))
}

go <- function(p, ch) {
  ctrl <- clarabel_control(verbose = FALSE, chordal_decomposition_enable = ch)
  t <- system.time(s <- clarabel(A = p$A, b = p$b, q = p$q, cones = p$cones, control = ctrl))
  list(t = unname(t[["elapsed"]]), status = s$status, obj = s$obj_val,
       x = s$x, z = s$z, iters = s$iterations)
}

blas <- extSoftVersion()[["BLAS"]]
cat("BLAS:", blas, "\n")

res <- list()
for (n in c(20L, 40L, 60L, 80L, 120L)) {
  p <- banded_sdp(n)
  res[[as.character(n)]] <- list(off = go(p, FALSE), on = go(p, TRUE))
  r <- res[[as.character(n)]]
  cat(sprintf("n=%3d  off %8.3fs (st %d)  on %7.3fs (st %d)\n",
              n, r$off$t, r$off$status, r$on$t, r$on$status))
}

saveRDS(list(blas = blas, res = res), out_path)
cat("saved:", out_path, "\n")
