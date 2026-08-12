## Sweep for the vignette table: same construction the vignette chunk uses
## (dense A, seed 7), so the live chunk and the table are the same experiment.
suppressMessages(library(clarabel))

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
  list(t = unname(t[["elapsed"]]), st = s$status, obj = s$obj_val, z = s$z, it = s$iterations)
}

nm <- names(clarabel::solver_status_descriptions())
cat(sprintf("%4s %8s %10s %10s %9s %-14s %-14s %11s %10s\n",
            "n", "psdlen", "off(s)", "on(s)", "speedup",
            "status_off", "status_on", "obj_absdiff", "max_dz"))
for (n in c(20L, 40L, 60L, 80L, 120L)) {
  p <- banded_sdp(n)
  o <- go(p, FALSE)
  c2 <- go(p, TRUE)
  cat(sprintf("%4d %8d %10.3f %10.3f %8.0fx %-14s %-14s %11.2e %10.2e\n",
              n, length(p$b), o$t, c2$t, o$t / c2$t,
              nm[o$st], nm[c2$st], abs(o$obj - c2$obj), max(abs(o$z - c2$z))))
}
