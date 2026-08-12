## Compare the two BLAS runs: speed, termination status, and how far the
## returned solutions differ.
S <- "/private/tmp/claude-503/-Users-naras-GitHub-clarabel-r/826c161b-9109-489f-b745-15dabc1b237a/scratchpad"
ref <- readRDS(file.path(S, "res_reference.rds"))
vec <- readRDS(file.path(S, "res_veclib.rds"))

st <- names(clarabel::solver_status_descriptions())
ns <- names(ref$res)

cat("reference:", basename(ref$blas), "\n")
cat("vecLib   :", basename(vec$blas), "\n\n")

cat("TIMING (seconds) and STATUS\n")
cat(sprintf("%4s | %9s %9s %7s | %9s %9s %7s\n",
            "n", "ref_off", "vec_off", "ratio", "ref_on", "vec_on", "ratio"))
for (n in ns) {
  a <- ref$res[[n]]; b <- vec$res[[n]]
  cat(sprintf("%4s | %9.3f %9.3f %6.2fx | %9.3f %9.3f %6.2fx\n",
              n, a$off$t, b$off$t, a$off$t / b$off$t,
              a$on$t, b$on$t, a$on$t / b$on$t))
}

cat("\nSTATUS (ref -> vecLib)\n")
for (n in ns) {
  a <- ref$res[[n]]; b <- vec$res[[n]]
  cat(sprintf("%4s | off: %-14s -> %-14s | on: %-14s -> %-14s\n",
              n, st[a$off$status], st[b$off$status],
              st[a$on$status], st[b$on$status]))
}

cat("\nSOLUTION DIFFERENCES between BLAS (same code, same data)\n")
cat(sprintf("%4s | %12s %12s %12s | %12s %12s %12s\n",
            "n", "off_dobj", "off_dx", "off_dz", "on_dobj", "on_dx", "on_dz"))
for (n in ns) {
  a <- ref$res[[n]]; b <- vec$res[[n]]
  cat(sprintf("%4s | %12.2e %12.2e %12.2e | %12.2e %12.2e %12.2e\n", n,
              abs(a$off$obj - b$off$obj), max(abs(a$off$x - b$off$x)), max(abs(a$off$z - b$off$z)),
              abs(a$on$obj  - b$on$obj),  max(abs(a$on$x  - b$on$x)),  max(abs(a$on$z  - b$on$z))))
}

cat("\nITERATIONS (ref / vecLib)\n")
for (n in ns) {
  a <- ref$res[[n]]; b <- vec$res[[n]]
  cat(sprintf("%4s | off: %2d / %2d | on: %2d / %2d\n",
              n, a$off$iters, b$off$iters, a$on$iters, b$on$iters))
}
