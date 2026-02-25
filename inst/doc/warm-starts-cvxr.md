# Clarabel Warm Starts for CVXR

## Overview

Starting with `clarabel` v0.11.9003, the R package exposes a
**persistent solver API** that allows CVXR to reuse the same Clarabel
solver object across multiple solves. When only the problem data
changes (but not the problem structure), this avoids reallocating and
refactoring the solver's internal KKT system, yielding significant
speedups in parametric programming workflows.

This document describes:

1. The new `clarabel` API for persistent solvers
2. How to integrate it into CVXR's `Clarabel_Solver` class
3. The warm start pattern, following the existing OSQP implementation

## 1. The `clarabel` Persistent Solver API

### Functions

| Function | Purpose |
|----------|---------|
| `clarabel_solver(A, b, q, P, cones, control, strict_cone_order)` | Create a persistent solver (same args as `clarabel()`) |
| `solver_solve(solver)` | Solve and return results (same format as `clarabel()`) |
| `solver_update(solver, P, q, A, b)` | Update problem data; `NULL` = leave unchanged |
| `solver_is_update_allowed(solver)` | Check if data updates are permitted |

### Basic usage

```r
library(clarabel)
library(Matrix)

P <- sparseMatrix(i = 1:2, j = 1:2, x = c(2, 1), dims = c(2, 2))
P <- as(P, "symmetricMatrix")
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 2)
q <- c(-2, -3)
cones <- list(l = 3L)

# Settings: presolve must be OFF for data updates
ctrl <- clarabel_control(presolve_enable = FALSE, verbose = FALSE)

# Create persistent solver
s <- clarabel_solver(A, b, q, P, cones, control = ctrl)

# First solve
sol1 <- solver_solve(s)
sol1$x  # [0.333, 1.667]

# Update linear objective and re-solve
solver_update(s, q = c(-4, -1))
sol2 <- solver_solve(s)
sol2$x  # [1.667, 0.333]
```

### Rules for data updates

- **Same sparsity pattern**: When updating `P` or `A`, pass a matrix
  with the same sparsity structure as the original. Only the nonzero
  values change; the row/column index arrays stay fixed.
- **Full vectors**: `q` and `b` are replaced entirely.
- **Selective updates**: Pass `NULL` (the default) for any component
  you don't want to change.
- **Settings constraints**: For updates to work, the solver must have
  been created with `presolve_enable = FALSE`,
  `chordal_decomposition_enable = FALSE`, and
  `input_sparse_dropzeros = FALSE`. Use
  `solver_is_update_allowed(solver)` to verify.

## 2. CVXR Integration

### 2.1 What changes in CVXR

Only one file needs modification:
`R/164_clarabel_conif.R` (or its source at
`rsrc_tree/reductions/solvers/conic_solvers/clarabel_conif.R`).

The `solve_via_data` method for `Clarabel_Solver` currently ignores
the `warm_start` parameter and creates a fresh solver on every call.
The change follows the same cache pattern used by
`OSQP_QP_Solver` in `R/168_osqp_qpif.R`.

### 2.2 Implementation

Replace the existing `solve_via_data` method with:

```r
method(solve_via_data, Clarabel_Solver) <- function(x, data, warm_start, verbose,
                                                    solver_opts, ...) {
  if (!requireNamespace("clarabel", quietly = TRUE)) {
    cli_abort("Package {.pkg clarabel} is required but not installed.")
  }

  dots <- list(...)
  solver_cache <- dots[["solver_cache"]]

  A <- data[[SD_A]]
  b <- data[[SD_B]]
  q <- data[[SD_C]]
  cones <- dims_to_solver_dict_clarabel(data[[SD_DIMS]])

  ## Build P (quadratic objective)
  nvars <- length(q)
  if (!is.null(data[[SD_P]])) {
    P <- Matrix::forceSymmetric(Matrix::triu(data[[SD_P]]), uplo = "U")
  } else {
    P <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(nvars, nvars))
  }

  ## Parse settings
  settings <- clarabel::clarabel_control()
  settings$verbose <- verbose
  for (opt_name in names(solver_opts)) {
    settings[[opt_name]] <- solver_opts[[opt_name]]
  }

  cache_key <- CLARABEL_SOLVER
  used_warm <- FALSE

  ## ── Warm path ──────────────────────────────────────────────────
  if (warm_start && !is.null(solver_cache) &&
      exists(cache_key, envir = solver_cache)) {
    cached <- get(cache_key, envir = solver_cache)
    old_solver <- cached$solver
    old_data   <- cached$data

    ## Dimension check: structure must match
    if (length(q) == length(old_data$q) &&
        nrow(A) == nrow(old_data$A) &&
        ncol(A) == ncol(old_data$A)) {

      ## Determine what changed and build update args
      new_P <- NULL
      new_q <- NULL
      new_A <- NULL
      new_b <- NULL

      if (!is.null(data[[SD_P]]) && !is.null(old_data$P)) {
        if (!identical(P@x, old_data$P@x)) new_P <- P
      }
      if (!identical(q, old_data$q)) new_q <- q
      if (!identical(A@x, old_data$A@x)) new_A <- A
      if (!identical(b, old_data$b)) new_b <- b

      ## Send incremental updates
      if (!is.null(new_P) || !is.null(new_q) ||
          !is.null(new_A) || !is.null(new_b)) {
        clarabel::solver_update(old_solver,
                                P = new_P, q = new_q,
                                A = new_A, b = new_b)
      }

      result <- clarabel::solver_solve(old_solver)
      used_warm <- TRUE
    }
  }

  ## ── Cold path ──────────────────────────────────────────────────
  if (!used_warm) {
    ## For warm-start-capable solver, disable features that block updates
    if (warm_start) {
      settings$presolve_enable <- FALSE
      settings$chordal_decomposition_enable <- FALSE
      settings$input_sparse_dropzeros <- FALSE
    }

    solver_obj <- clarabel::clarabel_solver(
      A = A, b = b, q = q, P = P, cones = cones, control = settings
    )
    result <- clarabel::solver_solve(solver_obj)
  }

  ## ── Cache for future warm starts ───────────────────────────────
  if (!is.null(solver_cache)) {
    ## Cache the solver object and data snapshot
    solver_to_cache <- if (used_warm) old_solver else solver_obj
    assign(cache_key,
           list(solver = solver_to_cache,
                data = list(P = P, q = q, A = A, b = b)),
           envir = solver_cache)
  }

  result
}
```

### 2.3 Key design decisions

**Why follow the OSQP pattern?**

The OSQP warm start in `R/168_osqp_qpif.R` is the most complete
reference implementation in CVXR. The Clarabel implementation mirrors
it:

| Aspect | OSQP | Clarabel |
|--------|------|----------|
| Cache key | `OSQP_SOLVER` | `CLARABEL_SOLVER` |
| Cached objects | `model`, `data`, `result` | `solver`, `data` |
| Dimension guard | `length(q)` + `nrow(A)` | `length(q)` + `nrow(A)` + `ncol(A)` |
| Incremental updates | `model@Update(Px=, Ax=, q=, l=, u=)` | `clarabel::solver_update(solver, P=, q=, A=, b=)` |
| Cold fallback | Creates new `osqp::osqp()` | Creates new `clarabel::clarabel_solver()` |

**Differences from OSQP:**

- Clarabel's `solver_update()` accepts full sparse matrices (the R
  wrapper extracts `@x` internally), so CVXR passes `P` and `A`
  directly rather than extracting `P@x`/`A@x` manually.
- No iterate warm-starting (no `x0`/`y0` passed). Clarabel
  reinitializes iterates on each `solve()` call but reuses the
  factored KKT system structure.
- When `warm_start = TRUE` on the cold path, the settings are
  automatically adjusted to disable presolve, chordal decomposition,
  and dropzeros — features that would block future updates.

**Why not cache the result?**

Unlike SCS (which passes `initial = list(x, y, s)` to warm-start
iterates), Clarabel does not accept initial iterates. The warm start
benefit comes entirely from reusing the solver's internal data
structures, so only the solver object and data snapshot need caching.

### 2.4 What stays the same

- `reduction_invert` — unchanged; the result format from
  `solver_solve()` is identical to `clarabel()`.
- `dims_to_solver_dict_clarabel` — unchanged.
- `clarabel_psd_format_mat_fn` — unchanged.
- All existing `psolve(..., solver = "CLARABEL")` calls without
  `warm_start = TRUE` behave identically.

## 3. Usage from CVXR

### Basic warm start

```r
library(CVXR)

x <- Variable(2)
P_data <- matrix(c(2, 0, 0, 1), 2, 2)
objective <- Minimize(0.5 * quad_form(x, P_data) - c(2, 3) %*% x)
constraints <- list(x >= 0, sum(x) <= 2)
prob <- Problem(objective, constraints)

# First solve — cold start, creates and caches solver
val1 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)

# Change a parameter and re-solve — warm start, reuses solver
prob2 <- Problem(
  Minimize(0.5 * quad_form(x, P_data) - c(4, 1) %*% x),
  constraints
)
val2 <- psolve(prob2, solver = "CLARABEL", warm_start = TRUE)
```

### Parametric programming loop

The primary use case: solving a sequence of related problems where
only the data changes (e.g., MPC, cross-validation, regularization
paths).

```r
library(CVXR)

n <- 50; m <- 100
A_data <- matrix(rnorm(m * n), m, n)
b_data <- rnorm(m)

x <- Variable(n)
lambda <- Parameter(pos = TRUE)
objective <- Minimize(sum_squares(A_data %*% x - b_data) + lambda * p_norm(x, 1))
prob <- Problem(objective)

# Regularization path
lambdas <- 10^seq(2, -2, length.out = 20)
results <- vector("list", length(lambdas))
for (i in seq_along(lambdas)) {
  value(lambda) <- lambdas[i]
  results[[i]] <- psolve(prob, solver = "CLARABEL", warm_start = TRUE, verbose = FALSE)
}
```

### When warm start helps most

- **Parametric sweeps**: Same structure, varying data.
- **Iterative algorithms**: ADMM outer loops, cutting-plane methods.
- **Real-time applications**: MPC with shifting horizon.

### When warm start doesn't help

- **Different problem structures**: Adding/removing constraints
  changes the sparsity pattern. The solver falls back to cold start
  automatically.
- **Very different data**: If the data changes drastically, the
  factored KKT system provides little benefit.

## 4. Testing

Add tests following the pattern in
`tests/testthat/test-osqp-warm-start.R`:

```r
test_that("Clarabel warm start produces same result as cold", {
  x <- Variable(2)
  obj <- Minimize(sum_squares(x) - 2 * x[1] - 3 * x[2])
  con <- list(x >= 0, sum(x) <= 2)
  prob <- Problem(obj, con)

  # Cold solve
  val_cold <- psolve(prob, solver = "CLARABEL", warm_start = FALSE)
  x_cold <- value(x)

  # Warm solve (first call is effectively cold, cached)
  val_warm1 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  # Second call reuses cached solver
  val_warm2 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)

  expect_equal(val_cold, val_warm1, tolerance = 1e-6)
  expect_equal(val_cold, val_warm2, tolerance = 1e-6)
})

test_that("Clarabel warm start with changed data", {
  x <- Variable(2)
  c_param <- Parameter(2)
  value(c_param) <- c(2, 3)
  obj <- Minimize(sum_squares(x) - t(c_param) %*% x)
  con <- list(x >= 0, sum(x) <= 2)
  prob <- Problem(obj, con)

  val1 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  x1 <- value(x)

  value(c_param) <- c(4, 1)
  val2 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  x2 <- value(x)

  # Solutions should differ

  expect_false(isTRUE(all.equal(x1, x2, tolerance = 1e-4)))
})

test_that("Clarabel warm start falls back on dimension change", {
  x <- Variable(2)
  prob1 <- Problem(Minimize(sum_squares(x)), list(x >= 0))
  psolve(prob1, solver = "CLARABEL", warm_start = TRUE)

  # Different problem structure — should fall back to cold
  y <- Variable(3)
  prob2 <- Problem(Minimize(sum_squares(y)), list(y >= 0))
  val <- psolve(prob2, solver = "CLARABEL", warm_start = TRUE)
  expect_true(is.finite(val))
})
```

## 5. Summary of changes

| Component | File | Change |
|-----------|------|--------|
| `clarabel` R package | `R/clarabel.R` | New: `clarabel_solver()`, `solver_solve()`, `solver_update()`, `solver_is_update_allowed()` |
| `clarabel` Rust | `src/rust/src/lib.rs` | New: `ClarabelSolver` struct with `new/solve/update_data/is_update_allowed` |
| CVXR | `R/164_clarabel_conif.R` | Modified: `solve_via_data` gains warm start cache logic |
| CVXR | `tests/testthat/` | New: `test-clarabel-warm-start.R` |
