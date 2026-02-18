# Control parameters with default values and types in parenthesis

Control parameters with default values and types in parenthesis

## Usage

``` r
clarabel_control(
  max_iter = 200L,
  time_limit = Inf,
  verbose = TRUE,
  max_step_fraction = 0.99,
  tol_gap_abs = 1e-08,
  tol_gap_rel = 1e-08,
  tol_feas = 1e-08,
  tol_infeas_abs = 1e-08,
  tol_infeas_rel = 1e-08,
  tol_ktratio = 1e-06,
  reduced_tol_gap_abs = 5e-05,
  reduced_tol_gap_rel = 5e-05,
  reduced_tol_feas = 1e-04,
  reduced_tol_infeas_abs = 5e-05,
  reduced_tol_infeas_rel = 5e-05,
  reduced_tol_ktratio = 1e-04,
  equilibrate_enable = TRUE,
  equilibrate_max_iter = 10L,
  equilibrate_min_scaling = 1e-04,
  equilibrate_max_scaling = 10000,
  linesearch_backtrack_step = 0.8,
  min_switch_step_length = 0.1,
  min_terminate_step_length = 1e-04,
  max_threads = 0L,
  direct_kkt_solver = TRUE,
  direct_solve_method = c("qdldl", "mkl", "cholmod"),
  static_regularization_enable = TRUE,
  static_regularization_constant = 1e-08,
  static_regularization_proportional = .Machine$double.eps * .Machine$double.eps,
  dynamic_regularization_enable = TRUE,
  dynamic_regularization_eps = 1e-13,
  dynamic_regularization_delta = 2e-07,
  iterative_refinement_enable = TRUE,
  iterative_refinement_reltol = 1e-13,
  iterative_refinement_abstol = 1e-12,
  iterative_refinement_max_iter = 10L,
  iterative_refinement_stop_ratio = 5,
  presolve_enable = TRUE,
  input_sparse_dropzeros = FALSE,
  chordal_decomposition_enable = FALSE,
  chordal_decomposition_merge_method = c("none", "parent_child", "clique_graph"),
  chordal_decomposition_compact = FALSE,
  chordal_decomposition_complete_dual = FALSE
)
```

## Arguments

- max_iter:

  maximum number of iterations (`200L`)

- time_limit:

  maximum run time (seconds) (`Inf`)

- verbose:

  verbose printing (`TRUE`)

- max_step_fraction:

  maximum interior point step length (`0.99`)

- tol_gap_abs:

  absolute duality gap tolerance (`1e-8`)

- tol_gap_rel:

  relative duality gap tolerance (`1e-8`)

- tol_feas:

  feasibility check tolerance (primal and dual) (`1e-8`)

- tol_infeas_abs:

  absolute infeasibility tolerance (primal and dual) (`1e-8`)

- tol_infeas_rel:

  relative infeasibility tolerance (primal and dual) (`1e-8`)

- tol_ktratio:

  KT tolerance (`1e-7`)

- reduced_tol_gap_abs:

  reduced absolute duality gap tolerance (`5e-5`)

- reduced_tol_gap_rel:

  reduced relative duality gap tolerance (`5e-5`)

- reduced_tol_feas:

  reduced feasibility check tolerance (primal and dual) (`1e-4`)

- reduced_tol_infeas_abs:

  reduced absolute infeasibility tolerance (primal and dual) (`5e-5`)

- reduced_tol_infeas_rel:

  reduced relative infeasibility tolerance (primal and dual) (`5e-5`)

- reduced_tol_ktratio:

  reduced KT tolerance (`1e-4`)

- equilibrate_enable:

  enable data equilibration pre-scaling (`TRUE`)

- equilibrate_max_iter:

  maximum equilibration scaling iterations (`10L`)

- equilibrate_min_scaling:

  minimum equilibration scaling allowed (`1e-4`)

- equilibrate_max_scaling:

  maximum equilibration scaling allowed (`1e+4`)

- linesearch_backtrack_step:

  linesearch backtracking (`0.8`)

- min_switch_step_length:

  minimum step size allowed for asymmetric cones with PrimalDual scaling
  (`1e-1`)

- min_terminate_step_length:

  minimum step size allowed for symmetric cones && asymmetric cones with
  Dual scaling (`1e-4`)

- max_threads:

  maximum solver threads for multithreaded KKT solvers, 0 lets the
  solver choose for itself (`0L`)

- direct_kkt_solver:

  use a direct linear solver method (required true) (`TRUE`)

- direct_solve_method:

  direct linear solver (`"qdldl"`, `"mkl"` or `"cholmod"`) (`"qdldl"`)

- static_regularization_enable:

  enable KKT static regularization (`TRUE`)

- static_regularization_constant:

  KKT static regularization parameter (`1e-8`)

- static_regularization_proportional:

  additional regularization parameter w.r.t. the maximum abs diagonal
  term (`.Machine.double_eps^2`)

- dynamic_regularization_enable:

  enable KKT dynamic regularization (`TRUE`)

- dynamic_regularization_eps:

  KKT dynamic regularization threshold (`1e-13`)

- dynamic_regularization_delta:

  KKT dynamic regularization shift (`2e-7`)

- iterative_refinement_enable:

  KKT solve with iterative refinement (`TRUE`)

- iterative_refinement_reltol:

  iterative refinement relative tolerance (`1e-12`)

- iterative_refinement_abstol:

  iterative refinement absolute tolerance (`1e-12`)

- iterative_refinement_max_iter:

  iterative refinement maximum iterations (`10L`)

- iterative_refinement_stop_ratio:

  iterative refinement stalling tolerance (`5.0`)

- presolve_enable:

  whether to enable presolvle (`TRUE`)

- input_sparse_dropzeros:

  explicitly drop structural zeros from sparse data inputs (`FALSE`);
  see details

- chordal_decomposition_enable:

  whether to enable chordal decomposition for SDPs (`FALSE`)

- chordal_decomposition_merge_method:

  chordal decomposition merge method, one of `'none'`, `'parent_child'`
  or `'clique_graph'`, for SDPs (`'none'`)

- chordal_decomposition_compact:

  a boolean flag for SDPs indicating whether to assemble decomposed
  system in *compact* form for SDPs (`FALSE`)

- chordal_decomposition_complete_dual:

  a boolean flag indicating complete PSD dual variables after
  decomposition for SDPs

## Value

a list containing the control parameters.

## Details

Setting `input_sparse_dropzeros` to `TRUE` will disable parametric
updating functionality. See documentation of 'dropzeros' in Rust
`struct CscMatrix` for dropping structural zeros before passing to the
solver.

## Examples

``` r
# Default control parameters
ctrl <- clarabel_control()
ctrl$max_iter
#> [1] 200
# Custom tolerances and quiet output
ctrl <- clarabel_control(verbose = FALSE, tol_gap_rel = 1e-7, max_iter = 100L)
```
