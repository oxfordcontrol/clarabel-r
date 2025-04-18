#' @title Interface to 'Clarabel', an interior point conic solver
#'
#' @description
#'
#' Clarabel solves linear programs (LPs), quadratic programs (QPs),
#'   second-order cone programs (SOCPs) and semidefinite programs
#'   (SDPs). It also solves problems with exponential and power cone
#'   constraints. The specific problem solved is:
#'
#' Minimize \deqn{\frac{1}{2}x^TPx + q^Tx} subject to \deqn{Ax + s =
#'   b} \deqn{s \in K} where \eqn{x \in R^n}, \eqn{s \in R^m}, \eqn{P
#'   = P^T} and nonnegative-definite, \eqn{q \in R^n}, \eqn{A \in
#'   R^{m\times n}}, and \eqn{b \in R^m}. The set \eqn{K} is a
#'   composition of convex cones.
#'
#' @param A a matrix of constraint coefficients.
#' @param b a numeric vector giving the primal constraints
#' @param q a numeric vector giving the primal objective
#' @param P a symmetric positive semidefinite matrix, default
#'   \code{NULL}
#' @param cones a named list giving the cone sizes, see \dQuote{Cone
#'   Parameters} below for specification
#' @param control a list giving specific control parameters to use in
#'   place of default values, with an empty list indicating the
#'   default control parameters. Specified parameters should be
#'   correctly named and typed to avoid Rust system panics as no
#'   sanitization is done for efficiency reasons
#' @param strict_cone_order a logical flag, default `TRUE` for forcing
#'   order of cones described below. If `FALSE` cones can be specified
#'   in any order and even repeated and directly passed to the solver
#'   without type and length checks
#' @return named list of solution vectors x, y, s and information
#'   about run
#' @seealso [clarabel_control()]
#' @export clarabel
#'
#' @details
#'
#' The order of the rows in matrix \eqn{A} has to correspond to the
#'   order given in the table \dQuote{Cone Parameters}, which means
#'   means rows corresponding to \emph{primal zero cones} should be
#'   first, rows corresponding to \emph{non-negative cones} second,
#'   rows corresponding to \emph{second-order cone} third, rows
#'   corresponding to \emph{positive semidefinite cones} fourth, rows
#'   corresponding to \emph{exponential cones} fifth and rows
#'   corresponding to \emph{power cones} at last.
#'
#' When the parameter `strict_cone_order` is `FALSE`, one can specify
#' the cones in any order and even repeat them in the order they
#' appear in the `A` matrix. See below.
#'
#' \subsection{Clarabel can solve}{ \enumerate{ \item linear programs
#' (LPs) \item second-order cone programs (SOCPs) \item exponential
#' cone programs (ECPs) \item power cone programs (PCPs) \item
#' problems with any combination of cones, defined by the parameters
#' listed in \dQuote{Cone Parameters} below } }
#'
#' \subsection{Cone Parameters}{
#' The table below shows the cone parameter specifications. Mathematical definitions are in the vignette.
#' \tabular{rllll}{
#'    \tab \bold{Parameter} \tab \bold{Type} \tab \bold{Length} \tab \bold{Description}                       \cr
#'    \tab \code{z}   \tab integer  \tab \eqn{1}       \tab number of primal zero cones (dual free cones),       \cr
#'    \tab            \tab          \tab               \tab which corresponds to the primal equality constraints \cr
#'    \tab \code{l}   \tab integer  \tab \eqn{1}       \tab number of linear cones (non-negative cones)          \cr
#'    \tab \code{q}   \tab integer  \tab \eqn{\ge 1}   \tab vector of second-order cone sizes                    \cr
#'    \tab \code{s}   \tab integer  \tab \eqn{\ge 1}   \tab vector of positive semidefinite cone sizes           \cr
#'    \tab \code{ep}  \tab integer  \tab \eqn{1}       \tab number of primal exponential cones                   \cr
#'    \tab \code{p}   \tab numeric  \tab \eqn{\ge 1}   \tab vector of primal power cone parameters               \cr
#'    \tab \code{gp}  \tab list     \tab \eqn{\ge 1}  \tab list of named lists of two items, `a` : a numeric vector of at least 2 exponent terms in the product summing to 1, and `n` : an integer dimension of generalized power cone parameters
#' } }
#'
#' When the parameter `strict_cone_order` is `FALSE`, one can specify
#' the cones in the order they appear in the `A` matrix. The `cones`
#' argument in such a case should be a named list with names matching
#' `^z*` indicating primal zero cones, `^l*` indicating linear cones,
#' and so on. For example, either of the following would be valid: `list(z =
#' 2L, l = 2L, q = 2L, z = 3L, q = 3L)`, or, `list(z1 =
#' 2L, l1 = 2L, q1 = 2L, zb = 3L, qx = 3L)`, indicating a zero
#' cone of size 2, followed by a linear cone of size 2, followed by a second-order
#' cone of size 2, followed by a zero cone of size 3, and finally a second-order
#' cone of size 3. Generalized power cones parameters have to specified as named lists, e.g., `list(z = 2L, gp1 = list(a = c(0.3, 0.7), n = 3L), gp2 = list(a = c(0.5, 0.5), n = 1L))`.
#'
#' _Note that when `strict_cone_order = FALSE`, types of cone parameters such as integers, reals have to be explicit since the parameters are directly passed to the Rust interface with no sanity checks.!_
#'
#' @examples
#' A <- matrix(c(1, 1), ncol = 1)
#' b <- c(1, 1)
#' obj <- 1
#' cone <- list(z = 2L)
#' control <- clarabel_control(tol_gap_rel = 1e-7, tol_gap_abs = 1e-7, max_iter = 100)
#' clarabel(A = A, b = b, q = obj, cones = cone, control = control)
#' 
#  ---------------------------------------------------------
clarabel <- function(A, b, q, P = NULL, cones, control = list(),
                     strict_cone_order = TRUE) {

  m <- length(b); n <- length(q);
  n_variables <- ncol(A)
  n_constraints <- nrow(A)
  
  if (m != n_constraints) stop("A and b incompatible dimensions.")
  if (n != n_variables) stop("A and q incompatible dimensions.")
  if (!is.null(P)) {
    if (n != ncol(P)) stop("P and q incompatible dimensions.")
    if (ncol(P) != nrow(P)) stop("P not square.")
  }

  # Sanitize control parameters
  control <- do.call(clarabel_control, control)
  if (strict_cone_order) {
    cones_and_nvars <- sanitize_cone_spec(cones)
    cones <- cones_and_nvars[["cones"]]
    nvars <- cones_and_nvars[["nvars"]]
  } else {
    nvars <- nvars(cones)
  }
  if (sum(nvars) != m) stop("Constraint dimensions inconsistent with size of cones.")
  
  ## TBD check box cone parameters, bsize > 0  & bl, bu have lengths bsize - 1

  if ( inherits(A, "dgCMatrix") ) {
    Ai <- A@i
    Ap <- A@p
    Ax <- A@x
  } else {
    csc <- make_csc_matrix(A)
    Ai <- csc[["matind"]]
    Ap <- csc[["matbeg"]]
    Ax <- csc[["values"]]
  }

  if (!is.null(P)) {
    if (inherits(P, "dsCMatrix") ) {
      Pi <- P@i
      Pp <- P@p
      Px <- P@x
    } else {
      csc  <- make_csc_symm_matrix(P)
      Pi <- csc[["matind"]]
      Pp <- csc[["matbeg"]]
      Px <- csc[["values"]]
    }
  } else {
    Pi <- integer(0)
    Pp <- integer(0)
    Px <- numeric(0)
  }

  .Call(savvy_clarabel_solve__impl, n_constraints, n_variables, Ai, Ap, Ax, b, q, Pi, Pp, Px, cones, control, PACKAGE = "clarabel")

  ## clarabel_solve(n_constraints, n_variables, Ai, Ap, Ax, b, q, Pi, Pp, Px, cones, control)
}

#' Control parameters with default values and types in parenthesis
#'
#' @param max_iter maximum number of iterations (`200L`)		 
#' @param time_limit maximum run time (seconds) (`Inf`)			 
#' @param verbose verbose printing (`TRUE`)		 
#' @param max_step_fraction maximum interior point step length (`0.99`)		 
#' @param tol_gap_abs absolute duality gap tolerance (`1e-8`)		 
#' @param tol_gap_rel relative duality gap tolerance (`1e-8`)		 
#' @param tol_feas feasibility check tolerance (primal and dual) (`1e-8`)		 
#' @param tol_infeas_abs absolute infeasibility tolerance (primal and dual) (`1e-8`)		 
#' @param tol_infeas_rel relative infeasibility tolerance (primal and dual) (`1e-8`)		 
#' @param tol_ktratio KT tolerance (`1e-7`)		 
#' @param reduced_tol_gap_abs reduced absolute duality gap tolerance (`5e-5`)		 
#' @param reduced_tol_gap_rel reduced relative duality gap tolerance (`5e-5`)		 
#' @param reduced_tol_feas reduced feasibility check tolerance (primal and dual) (`1e-4`)		 
#' @param reduced_tol_infeas_abs reduced absolute infeasibility tolerance (primal and dual) (`5e-5`)		 
#' @param reduced_tol_infeas_rel reduced relative infeasibility tolerance (primal and dual) (`5e-5`)		 
#' @param reduced_tol_ktratio reduced KT tolerance (`1e-4`)		 
#' @param equilibrate_enable enable data equilibration pre-scaling (`TRUE`)		 
#' @param equilibrate_max_iter maximum equilibration scaling iterations (`10L`)			 
#' @param equilibrate_min_scaling minimum equilibration scaling allowed (`1e-4`)		 
#' @param equilibrate_max_scaling maximum equilibration scaling allowed (`1e+4`)		 
#' @param linesearch_backtrack_step linesearch backtracking (`0.8`)			 
#' @param min_switch_step_length minimum step size allowed for asymmetric cones with PrimalDual scaling (`1e-1`)		 
#' @param min_terminate_step_length minimum step size allowed for symmetric cones && asymmetric cones with Dual scaling (`1e-4`)		 
#' @param max_threads maximum solver threads for multithreaded KKT solvers, 0 lets the solver choose for itself (`0L`)
#' @param direct_kkt_solver use a direct linear solver method (required true) (`TRUE`)		 
#' @param direct_solve_method direct linear solver (`"qdldl"`, `"mkl"` or `"cholmod"`) (`"qdldl"`)		 
#' @param static_regularization_enable enable KKT static regularization (`TRUE`)		 
#' @param static_regularization_constant KKT static regularization parameter (`1e-8`)		 
#' @param static_regularization_proportional additional regularization parameter w.r.t. the maximum abs diagonal term (`.Machine.double_eps^2`) 
#' @param dynamic_regularization_enable enable KKT dynamic regularization (`TRUE`)		 
#' @param dynamic_regularization_eps KKT dynamic regularization threshold (`1e-13`)		 
#' @param dynamic_regularization_delta KKT dynamic regularization shift (`2e-7`)		 
#' @param iterative_refinement_enable KKT solve with iterative refinement (`TRUE`)		 
#' @param iterative_refinement_reltol iterative refinement relative tolerance (`1e-12`)		 
#' @param iterative_refinement_abstol iterative refinement absolute tolerance (`1e-12`)		 
#' @param iterative_refinement_max_iter iterative refinement maximum iterations (`10L`)			 
#' @param iterative_refinement_stop_ratio iterative refinement stalling tolerance (`5.0`)
#' @param presolve_enable whether to enable presolvle (`TRUE`)
#' @param chordal_decomposition_enable whether to enable chordal decomposition for SDPs (`FALSE`)
#' @param chordal_decomposition_merge_method chordal decomposition merge method, one of `'none'`, `'parent_child'` or `'clique_graph'`, for SDPs (`'none'`)
#' @param chordal_decomposition_compact a boolean flag for SDPs indicating whether to assemble decomposed system in _compact_ form for SDPs (`FALSE`)
#' @param chordal_decomposition_complete_dual a boolean flag indicating complete PSD dual variables after decomposition for SDPs
#' @return a list containing the control parameters.
#' @export clarabel_control
clarabel_control <- function(
                             ## Main algorithm settings
                             max_iter = 200L,
                             time_limit = Inf,
                             verbose = TRUE,
                             max_step_fraction = 0.99,
                             ## Full accuracy settings
                             tol_gap_abs = 1e-8,
                             tol_gap_rel = 1e-8,
                             tol_feas = 1e-8,
                             tol_infeas_abs = 1e-8,
                             tol_infeas_rel = 1e-8,
                             tol_ktratio = 1e-6,
                             ## Reduced accuracy settings
                             reduced_tol_gap_abs = 5e-5,
                             reduced_tol_gap_rel = 5e-5,
                             reduced_tol_feas = 1e-4,
                             reduced_tol_infeas_abs = 5e-5,
                             reduced_tol_infeas_rel = 5e-5,
                             reduced_tol_ktratio = 1e-4,
                             ## data equilibration settings
                             equilibrate_enable = TRUE,
                             equilibrate_max_iter = 10L,
                             equilibrate_min_scaling = 1e-4,
                             equilibrate_max_scaling = 1e4,
                             ## Step size settings
                             linesearch_backtrack_step = 0.8,
                             min_switch_step_length = 1e-1,
                             min_terminate_step_length = 1e-4,
                             ## maximum solver threads for multithreaded KKT solvers
                             ## choosing 0 lets the solver choose for itself
                             max_threads = 0L,
                             ## Linear solver settings
                             direct_kkt_solver = TRUE,
                             direct_solve_method = c("qdldl", "mkl", "cholmod"),
                             ## static regularization parameters
                             static_regularization_enable = TRUE,
                             static_regularization_constant = 1e-8,
                             static_regularization_proportional = .Machine$double.eps * .Machine$double.eps,
                             ## dynamic regularization parameters
                             dynamic_regularization_enable = TRUE,
                             dynamic_regularization_eps = 1e-13,
                             dynamic_regularization_delta = 2e-7,
                             iterative_refinement_enable = TRUE,
                             iterative_refinement_reltol = 1e-13,
                             iterative_refinement_abstol = 1e-12,
                             iterative_refinement_max_iter = 10L,
                             iterative_refinement_stop_ratio = 5.0,
                             presolve_enable = TRUE,
                             chordal_decomposition_enable = FALSE,
                             chordal_decomposition_merge_method = c('none', 'parent_child', 'clique_graph'),
                             chordal_decomposition_compact = FALSE,
                             chordal_decomposition_complete_dual = FALSE
                             ) {

  params <- as.list(environment())

  ## Match string args to make sure it is kosher
  direct_solve_method <- match.arg(direct_solve_method)
  params$direct_solve_method <- direct_solve_method
  chordal_decomposition_merge_method <- match.arg(chordal_decomposition_merge_method)
  params$chordal_decomposition_merge_method <- chordal_decomposition_merge_method
  
  ## Rust has strict type and length checks, so try to avoid panics
  bool_params <- c("verbose", "equilibrate_enable", "direct_kkt_solver", "static_regularization_enable",
                   "dynamic_regularization_enable", "iterative_refinement_enable", "presolve_enable",
                   "chordal_decomposition_enable", "chordal_decomposition_compact",
                   "chordal_decomposition_complete_dual")

  int_params <- c("max_iter", "max_threads", "equilibrate_max_iter", "iterative_refinement_max_iter")

  string_params <- c("direct_solve_method", "chordal_decomposition_merge_method") # Might need to uncomment character coercion below, if length > 1
  
  if (any(sapply(params, length) != 1L)) stop("clarabel_control: arguments should be scalars!")
  if (any(unlist(params[int_params]) < 0)) stop("clarabel_control: integer arguments should be >= 0!")
 
  ## The rest
  float_params <- setdiff(names(params), c(bool_params, int_params, string_params))

  for (x in bool_params) {
    params[[x]] <- as.logical(params[[x]])
  }
  for (x in int_params) {
    params[[x]] <- as.integer(params[[x]])
  }
  ## Not needed since match.arg takes care of this for the single string param
  ## for (x in string_params) {
  ##   params[[x]] <- as.character(params[[x]])
  ## }
  for (x in float_params) {
    params[[x]] <- as.numeric(params[[x]])
  }
  params
}

#' Return the solver status description as a named character vector
#' @return a named list of solver status descriptions, in order of status codes returned by the solver
#' @examples
#' solver_status_descriptions()[2] ## for solved problem
#' solver_status_descriptions()[8] ## for max iterations limit reached
#' @export
solver_status_descriptions <- function() {
  c(Unsolved = "Problem is not solved (solver hasn't run).",
    Solved = "Solver terminated with a solution.",
    PrimalInfeasible = "Problem is primal infeasible.  Solution returned is a certificate of primal infeasibility.",
    DualInfeasible = "Problem is dual infeasible.  Solution returned is a certificate of dual infeasibility.",
    AlmostSolved = "Solver terminated with a solution (reduced accuracy)",
    AlmostPrimalInfeasible = "Problem is primal infeasible.  Solution returned is a certificate of primal infeasibility (reduced accuracy)",
    AlmostDualInfeasible = "Problem is dual infeasible.  Solution returned is a certificate of dual infeasibility (reduced accuracy)",
    MaxIterations = "Iteration limit reached before solution or infeasibility certificate found",
    MaxTime = "Time limit reached before solution or infeasibility certificate found",
    NumericalError = "Solver terminated with a numerical error",
    InsufficientProgress = "Solver terminated due to lack of progress"
    )
}

### Sanitize cone specifications
### @param cone_spec a list of cone specifications
### @return a named list of sanitized cone specifications and the number of variables for each cone (`nvars`)
sanitize_cone_spec <- function(cone_spec) {
  cone_names <- names(cone_spec)
  
  ## Simple sanity checks
  if ((nc <- length(cone_names)) == 0L) {
    stop("sanitize_cone_spec: no cone parameters specified")    
  } 
  if (length(intersect(cone_names, c("z", "l", "q", "s", "ep", "p"))) != nc) {
    stop("sanitize_cone_spec: repeated cone parameters or unknown cone parameters specified")
  }
  
  ## Check lengths as noted cone parameters table for ?clarabel
  ## First, scalars
  z <- as.integer(cone_spec[["z"]]); zl <- length(z); nvar_z <- sum(z);
  l <- as.integer(cone_spec[["l"]]); ll <- length(l); nvar_l <- sum(l);
  ep <- as.integer(cone_spec[["ep"]]); epl <- length(ep); nvar_ep <- sum(ep) * 3 ## 3 variables per exp cone
  if (zl > 1 || ll > 1 || epl > 1) {
    stop("sanitize_cone_spec: z, l, ep should be scalars")
  }
  if (any(c(z, l, ep) < 0L)) {
    stop("sanitize_cone_spec: z, l, ep should be scalars > 0")
  }
  
  ## Now the others
  ## SOC 
  q <- as.integer(cone_spec[["q"]]); ql <- length(q); nvar_q <- sum(q);
  if (any(q <= 0L)) stop("sanitize_cone_spec: SOC dimensions should be > 0")
  q <- as.list(q); names(q) <- rep("q", ql);
  
  ## PSD 
  s <- as.integer(cone_spec[["s"]]); sl <- length(s); nvar_s <- if (sl > 0) sum(sapply(s, triangular_number)) else 0; ## triangular number
  if (any(s <= 0L)) stop("sanitize_cone_spec: PSD dimensions should be > 0")
  s <- as.list(s); names(s) <- rep("s", sl);
  
  ## Power Cone
  p <- as.numeric(cone_spec[["p"]]); pl <- length(p); nvar_p <- pl * 3; ## 3 variables per power cone
  if (any(p <= 0)) stop("sanitize_cone_spec: Power cone parameter should be > 0")
  p <- as.list(p); names(p) <- rep("p", pl);  

  ## Power Cone
  gp <- cone_spec[["gp"]];  gpl <- length(gp); result <- sanitize_gp_params_and_get_nvars(gp);
  gp <- result$sanitized_gp_params; names(gp) <- rep("gp", gpl);  nvar_gp <- result$nvar;
  
  cones <- c(
    if (zl > 0) list(z = z),
    if (ll > 0) list(l = l),
    if (ql > 0) q,
    if (sl > 0) s,
    if (epl > 0) list(ep = ep),
    if (pl > 0) p,
    if (gpl > 0) gp
  )

  nvars <- c(z = nvar_z,
             l = nvar_l,
             q = nvar_q,
             s = nvar_s,
             ep = nvar_ep,
             p = nvar_p,
             gp = nvar_gp)
  list(cones = cones, nvars = nvars)
}

### Return the number of variables used for each type of cone based on cone specification
### @param cones the cone specifications, expecting strict_cone_order to be FALSE
### @return a named vector of number of variables per each type of cone
nvars <- function(cones) {
  cone_names <- names(cones)
  ## separate out gp cones from others.
  gp_indices <- grep("^gp", cone_names)
  if (length(gp_indices) > 0) {
    gp_cones <- cones[gp_indices]
    nvar_gp <- sum(sapply(gp_cones, function(x) length(x[["a"]]) + x[["n"]]))
  } else {
    nvar_gp <- 0L
  }
  other_indices <- grep("^gp", cone_names, invert = TRUE)
  cones <- unlist(cones[other_indices])
  cone_names <- names(cones)
  nvar_z <- if (length(matched <- grep("^z", cone_names)) > 0) sum(cones[matched]) else 0L
  nvar_l <- if (length(matched <- grep("^l", cone_names)) > 0) sum(cones[matched]) else 0L
  nvar_q <- if (length(matched <- grep("^q", cone_names)) > 0) sum(cones[matched]) else 0L
  nvar_s <- if (length(matched <- grep("^s", cone_names)) > 0) sum(sapply(cones[matched], triangular_number)) else 0L
  nvar_ep <- if (length(matched <- grep("^ep", cone_names)) > 0) 3 * sum(cones[matched]) else 0L
  nvar_p <- if (length(matched <- grep("^p", cone_names)) > 0) 3 * length(cones[matched]) else 0L
  c(z = nvar_z,
    l = nvar_l,
    q = nvar_q,
    s = nvar_s,
    ep = nvar_ep,
    p = nvar_p,
    gp = nvar_gp)
}
  
## Return the n-th triangular number
triangular_number <- function(n) {
  n * (n + 1) / 2
}

## Check Generalized Power Cone Params
sanitize_gp_params_and_get_nvars <- function(gp) {
  # Generalized power cone α.len() + dim
  if (length(gp) == 0) {
    list(sanitized_gp_params = list(), nvar = 0L)
  } else  {
    sanitized_gp_params <- 
      lapply(gp, function(x) {
        par_names <- sort(names(x))
        if (length(x) != 2L || !identical(par_names, c("a", "n"))) {
          stop("Generalized power cone param list should be a list of two elements named 'a' and 'n'")
        }
        exps <- x[["a"]]
        if (length(exps) < 2L || any(exps <= 0) || any(exps >= 1) || abs(sum(exps) - 1.0) > 0.0) {
          stop("Improper Generalized power cone exponents!")
        }
        n <- x[["n"]]
        if (length(n) != 1L || n <= 0) stop("Improper Generalized power cone dimension!")
        list(a = as.numeric(exps), n = as.integer(n))
      })
    list(sanitized_gp_params = sanitized_gp_params, nvar = sum(sapply(sanitized_gp_params, function(x) length(x[["a"]]) + x[["n"]])))
  }
}


