#' @title Clarabel - Interior point conic solver
#'
#' @description Solves convex cone programs using an interior point
#'   method. The specific problem solved is: Minimize
#'   \deqn{\frac{1}{2}x^TPx + q^Tx} subject to \deqn{Ax + s = b}
#'   \deqn{s \in K} where \eqn{x \in R^n}, \eqn{s \in R^m}, \eqn{P =
#'   P^T} and nonnegative-definite, \eqn{q \in R^n}, \eqn{A \in
#'   R^{m\times n}}, and \eqn{b \in R^m}. The set \eqn{K} is a
#'   composition of convex cones.
#' @param A a matrix of constraint coefficients.
#' @param b a numeric vector giving the primal constraints
#' @param q a numeric vector giving the primal objective
#' @param P a symmetric positive semidefinite matrix, default
#'   \code{NULL}
#' @param cones a named list giving the cone sizes, see \dQuote{Cone
#'   Parameters} below for specification
#' @param control a list giving specific control parameters to use in place of default values, with an empty list indicating the default control parameters. Specified parameters should be correctly named and typed to avoid Rust system panics as no sanitization is done for efficiency reasons
#' @param strict_cone_order a logical flag, default `TRUE` for forcing order of cones described below. If `FALSE` cones can be specified in any order and even repeated and directly passed to the solver without type and length checks
#' @return named list of solution vectors x, y, s and information about run
#' @seealso [clarabel_control()]
#' @export clarabel
#'
#' @details
#'
#' The order of the rows in matrix \eqn{A} has to correspond to the
#' order given in the table \dQuote{Cone Parameters}, which means
#' means rows corresponding to \emph{primal zero cones} should be
#' first, rows corresponding to \emph{non-negative cones} second, rows
#' corresponding to \emph{second-order cone} third, rows corresponding to
#' \emph{exponential cones} fourth and rows corresponding to
#' \emph{power cones} at last.
#'
#' When the parameter `strict_cone_order` is `FALSE`, one can specify the cones in any order and even repeat them in the order they appear in the `A` matrix. See below.
#'
#' \subsection{Clarabel can solve}{ \enumerate{ \item linear programs
#' (LPs) \item second-order cone programs (SOCPs) \item exponential
#' cone programs (ECPs) \item power cone programs (PCPs) \item
#' problems with any combination of cones, defined by the parameters
#' listed in \dQuote{Cone Parameters} below } }
#'
#' \subsection{Cone Parameters}{
#' The table below shows the cone parameter specifications
#' \tabular{rllll}{
#'    \tab \bold{Parameter} \tab \bold{Type} \tab \bold{Length} \tab \bold{Description}                       \cr
#'    \tab \code{z}         \tab integer     \tab \eqn{1}       \tab number of primal zero cones (dual free cones),       \cr
#'    \tab                  \tab             \tab               \tab which corresponds to the primal equality constraints \cr
#'    \tab \code{l}         \tab integer     \tab \eqn{1}       \tab number of linear cones (non-negative cones)          \cr
#'    \tab \code{q}         \tab integer     \tab \eqn{\geq1}   \tab vector of second-order cone sizes                    \cr
#'    \tab \code{ep}        \tab integer     \tab \eqn{1}       \tab number of primal exponential cones                   \cr
#'    \tab \code{p}         \tab numeric     \tab \eqn{\geq1}   \tab vector of primal power cone parameters
#' } }
#'
#' When the parameter `strict_cone_order` is `FALSE`, one can specify
#' the cones in the order they appear in the `A` matrix. The `cones`
#' argument in such a case should be a named list with names matching
#' `^z*` indicating primal zero cones, `^l*` indicating linear cones,
#' and so on. For example, the following would be valid: `list(z1 =
#' 2L, l1 = 2L, q1 = 2L, zb = 3L, qx = 3L)`, indicating three zero
#' cones, followed by two linear cones, followed by two second-order
#' cones, followed by two zero cones, and finally 3 second-order
#' cones.
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
clarabel <- function(A, b, q, P = NULL, cones, control = list(), strict_cone_order = TRUE) {

  ## TBD check box cone parameters, bsize > 0  & bl, bu have lengths bsize - 1

  n_variables <- NCOL(A)
  n_constraints <- NROW(A)

  ## if (n_variables <= 0) {
  ##   stop("No variables")
  ## }

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

  ## data <- list(m = n_constraints, n = n_variables, c = obj,
  ##              Ai = Ai, Ap = Ap, Ax = Ax, b = b)

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

  # Sanitize control parameters
  control <- do.call(clarabel_control, control)
  if (strict_cone_order) cones <- sanitize_cone_spec(cones)
  clarabel_solve(n_constraints, n_variables, Ai, Ap, Ax, b, q, Pi, Pp, Px, cones, control)
}

#' Clarabel control parameters with default values and types in parenthesis
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
                             presolve_enable = TRUE) {

  params <- as.list(environment())
  ## Match string arg to make sure it is kosher
  direct_solve_method <- match.arg(direct_solve_method)
  params$direct_solve_method <- direct_solve_method
  
  ## Rust has strict type and length checks, so try to avoid panics
  bool_params <- c("verbose", "equilibrate_enable", "direct_kkt_solver", "static_regularization_enable",
                 "dynamic_regularization_enable", "iterative_refinement_enable", "presolve_enable")

  int_params <- c("max_iter", "equilibrate_max_iter", "iterative_refinement_max_iter")

  string_params <- "direct_solve_method"
  
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
  for (x in string_params) {
    params[[x]] <- as.character(params[[x]])
  }
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
### @return sanitized cone specifications
sanitize_cone_spec <- function(cone_spec) {
  result <- list()
  cone_names <- names(cone_spec)
  if (length(intersect(cone_names, c("z", "l", "q", "ep", "p"))) != length(cone_spec)) {
    stop("sanitize_cone_spec: unspecified cone parameters")
  }
  
  for (name in cone_names) {
    value <- cone_spec[[name]]    
    if (grepl("^z|^l|^ep", name)) {
      if (length(value) != 1L || value <= 0) {
        stop("sanitize_cone_spec: z, l, ep should be scalars > 0")
      }
      result[[name]] <- as.integer(value)
    } else if (grepl("^q", name)) {
      value <- as.integer(value)
      qnames <- paste0("q", seq_along(value))
      for (i in seq_along(value)) {
        result[[qnames[i]]] <- value[i]
      }
    } else { ## must be "p", power cone
      if (any(value <= 0.0 | value >= 1.0)) {
        stop("sanitize_cone_spec: Power cone parameter should be in (0, 1)")
      }
      pnames <- paste0("p", seq_along(value))
      for (i in seq_along(value)) {
        result[[pnames[i]]] <- value[i]
      }
    }
  }
  result
}
