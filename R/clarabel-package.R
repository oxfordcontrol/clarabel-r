#' Interface to Clarabel solver implemented in Rust.
#'
#' @description Clarabel is a versatile interior point solver for convex programs using a new homogeneous embedding. It solves solves linear programs (LPs), quadratic programs (QPs), second-order cone programs (SOCPs), and problems with exponential and power cone constraints. For quadratic objectives, unlike interior point solvers based on the standard homogeneous self-dual embedding (HSDE) model, Clarabel handles quadratic objective without requiring any epigraphical reformulation of its objective function. It can therefore be significantly faster than other HSDE-based solvers for problems with quadratic objective functions. Infeasible problems are detected using a homogeneous embedding technique. See <https://clarabel.org/stable/>.
#'
#' @name clarabel-package
#' @useDynLib clarabel
#' @author Balasubramanian Narasimhan, Paul Goulart, Yuwen Chen
#' @keywords internal
"_PACKAGE"
