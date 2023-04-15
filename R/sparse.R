##
## copied from the Rsymphony package from
## Kurt Hornik and Stefan Theussl and Reinhard Harter
##
## Simple functions for converting "matrix" type objects into the
## sparse "column major order" (CSC, modulo offsets) format used by
## SYMPHONY.

## matind: vector of the row indices corresponding to each entry of
##   value
## values: vector of the values of nonzero entries of the constraint
##   matrix in column order.

make_csc_matrix <- function(x) UseMethod("make_csc_matrix")

make_csc_matrix.matrix <- function(x) {
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")

    ind <- which(x != 0, arr.ind = TRUE)
    list(matbeg = c(0L, cumsum(tabulate(ind[, 2L], ncol(x)))),
         matind = ind[, 1] - 1L,
         values = x[ind])
}

make_csc_matrix.simple_triplet_matrix <- function(x) {
    if(!inherits(x, "simple_triplet_matrix"))
        stop("Argument 'x' must be of class 'simple_triplet_matrix'.")

    ## The matrix method assumes that indices for non-zero entries are
    ## in row-major order, but the simple_triplet_matrix() constructor
    ## currently does not canonicalize accordingly ...
    ind <- order(x$j, x$i)
    list(matbeg = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
         matind = x$i[ind] - 1L,
         values = x$v[ind])
}

## Added by @bnaras for sparse symmetric matrix.

make_csc_symm_matrix <- function(x) UseMethod("make_csc_symm_matrix")

make_csc_symm_matrix.matrix  <- function(m) {
  ind <- which(m!=0 ,arr.ind = TRUE)
  ind  <- ind[ind[, 1] <= ind[, 2], ] ## keep upper part only
  values  <- m[ind]
  x  <- list(i = ind[, 1], j = ind[, 2])  ## triplet form
  ind  <- order(x$j, x$i)  ## may not be needed
  list(matbeg = c(0L, cumsum(tabulate(x$j[ind], ncol(m)))),
       matind = x$i[ind] - 1L,
       values = values)
}

make_csc_symm_matrix.simple_triplet_matrix  <- function(m) {
  ind <- which(m$i <= m$j)
  x <- list(i = m$i[ind] + 1L, j = m$j[ind] + 1L) ##make it 1-based
  values <- m$v[ind]
  ind  <- order(x$j, x$i)  ## may not be needed
  list(matbeg = c(0L, cumsum(tabulate(x$j[ind], m$ncol))),
       matind = x$i[ind] - 1L,
       values = values)
}

