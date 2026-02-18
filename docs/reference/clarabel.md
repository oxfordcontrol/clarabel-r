# Interface to 'Clarabel', an interior point conic solver

Clarabel solves linear programs (LPs), quadratic programs (QPs),
second-order cone programs (SOCPs) and semidefinite programs (SDPs). It
also solves problems with exponential and power cone constraints. The
specific problem solved is:

Minimize \$\$\frac{1}{2}x^TPx + q^Tx\$\$ subject to \$\$Ax + s = b\$\$
\$\$s \in K\$\$ where \\x \in R^n\\, \\s \in R^m\\, \\P = P^T\\ and
nonnegative-definite, \\q \in R^n\\, \\A \in R^{m\times n}\\, and \\b
\in R^m\\. The set \\K\\ is a composition of convex cones.

## Usage

``` r
clarabel(A, b, q, P = NULL, cones, control = list(), strict_cone_order = TRUE)
```

## Arguments

- A:

  a matrix of constraint coefficients.

- b:

  a numeric vector giving the primal constraints

- q:

  a numeric vector giving the primal objective

- P:

  a symmetric positive semidefinite matrix, default `NULL`

- cones:

  a named list giving the cone sizes, see “Cone Parameters” below for
  specification

- control:

  a list giving specific control parameters to use in place of default
  values, with an empty list indicating the default control parameters.
  Specified parameters should be correctly named and typed to avoid Rust
  system panics as no sanitization is done for efficiency reasons

- strict_cone_order:

  a logical flag, default `TRUE` for forcing order of cones described
  below. If `FALSE` cones can be specified in any order and even
  repeated and directly passed to the solver without type and length
  checks

## Value

named list of solution vectors x, y, s and information about run

## Details

The order of the rows in matrix \\A\\ has to correspond to the order
given in the table “Cone Parameters”, which means means rows
corresponding to *primal zero cones* should be first, rows corresponding
to *non-negative cones* second, rows corresponding to *second-order
cone* third, rows corresponding to *positive semidefinite cones* fourth,
rows corresponding to *exponential cones* fifth and rows corresponding
to *power cones* at last.

When the parameter `strict_cone_order` is `FALSE`, one can specify the
cones in any order and even repeat them in the order they appear in the
`A` matrix. See below.

### Clarabel can solve

1.  linear programs (LPs)

2.  second-order cone programs (SOCPs)

3.  exponential cone programs (ECPs)

4.  power cone programs (PCPs)

5.  problems with any combination of cones, defined by the parameters
    listed in “Cone Parameters” below

### Cone Parameters

The table below shows the cone parameter specifications. Mathematical
definitions are in the vignette.

|  |  |  |  |  |
|----|----|----|----|----|
|  | **Parameter** | **Type** | **Length** | **Description** |
|  | `z` | integer | \\1\\ | number of primal zero cones (dual free cones), |
|  |  |  |  | which corresponds to the primal equality constraints |
|  | `l` | integer | \\1\\ | number of linear cones (non-negative cones) |
|  | `q` | integer | \\\ge 1\\ | vector of second-order cone sizes |
|  | `s` | integer | \\\ge 1\\ | vector of positive semidefinite cone sizes |
|  | `ep` | integer | \\1\\ | number of primal exponential cones |
|  | `p` | numeric | \\\ge 1\\ | vector of primal power cone parameters |
|  | `gp` | list | \\\ge 1\\ | list of named lists of two items, `a` : a numeric vector of at least 2 exponent terms in the product summing to 1, and `n` : an integer dimension of generalized power cone parameters |

When the parameter `strict_cone_order` is `FALSE`, one can specify the
cones in the order they appear in the `A` matrix. The `cones` argument
in such a case should be a named list with names matching `^z*`
indicating primal zero cones, `^l*` indicating linear cones, and so on.
For example, either of the following would be valid:
`list(z = 2L, l = 2L, q = 2L, z = 3L, q = 3L)`, or,
`list(z1 = 2L, l1 = 2L, q1 = 2L, zb = 3L, qx = 3L)`, indicating a zero
cone of size 2, followed by a linear cone of size 2, followed by a
second-order cone of size 2, followed by a zero cone of size 3, and
finally a second-order cone of size 3. Generalized power cones
parameters have to specified as named lists, e.g.,
`list(z = 2L, gp1 = list(a = c(0.3, 0.7), n = 3L), gp2 = list(a = c(0.5, 0.5), n = 1L))`.

*Note that when `strict_cone_order = FALSE`, types of cone parameters
such as integers, reals have to be explicit since the parameters are
directly passed to the Rust interface with no sanity checks!*

## See also

[`clarabel_control()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_control.md)

## Examples

``` r
A <- matrix(c(1, 1), ncol = 1)
b <- c(1, 1)
obj <- 1
cone <- list(z = 2L)
control <- clarabel_control(tol_gap_rel = 1e-7, tol_gap_abs = 1e-7, max_iter = 100)
clarabel(A = A, b = b, q = obj, cones = cone, control = control)
#> $x
#> [1] 1
#> 
#> $z
#> [1] -0.5 -0.5
#> 
#> $s
#> [1] 0 0
#> 
#> $obj_val
#> [1] 1
#> 
#> $obj_val_dual
#> [1] 1
#> 
#> $status
#> [1] 2
#> 
#> $solve_time
#> [1] 0.000189584
#> 
#> $iterations
#> [1] 0
#> 
#> $r_prim
#> [1] 0
#> 
#> $r_dual
#> [1] 0
#> 
#> $info
#> $info$μ
#> [1] 1
#> 
#> $info$sigma
#> [1] 1
#> 
#> $info$step_length
#> [1] 0
#> 
#> $info$cost_primal
#> [1] 1
#> 
#> $info$cost_dual
#> [1] 1
#> 
#> $info$res_primal
#> [1] 0
#> 
#> $info$res_dual
#> [1] 0
#> 
#> $info$res_primal_inf
#> [1] 1
#> 
#> $info$res_dual_inf
#> [1] 1.414214
#> 
#> $info$gap_abs
#> [1] 0
#> 
#> $info$gap_rel
#> [1] 0
#> 
#> $info$ktratio
#> [1] 1
#> 
#> $info$solve_time
#> [1] 0.000189584
#> 
#> $info$iterations
#> [1] 0
#> 
#> $info$status
#> [1] 2
#> 
#> 
```
