# Clarabel Solver Examples

## Introduction

The first two examples are from the original [Clarabel
documentation](https://clarabel.org) and the third is from
[SCS](https://www.cvxgrp.org/scs/).

## 1. Basic Quadratic Program Example

Suppose that we want to solve the following 2-dimensional quadratic
programming problem:

``` math
\begin{array}{ll} \text{minimize} &  3x_1^2 + 2x_2^2 - x_1 - 4x_2\\
\text{subject to} &  -1 \leq x \leq 1, ~ x_1 = 2x_2
\end{array}
```

We will show how to solve this problem using Clarabel in R.

The first step is to put the problem data into the standard form
expected by the solver.

### 1.1. Objective function

The Clarabel solver’s default configuration expects problem data in the
form $`\frac{1}{2}x^\top P x + q^\top x`$.\
We therefore define the objective function data as

``` math
P = 2 \cdot \begin{bmatrix} 3 & 0 \\ 0 & 2\end{bmatrix}
\mbox{ and }
q = \begin{bmatrix} -1 \\ -4\end{bmatrix}.
```

### 1.2. Constraints

The solver’s default configuration expects constraints in the form
$`Ax + s = b`$, where $`s \in \mathcal{K}`$ for some composite cone
$`\mathcal{K}`$. We have 1 equality constraint and 4 inequalities, so we
require the first element of $`s`$ to be zero (i.e. the first constraint
will correspond to the equality) and all other elements $`s_i \ge 0`$.
Our cone constraint on $`s`$ is therefore

``` math
s \in \mathcal K = \{0\}^1 \times \mathbb{R}^4_{\ge 0}.
```

Define the constraint data as

``` math
A =
\begin{bmatrix} 1 & -2 \\ 1 & 0 \\ 0 & 1 \\ -1 & 0 \\ 0 & -1\end{bmatrix}
\mbox{ and }
b=\begin{bmatrix} 0 \\ 1 \\ 1 \\ 1 \\ 1 \end{bmatrix}.
```

Note that Clarabel expects inputs in Compressed Sparse Column (CSC)
format for both $`P`$ and $`A`$ and will try to convert them if not so.

### 1.3. Solution

`P`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(``2`` ``*`` `[`c`](https://rdrr.io/r/base/c.html)`(``3``, ``0``, ``0``, ``2``)``, nrow ``=`` ``2``, ncol ``=`` ``2``, sparse ``=`` ``TRUE``)`` ``P`` ``<-`` ``as``(``P``, ``"symmetricMatrix"``)`` ``# P needs to be a symmetric matrix`` ``q`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``-``1``, ``-``4``)`` ``A`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``1``, ``0``, ``-``1``, ``0``, ``-``2``, ``0``, ``1``, ``0``, ``-``1``)``, ncol ``=`` ``2``, sparse ``=`` ``TRUE``)`` ``b`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``1``, ``1``, ``1``, ``1``)`` ``cones`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``z ``=`` ``1L``, l ``=`` ``4L``)`` ``## 1 equality and 4 inequalities, in order`` ``s`` ``<-`` `[`clarabel`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md)`(``A ``=`` ``A``, b ``=`` ``b``, q ``=`` ``q``, P ``=`` ``P``, cones ``=`` ``cones``)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution status, description: = (%d, %s)\n"``,`` `` ``s``$``status``, `[`solver_status_descriptions`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_status_descriptions.md)`(``)``[``s``$``status``]``)``)`` ``#> Solution status, description: = (2, Solver terminated with a solution.)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution: (x1, x2) = (%f, %f)\n"``, ``s``$``x``[``1``]``, ``s``$``x``[``2``]``)``)`` ``#> Solution: (x1, x2) = (0.428571, 0.214286)`

## 2. Basic Second-order Cone Programming Example

We want to solve the following 2-dimensional optimization problem:

``` math
\begin{array}{ll} \text{minimize} & x_2^2\\[2ex]
\text{subject to} &  \left\|\begin{pmatrix} 2x_1 \\ x_2 \end{pmatrix}
- \begin{pmatrix} 2 \\ 2 \end{pmatrix}\right\|_2 \le 1
\end{array}
```

### 2.1. Objective function

The Clarabel solver’s default configuration expects problem data in the
form $`\frac{1}{2}x^\top P x + q^\top x`$.\
We therefore define the objective function data as

``` math
P = 2 \cdot \begin{bmatrix} 0 & 0 \\ 0 & 1\end{bmatrix}
\mbox{ and }
q = \begin{bmatrix} 0 \\ 0\end{bmatrix}.
```

### 2.2. Constraints

The solver’s default configuration expects constraints in the form
$`Ax + s = b`$, where $`s \in \mathcal{K}`$ for some composite cone
$`\mathcal{K}`$. We have a single constraint on the 2-norm of a vector,
so we rewrite

``` math
\left\|\begin{pmatrix} 2x_1 \\ x_2 \end{pmatrix} - \begin{pmatrix} 2 \\ 2 \end{pmatrix}\right\|_2 \le 1
\quad \Longleftrightarrow \quad
\begin{pmatrix} 1 \\ 2x_1 - 2\\ x_2 - 2 \end{pmatrix} \in \mathcal{K}_{SOC}
```
which puts our constraint in the form $`b - Ax \in \mathcal{K}_{SOC}`$.

### 2.3. Solution

`P`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(``2`` ``*`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``0``, ``0``, ``1``)``, nrow ``=`` ``2``, ncol ``=`` ``2``, sparse ``=`` ``TRUE``)`` ``P`` ``<-`` ``as``(``P``, ``"symmetricMatrix"``)`` ``# P needs to be a symmetric matrix`` ``q`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``0``)`` ``A`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``-``2.0``, ``0``, ``0``, ``0``, ``1.0``)``, nrow ``=`` ``3``, ncol ``=`` ``2``, sparse ``=`` ``TRUE``)`` ``b`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``-``2``, ``-``2``)`` ``cones`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``q ``=`` ``3L``)`` ``s`` ``<-`` `[`clarabel`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md)`(``A ``=`` ``A``, b ``=`` ``b``, q ``=`` ``q``, P ``=`` ``P``, cones ``=`` ``cones``)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution status, description: = (%d, %s)\n"``,`` `` ``s``$``status``, `[`solver_status_descriptions`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_status_descriptions.md)`(``)``[``s``$``status``]``)``)`` ``#> Solution status, description: = (2, Solver terminated with a solution.)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution (x1, x2) = (%f, %f)\n"``, ``s``$``x``[``1``]``, ``s``$``x``[``2``]``)``)`` ``#> Solution (x1, x2) = (1.000000, -1.000000)`

## 3. Semidefinite Cone Programming

Semidefinite problems are the only ones here that use dense linear
algebra, which they get from your R installation’s BLAS and LAPACK
rather than from this package. That could affect both speed and results;
see Section 7.

Semidefinite cones have to be specified in a particular form. We borrow
from the [documentation for the SCS
solver](https://www.cvxgrp.org/scs/api/cones.html#sdcone) which has
similar calling conventions.

The symmetric positive semidefinite cone of matrices is the set

``` math
\{S \in \mathbf{R}^{k \times k} \mid  S = S^\top,  x^\top S x \geq 0 \
\forall x \in \mathbb{R}^k \}
```

and for short, we use $`S \succeq 0`$ to denote membership. Clarabel
vectorizes this cone in a special way which we detail here.

Clarabel assumes that the input data corresponding to semidefinite cones
have been vectorized by scaling the off-diagonal entries by $`\sqrt{2}`$
and stacking the upper triangular elements column-wise.
([SCS](https://www.cvxgrp.org/scs/index.html) uses the lower triangular
elements.) For a $`k \times k`$ matrix variable (or data matrix) this
operation would create a vector of length $`k(k+1)/2`$. Scaling by
$`\sqrt{2}`$ is required to preserve the inner-product.

*This must be done for the rows of both $`A`$ and $`b`$ that correspond
to semidefinite cones and must be done independently for each
semidefinite cone.*

More explicitly, we want to express $`\text{Trace}(Y S)`$ as
$`\text{vec}(Y)^\top \text{vec}(S)`$, where the $`\text{vec}`$ operation
takes the (assumed to be symmetric) $`k \times k`$ matrix

``` math
\begin{aligned}
S =  \begin{bmatrix}
        S_{11} & S_{12} & \ldots & S_{1k}  \\
        S_{21} & S_{22} & \ldots & S_{2k}  \\
        \vdots & \vdots & \ddots & \vdots  \\
        S_{k1} & S_{k2} & \ldots & S_{kk}  \\
      \end{bmatrix}
\end{aligned}
```

and produces a vector consisting of the upper triangular elements scaled
and arranged as

``` math
\text{vec}(S) = (S_{11}, \sqrt{2} S_{12}, S_{22}, \sqrt{2}S_{13}, \ldots, \sqrt{2}
S_{1k}, \sqrt{2}S_{2k}, \sqrt{2}S_{3k}, \dots, \sqrt{2}S_{k-1,k}, S_{kk}) \in
\mathbb{R}^{k(k+1)/2}.
```

To recover the matrix solution this operation must be inverted on the
components of the vectors returned by Clarabel corresponding to each
semidefinite cone. That is, the off-diagonal entries must be scaled by
$`1/\sqrt{2}`$ and the upper triangular entries are filled in by copying
the values of lower triangular entries. Explicitly, the inverse
operation takes vector $`s \in
\mathbb{R}^{k(k+1)/2}`$ and produces the matrix

``` math
\begin{aligned}
\text{mat}(s) =  \begin{bmatrix}
s_{1} & s_{2} / \sqrt{2} & \ldots & s_{k(k-1)/2-1} / \sqrt{2}  \\
s_{2} / \sqrt{2} & s_{3} & \ldots & s_{k(k-1)/2} / \sqrt{2}  \\
\vdots & \vdots & \ddots & \vdots  \\
s_{k(k-1)/2-1} / \sqrt{2} & s_{k(k-1)/2} / \sqrt{2} & \ldots & s_{k(k+1) / 2}  \\
\end{bmatrix}
\in \mathbb{R}^{k \times k}.
\end{aligned}
```

So the cone definition that Clarabel uses is

``` math
\mathcal{S}_+^k = \{ \text{vec}(S) \mid S \succeq 0\} = \{s \in
\mathbb{R}^{k(k+1)/2} \mid \text{mat}(s) \succeq 0 \}.
```

Below are two functions to implement both $`\text{vec}`$ and
$`\text{mat}`$.

`#' Return an vectorization of symmetric matrix using the upper triangular part,`` ``#' still in column order.`` ``#' @param S a symmetric matrix`` ``#' @return vector of values`` ``vec`` ``<-`` ``function``(``S``)`` ``{`` `` ``n`` ``<-`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``S``)`` `` ``sqrt2`` ``<-`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``2.0``)`` `` ``upper_tri`` ``<-`` `[`upper.tri`](https://rdrr.io/r/base/lower.tri.html)`(``S``, diag ``=`` ``FALSE``)`` `` ``S``[``upper_tri``]`` ``<-`` ``S``[``upper_tri``]`` ``*`` ``sqrt2`` `` ``S``[`[`upper.tri`](https://rdrr.io/r/base/lower.tri.html)`(``S``, diag ``=`` ``TRUE``)``]`` ``}`` `` ``#' Return the symmetric matrix from the [vec] vectorization`` ``#' @param v a vector`` ``#' @return a symmetric matrix`` ``mat`` ``<-`` ``function``(``v``)`` ``{`` `` ``n`` ``<-`` ``(`[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``8`` ``*`` `[`length`](https://rdrr.io/r/base/length.html)`(``v``)`` ``+`` ``1``)`` ``-`` ``1``)`` ``/`` ``2`` `` ``sqrt2`` ``<-`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``2.0``)`` `` ``S`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0``, ``n``, ``n``)`` `` ``upper_tri`` ``<-`` `[`upper.tri`](https://rdrr.io/r/base/lower.tri.html)`(``S``, diag ``=`` ``TRUE``)`` `` ``S``[``upper_tri``]`` ``<-`` ``v`` ``/`` ``sqrt2`` `` ``S`` ``<-`` ``S`` ``+`` `[`t`](https://rdrr.io/r/base/t.html)`(``S``)`` `` `[`diag`](https://rdrr.io/r/base/diag.html)`(``S``)`` ``<-`` `[`diag`](https://rdrr.io/r/base/diag.html)`(``S``)`` ``/`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``2``)`` `` ``S`` ``}`

### 3.1. Example

Consider the problem:

``` math
\begin{array}{ll} \text{minimize} &  x_1 - x_2 + x_3\\
\text{subject to} & B_1 - A_{11}x_1 - A_{12}x_2 - A_{13}x_3 \succeq 0\\
                  & B_2 - A_{21}x_1 - A_{22}x_2 - A_{23}x_3 \succeq 0
\end{array}
```
where
``` math
A_{11}=\begin{bmatrix} -7 & -11 \\ -11 & 3\end{bmatrix}\mbox{, }
A_{12}=\begin{bmatrix} 7 & -18 \\ -18 & 8\end{bmatrix}\mbox{, }
A_{13}=\begin{bmatrix} -2 & -8 \\ -8 & 1\end{bmatrix},
```

``` math
A_{21}=\begin{bmatrix} -21 & -11 & 0\\ -11 & 10 & 8\\ 0 & 8 &
5\end{bmatrix}\mbox{, }
A_{22}=\begin{bmatrix} 0 & 10 & 16\\ 10 & -10 & -10\\ 16 & -10 &
3\end{bmatrix}\mbox{, }
A_{23}=\begin{bmatrix} -5 & 2 & -17\\ 2 & -6 & 8\\ -17 & 8 &
6\end{bmatrix},
```

and

``` math
B_1=\begin{bmatrix} 33 & -9 \\ -9 & 26\end{bmatrix}\mbox{, }
B_2=\begin{bmatrix} 14 & 9 & 40\\ 9 & 91 & 10\\ 40 & 10 &
15\end{bmatrix}.
```

The constraints involve symmetric positive semidefinite cones over
variables $`x \in \mathbb{R}^n`$ and $`S \in
\mathbb{R}^{k \times k}`$

``` math
B - \sum_{i=1}^n \mathcal{A}_i x_i = S \succeq 0
```

where data $`B, \mathcal{A}_1, \ldots, \mathcal{A}_n \in \mathbb{R}^{k
\times k}`$ are symmetric. We can write this in the canonical form over
a new variable $`s \in \mathcal{S}_+^k`$:

``` math
\begin{aligned}
\begin{align}
s &= \text{vec}(S)\\
&= \text{vec}(B - \sum_{i=1}^n \mathcal{A}_i x_i) \\
&= \text{vec}(B) - \sum_{i=1}^n \text{vec}(\mathcal{A}_i) x_i \\
&= b - Ax
\end{align}
\end{aligned}
```

using the fact that $`\text{vec}`$ is linear, where $`b =
\text{vec}(B)`$ and

``` math
A =
\begin{bmatrix}
\text{vec}(\mathcal{A}_1) & \text{vec}(\mathcal{A}_2) & \cdots & \text{vec}(\mathcal{A}_n)
\end{bmatrix}
```

i.e., the vectors $`\text{vec}(\mathcal{A}_i)`$ stacked columnwise. This
is in a form that we can input into Clarabel. To recover the matrix
solution from the optimal solution returned by Clarabel, we simply use
$`S^\star =
\text{mat}(s^\star)`$.

We have two such constraints and will therefore construct the vectors
for both cones in order and specify the appropriate dimensions, 2 and 3,
respectively.

`q`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``-``1``, ``1``)`` ``# objective: x_1 - x2 + x_3`` ``A11`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``7``, ``-``11``, ``-``11``, ``3``)``, nrow ``=`` ``2``)`` ``A12`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``7``, ``-``18``, ``-``18``, ``8``)``, nrow ``=`` ``2``)`` ``A13`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``2``, ``-``8``, ``-``8``, ``1``)``, nrow ``=`` ``2``)`` `` ``A21`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``21``, ``-``11``, ``0``, ``-``11``, ``10``, ``8``, ``0``, ``8``, ``5``)``, nrow ``=`` ``3``)`` ``A22`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``10``, ``16``, ``10``, ``-``10``, ``-``10``, ``16``, ``-``10``, ``3``)``, nrow ``=`` ``3``)`` ``A23`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``5``, ``2``, ``-``17``, ``2``, ``-``6``, ``8``, ``-``17``, ``8``, ``6``)``, nrow ``=`` ``3``)`` `` ``B1`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``33``, ``-``9``, ``-``9``, ``26``)``, nrow ``=`` ``2``)`` ``B2`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``14``, ``9``, ``40``, ``9``, ``91``, ``10``, ``40``, ``10``, ``15``)``, nrow ``=`` ``3``)`` `` ``A`` ``<-`` `[`rbind`](https://rdrr.io/r/base/cbind.html)`(`` `` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``vec``(``A11``)``, ``vec``(``A12``)``, ``vec``(``A13``)``)``, ``# first psd constraint`` `` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``vec``(``A21``)``, ``vec``(``A22``)``, ``vec``(``A23``)``)`` ``# second psd constraint`` ``)`` ``b`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``vec``(``B1``)``, ``vec``(``B2``)``)`` ``# stack both psd constraints`` ``cones`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``s ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``3``)``)`` ``# cone dimensions`` ``s`` ``<-`` `[`clarabel`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md)`(``A ``=`` ``A``, b ``=`` ``b``, q ``=`` ``q``, cones ``=`` ``cones``)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution status, description: = (%d, %s)\n"``,`` `` ``s``$``status``, `[`solver_status_descriptions`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_status_descriptions.md)`(``)``[``s``$``status``]``)``)`` ``#> Solution status, description: = (2, Solver terminated with a solution.)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution (x1, x2, x3) = (%f, %f, %f)\n"``, ``s``$``x``[``1``]``, ``s``$``x``[``2``]``, ``s``$``x``[``3``]``)``)`` ``#> Solution (x1, x2, x3) = (-0.367746, 1.898333, -0.887466)`

### 3.2. Chordal Decomposition of Sparse Semidefinite Programs

The cost of a semidefinite program grows very quickly with the side
length of its semidefinite cones. When the aggregate sparsity pattern of
those constraints is *chordal*, however, a large cone can be replaced by
several small ones that overlap only on their shared entries, and the
solver works with the small ones instead. This is enabled by default
through `chordal_decomposition_enable`.

The following builds a linear matrix inequality
$`F_0 + \sum_i x_i F_i \succeq 0`$ in which every data matrix is banded
with half-bandwidth 2. A banded pattern is chordal, with cliques of size
3, so an $`n \times n`$ cone decomposes into roughly $`n`$ tiny ones.

`banded_sdp`` ``<-`` ``function``(``n``, ``w`` ``=`` ``2``, ``m`` ``=`` ``5``, ``seed`` ``=`` ``7``)`` ``{`` `` ``band_sym`` ``<-`` ``function``(``s``)`` ``{`` `` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``s``)`` `` ``M`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0``, ``n``, ``n``)`` `` ``for`` ``(``i`` ``in`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n``)``)`` `` ``for`` ``(``j`` ``in`` ``i``:`[`min`](https://rdrr.io/r/base/Extremes.html)`(``i`` ``+`` ``w``, ``n``)``)`` ``M``[``i``, ``j``]`` ``<-`` ``M``[``j``, ``i``]`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``1``)`` `` ``M`` `` ``}`` `` ``Fs`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``m``)``, ``function``(``k``)`` ``band_sym``(``1000`` ``*`` ``seed`` ``+`` ``k``)``)`` `` ``F0`` ``<-`` ``band_sym``(``1000`` ``*`` ``seed`` ``+`` ``m`` ``+`` ``1``)`` `` ``## shift F0 so that x = 0 is strictly feasible`` `` ``F0`` ``<-`` ``F0`` ``+`` `[`diag`](https://rdrr.io/r/base/diag.html)`(``n``)`` ``*`` ``(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(`[`min`](https://rdrr.io/r/base/Extremes.html)`(`[`eigen`](https://rdrr.io/r/base/eigen.html)`(``F0``, only.values ``=`` ``TRUE``)``$``values``)``)`` ``+`` ``1``)`` `` `[`list`](https://rdrr.io/r/base/list.html)`(``A ``=`` ``-`[`do.call`](https://rdrr.io/r/base/do.call.html)`(``cbind``, `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``Fs``, ``vec``)``)``,`` `` b ``=`` ``vec``(``F0``)``,`` `` q ``=`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``Fs``, ``function``(``Fi``)`` `[`sum`](https://rdrr.io/r/base/sum.html)`(`[`diag`](https://rdrr.io/r/base/diag.html)`(``Fi``)``)``, `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``1``)``)``,`` `` cones ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``s ``=`` `[`as.integer`](https://rdrr.io/r/base/integer.html)`(``n``)``)``)`` ``}`` `` ``solve_timed`` ``<-`` ``function``(``p``, ``chordal``)`` ``{`` `` ``ctrl`` ``<-`` `[`clarabel_control`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_control.md)`(``verbose ``=`` ``FALSE``,`` `` chordal_decomposition_enable ``=`` ``chordal``)`` `` ``elapsed`` ``<-`` `[`system.time`](https://rdrr.io/r/base/system.time.html)`(`` `` ``s`` ``<-`` `[`clarabel`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md)`(``A ``=`` ``p``$``A``, b ``=`` ``p``$``b``, q ``=`` ``p``$``q``, cones ``=`` ``p``$``cones``, control ``=`` ``ctrl``)`` `` ``)``[[``"elapsed"``]``]`` `` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``chordal ``=`` ``chordal``,`` `` seconds ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``elapsed``)``,`` `` status ``=`` `[`names`](https://rdrr.io/r/base/names.html)`(`[`solver_status_descriptions`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_status_descriptions.md)`(``)``)``[``s``$``status``]``,`` `` objective ``=`` ``s``$``obj_val``)`` ``}`` `` ``p40`` ``<-`` ``banded_sdp``(``n ``=`` ``40``)`` `[`rbind`](https://rdrr.io/r/base/cbind.html)`(``solve_timed``(``p40``, ``FALSE``)``, ``solve_timed``(``p40``, ``TRUE``)``)`` ``#> chordal seconds status objective`` ``#> 1 FALSE 0.473 AlmostSolved -6.603129`` ``#> 2 TRUE 0.003 Solved -6.603129`

Two things are worth reading off that comparison. The decomposed solve
is faster by more than two orders of magnitude, and the two solves do
not terminate at the same accuracy: on this instance the undecomposed
solve stops at `AlmostSolved` while the decomposed one reaches `Solved`.
That is not a general rule — which side stops early depends on the
instance, and it goes both ways — but it is a reminder that
decomposition produces a genuinely different numerical problem, and that
the status is worth checking rather than assuming.

The gap widens sharply with the size of the cone. Measured separately,
since the larger rows are too slow to run while this vignette is built,
on R’s reference BLAS (Section 7):

| $`n`$ | cone length | decomposition off |     on | speedup | status (off / on)     |
|------:|------------:|------------------:|-------:|--------:|:----------------------|
|    20 |         210 |            0.013s | 0.001s |     13x | Solved / Solved       |
|    40 |         820 |            0.420s | 0.002s |    210x | AlmostSolved / Solved |
|    60 |        1830 |            4.393s | 0.004s |   1098x | Solved / Solved       |
|    80 |        3240 |           29.855s | 0.004s |   7464x | Solved / Solved       |
|   120 |        7260 |          424.978s | 0.008s |  53122x | Solved / AlmostSolved |

The undecomposed column grows roughly like $`n^6`$; the decomposed one
is close to linear in $`n`$. Seven minutes becomes eight milliseconds.
Note also that the reduced-accuracy termination lands on the decomposed
side at $`n = 120`$ and on the undecomposed side at $`n = 40`$, which is
what is meant above by it going both ways.

Two consequences are worth knowing. Dual variables reported for
semidefinite constraints are not in general the same as those obtained
without decomposition, because the decomposition is reversed and the
dual completed afterwards, and a positive semidefinite completion is not
unique. And a problem that was actually decomposed cannot have its data
updated in place, which matters for the warm starts described next. Pass
`chordal_decomposition_enable = FALSE` to recover the previous behavior
in either case.

## 4. Updating Problem Data (Warm Starts)

When solving a sequence of related problems that share the same sparsity
structure, it is more efficient to create a persistent solver and update
only the data that changes between solves. This avoids rebuilding the
solver’s internal data structures each time.

Updates are refused when the solver has altered the structure of the
problem, since the stored factorization no longer corresponds to the
data you would be updating. Three settings can do that:
`presolve_enable`, `input_sparse_dropzeros`, and
`chordal_decomposition_enable`.

What matters is whether the transformation *actually happened*, not
whether it was permitted. Presolve only rewrites the problem if it finds
something to remove, and chordal decomposition only engages for
semidefinite constraints whose sparsity pattern is decomposable. Linear,
quadratic and second-order cone problems are never chordally decomposed,
so leaving `chordal_decomposition_enable` at its default costs them
nothing. The example below is a semidefinite program, but its two blocks
are small and dense, hence not decomposable, so its updates are allowed
even with chordal decomposition enabled.

Rather than reasoning about which transformation might apply, ask the
solver:
[`solver_is_update_allowed()`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_is_update_allowed.md)
reports whether the instance you actually built can be updated.

### 4.1. Example: SDP with a changed constraint

We revisit the semidefinite program from Section 3, but now suppose the
constraint matrix $`B_1`$ changes from
$`\begin{bmatrix} 33 & -9 \\ -9 & 26 \end{bmatrix}`$ to
$`\begin{bmatrix} 40 & -12 \\ -12 & 30 \end{bmatrix}`$.

First, we set up and solve the original problem using a persistent
solver.

`q_sdp`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``-``1``, ``1``)`` ``A11`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``7``, ``-``11``, ``-``11``, ``3``)``, nrow ``=`` ``2``)`` ``A12`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``7``, ``-``18``, ``-``18``, ``8``)``, nrow ``=`` ``2``)`` ``A13`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``2``, ``-``8``, ``-``8``, ``1``)``, nrow ``=`` ``2``)`` `` ``A21`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``21``, ``-``11``, ``0``, ``-``11``, ``10``, ``8``, ``0``, ``8``, ``5``)``, nrow ``=`` ``3``)`` ``A22`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``10``, ``16``, ``10``, ``-``10``, ``-``10``, ``16``, ``-``10``, ``3``)``, nrow ``=`` ``3``)`` ``A23`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``-``5``, ``2``, ``-``17``, ``2``, ``-``6``, ``8``, ``-``17``, ``8``, ``6``)``, nrow ``=`` ``3``)`` `` ``B1`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``33``, ``-``9``, ``-``9``, ``26``)``, nrow ``=`` ``2``)`` ``B2`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``14``, ``9``, ``40``, ``9``, ``91``, ``10``, ``40``, ``10``, ``15``)``, nrow ``=`` ``3``)`` `` ``A_sdp`` ``<-`` `[`rbind`](https://rdrr.io/r/base/cbind.html)`(`` `` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``vec``(``A11``)``, ``vec``(``A12``)``, ``vec``(``A13``)``)``,`` `` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``vec``(``A21``)``, ``vec``(``A22``)``, ``vec``(``A23``)``)`` ``)`` ``b_sdp`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``vec``(``B1``)``, ``vec``(``B2``)``)`` ``cones_sdp`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``s ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``3``)``)`` `` ``## Create a persistent solver with updates enabled`` ``ctrl`` ``<-`` `[`clarabel_control`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_control.md)`(``presolve_enable ``=`` ``FALSE``, verbose ``=`` ``FALSE``)`` ``solver`` ``<-`` `[`clarabel_solver`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.md)`(``A ``=`` ``A_sdp``, b ``=`` ``b_sdp``, q ``=`` ``q_sdp``, cones ``=`` ``cones_sdp``, control ``=`` ``ctrl``)`` `[`solver_is_update_allowed`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_is_update_allowed.md)`(``solver``)`` ``# should be TRUE`` ``#> [1] TRUE`` ``sol1`` ``<-`` `[`solver_solve`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_solve.md)`(``solver``)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Original solution (x1, x2, x3) = (%f, %f, %f)\n"``,`` `` ``sol1``$``x``[``1``]``, ``sol1``$``x``[``2``]``, ``sol1``$``x``[``3``]``)``)`` ``#> Original solution (x1, x2, x3) = (-0.367746, 1.898333, -0.887466)`

Now update the constraint with the new $`B_1`$ and re-solve.

`B1_new`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``40``, ``-``12``, ``-``12``, ``30``)``, nrow ``=`` ``2``)`` ``b_sdp_new`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``vec``(``B1_new``)``, ``vec``(``B2``)``)`` ``# only B1 changes`` `[`solver_update`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_update.md)`(``solver``, b ``=`` ``b_sdp_new``)`` ``sol2`` ``<-`` `[`solver_solve`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_solve.md)`(``solver``)`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Updated solution (x1, x2, x3) = (%f, %f, %f)\n"``,`` `` ``sol2``$``x``[``1``]``, ``sol2``$``x``[``2``]``, ``sol2``$``x``[``3``]``)``)`` ``#> Updated solution (x1, x2, x3) = (-0.404810, 2.166303, -0.649288)`

The solver reuses its internal factorization, making the second solve
faster than constructing a new solver from scratch.

## 5. Cone Specifications

The following cones can be specified in Clarabel.

| Parameter | Type | Length | Description | Definition (per parameter element) |
|:---|:---|:---|:---|:---|
| z | integer | 1 | Number of primal zero cones (dual free cones), which corresponds to the primal equality constraints | $`\{ 0 \}^{z}`$ |
| l | integer | 1 | Number of linear cones (non-negative cones) | $`\{ x \in \mathbb{R}^{l} : x_i \ge 0, \forall i=1,\dots,l \}`$ |
| q | integer | \>= 1 | Vector of second-order cone sizes | $`\{ (t,x) \in \mathbb{R}^{q}  :  \lVert x\rVert_2  \leq t \}`$ |
| s | integer | \>= 1 | Vector of positive semidefinite cone sizes | Upper triangular part of the positive semidefinite cone $`S^s_+`$. The elements $`x`$ of this cone represent the columnwise stacking of the upper triangular part of a positive semidefinite matrix $`X \in S^s_+`$, so that $`x \in R^d`$ with $`d = s(s+1)/2`$ |
| ep | integer | 1 | Number of primal exponential cones | $`\{(x, y, z) : y > 0,~~ ye^{x/y} \le z \}`$ |
| p | numeric | \>= 1 | Vector of primal power cone parameters | $`\{(x, y, z) : x^p y^{(1-p)} \ge  \lVert z\rVert,~ (x,y) \ge 0 \}`$ with $`p \in (0,1)`$ |
| gp | list | \>= 1 | List of named lists of two items, `a` : the numeric vector of at least 2 exponent terms, and `n` : an integer dimension of generalized power cone parameters | $`\{(x, y) \in R^{len(a)} \times R^n : \prod\limits_{a_i \in a} x_i^{a_i} \ge \lVert y\rVert_2,~ x \ge 0 \}`$ with $`a_i \in (0,1)`$ and $`\sum a_i = 1`$ |

Generalized power cone parameters are specified as list of two-item
lists, with component named $`a`$ denoting the exponents and the named
component $`n`$ denoting the dimension.

One can specify cones in any order if `strict_cone_order` is set to
`FALSE` in the call to
[`clarabel()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md)
but one has to ensure that parameter types are strictly specified for
the values, e.g. `5L` for integers, `0.` for reals etc.

## 6. Control parameters

Clarabel has a number of parameters that control its behavior, including
verbosity, time limits, and tolerances; see help on
[`clarabel_control()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_control.md).
As an example, in the last problem, we can reduce the number of
iterations.

`P`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(``2`` ``*`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``0``, ``0``, ``1``)``, nrow ``=`` ``2``, ncol ``=`` ``2``, sparse ``=`` ``TRUE``)`` ``P`` ``<-`` ``as``(``P``, ``"symmetricMatrix"``)`` ``# P needs to be a symmetric matrix`` ``q`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``0``)`` ``A`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``-``2.0``, ``0``, ``0``, ``0``, ``1.0``)``, nrow ``=`` ``3``, ncol ``=`` ``2``, sparse ``=`` ``TRUE``)`` ``b`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``-``2``, ``-``2``)`` ``cones`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``q ``=`` ``3L``)`` ``s`` ``<-`` `[`clarabel`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.md)`(``A ``=`` ``A``, b ``=`` ``b``, q ``=`` ``q``, P ``=`` ``P``, cones ``=`` ``cones``,`` `` control ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``max_iter ``=`` ``3``)``)`` ``## Reduced number of iterations`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution status, description: = (%d, %s)\n"``,`` `` ``s``$``status``, `[`solver_status_descriptions`](https://oxfordcontrol.github.io/clarabel-r/reference/solver_status_descriptions.md)`(``)``[``s``$``status``]``)``)`` ``#> Solution status, description: = (5, Solver terminated with a solution (reduced accuracy))`` `[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Solution (x1, x2) = (%f, %f)\n"``, ``s``$``x``[``1``]``, ``s``$``x``[``2``]``)``)`` ``#> Solution (x1, x2) = (1.000000, -0.999998)`

Note the different status, which should always be checked in code.

## 7. A Note on BLAS and LAPACK

This package links whatever BLAS and LAPACK your R installation uses;
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) reports
which one is active.

On macOS the CRAN binary ships both R’s reference BLAS and a vecLib
(Accelerate) build, and defaults to the reference one. Switching is a
symlink:

``` sh
cd /Library/Frameworks/R.framework/Resources/lib
ln -sf libRblas.vecLib.dylib libRblas.dylib   # vecLib
ln -sf libRblas.0.dylib      libRblas.dylib   # reference: the default
```

Repeating the Section 3.2 sweep under both libraries on a ARM Mac, moved
solve times by under 1% at every size, in both directions:

| $`n`$ | reference, off | vecLib, off | reference, on | vecLib, on |
|------:|---------------:|------------:|--------------:|-----------:|
|    20 |         0.015s |      0.016s |        0.001s |     0.002s |
|    40 |         0.429s |      0.444s |        0.003s |     0.003s |
|    60 |         4.475s |      4.576s |        0.003s |     0.004s |
|    80 |        30.078s |     30.535s |        0.005s |     0.005s |
|   120 |       426.752s |    429.727s |        0.008s |     0.009s |

For these problems not much time is not spent in the BLAS and LAPACK
(`dgemm`, `dsyrk`, `dpotrf`, `dsyevr`) routines: the dense operations
act on the $`n \times n`$ cone, while the cost is dominated by
factorizing a KKT system whose size follows the vectorized cone length
$`n(n+1)/2`$, which is 7260 rows at $`n = 120`$.

The results do differ between the two libraries. At $`n = 40`$ the
decomposed dual differed by $`1.3 \times 10^{-1}`$ and the termination
status moved from `Solved` to `AlmostSolved`. (See the [R for macOS
FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html) and [R
Installation and
Administration](https://cran.r-project.org/doc/manuals/r-release/R-admin.html)
on how Accelerate differs from the reference BLAS.)

The Section 3.2 timings were measured on the reference BLAS.
