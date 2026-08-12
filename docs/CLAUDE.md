# NA

## Clarabel-r

The main goal here is to ensure that this package has full support for
warm starts for `CVXR`. So in order, the tasks are:

1.  Robustify this rust-based R package by examining any bad choices.
2.  Note that we use the latest [`savvy`
    crate](file:///Users/naras/GitHub/savvy). This repo is local and you
    are free to read it fully.
3.  We also need `BLAS` so the `r-src` crate provides R’s BLAS.
4.  Provide facilities for warm starts especially in this [`CVXR`
    project](https://github.com/bnaras/new_design/) which uses
    `clarabel`.

## Tasks

1.  Add simple examples to functions where not present.
2.  The vignette should add a section on updating. The Semi-definite
    program with a changed matrix would be a good example. It can go
    right after it and before `Control Parameters`.
3.  Update `NEWS.md`.

## New Tasks

1.  I thought I’d used the values from
    `~/GitHub/Clarabel.rs/src/solver/implementations/default/settings.rs`
    to populate our default settings. However, some seem to be off in
    `R/clarabel.R` for
    [`clarabel_control()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_control.md).
    We need to align in a dev branch of this repo.
