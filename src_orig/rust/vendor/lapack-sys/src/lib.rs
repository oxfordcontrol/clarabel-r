//! Bindings to [LAPACK] \(Fortran).
//!
//! ## [Architecture]
//!
//! [architecture]: https://blas-lapack-rs.github.io/architecture
//! [lapack]: https://en.wikipedia.org/wiki/LAPACK

#![allow(non_camel_case_types)]
#![no_std]

extern crate libc;

use libc::{c_char, c_int, size_t};

/// A complex number with 64-bit parts.
pub type c_double_complex = __BindgenComplex<f64>;

/// A complex number with 32-bit parts.
pub type c_float_complex = __BindgenComplex<f32>;

pub type lapack_complex_double = c_double_complex;
pub type lapack_complex_float = c_float_complex;
pub type lapack_int = i32;
pub type lapack_logical = lapack_int;

include!("lapack.rs");
