/// Copyright (c) 1996-2015 Timothy A. Davis, Patrick R. Amestoy and Iain S. Duff.
/// Copyright (c) 2011-2021 Richard Lincoln.
/// All Rights Reserved.
extern crate num_traits;

mod aat;
mod amd;
mod amd_1;
mod amd_2;
mod control;
mod dump;
mod info;
mod internal;
mod post_tree;
mod postorder;
mod preprocess;
mod valid;

use num_traits::{NumAssignOps, PrimInt};

use aat::aat;
pub use crate::amd::*;
use amd_1::amd_1;
pub use control::control;
pub use info::info;
use internal::*;
use preprocess::preprocess;
use std::cmp::max;
use std::fmt::Display;
use valid::valid;

pub fn order<I: PrimInt + NumAssignOps + Display>(
    n: I,
    a_p: &[I],
    a_i: &[I],
    control: &Control,
) -> Result<(Vec<I>, Vec<I>, Info), Status> {
    let un = n.to_usize().unwrap();

    let mut info = Info::new(un);

    if n == I::zero() {
        let p = Vec::new();
        let p_inv = Vec::new();
        return Ok((p, p_inv, info));
    }
    let nz = a_p[un].to_usize().unwrap();
    info.nz = nz;

    // Check the input matrix.
    let status = valid(n, n, a_p, a_i);
    if status == Status::Invalid {
        return Err(Status::Invalid);
    }

    let (c_p, c_i) = if status == Status::OkButJumbled {
        // Sort the input matrix and remove duplicate entries.
        debug1_println!("Matrix is jumbled");
        preprocess(n, a_p, a_i) // R = A'.
    } else {
        // Order the input matrix as-is. No need to compute R = A' first.
        (a_p.to_vec(), a_i.to_vec())
    };

    // Determine the symmetry and count off-diagonal nonzeros in A+A'.
    let (nzaat, mut len) = aat(n, &c_p, &c_i, &mut info);
    debug1_print!("nzaat: {}\n", nzaat);
    debug_assert!((max(nz - un, 0) <= nzaat) && (nzaat <= 2 * nz));

    let iwlen = nzaat + (nzaat / 5) + un; // Space for matrix + elbow room.
    debug1_print!("iwlen {}\n", iwlen);

    // Order the matrix.
    let (p, p_inv) = amd_1(n, &c_p, &c_i, &mut len, iwlen, control, &mut info);

    info.status = status;

    Ok((p, p_inv, info))
}
