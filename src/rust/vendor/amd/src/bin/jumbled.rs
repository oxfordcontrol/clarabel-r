extern crate amd;

// A simple test that illustrates the use of the interface to AMD.
//
// Identical to order.rs, except that it operates on an input matrix
// that has unsorted columns and duplicate entries.
fn main() {
    // The symmetric can_24 Harwell/Boeing matrix (jumbled, and not symmetric).
    // Since AMD operates on A+A', only A(i,j) or A(j,i) need to be specified,
    // or both.  The diagonal entries are optional (some are missing).
    // There are many duplicate entries, which must be removed. */
    let n: usize = 24;

    let a_p = vec![
        0, 9, 14, 20, 28, 33, 37, 44, 53, 58, 63, 63, 66, 69, 72, 75, 78, 82, 86, 91, 97, 101, 112,
        112, 116,
    ];

    let a_i = vec![
        0, 17, 18, 21, 5, 12, 5, 0, 13, // column: 0
        14, 1, 8, 13, 17, // column: 1
        2, 20, 11, 6, 11, 22, // column: 2
        3, 3, 10, 7, 18, 18, 15, 19, // column: 3
        7, 9, 15, 14, 16, // column: 4
        5, 13, 6, 17, // column: 5
        5, 0, 11, 6, 12, 6, 23, // column: 6
        3, 4, 9, 7, 14, 16, 15, 17, 18, // column: 7
        1, 9, 14, 14, 14, // column: 8
        7, 13, 8, 1, 17, // column: 9
        // column: 10
        2, 12, 23, // column: 11
        5, 11, 12, // column: 12
        0, 13, 17, // column: 13
        1, 9, 14, // column: 14
        3, 15, 16, // column: 15
        16, 4, 4, 15, // column: 16
        13, 17, 19, 17, // column: 17
        15, 17, 19, 9, 10, // column: 18
        17, 19, 20, 0, 6, 10, // column: 19
        22, 10, 20, 21, // column: 20
        6, 2, 10, 19, 20, 11, 21, 22, 22, 22, 22, // column: 21
        // column: 22
        12, 11, 12, 23, // column: 23
    ];

    let mut p_inv = vec![0; 24];
    let control = amd::Control::default();
    let mut a = [[""; 24]; 24];

    println!("AMD demo, with a jumbled version of the 24-by-24");
    println!("Harwell/Boeing matrix, can_24:");

    amd::control(&control);

    // Print the input matrix.
    let nz = a_p[n];
    println!(
        "\nJumbled input matrix:  {}-by-{}, with {} entries.
    Note that for a symmetric matrix such as this one, only the
    strictly lower or upper triangular parts would need to be
    passed to AMD, since AMD computes the ordering of A+A'.  The
    diagonal entries are also not needed, since AMD ignores them.
    This version of the matrix has jumbled columns and duplicate
    row indices.",
        n, n, nz
    );
    for j in 0..n {
        print!(
            "\nColumn: {}, number of entries: {}, with row indices in
 Ai [{} ... {}]:
    row indices:",
            j,
            a_p[j + 1] - a_p[j],
            a_p[j],
            a_p[j + 1] - 1
        );
        for pj in a_p[j]..a_p[j + 1] {
            let i = a_i[pj as usize];
            print!(" {}", i);
        }
        println!();
    }

    // Print a character plot of the input matrix. This is only reasonable
    // because the matrix is small.
    println!("\nPlot of (jumbled) input matrix pattern:");
    for j in 0..n {
        for i in 0..n {
            a[i][j] = ".";
        }
        for pj in a_p[j]..a_p[j + 1] {
            let i = a_i[pj as usize] as usize;
            a[i][j] = "X";
        }
    }
    print!("   ");
    for j in 0..n {
        print!(" {}", j % 10);
    }
    println!();
    for i in 0..n {
        print!("{}: ", i);
        for j in 0..n {
            print!(" {}", a[i][j]);
        }
        println!();
    }

    // Print a character plot of the matrix A+A'.
    println!("\nPlot of symmetric matrix to be ordered by amd::order:");
    for j in 0..n {
        for i in 0..n {
            a[i][j] = ".";
        }
    }
    for j in 0..n {
        a[j][j] = "X";
        for pj in a_p[j]..a_p[j + 1] {
            let i = a_i[pj as usize] as usize;
            a[i][j] = "X";
            a[j][i] = "X";
        }
    }
    print!("   ");
    for j in 0..n {
        print!(" {}", j % 10);
    }
    println!();
    for i in 0..n {
        print!("{}: ", i);
        for j in 0..n {
            print!(" {}", a[i][j]);
        }
        println!();
    }

    // Order the matrix.
    let (p, _p_inv, info) = amd::order(n, &a_p, &a_i, &control).unwrap();
    println!(
        "return value from amd::order: {:?} (should be {:?})",
        info.status,
        amd::Status::OkButJumbled
    );

    // Print the statistics.
    amd::info(&info);

    if info.status != amd::Status::OkButJumbled {
        println!("AMD failed");
        return;
    }

    // Print the permutation vector, P, and compute the inverse permutation.
    println!("Permutation vector:");
    for k in 0..n {
        // Row/column j is the kth row/column in the permuted matrix.
        let j = p[k];
        p_inv[j as usize] = k as i32;
        print!(" {}", j);
    }
    println!();
    println!();

    println!("Inverse permutation vector:");
    for j in 0..n {
        let k = p_inv[j];
        print!(" {}", k);
    }
    println!();
    println!();

    // Print a character plot of the permuted matrix.
    println!("\nPlot of (symmetrized) permuted matrix pattern:");
    for j in 0..n {
        for i in 0..n {
            a[i][j] = ".";
        }
    }
    for jnew in 0..n {
        let j = p[jnew] as usize;
        a[jnew][jnew] = "X";
        for pi in a_p[j]..a_p[j + 1] {
            let i = a_i[pi as usize];
            let inew = p_inv[i as usize] as usize;
            a[inew][jnew] = "X";
            a[jnew][inew] = "X";
        }
    }
    print!("   ");
    for j in 0..n {
        print!(" {}", j % 10);
    }
    println!();
    for i in 0..n {
        print!("{}: ", i);
        for j in 0..n {
            print!(" {}", a[i][j]);
        }
        println!();
    }
}
