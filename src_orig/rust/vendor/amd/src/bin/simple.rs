extern crate amd;

fn main() {
    let n = 5;
    let a_p = vec![0, 2, 6, 10, 12, 14];
    let a_i = vec![
        0, 1, // 1st column
        0, 1, 2, 4, // 2nd column
        1, 2, 3, 4, // 3rd column
        2, 3, // 4th column
        1, 4, // 5th column
    ];
    let control = amd::Control::default();

    let (p, p_inv, info) = amd::order::<i32>(n, &a_p, &a_i, &control).unwrap();

    println!("P = {:?}", p);
    // Output:
    //   P = [0, 3, 2, 4, 1]

    println!("PInv = {:?}", p_inv);
    // Output:
    //   PInv = [0, 4, 2, 1, 3]

    amd::info(&info);
}
