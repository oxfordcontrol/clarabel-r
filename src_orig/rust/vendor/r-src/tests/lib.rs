extern crate libc;
extern crate r_src;

use libc::{c_double, c_int};

extern "C" {
    pub fn drotg_(a: *mut c_double, b: *mut c_double, c: *mut c_double, s: *mut c_double);
    pub fn dgesv_(
        n: *mut c_int,
        nrhs: *mut c_int,
        a: *mut c_double,
        lda: *mut c_int,
        ipiv: *mut c_int,
        b: *mut c_double,
        ldb: *mut c_int,
        info: *mut c_int,
    );
}

#[test]
fn blas_link() {
    unsafe {
        let mut a: f64 = 0.0;
        let mut b: f64 = 0.0;
        let mut c: f64 = 42.0;
        let mut d: f64 = 42.0;
        drotg_(
            &mut a as *mut _,
            &mut b as *mut _,
            &mut c as *mut _,
            &mut d as *mut _,
        );
        assert!(c == 1.0);
        assert!(d == 0.0);
    }
}

#[test]
fn lapack_test() {
    unsafe {
        let mut n = 3;
        let mut nrhs = 1;
        let mut lda = 3;
        let mut ldb = 3;
        let mut ipiv = vec![1, 2, 3];
        let mut info = 0;

        let mut a = vec![
            6.0f64, -4.0f64, 1.0f64, -4.0f64, 6.0f64, -4.0f64, 1.0f64, -4.0f64, 6.0f64,
        ];

        let mut b = vec![4.0f64, 4.0f64, -1f64];

        dgesv_(
            &mut n as *mut _,
            &mut nrhs as *mut _,
            a.as_mut_ptr(),
            &mut lda as *mut _,
            ipiv.as_mut_ptr(),
            b.as_mut_ptr(),
            &mut ldb as *mut _,
            &mut info as *mut _,
        );

        assert!(info == 0);
        for (one, another) in b.iter().zip(&[3.0, 4.0, 2.0]) {
            assert!((one - another).abs() < 1e-14);
        }
    }
}
