use faer::{assert_matrix_eq, dbgf, mat, prelude::*, Side};
use faer::complex_native::c64;
use nalgebra:: DMatrix;
use num_complex::Complex64;
use std::time::Instant;

fn main() {
    let matrix = mat![
        [2.28583329, -0.90628668, -1.71493024],
        [-0.90628668, 4.00729077, 2.17332502],
        [-1.71493024, 2.17332502, 1.97196187]
    ];

    let chol = match matrix.cholesky(Side::Lower) {
        Ok(chol) => chol,
        Err(_) => panic!("Non positive definite matrix"),
    };

    let rhs = mat![
        [-0.29945184, -0.5228196],
        [0.84136125, -1.15768694],
        [1.25678304, -0.46203532]
    ];

    let sol = chol.solve(&rhs);
    assert_matrix_eq!(rhs, &matrix * &sol, comp = abs, tol = 1e-10);

    let mut shah_matrix = mat![
        [c64::new(0.0, 0.0) , c64::new(0.0, -1.0) ],
        [c64::new(0.0, 1.0) , c64::new(0.0, 0.0)]
    ];

    let mut shah_matrix2 = shah_matrix.clone();
    shah_matrix2.write(1,1, c64::new(1.0, 2.0));
    println!("{}", shah_matrix2.read(1,1));

    let temp = shah_matrix * shah_matrix2;
    print!("{:?}", temp);

    dbg!(temp.sum());

    println!("");
    println!("nalgbera now");

    let max : usize = 1000;

    let mut temp_new_nalg = DMatrix::<Complex64>::zeros( max, max);

    let now = Instant::now();
    for i in 0..temp_new_nalg.nrows(){
        for j in 0..temp_new_nalg.ncols() {
            temp_new_nalg[(i,j)] = Complex64::new( i as f64 + j as f64, i as f64 - j as f64);
        }
    }
    let evd = temp_new_nalg.symmetric_eigenvalues();
    println!("{:?}", now.elapsed());
    
    println!("Faers Now");
    let now = Instant::now();
    let mut temp_new = mat::Mat::<c64>::zeros(max,max);
    for i in 0..temp_new.nrows() {
        for j in 0..temp_new.ncols() {
            temp_new.write(i, j , c64::new(i as f64 + j as f64,i as f64 - j as f64));
        }
    }
    let evd2 = temp_new.selfadjoint_eigenvalues(Side::Lower);
    println!("{:?}", now.elapsed());
    std::thread::sleep(std::time::Duration::from_secs(10));
    println!("{:?}, {:?}", evd, evd2);
}