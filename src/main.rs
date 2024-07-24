
use faer::{assert_matrix_eq, dbgf, linalg, mat, prelude::*, Side};
use faer::complex_native::c64;
use faer::linalg::zip;
use faer::col::*;
use num::*;
const N_STATES : usize = 7;

const N_FINER : usize = 9;
const N_COARSE : usize = 5;
use faer::iter;
fn main () {

    let myvec = vec![0,1];
    println!("{:?}", myvec );
    let max_p = 6.0;
    let mut a = mat::Mat::<c64>::identity(2, 1);
    a.write(0, 0, c64::new(0.0,-1.0));
    a.write(1, 0, c64::new(0.0, 2.0));

    println!("{:?}", a.adjoint());
    let q=0.0;
    println!("{:?}", a.adjoint()*&a);

    let fam = mat::Mat::<c64>::from_fn( N_STATES, 1, |i,_j| -> c64 {c64::new((i as f64 * 2.0 - max_p + q).powi(2), 0.0)} );

    let mut h0new = mat::Mat::<c64>::identity(N_STATES, N_STATES)*fam;

    h0new[(1,0)] = c64::new(1.0,2.0);
    println!("{:?}", h0new);

    let out = mat![ [c64::from(2.0), c64::from(1.0)]];
    let matr: Mat<f64> = zipped!(&out).map( |unzipped!( x)| (*x).norm_sqr() );
    println!("{:?}", matr);

    let mut matrix = mat! [[c64::from(2.0), c64::from(1.0)], [c64::from(3.0),c64::from( 4.0)]];

    let mut matrix2 = mat::Mat::<c64>::with_capacity(2, 2);
    unsafe { matrix2.set_dims(2, 2);}
    for i in 0..2 {
        for j in 0..2 {
            matrix2[(i,j)]= matrix[(i,j)];
        }
    }
    println!("{:?}", &matrix.row_capacity());

    let wow = matrix.as_ptr_mut();
    let wow2 = matrix2.as_ptr_mut();
    unsafe { 
        for i in 0..matrix.row_capacity()*matrix.col_capacity()  {
            println!("{:?}", *(wow.wrapping_add(i))); 
        }
    }
    println!("Now the ohter one");
    unsafe { 
        for i in 0..matrix2.row_capacity()*matrix2.col_capacity()  {
            println!("{:?}", *(wow2.wrapping_add(i))); 
        }
    }

    let psi = mat![ [ c64::from(2.0)],[ c64::from(3.0)]];
    let target = mat! [[c64::from(3.0)], [c64::from(-2.0)]];

    println!("{:?}", target.adjoint().shape());
    println!("{:?}", &psi.shape());
    println!("{:?}", (target.adjoint()*psi)[(0,0)]);



    let delta_p = 0.01;
    let max_coarse = N_COARSE as f64 - 1.0;
    let mut ones =  mat::Mat::<f64>::zeros(1,N_FINER); ones.fill(1.0) ;
    let momentum_states = mat::Mat::<f64>::from_fn( N_COARSE, 1, |i,_j|    2.0 * i as f64 - max_coarse) * ones.clone();
    
    println!("Momentum States are {:?}", momentum_states);

    ones = mat::Mat::<f64>::zeros(N_COARSE, 1); ones.fill(1.0) ;
    let finer = ones * mat::Mat::<f64>::from_fn( 1, N_FINER, |_i,j|    delta_p* (j as f64 - (N_FINER as f64 - 1.0)/2.0 ) ) ;
    println!("The finer momentum is {:?}", finer);

    let mut p = momentum_states + finer;
    println!("Our convolved momentum States are {:?}", p);

    let h0  = mat::Mat::<f64>::from_fn( p.nrows(), p.ncols(), |i,j| p[(i,j)].powi(2));
    println!("{:?}", h0);

    let mut h2 = mat::Mat::<f64>::zeros( N_COARSE , N_COARSE );

    let depth = 10.0;
    for i in 0..N_COARSE {
        if i + 1 < N_COARSE  {

            h2.write(i, i + 1,  (depth / 4.0) );
            h2.write(i + 1, i ,  (depth / 4.0) );

        }
    }

    println!("Operating h2 on momentum{:?}", h0+ h2 * p);

/* 


   // in place mutation
let mut mat: Mat<f64> = todo!();
zipped!(&mut mat).for_each(|unzipped!(mut x)| *x = 2.0 * (*x) + 1.0);

// map + collect into a new matrix
let mat: Mat<f64> = todo!();
let new_mat = zipped!(&mat).map(|unzipped!(x)| (*x) * (*x) + 1.0);

*/


}