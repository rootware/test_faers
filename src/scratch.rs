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

    let max : usize = 100;

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
    std::thread::sleep(std::time::Duration::from_secs(1));
    println!("{:?}, {:?}", evd, evd2);


    println!("testing code from discord");

    fn fd2(z: &c64) -> c64 {
        (1.0 - (1.0 / z.powf(2.0) + z.powf(3.0))).sqrt()
    }

    fn fd3(z: &c64) -> c64 {
        (c64::new(1.0, 0.0) - (c64::new(1.0,0.0)/ z.powf(2.0) + z.powf(3.0))).sqrt()
    }

    let mut f2 =c64::new(0.06666666666666665, 0.0);
    println!("{:?}", fd2(&f2));
    let mut f2 =c64::new(0.06666666666666665, 0.0);
    println!("{:?}", fd3(&f2));

    let f3 = c64::new(0.06666666666666665, 0.0);

    let result_1 = 1.0 / f3.powf(2.0) + f3.powf(3.0);
    let result_2 = c64::new(1.0, 0.0)/f3.powf(2.0) + f3.powf(3.0);
    println!( "Firstly, {:?}\t{:?}", result_1, result_2);

    let result_1 = 1.0 - 1.0 / f3.powf(2.0) + f3.powf(3.0);
    let result_2 = c64::new(1.0, 0.0) - c64::new(1.0, 0.0)/f3.powf(2.0) + f3.powf(3.0);
    println!( "Secondly, {:?}\t{:?}", result_1, result_2);

    let result_1 = 1.0 -  ( 1.0 / f3.powf(2.0) + f3.powf(3.0) ); 
    let result_2 = c64::new(1.0, 0.0) - (c64::new(1.0, 0.0)/f3.powf(2.0) + f3.powf(3.0) );
    println!( "What if I add brackets, {:?}\t{:?}", result_1, result_2);

    let result_1 = (1.0 -   1.0 / f3.powf(2.0) ) + f3.powf(3.0) ; 
    let result_2 =  (c64::new(1.0, 0.0) - c64::new(1.0, 0.0)/f3.powf(2.0) ) + f3.powf(3.0) ;
    println!( "What if I add brackets differently, {:?}\t{:?}", result_1, result_2);


    let result_1 = (1.0 - 1.0 / f3.powf(2.0) + f3.powf(3.0) ).sqrt();
    let result_2 = (c64::new(1.0, 0.0) - c64::new(1.0, 0.0)/f3.powf(2.0) + f3.powf(3.0)).sqrt();
    println!( "Thirdly, {:?}\t{:?}", result_1, result_2);

    println!("{:?}", 1.0 / Complex64::new(1e200, 0.0));
    println!("{:?}", Complex64::new(1.0/1e200, 0.0));
    println!("{:?}", std::f32::MAX);

    println!("With floats");
    println!("Add 2nd and 3rd first {:?}", 1.0 - (1.0 / f3.powf(2.0)+ f3.powf(3.0) ) +c64::new(num::zero::<f64>(), num::zero::<f64>() ) );
    println!("Add 1st and 2nd first {:?}", 1.0 -  1.0 / f3.powf(2.0) + f3.powf(3.0));
    println!("{:?}", 1.0 / f3.powf(2.0)+ f3.powf(3.0) );
    println!("With c64");
    println!("Add 2nd and 3rd first {:?}", c64::new(1.0,0.0) - (1.0 / f3.powf(2.0)+ f3.powf(3.0) ) );
    println!("Add 1st and 2nd first {:?}", c64::new(1.0,0.0) -  1.0 / f3.powf(2.0) + f3.powf(3.0));

    println!("{:?}", num::zero::<f64>());

    let test = c64::new(1.0, 0.0);
    println!("{:?}", 0.0 - test);
    println!("{:?}", c64::new(0.0, 0.0) - test);

    println!("{:?}",  0.0 - c64::new(0.0, 0.0).im() );


    let a =  mat![
        [c64::new(0.0, 0.0) , c64::new(0.0, -1.0) ],
        [c64::new(0.0, 1.0) , c64::new(0.0, 0.0)]
    ];

    let b =  mat![
        [c64::new(1.0, 0.0) , c64::new(0.0, 0.0) ],
        [c64::new(0.0, 0.0) , c64::new(1.0, 0.0)]
    ];

    println!("{:?}", a.adjoint() );

    println!("{:?}", a - a.clone().adjoint() );