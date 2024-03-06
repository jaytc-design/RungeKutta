use ndarray::prelude::*;
use ndarray::Array;

fn main() {
    const GRAVITY: f64 = 9.80665;
    const LENGTH1: f64 = 1.0;
    const LENGTH2: f64 = 1.0;
    const MASS1: f64 = 2.0;
    const MASS2: f64 = 2.0;

    for i in 
        rk4(|_t, y| { 
            let theta1 = y[0];
            let theta2 = y[1];
            let omega1 = y[2];
            let omega2 = y[3];


            let mut omega1_prime = -GRAVITY * (2.0 * MASS1 + MASS2) * (theta1).sin() - MASS2 * GRAVITY * (theta1 - 2.0 * theta2).sin() - 2.0 * (theta1 - theta2).sin() * MASS2 * (omega2 * omega2 * LENGTH2 + omega1 * omega1 * LENGTH1 * (theta1 - theta2).cos());
            let mut omega2_prime = 2.0 * (theta1 - theta2).sin() * (omega1 * omega1 * LENGTH1 * (MASS1 + MASS2) + GRAVITY * (MASS1 + MASS2) * (theta1).cos()) + omega2 * omega2 * LENGTH2 * MASS2 * (theta1 - theta2).cos();
            omega1_prime /= LENGTH1 * (2.0 * MASS1 + MASS2 - MASS2 * (2.0 * theta1 - 2.0 * theta2).cos());
            omega2_prime /= LENGTH2 * (2.0 * MASS1 + MASS2 - MASS2 * (2.0 * theta1 - 2.0 * theta2).cos());

            // < theta1' ,theta2', omega1', omega2' > (feed into runge kutta) 
            Array::from(vec![ // CHECK USE OF OMEGAS
                omega1, 
                omega2, 
                omega1_prime,
                omega2_prime
            ])
        }, 
            arr1(&[0.994837673637, 0.436332312999, 0.0, 0.0]), // < theta1, theta2, omega1, omega2 >
            0.0, 
            10.0, 
            10000.0,
        )
    {
        println!("theta1: {}\ntheta2: {}\nomega1: {}\nomega2: {}\n------------", i[0], i[1], i[2], i[3]);
    }
}

fn rk4<T, F>(f: F , initial: Array<f64, T>, a: f64, b: f64, n: f64) -> Vec<Array<f64, T>>
where 
    F: Fn(f64, &Array<f64, T>) -> Array<f64, T>, // mutating closure
    T: Dimension // n-dimensions
{
    let h = (b-a)/n; // Step size

    let mut v: Vec<Array<f64, T>> = Vec::new(); // Record of n-dimension arrays
    v.push(initial.clone());

    let mut t: f64 = a; // Time
    let mut y: Array<f64, T> = initial.clone(); // Initial value to be changed

    while t<b-h
    {
        let m1=f(t, &y);
        let m2=f(t+h/2.0, &(&y+&m1*h/2.0));
        let m3=f(t+h/2.0,&(&y+&m2*h/2.0));
        let m4=f(t+h,&(&y+&m3*h));
        let m=(m1+2.0*m2+2.0*&m3+m4)/6.0;

        y = y+h*m;

        v.push(y.clone());

        t+=h;
    }

    v
}