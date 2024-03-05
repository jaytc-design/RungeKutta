use ndarray::prelude::*;
use ndarray::Array;

fn main() {
    // println!("Enter x0,y0,xn,h:");
    // let (x0, y0, xn, h) = (inp(), inp(), inp(), inp());
    static g: f64 = 9.8;

    let (const L1: f64, const L2: f64, const theta1: f64, const theta2: f64) = (1.0, 1.0, 0.0, 0.0);
    
    for i in 
        rk4(|_t, y| y.clone(), 
            arr1(&[0.0, 1.0]), 
            0.0, 
            10.0, 
            100.0
        )
    {
        println!("{:?}", i);
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
