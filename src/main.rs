use ndarray::prelude::*;
use ndarray::Array;

use std::io::{self, Write};
use std::fs::File;

mod rungekutta;
use rungekutta::rk4;

#[allow(non_snake_case)]
fn main() {
    // const GRAVITY: f64 = 9.80665;
    // const LENGTH1: f64 = 1.0;
    // const LENGTH2: f64 = 1.0;
    // const MASS1: f64 = 2.0;
    // const MASS2: f64 = 2.0;

    let mut f = File::create("output.txt").expect("File can't be created!");
    
    let GRAVITY: f64 = input("Gravity: ");
    let LENGTH1 = input("Length 1: ");
    let LENGTH2 = input("Length 2: ");
    let MASS1 = input("Mass 1: ");
    let MASS2 = input("Mass 2: ");

    for (i, v) in 
        rk4(
            |_t, y| { 
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
        arr1(&[
            input("Angle 1: "),
            input("Angle 2: "), 
            input("Angle Velocity 1: "), 
            input("Angle Velocity 2: ")
        ]), // < theta1, theta2, omega1, omega2 >
        input("Start Time: "), 
        input("End Time: "), 
        input("Steps: "),
        ).iter().enumerate()
    {
        // println!("theta1: {}\ntheta2: {}\nomega1: {}\nomega2: {}\n------------", i[0], i[1], i[2], i[3]);
        writeln!(f, "------STEP {}------\ntheta1: {}\ntheta2: {}\nomega1: {}\nomega2: {}\nproduct: {}", i, v[0], v[1], v[2], v[3], v[0] * v[1]).unwrap();
    }

    drop(f);
}

fn input(prompt: &str) -> f64 {
    println!("{}", prompt);

    let mut input = String::new();

    loop {
        match io::stdin().read_line(&mut input) {
            Ok(_) => {
                match input.trim().parse::<f64>() {
                    Ok(num) => return num,
                    Err(_) => {
                        println!("Please enter a valid number...");
                        input.clear();
                    }
                }
            },
            Err(_) => {
                println!("Please enter a valid number...");
                input.clear();
            }
        }

    }
    


}
