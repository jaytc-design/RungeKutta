use std::io;

fn main() {
    // Take in user input
    let (mut m1, mut m2, mut m3, mut m4, mut m): (f64, f64, f64, f64, f64);

    println!("Enter x0,y0,xn,h:");
    let (x0, y0, xn, h) = (inp(), inp(), inp(), inp());

    let (mut x, mut y) = (x0, y0);

    while x<xn
    {
        m1=f(x0,y0);
        m2=f(x0+h/2.0,y0+m1*h/2.0);
        m3=f(x0+h/2.0,y0+m2*h/2.0);
        m4=f(x0+h,y0+m3*h);
        m=(m1+2.0*m2+2.0*m3+m4)/6.0;
        y+=m*h;
        x+=h;
        println!("{}\n{}\n",x,y);
    }
}

fn f(x: f64, y: f64) -> f64 {
    (x-y)/(x+y)
}

fn inp() -> f64 {
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    input.trim().parse().unwrap()
}
