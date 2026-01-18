mod numbers;
mod checker;

pub use numbers::{Real, Complex, cross_product};

#[cfg(test)]
pub use numbers::{f, cf};

fn main() {
    println!("Hello, world!");
}

