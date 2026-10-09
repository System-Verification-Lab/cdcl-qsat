#![allow(unused)]

mod numbers;
mod checker;
mod solver;
mod projectors;
mod utils;

use numbers::{Real, Complex, cross_product, cf};
use numbers::{ComplexFloat, Float};
use crate::utils::benchmark;

fn main() {
    benchmark();
}
